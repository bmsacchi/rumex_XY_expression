#!/usr/bin/env python3
"""Simulate paired-end reads from gametolog transcripts and measure kallisto assignment."""

import argparse
import csv
import gzip
import shutil
import subprocess
import tempfile
from pathlib import Path


def normalize_id(identifier):
    return identifier[:-3] if identifier.endswith("-RA") else identifier


def read_fasta(path):
    sequences = {}
    identifier = None
    chunks = []

    def save_record():
        if identifier is None:
            return
        normalized = normalize_id(identifier)
        if normalized in sequences:
            raise ValueError(
                f"Duplicate FASTA IDs after removing -RA suffix: {normalized}"
            )
        sequences[normalized] = (identifier, "".join(chunks))

    with open(path) as fasta:
        for line in fasta:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                save_record()
                identifier = line[1:].split()[0]
                chunks = []
            else:
                if identifier is None:
                    raise ValueError(f"Sequence found before FASTA header in {path}")
                chunks.append(line)
    save_record()
    return sequences


def read_pairs(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    pairs = []
    with opener(path, "rt", newline="") as table:
        reader = csv.DictReader(table)
        required = {"id_tx_mat", "id_tx_pat"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(
                f"{path} must contain columns id_tx_mat and id_tx_pat"
            )
        seen = set()
        for row in reader:
            mat = normalize_id(row["id_tx_mat"])
            pat = normalize_id(row["id_tx_pat"])
            if not mat or not pat or (mat, pat) in seen:
                continue
            seen.add((mat, pat))
            pairs.append((mat, pat))
    return pairs


def count_fastq_pairs(path):
    with open(path) as fastq:
        lines = sum(1 for _ in fastq)
    if lines % 4:
        raise ValueError(f"Malformed FASTQ file: {path}")
    return lines // 4


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Simulate each maternal/paternal gametolog transcript separately, "
            "quantify against the full kallisto index, and summarize cross-assignment."
        )
    )
    parser.add_argument("--transcripts", required=True, help="tx_nomatpar transcript FASTA")
    parser.add_argument("--index", required=True, help="Same kallisto index used for real reads")
    parser.add_argument("--pairs", required=True, help="pgMatOrths CSV or CSV.GZ")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--read-length", type=int, default=150)
    parser.add_argument("--fragment-mean", type=float, default=300)
    parser.add_argument("--fragment-sd", type=float, default=50)
    parser.add_argument(
        "--coverage",
        type=float,
        default=100,
        help="ART fold coverage per source transcript; increase for stable estimates",
    )
    parser.add_argument(
        "--art-profile",
        default="HS25",
        help="ART sequencing profile; choose one matching the real sequencer",
    )
    parser.add_argument(
        "--limit-pairs",
        type=int,
        help="Optional limit for a small test run; pairs are taken in table order",
    )
    parser.add_argument(
        "--bias",
        action="store_true",
        help="Pass kallisto --bias (leave off to match the current Snakefile)",
    )
    args = parser.parse_args()

    transcripts_path = Path(args.transcripts).resolve()
    index_path = Path(args.index).resolve()
    pairs_path = Path(args.pairs).resolve()
    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    work_dir = out_dir / "work"
    work_dir.mkdir(exist_ok=True)
    output_path = out_dir / "gametolog_simulation_summary.csv"

    for executable in ("art_illumina", "kallisto"):
        if shutil.which(executable) is None:
            raise RuntimeError(f"Required executable not found on PATH: {executable}")
    for path in (transcripts_path, index_path, pairs_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.read_length <= 0 or args.coverage <= 0:
        raise ValueError("read length and coverage must be greater than zero")
    if args.fragment_mean <= 0 or args.fragment_sd < 0:
        raise ValueError("fragment mean must be positive and SD cannot be negative")

    sequences = read_fasta(transcripts_path)
    pairs = read_pairs(pairs_path)
    if args.limit_pairs is not None:
        pairs = pairs[: args.limit_pairs]
    if not pairs:
        raise ValueError("No gametolog pairs found in the pair table")

    missing = sorted(
        {transcript for pair in pairs for transcript in pair}
        - sequences.keys()
    )
    if missing:
        preview = ", ".join(missing[:10])
        raise ValueError(
            f"{len(missing)} pair-table transcript IDs are absent from the FASTA "
            f"(examples: {preview})"
        )

    columns = [
        "pair_id",
        "source_sex",
        "source_transcript",
        "counterpart_transcript",
        "simulated_read_pairs",
        "estimated_source",
        "estimated_counterpart",
        "estimated_other_targets",
        "fraction_assigned_to_pair",
        "fraction_of_pair_assigned_to_counterpart",
    ]
    with output_path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=columns)
        writer.writeheader()

        for pair_number, (mat, pat) in enumerate(pairs, start=1):
            pair_id = f"{mat}__{pat}"
            for source_sex, source, counterpart in (
                ("mat", mat, pat),
                ("pat", pat, mat),
            ):
                original_id, sequence = sequences[source]
                with tempfile.TemporaryDirectory(
                    prefix="sim_", dir=work_dir
                ) as temporary:
                    temporary = Path(temporary)
                    source_fasta = temporary / "source.fa"
                    source_fasta.write_text(f">{original_id}\n{sequence}\n")
                    prefix = temporary / "reads"
                    subprocess.run(
                        [
                            "art_illumina",
                            "-ss",
                            args.art_profile,
                            "-i",
                            str(source_fasta),
                            "-p",
                            "-l",
                            str(args.read_length),
                            "-f",
                            str(args.coverage),
                            "-m",
                            str(args.fragment_mean),
                            "-s",
                            str(args.fragment_sd),
                            "-o",
                            str(prefix),
                        ],
                        check=True,
                    )
                    read1 = Path(f"{prefix}1.fq")
                    read2 = Path(f"{prefix}2.fq")
                    simulated_pairs = count_fastq_pairs(read1)
                    quant_dir = temporary / "quant"
                    command = [
                        "kallisto",
                        "quant",
                        "-i",
                        str(index_path),
                        "-o",
                        str(quant_dir),
                        "-t",
                        "1",
                        str(read1),
                        str(read2),
                    ]
                    if args.bias:
                        command.insert(2, "--bias")
                    subprocess.run(command, check=True)

                    estimates = {}
                    with (quant_dir / "abundance.tsv").open() as abundance:
                        reader = csv.DictReader(abundance, delimiter="\t")
                        for row in reader:
                            target = normalize_id(row["target_id"])
                            estimates[target] = (
                                estimates.get(target, 0.0)
                                + float(row["est_counts"])
                            )

                    source_count = estimates.get(source, 0.0)
                    counterpart_count = estimates.get(counterpart, 0.0)
                    all_estimates = sum(estimates.values())
                    pair_estimates = source_count + counterpart_count
                    writer.writerow(
                        {
                            "pair_id": pair_id,
                            "source_sex": source_sex,
                            "source_transcript": source,
                            "counterpart_transcript": counterpart,
                            "simulated_read_pairs": simulated_pairs,
                            "estimated_source": source_count,
                            "estimated_counterpart": counterpart_count,
                            "estimated_other_targets": max(
                                0.0, all_estimates - pair_estimates
                            ),
                            "fraction_assigned_to_pair": (
                                pair_estimates / all_estimates
                                if all_estimates > 0
                                else ""
                            ),
                            "fraction_of_pair_assigned_to_counterpart": (
                                counterpart_count / pair_estimates
                                if pair_estimates > 0
                                else ""
                            ),
                        }
                    )
                    output.flush()
                print(
                    f"[{pair_number}/{len(pairs)}] "
                    f"{source_sex} source {source}",
                    flush=True,
                )

    print(f"Summary written to {output_path}")


if __name__ == "__main__":
    main()
