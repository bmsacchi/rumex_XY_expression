import csv
from Bio import SeqIO

# Step 1: Parse the CSV file to extract the pgID and gene IDs
pangene_dict = {}
with open("../orths.csv", mode="r") as infile:
    reader = csv.DictReader(infile)
    for row in reader:
        pgID = row["pgID"]
        pangene_dict[pgID] = {
            "TX_maternal": row["id_tx_mat"],
            "TX_paternal": row["id_tx_pat"],
            "buc": row["id_bucephalophorus"],
            "Rsag": row["id_sagittatus"],
        }

# Step 2: Parse the fasta files for each species to retrieve sequences
def get_gene_sequences(fasta_file):
    gene_sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        gene_sequences[record.id] = str(record.seq)
    return gene_sequences

tx_maternal_seqs = get_gene_sequences("../nuc_fastas/transcript_fastas/TX_maternal.all.maker.transcripts.fasta")
tx_paternal_seqs = get_gene_sequences("../nuc_fastas/transcript_fastas/TX_paternal.all.maker.transcripts.fasta")
buc_seqs = get_gene_sequences("../nuc_fastas/transcript_fastas/buc.all.maker.transcripts.fasta")
rsag_seqs = get_gene_sequences("../nuc_fastas/transcript_fastas/sagittatus.transcripts.fasta")

# Step 3: Create a new fasta file for each pgID and write the corresponding sequences
for pgID, genes in pangene_dict.items():
    filename = f"../pgIDfastas/pgID{pgID}.fasta"
    with open(filename, "w") as outfile:
        for species, gene_id in genes.items():
            seq = None
            if species == "TX_maternal":
                seq = tx_maternal_seqs.get(gene_id)
            elif species == "TX_paternal":
                seq = tx_paternal_seqs.get(gene_id)
            elif species == "buc":
                seq = buc_seqs.get(gene_id)
            elif species == "Rsag":
                seq = rsag_seqs.get(gene_id)
            if seq:
                outfile.write(f">{gene_id}\n{seq}\n")
                print(f"Written {gene_id} to {filename}")

            else:
                print(f"Gene {gene_id} not found in {species} fasta file.")

print("Fasta files created successfully.")
