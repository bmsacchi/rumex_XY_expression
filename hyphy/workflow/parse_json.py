import json
import pandas as pd
import os
import re

file_path = snakemake.input[0]
output_file = snakemake.output[0]

with open(file_path, 'r') as file:
    data = json.load(file)

file_name = os.path.basename(file_path)

patterns = [r'Node\d+$',r'^Rsag_hap1_\d+_RA_1$',r'^TX_maternal_\d+_RA_1$',r'^TX_paternal_\d+_RA_1$',r'^buc_\d+_RA_1$']

branch_attributes = data.get("branch attributes", {}).get("0", {})

nodes_of_interest = []

for pattern in patterns:
    nodes_of_interest.extend([node for node in branch_attributes.keys() if re.match(pattern, node)])

rows = []

for node in nodes_of_interest:
    branch_info = branch_attributes[node]
    confidence_intervals = branch_info.get("Confidence Intervals", {})

    row = {
        "OG": file_name,  # Use the file name for the OG column
        "branch": node,
        "dN": branch_attributes[node].get("dN"),
        "dS": branch_attributes[node].get("dS"),
        "LB": confidence_intervals.get("LB"),
        "MLE": confidence_intervals.get("MLE"),
        "UB": confidence_intervals.get("UB")
    }
    rows.append(row)

df = pd.DataFrame(rows, columns=["OG", "branch", "dN", "dS", "LB", "MLE", "UB"])


# Save the DataFrame to a CSV file if needed
df.to_csv(output_file, index=False)

