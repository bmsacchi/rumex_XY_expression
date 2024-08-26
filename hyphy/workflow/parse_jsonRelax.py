import json
import pandas as pd
import os
import re

file_path = snakemake.input[0]
output_file = snakemake.output[0]

with open(file_path, 'r') as file:
    data = json.load(file)

file_name = os.path.basename(file_path)

# Get the LRT, p-value, and relaxation or intensification parameter
lrt = data.get("test results", {}).get("LRT")
p_value = data.get("test results", {}).get("p-value")
relaxation_parameter = data.get("test results", {}).get("relaxation or intensification parameter")


# Extract nodes of interest from the "tested" section using the same patterns
#patterns = [r'Node\d+$',r'^Rsag_hap1_\d+_RA_1$',r'^TX_maternal_\d+_RA_1$',r'^TX_paternal_\d+_RA_1$',r'^buc_\d+_RA_1$']

patterns = [r'Node\d+$', r'^Rsag_hap1_\d+_RA_1$', r'^TX_maternal_\d+_RA_1$', r'^TX_paternal_\d+_RA_1$', r'^buc_\d+_RA_1$']
tested_nodes = data.get("tested", {}).get("0", {})


branch_attributes = data.get("branch attributes", {}).get("0", {})

nodes_of_interest = []

for pattern in patterns:
    nodes_of_interest.extend([node for node in branch_attributes.keys() if re.match(pattern, node)])

rows = []
# get values
for node in nodes_of_interest:
    row = {
        "OG": file_name,  # Use the file name for the OG column
        "branch": node,
        "LRT": lrt,
        "p-value": p_value,
        "relaxation or intensification parameter": relaxation_parameter
    }
    rows.append(row)


# dataframe w cols
#df = pd.DataFrame(rows, columns=["OG", "branch", "dN", "dS", "LB", "MLE", "UB"])
df = pd.DataFrame(rows, columns=["OG", "branch", "LRT", "p-value", "relaxation or intensification parameter"])


# Save the DataFrame to a CSV file if needed
df.to_csv(output_file, index=False)

