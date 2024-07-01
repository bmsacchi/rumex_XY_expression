import pandas as pd
import sys

# Reading the input file with specified data types
input_file = sys.argv[1]
pgMat_tx = pd.read_csv(input_file, delimiter='\t', dtype={'pgID': str, 'genome': str, 'id': str, 'chr': str, 'start': int, 'end': int, 'flag': str}, low_memory=False)

# Filtering for PASS
pgMat_PASS_ks = (pgMat_tx[pgMat_tx['flag'] == 'PASS']
                 .loc[:, ['pgID', 'genome', 'id', 'chr', 'start', 'end', 'flag']]
                 .groupby(['pgID', 'genome', 'flag'])
                 .filter(lambda x: len(x) == 1)
                 .drop_duplicates()
                 .pivot(index='pgID', columns='genome', values=['id', 'chr', 'start', 'end', 'flag'])
                 .reset_index())

# Flattening the MultiIndex columns
pgMat_PASS_ks.columns = ['_'.join(col).strip() if type(col) is tuple else col for col in pgMat_PASS_ks.columns.values]
pgMat_PASS_ks['pairID'] = pgMat_PASS_ks['id_tx_mat'].astype(str) + '_' + pgMat_PASS_ks['id_tx_pat'].astype(str)
pgMat_PASS_ks = pgMat_PASS_ks.dropna().drop_duplicates()

# Selecting pgID
pgIDs_PASS = pgMat_PASS_ks[['pgID_']]

# Processing the second part
pgMat_NSOrths_ks = (pgMat_tx[pgMat_tx['flag'] != 'array']
                    .loc[:, ['pgID', 'genome', 'id', 'chr', 'start', 'end', 'flag']]
                    .groupby(['pgID', 'genome'])
                    .filter(lambda x: len(x) == 1)
                    .loc[~pgMat_tx['pgID'].isin(pgIDs_PASS['pgID_'])]
                    .drop_duplicates()
                    .pivot(index='pgID', columns='genome', values=['id', 'chr', 'start', 'end', 'flag'])
                    .reset_index())

# Flattening the MultiIndex columns
pgMat_NSOrths_ks.columns = ['_'.join(col).strip() if type(col) is tuple else col for col in pgMat_NSOrths_ks.columns.values]
pgMat_NSOrths_ks = (pgMat_NSOrths_ks
                    .dropna()
                    .query('chr_tx_mat == "X" & chr_tx_pat == "Y"'))
pgMat_NSOrths_ks['pairID'] = pgMat_NSOrths_ks['id_tx_mat'].astype(str) + '_' + pgMat_NSOrths_ks['id_tx_pat'].astype(str)

# Combining the results
pgMatOrths_ks = pd.concat([pgMat_PASS_ks, pgMat_NSOrths_ks], ignore_index=True)
pgMatOrths_ks = pgMatOrths_ks.loc[:, ['pgID_', 'id_tx_mat', 'id_tx_pat']]

# Writing the output to a gzipped CSV file
pgMatOrths_ks.to_csv("orthsUpdatedAnnos_py.csv", index=False)

