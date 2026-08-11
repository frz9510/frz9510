# Fix and simplify patient haploblock export section in your Python script
# This ensures the haploblock .tsv output matches the recombination tree/haplo graph

import pandas as pd
import numpy as np
import argparse
import re
import os

# --- ARGUMENTS ---
parser = argparse.ArgumentParser()
parser.add_argument("--hpt_file", required=True)
parser.add_argument("--blocks_file", required=True)
parser.add_argument("--patients_file", required=True)
parser.add_argument("--out_dir", required=True)
parser.add_argument("--prefix", default="")
args = parser.parse_args()

# --- LOAD ---
patients = pd.read_csv(args.patients_file, header=None)[0].tolist()
hpt = pd.read_csv(args.hpt_file, sep='\t', header=None, skiprows=5)
pos = pd.read_csv(args.hpt_file, sep='\t', header=None, skiprows=1, nrows=1)
hpt.columns = ['H17N'] + pos.iloc[0, 1:-1].tolist() + ['sample']
hpt = hpt[hpt['sample'].isin(patients)].copy()
print (hpt)

# --- BLOCKS ---
def num2acgt(seq_list):
    d = {'1': 'A', '2': 'C', '3': 'G', '4': 'T'}
    return [''.join([d.get(c, '0') for c in seq]) for seq in seq_list]

def read_blocks(path):
    blocks = []
    with open(path) as f:
        markers, blockname, seqs = [], None, []
        for line in f:
            if line.startswith("BLOCK"):
                if blockname:
                    blocks.append({'blockname': blockname, 'markers': markers, 'seq': num2acgt(seqs)})
                parts = line.split()
                blockname = parts[1]
                markers = [int(x) for x in parts[3:]]
                seqs = []
            elif re.match("^[0-4]", line):
                seqs.append(line.strip().split()[0])
        blocks.append({'blockname': blockname, 'markers': markers, 'seq': num2acgt(seqs)})
    return pd.DataFrame(blocks)

blocks = read_blocks(args.blocks_file)
print (blocks)


# --- HAPLOBLOCK ASSIGNMENT ---
rows = []
for _, row in hpt.iterrows():
    haplo = row['H17N']
    patient = row['sample']
    snps = row[1:-1]  # Exclude 'haplo' and 'ind'
    # print (snps)
    # print (snps.iloc[2])
    # snp_dict = snps.to_dict()
    # print (snp_dict)
    for block in blocks.itertuples():
        # print (block.markers)
        try:
            # seq = ''.join([str(snp_dict.get(p, '0')) for p in block.markers])
            seq = ''.join([str(snps.iloc[p-1]) for p in block.markers])

        except:
            continue
        color = '99'
        for i, known in enumerate(block.seq):
            # print (i)
            # print (known)
            # print (seq)
            if known == seq:
                color = str(i)
                break
        rows.append({"patient": patient, "haplo": haplo, "start": min(block.markers), "end": max(block.markers), "color": color})
        # break

# --- SAVE ---
os.makedirs(args.out_dir, exist_ok=True)
outfile = os.path.join(args.out_dir, f"{args.prefix}haploblock_segments.tsv")
pd.DataFrame(rows).to_csv(outfile, sep='\t', index=False)
print(f"Saved: {outfile}")
