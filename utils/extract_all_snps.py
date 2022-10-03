import pandas as pd
import numpy as np

import sys

processedfiles = sys.argv[1]
outfile = sys.argv[2]


f = open(processedfiles, 'r')

lines = f.readlines()
num_files = len(lines)

count = 0
assert(len(lines) >= 2)
df_all_pos = pd.read_csv(lines[0].strip(), sep=" ")[["pos"]]
print(df_all_pos.shape)
for idx in range(1, len(lines)):
    count += 1
    df_next = pd.read_csv(lines[idx].strip(), sep = " ")[["pos"]]
    print(df_next.shape)
    df_all_pos = df_all_pos.merge(df_next, how="outer", on="pos")
    print(df_all_pos)

print(df_all_pos.shape)
#df_all.sort_values(by="pos", inplace=True)

f.close()
f = open(processedfiles, 'r')
lines = f.readlines()
print(len(lines))

for idx in range(len(lines)):
    df = pd.read_csv(lines[idx].strip(), sep=" ")[["rsid","pos"]]
    print(df)
    df_all_pos = df_all_pos.merge(df, on="pos", how="outer")
    col_label = "study" + str(idx)
    df_all_pos.rename(columns={"rsid":col_label}, inplace=True)

df_all_pos.sort_values("pos", inplace=True, ignore_index=True)
print(df_all_pos)

def mk_idx(col):
    counter = 0
    newcol = col.copy()
    print(newcol)
    for i in range(len(col)):
        print("this")
        print(type(newcol[i]))
        if newcol[i] == pd.NA:
            newcol[i] = -1
        else:
            newcol[i] = counter
            counter += 1
    return newcol

print("num files", num_files)
for i in range(num_files):
    col_label = "study" + str(i)
    print("fixing", col_label)
    na_idxs = df_all_pos[col_label].isna()
    not_na = df_all_pos[col_label].notna()
    idx = np.arange(np.sum(not_na))
    df_all_pos[col_label][not_na] = idx
    df_all_pos[col_label][na_idxs] = -1

print(df_all_pos)
df_all_pos.to_csv(outfile, sep=',', header=False, index=False)
