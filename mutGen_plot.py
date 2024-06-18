#!/bin/python3    

# Written by Elton J. F. Chaves
# -----------------------------
import pandas as pd
import os
import glob
import sys

# ---
arg1 = sys.argv[1]
arg2 = sys.argv[2]

# ---
df             = pd.DataFrame(columns=['pdb', 'mutation', 'dGpred_wt', 'dGpred_mt', 'ddGpred'])
pdb_list       = glob.glob(f'{arg1}/*')
pdbs           = []
mutants        = []
dGpred_wt_list = []
dGpred_mt_list = []

# ---
for item in pdb_list:
    pdb = os.path.basename(item)

    if item.__contains__('_wt_'):
        condition1 = os.path.isfile(f'{item}/dG_pred.csv')
    else:
        condition1 = False
        
    if condition1:
        dGpred_wt = pd.read_csv(f'{item}/dG_pred.csv')['dG_pred'][0]
        dGpred_wt_list.append(dGpred_wt)        
    else:
        dGpred_mt = pd.read_csv(f'{item}/dG_pred.csv')['dG_pred'][0]
        mutchain  = os.path.basename(item).split('_')[1]
        mutinfos  = os.path.basename(item).split('_')[2]
        mutant    = mutinfos[0] + mutchain + mutinfos[1:]
        dGpred_mt_list.append(dGpred_mt)
        mutants.append(mutant)
        pdbs.append(pdb)

# ---
df['pdb'] = pdbs
df['dGpred_wt'] = dGpred_wt
df['dGpred_mt'] = dGpred_mt_list
df['mutation']  = mutants
df['ddGpred']   = df['dGpred_mt'] - df['dGpred_wt']
df.to_csv('data.csv', index=False)