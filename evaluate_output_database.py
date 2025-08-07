import pandas as pd
import numpy as np
from tools import *
import sys
import IMP
import IMP.pmi
import IMP.pmi.topology
import os
import math

def read_sequences(sequence_file):
    sequences = IMP.pmi.topology.Sequences(sequence_file)
    seqs = {}
    for s in sequences.sequences.keys():
        s_cs = s.split("|")[-1].strip()
        if "Chains" in s_cs:
            chains = s_cs.split()[1].split(",")
        else:
            chains = s_cs.split()
        for c in chains:
            seqs[c]=sequences[s]
    return seqs

def tint(x, n=3):
    if not math.isnan(x):
    	return int(x*10**3)/10**3
    else:
        return x

priorh = {}
priors = {}
for r in residues:
    priorh[r] = np.log(h_residue_propensities[r])
    priors[r] = np.log(s_residue_propensities[r])

dbh_file = sys.argv[1]
out_file = sys.argv[2]

# Open output file. Line buffered, so we actually see the printout before the system
# buffering kicks in
of = open(out_file, "w", buffering=1)

# 1 - open scoring function database
dtypes = {'EMDB': 'object'}
dfs = pd.read_csv(dbh_file, delim_whitespace=True, dtype=dtypes)
dfs.resid = dfs.resid.astype(int)
# 2 - Get list of EMDBs
emdbs = get_emdbs(dfs)

print("Total number of systems:", len(emdbs))

# 3 - Cycle through EMDBs
for emdb in emdbs:
    sedfs = get_emdb_SEs_from_db(dfs, emdb)
#    print(emdb)
    database_home_ = "."
    fasta_file = f"{emdb}.fasta"
#    seqs = read_sequences(os.path.join(database_home, emdb, "0system", "native.fasta"))
    seqs = read_sequences(os.path.join(database_home_, emdb, "0system", fasta_file))
#    print(seqs)

    if len(sedfs) == 0:
        print("EMDB", emdb, "has no structure elements")
        continue

    seids = list(sedfs.keys())
    resolution = float(sedfs[seids[0]]["resolution"].values[0])
    print(resolution)
    for seid in seids:
        print(seid)
        sedf = sedfs[seid]
#        se_ccc = sedf["se_ccc"].values[0]
        length = int(seid.split("_")[-1])
        ss = get_se_ss(sedf)
        seq = get_se_sequence(sedf)
        print(seq)
        all_score_dict = score_sesf_over_against_sequences(sedf, seqs, log=False)
        if ss=="H":
            prior = priorh
        elif ss=="S":
            prior = priors
        prior_score_dict = score_prior_over_against_sequences(sedf, seqs, prior)
        actual_prior_score = score_prior_sequence(prior, seq)
        all_scores = list(NestedDictValues(all_score_dict))
        all_prior_scores = list(NestedDictValues(prior_score_dict))
        actual_score = score_sequence(sedf, seq, log=False)
        act_rank = sum(np.abs(all_scores) < actual_score)
        act_pct = act_rank / float(len(all_scores))
        prior_rank = sum(np.abs(all_prior_scores) < actual_prior_score)
        prior_pct = prior_rank / float(len(all_prior_scores))
#        outstring = catstring([emdb, resolution, seid, se_ccc, "|", act_rank, tint(act_pct,3), tint(actual_score,2), tint(min(all_scores)), tint(np.average(all_scores)), tint(max(all_scores)), "|", prior_rank, tint(prior_pct,3), tint(actual_prior_score,2), tint(min(all_prior_scores)), tint(np.average(all_prior_scores)), tint(max(all_prior_scores)), "|"] + seq)
        outstring = catstring([emdb, resolution, seid, "|", act_rank, tint(act_pct,3), tint(actual_score,2), tint(min(all_scores)), tint(np.average(all_scores)), tint(max(all_scores)), "|", prior_rank, tint(prior_pct,3), tint(actual_prior_score,2), tint(min(all_prior_scores)), tint(np.average(all_prior_scores)), tint(max(all_prior_scores)), "|"] + seq)
        of.write(outstring+"\n")
        #print(outstring)

