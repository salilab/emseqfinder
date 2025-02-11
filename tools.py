import pandas as pd
import numpy as np
import os

residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS',
 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
features = ["EMDB", "resolution", "chain", "resid", "resname", "ss", "se_ccc", "res_ccc"]
all_columns = features+residues
one_to_three = { "A":"Ala", "R":"Arg", "N":"Asn", "D":"Asp",
        "C":"Cys", "E":"Glu", "Q":"Gln", "G":"Gly",
        "H":"His", "I":"Ile", "L":"Leu", "K":"Lys",
        "M":"Met", "F":"Phe", "P":"Pro", "S":"Ser",
              "T":"Thr", "W":"Trp", "Y":"Tyr", "V":"Val"}

emdb_header_home = "/wynton/group/sali/saltzberg/database/emdb_headers"

h_residue_propensities={'ALA': 0.09743660383807928, 'ARG': 0.052206653248949665, 'ASN': 0.033516232631581634, 'ASP': 0.03697014685512692, 'CYS': 0.014795201253000493, 'GLN': 0.04344432285169852, 'GLU': 0.06778753070641984, 'GLY': 0.03754920042140815, 'HIS': 0.01862409730140631, 'ILE': 0.07440455692198593, 'LEU': 0.14053655562613038, 'LYS': 0.05310967070473625, 'MET': 0.030689839012904986, 'PHE': 0.05122200709659481, 'PRO': 0.017251714840352126, 'SER': 0.05328313168494384, 'THR': 0.04712016509403881, 'TRP': 0.016744086383568144, 'TYR': 0.038248146135774035, 'VAL': 0.07506013739129991}
s_residue_propensities={'ALA': 0.06003398703023674, 'ARG': 0.044378467067739666, 'ASN': 0.03278576451285257, 'ASP': 0.03541292288460036, 'CYS': 0.018204547230252364, 'GLN': 0.03322525197281038, 'GLU': 0.04347995937182592, 'GLY': 0.04848035002734589, 'HIS': 0.017638096726306743, 'ILE': 0.09542737713883898, 'LEU': 0.10406086413001016, 'LYS': 0.04244472224392531, 'MET': 0.02161301664192515, 'PHE': 0.05539495273068208, 'PRO': 0.019483944058129542, 'SER': 0.06141104773810454, 'THR': 0.0748593640128135, 'TRP': 0.015216032502539261, 'TYR': 0.045706695835612154, 'VAL': 0.13074263614344872}


#database_home = "/wynton/group/sali/saltzberg/database/em_residue_densities"
#database_home = "/wynton/group/sali/saltzberg/database/em_residue_densities_lr"
#database_home = "/wynton/group/sali/dibyendu/database/em_residue_density/EMDB"
database_home="/wynton/group/sali/saltzberg/database/em_residue_density/EMDB"

def read_sequences(sequence_file):
    import IMP
    import IMP.pmi
    import IMP.pmi.topology
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
    return int(x*10**3)/10**3


def process_input_database(db, reslow=0.0, reshigh=10.0, normalize=True, threshold=True, nvoxels=2744):

    df = import_voxel_file(db)
    print("Imported voxel file of length", len(df))

    # Remove unknown residues
    unknown_resis = set(df["resname"].values)-set(residues)
    df = df[df["resname"]!="UNK"]

    all_voxel_cols = df.columns[-1*nvoxels:]
    # Get all EMDBs
    emdbs = df.EMDB.unique()

    for emdb in emdbs:
        print(emdb)
        # Get threshold and min/max stats
        xml = os.path.join(emdb_header_home, "headers", emdb+".xml")
        thresh = get_threshold_from_xml(xml)
        (vmin,vmax,vavg,vstd) = get_statistics_from_xml(xml)
        if thresh is None or vavg is None or vstd==0:
            index_names = df[ df['EMDB'] == emdb ].index
            df.drop(index_names, inplace = True)
            continue

        # Threshold voxels at user-defined thresholding value
        df_vals = df.EMDB[all_voxel_cols]

        df_vals.loc[df_vals<float(thresh)] = 0

        # Now scale voxels by mean and std of map
        #df.loc[(df.EMDB==emdb)][all_voxel_cols] = df.loc[(df.EMDB==emdb)][all_voxel_cols] / float(vstd) - float(vagv)
        df_vals = df_vals / float(vstd) - float(vavg)
        for v in all_voxel_cols:
            df.EMDB[v] = df_vals[v]

    df_x = df[(df["resolution"]<reshigh) & (df["resolution"]>reslow)]

    print("Filtered database #entries:", len(df_x))
    df_h = df_x[df_x["ss"]=="H"]
    df_s = df_x[df_x["ss"]=="S"]
    df_h = df_h.dropna()
    df_s = df_s.dropna()

    return df_h, df_s

def import_voxel_file(fname):
    f = open(fname, "r")
    columns = ["EMDB", "resolution", "chain", "resid", "resname", "ss", "se_ccc", "res_ccc"]
    data = []
    for line in f.readlines():
        try:
            int(line[0])
            dlist = []
            fields = line.split()
            for i in range(len(fields)):
                if i == 3:
                    dlist.append(int(fields[i]))
                elif i in [0,2,4,5]:
                    dlist.append(fields[i])
                else:
                    dlist.append(float(fields[i]))
            data.append(dlist)
        except:
            continue

    for i in range(len(data[0])-len(columns)):
        columns.append("v"+str(i))
    database = pd.DataFrame(data, columns=columns)
    return database

def get_xml_line(xml, field):
    f = open(xml, "r")

    for l in f.readlines():
        if field in l:
            f.close()
            return l
    f.close()
    return None

def get_feature_points(fdict, feature):
    cdict = features["Coordinates"]

    coords = []

    for v in fdict[feature]:
        coords.append(cdict[v])

    return coords

def get_threshold_from_xml(xml):
    l = get_xml_line(xml, "contourLevel")
    if l is None:
        return l
    return float(l.strip().split(">")[1].split("<")[0])

def get_statistics_from_xml(xml):
    # Returns voxel minimum, maximum, average and SD as a tuple
    l = get_xml_line(xml, "<minimum>")
    minimum = float(l.strip().split(">")[1].split("<")[0])
    l = get_xml_line(xml, "<maximum>")
    maximum = float(l.strip().split(">")[1].split("<")[0])
    l = get_xml_line(xml, "<average>")
    average = float(l.strip().split(">")[1].split("<")[0])
    l = get_xml_line(xml, "<std>")
    std = float(l.strip().split(">")[1].split("<")[0])
    return (minimum, maximum, average, std)


def make_simulated_scores(sequence, sens=0.05, acc=0.5):
    pass

def get_auc(pcts, binwidth=0.05):
    histbins = np.arange(0,1 + binwidth,binwidth)
    
    hh, hb = np.histogram(pcts, bins=histbins, density=True)
    print(hh)    
    auch = 0
    sumh = 0
    for i in range(len(hb)-1):
        sumh+=hh[i]/len(histbins)
        auch+=(sumh-hb[i+1])*binwidth
        
    return auch

def score_prior_over_against_sequences(sedf, sequences, prior, unavailable=None):
    se_len = len(sedf)
    scores = {}

    # Loop over all sequence chains
    for c in sequences.keys():
        seq = sequences[c]
        scores[c] = {}

        # Loop over all starting residues
        for r in range(len(seq)-se_len):
            try:
                segment = [one_to_three[s].upper() for s in seq[r:r+se_len]]
                score = score_prior_sequence(prior, segment)
                scores[c][r+1] = score
            except:
                continue
        if unavailable is not None:
            for r in range(len(seq)-se_len, len(seq)):
                scores[c][r+1] = unavailable
    return scores

def score_sesf_over_against_sequences(sedf, sequences, count=False, unavailable=None, log=True):
    se_len = len(sedf)
    scores = {}
    skeys = list(sequences.keys())
    # Loop over all sequence chains 
    for c in range(len(skeys)):
        if count:
            if c%100==0:
                print("Doing sequence", c, "out of", len(skeys))
        sk = skeys[c]
        seq = sequences[sk]
        scores[sk] = {}

        # Loop over all starting residues
        for r in range(len(seq)-se_len):
            try:
                segment = [one_to_three[s].upper() for s in seq[r:r+se_len]]
                score = score_sequence(sedf, segment, log)
                scores[sk][r+1] = score
            except:
                continue   
        if unavailable is not None:
            for r in range(len(seq)-se_len, len(seq)):
                scores[sk][r+1] = unavailable
    return scores

def NestedDictValues(d):
  for v in d.values():
    if isinstance(v, dict):
      yield from NestedDictValues(v)
    else:
      yield v

def catstring(stuff, delimiter=" "):
    outstring = ""
    for s in stuff:
        outstring+=str(s)+delimiter

    return outstring[0:-1]

def get_se_ss(df_se):
    ss = np.unique(df_se["ss"].values)
    if len(ss) > 1:
        print("Multiple secondary structures in this element!  That's bad!!", df_se["emdb"].values[0])
        #print(df_se["resid"])
    return ss[0]

def get_se_sequence(sedf):
    seq = []
    for i in range(len(sedf)):
        seq.append(sedf.iloc[i]["resname"])
    return seq

def score_prior_sequence(prior, seq):
    score = 0
    for i in range(len(seq)):
        score+=-1*prior[seq[i]]
    return score

def score_sequence(sedf, seq, log=True):
    if len(seq) != len(sedf):
        raise Exception("Input sequence and structure element must be the same size")
    score = 0
    for i in range(len(seq)):
        if not log:
#            print("seq[i] : ", seq[i])		
#            score+=-1*np.log(sedf.iloc[i][seq[i]])
            score+=-1*np.log(sedf.iloc[i]["p_" + seq[i]])
        else:
#            score+=-1*sedf.iloc[i][seq[i]]
            score+=-1*sedf.iloc[i]["p_" + seq[i]]
    return score

def get_emdbs(df):
    return list(df["EMDB"].unique())

def get_correct_scores(df):
    correct_scores = {}
    for i, row in df.iterrows():
        correct_scores[i] = row[row["resname"]]
    return correct_scores

def get_emdb_SEs_from_db(df, emdb):
    # Finds StructureElements (contiguous sequence elements) in a database
    # Returns a list of DFs with those entries in it
    df_emdb = df[df["EMDB"]==emdb]
    if len(df_emdb) == 0:
        print("No entries for EMDB", emdb, "found")

    chains = np.unique(df_emdb["chain"].values)

    # Dictionary of structure elements
    ses = {}
    for c in chains:
        df_chain = df_emdb[df_emdb["chain"]==c]
        sorted_resids = np.sort(df_chain["resid"].values)
        rid_0 = sorted_resids[0]
        ss = df_chain[df_chain["resid"]==rid_0]["ss"].values[0]
        rid_1 = rid_0

        for i in range(1,len(sorted_resids)):
            ridx = sorted_resids[i]
            this_ss = df_chain[df_chain["resid"]==ridx]["ss"].values[0]
            if ridx - rid_1 > 1 or this_ss != ss:
                if rid_1-rid_0 > 3:
                    seid = ss+"_"+c+"_"+str(int(rid_0))+"_"+str(int(rid_1-rid_0)+1)

                    idx = np.where((df_chain['resid']>=rid_0) & (df_chain['resid']<= rid_1))
                    if len(idx[0]) == rid_1-rid_0+1:
                        ses[seid] = df_chain.iloc[idx[0]]
                rid_0 = ridx
                rid_1 = rid_0
                ss = this_ss
            else:
                rid_1 = sorted_resids[i]

# for last seid
        if rid_1-rid_0 > 3:
            seid = ss+"_"+c+"_"+str(int(rid_0))+"_"+str(int(rid_1-rid_0)+1)
            idx = np.where((df_chain['resid']>=rid_0) & (df_chain['resid']<= rid_1))
            if len(idx[0]) == rid_1-rid_0+1:
                ses[seid] = df_chain.iloc[idx[0]]

    return ses
