import pandas as pd
import sys
import time
import numpy as np
import os

# usage: python convert_MLDB_topkl.py 10k.dat


def import_voxel_file_pickle(fname, other_columns, n_voxels):
    print("Extracting data from "+fname, flush=True)
    columns = [oc[0] for oc in other_columns]
    dtype = [oc[1] for oc in other_columns]
    for i in range(n_voxels):
        dtype.append("f")
        columns.append("v"+str(i))

    dtypes = {}
    for d in range(len(dtype)):
        if dtype[d] == "f":
            dtypes[columns[d]] = np.float32
        elif dtype[d] == "s":
            dtypes[columns[d]] = str
        elif dtype[d] == "i":
            dtypes[columns[d]] = np.int16
    database = pd.read_csv(fname, sep=" ",  names=columns,
                           na_values=["None"], dtype=dtypes)
    database.drop(database[database["resname"] == "UNK"].index, inplace=True)
    return database


fname = sys.argv[1]
if sys.argv[2]:
    basename = sys.argv[2]
else:
    basename = os.path.basename(fname).split('.')[0]

n_voxels = 2744
# for parts fitting we are using this
other_columns = [("EMDB", "s"), ("resolution", "f"), ("pdbname", "s"),
                 ("chain", "s"), ("resid", "i"), ("resname", "s"), ("ss", "s")]

# read_pickle
start = time.time()
f_name = basename + '.pkl'
print(f_name)
df = import_voxel_file_pickle(fname, other_columns, n_voxels)
df.to_pickle(f_name)
end = time.time()
print('for pickle convert time', end-start)
