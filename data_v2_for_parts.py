import sys
import pandas as pd
import numpy as np
import sklearn
from sklearn.preprocessing import LabelBinarizer
import time
# This script contains all important scripts for preparing input file and running the ML model

residues = ["GLY","ALA", "THR", "HIS", "SER", "ASP", "ASN", "GLU", "GLN", "VAL", "LEU", "ILE", "LYS", "ARG", "TRP", "TYR", "PHE", "CYS", "MET", "PRO"]

def split_df_random(df_in, x=0.8,y=0.1,z=0.1):
    shuffled = df_in.sample(frac=1)
    n_arr = len(df_in)
    nx = int(x*n_arr)
    ny = int(y*n_arr)
    return shuffled.iloc[0:nx], shuffled.iloc[nx:nx+ny], shuffled.iloc[nx+ny:]

def import_voxel_file_pandas(fname, other_columns, n_voxels, read_rows='all'):
    print("Extracting data from "+fname, flush=True)
    data = []
    columns = [oc[0] for oc in other_columns]
    dtype = [oc[1] for oc in other_columns]
    for i in range(n_voxels):
        dtype.append("float32")
        columns.append("v"+str(i))
    dtypes = {}
    for col in columns:
        dtypes[col] = dtype[columns.index(col)]
    if read_rows == 'all':
        database = pd.read_csv(fname, sep=" ", names=columns, na_values=["None"], dtype=dtypes)
    elif type(read_rows) == int or type(read_rows) == float:
        database = pd.read_csv(fname, sep=" ", names=columns, nrows=read_rows,na_values=["None"], dtype=dtypes)
    else:
        print("Processing all lines in the file", fname)
        database = pd.read_csv(fname, sep=" ", names=columns, na_values=["None"], dtype=dtypes)
    #print("Extracted "+str(len(data))+" lines in "+str(time.time()-t0), flush=True)
    t0 = time.time()
    print("Extracted "+str(len(database))+" lines in "+str(time.time()-t0), flush=True)
    #database = pd.DataFrame(data, columns=columns)

    # Drop all instances of "UNK"
    database.drop(database[database["resname"]=="UNK"].index, inplace=True)
    print("Processed database in "+str(time.time()-t0), flush=True)

    return database

def import_voxel_file_pickle(fname, other_columns, n_voxels):
    print("Extracting data from "+fname, flush=True)
    data = []
    columns = [oc[0] for oc in other_columns]
    dtype = [oc[1] for oc in other_columns]
    for i in range(n_voxels):
        dtype.append("f")
        columns.append("v"+str(i))

    dtypes = {}
    for d in range(len(dtype)):
        if dtype[d]=="f":
            dtypes[columns[d]]=np.float32
        elif dtype[d]=="s":
            dtypes[columns[d]]=str
        elif dtype[d]=="i":
            dtypes[columns[d]]=np.int16
    t0 = time.time()
    database = pd.read_pickle(fname)
#    print("Database size before datatypes change (MB):", str(sys.getsizeof(database)/1024/1024), len(database))
#    print("changing data types now ...", flush=True)
#    database = database.astype(dtypes)
    #print("Extracted "+str(len(data))+" lines in "+str(time.time()-t0), flush=True)
    print("Extracted "+str(len(database))+" lines in "+str(time.time()-t0), flush=True)
    t0 = time.time()
    #database = pd.DataFrame(data, columns=columns)
    print("Database size (MB):", str(sys.getsizeof(database)/1024/1024), len(database))
    
    print("Processed database in "+str(time.time()-t0), flush=True)
    return database


def import_voxel_file(fname, other_columns, n_voxels):
    print("Extracting data from "+fname, flush=True)
    t0 = time.time() 
    #f = open(fname, "r")
    data = []
    columns = [oc[0] for oc in other_columns]
    dtype = [oc[1] for oc in other_columns]
    for i in range(n_voxels):
        dtype.append("f")
        columns.append("v"+str(i))
    dtypes = {}
    for d in range(len(dtypes)):
        if dtype[d]=="f":
            dtypes[columns[d]]=np.float32
        elif dtype[d]=="s":
            dtypes[columns[d]]=str
        elif dtype[d]=="i":
            dtypes[columns[d]]=np.int16

    k = 0
    l = 0
    for line in open(fname, "r"):
        k+=1
        if k%100000==0:
            print("-- Line "+str(k)+" entries: "+str(l)+" in "+str(time.time()-t0), flush=True)
        try:
            int(line[0])
        except:
            print(k, line[0], line.split()[0:8], flush=True)
            continue
        dlist = []
        fields = line.split()
        for i in range(len(fields)):
            if dtype[i] == "i":
                dlist.append(int(fields[i]))
            elif dtype[i] == "s":
                dlist.append(fields[i])
            else:
                if fields[i] == "None":
                   dlist.append(None)
                else:
                    try:
                        dlist.append(float(fields[i]))
                    except:
                        dlist.append(None)
            data.append(dlist)
        l+=1
    print("Extracted "+str(len(data))+" lines in "+str(time.time()-t0), flush=True)
    t0 = time.time()
    database = pd.DataFrame(data, columns=columns)
    
    # Drop all instances of "UNK"
    database.drop(database[database["resname"]=="UNK"].index, inplace=True)
    print("Processed database in "+str(time.time()-t0), flush=True)

    return database

def reshape_df(df, dims):
    outlist = []
    for i, row in df.iterrows():
        outlist.append(np.reshape(row.values, dims))
    return np.array(outlist)

def split_image_and_other_features(df, image_prefix="v", other_columns=["resolution"]):
    # Split a dataframe into two dataframes: one with image information, one with other features
    # other_columns :: column labels to extract and return in other_features - these should be ones used for learning
    vcols = [c for c in df.columns if c[0]==image_prefix]
    print("SIOF",other_columns)
    print(df.columns)
    image_features = df[vcols]
    other_features = df[other_columns]
    return image_features, other_features


class ResidueVoxelDataset():
    def __init__(self, infile, classes=residues, target=None, voxel_dim=(14,14,14), 
                 other_columns=[("EMDB", "object"),
                                ("resolution", "float32"),
                                ("pdbname","object"),
                                ("chain", "object"),
                                ("resid", "int16"),
                                ("resname", "object"),
                                ("ss", "object")],
                 train_features = ["resolution", "H", "S"], binarize_ss=True):
        
        self.infile = infile
        self.target = target
        self.voxel_dim = voxel_dim
        self.other_columns = other_columns
        self.n_voxels = voxel_dim[0]*voxel_dim[1]*voxel_dim[2]
        self.train_features = train_features
        self.data_df = import_voxel_file_pickle(self.infile, other_columns=other_columns, n_voxels = self.n_voxels)
        self.pred_classes = classes
        self.zb = self.get_binarizer()

    def get_voxel_columns(self):
        # Return a numpy array of the voxel values only
        vcols = []
        for i in range(self.n_voxels):
            vcols.append("v"+str(i))
        return vcols

    def drop_non_residues(self):
        self.data_df = self.data_df[self.data_df[self.target] in residues]

    def binarize_ss(self, df, add_to_training_features=True):
#        ss_cats = self.data_df["ss"].unique()
        ss_cats = ["H", "S"]
        for c in ss_cats:
            bc = self.data_df["ss"]==c
#            print(bc)
            self.data_df[c] = bc.astype(int)
            
        if add_to_training_features:
            self.train_features = [y for x in [self.train_features, ss_cats] for y in x]
#        print(self.data_df[c])
        return ss_cats

    def split_features_and_reshape(self, df):
        # First, split the image voxels from the other features
        ximages, xother = split_image_and_other_features(df, other_columns=self.train_features)

        # Reshape the images to a 14x14x14 box
        ximages_reshape = reshape_df(ximages, self.voxel_dim)

        if self.target is not None:
            y = df[self.target]
        else:
            y = None 

        return {"image_features": ximages_reshape,
                "other_features": xother, 
                "target": y}
    
    def get_train_test_val_sets(self, train=0.75, test=0.15, val=0.1):
        df_train, df_test, df_val = split_df_random(self.data_df, train, test, val)

        train_test_val = [] 
        for dft in [df_train, df_test, df_val]:
            if len(dft) > 0:
                train_test_val.append(self.split_features_and_reshape(dft))
            else:
                train_test_val.append(None)
        
        return train_test_val

    def get_binarizer(self):
        return LabelBinarizer().fit(self.pred_classes) 

    def get_weights_by_class(self, ser):
        class_weights = {}
        for i in range(len(self.pred_classes)):
            res = self.pred_classes[i]
            res_frac = sum(x==res for x in ser)/len(ser)
            class_weights[i] = (1/len(self.pred_classes)) / res_frac
        
        return class_weights
       

