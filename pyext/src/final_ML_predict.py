from tensorflow import keras
import tensorflow as tf
import pandas as pd
import sys
import os
import time
from .data_v2_for_parts import (reshape_df, split_image_and_other_features,
                                ResidueVoxelDataset)

# Non-voxel columns in the input dataframe
other_columns = ["EMDB", "resolution", "pdbname", "chain", "resid",
                 "resname", "ss"]
data_feats = ["EMDB", "resolution", "pdbname", "chain", "resid",
              "resname", "ss"]
other_features = ["resolution"]

# Learning parameters
param_BATCH_SIZE = 64
param_MAX_EPOCHS = 100
param_VAL_SPLIT = 0.25
param_LOG_LRATE = -4  # -5  # -6 --> -1 (DIB : changed from -4 to -5)
# 1--> MAX_EPOCHS, but should be much less than MAX_EPOCHS
param_PATIENCE = 10  # 5
param_DROPOUT = 0.0  # 0.1  # 0-->1


def send_chunk(df, counter, chunk_size=20000):
    df_chunk = df[counter:counter+chunk_size]
    return df_chunk


param_LRATE = 2*10**param_LOG_LRATE

# Label prefix for output/plot files
prefix = ("test_nn_"+str(param_BATCH_SIZE)+"_"+str(param_LOG_LRATE)
          + "_" + str(param_DROPOUT)+"_")

# Set dynamic memory growth
gpus = tf.config.experimental.list_physical_devices('GPU')
for gpu in gpus:
    tf.config.experimental.set_memory_growth(gpu, True)


t0 = time.time()
###############################
# Import input data
spl_input_file = sys.argv[1]
print("Begin", flush=True)
indat = ResidueVoxelDataset(
    spl_input_file, target="resname", train_features=["resolution"])
print("Time to open data:"+str(time.time()-t0), flush=True)
print("Database Size:", len(indat.data_df))

zb = indat.get_binarizer()
pred_classes = ["p_"+c for c in zb.classes_]   # Make column labels

print(indat.binarize_ss(indat.data_df))

##################################
# Seup multi-GPU strategy
strategy = tf.distribute.MirroredStrategy()
print("Strategy:")
print('Number of devices: {}'.format(strategy.num_replicas_in_sync))

with strategy.scope():
    # Load checkpoint:
    from . import get_data_path
    checkpoint_path = get_data_path(
        "finalmodel_no_amino_weights/finalmodel.tfl")
#    checkpoint_path = None
    test_mode = False
    evaluate_mode = True
    if checkpoint_path is not None and test_mode:
        # Load model:
        print("Loading Model .....")
        model = keras.models.load_model(checkpoint_path)
        model.summary()
        train_data, test_data, val_data = indat.get_train_test_val_sets(
            0.01, 0.98, 0.01)
        score = model.evaluate(
            [test_data["image_features"], test_data["other_features"]],
            zb.transform(test_data["target"]))
        print("test results: ", score)

    elif checkpoint_path is not None and evaluate_mode:
        # Load model:
        print("Loading Model .....")
        model = keras.models.load_model(checkpoint_path)
        model.summary()
        df_pred = pd.DataFrame(columns=pred_classes)
        data_set_size = indat.data_df.shape[0]
        chunk_size = int(sys.argv[2])
        for row_num in range(0, data_set_size, chunk_size):
            print(indat.data_df["H"].head(5))
            print(indat.data_df["S"].head(5))
            Ximageall, Ximageother = split_image_and_other_features(
                indat.data_df[row_num:row_num+chunk_size],
                other_columns=["resolution", "H", "S"])
            Xinall = reshape_df(Ximageall, (14, 14, 14))
            probsall = model.predict([Xinall, Ximageother])
            df_pred_temp = pd.DataFrame(data=probsall, columns=pred_classes)
            print("shape of the chunk of the current prediction database is: ",
                  df_pred_temp.shape)
            df_pred = df_pred.append(df_pred_temp)
            print("shape of the current prediction database is: ",
                  df_pred.shape)

    else:
        print("No model to load, please train your model first")
        sys.exit()

# df_pred = pd.DataFrame(data=probsall, columns=pred_classes)
df_oc = indat.data_df[data_feats]
print("---")
# print(len(Xinall), len(df_oc), len(indat.data_df))
print("OC columns", len(df_oc.columns), len(data_feats))


def get_correct_probs(df, probs):
    prob_list = []
    resnames = ["p_"+res for res in df["resname"].values]
    i = 0
    for j, row in probs.iterrows():
        # print(j, i, resnames[i], row[resnames[i]])
        prob_list.append(row[resnames[i]])
        i += 1
    return pd.Series(prob_list)


# Use the same prefix as the input .pkl
base_name = os.path.splitext(sys.argv[1])[0]  # e.g., "3j5r_ML_side"
output_file = f"{base_name}_ML_prob.dat"
df_oc["c_prob"] = get_correct_probs(indat.data_df, df_pred)
# df_oc["vavg"] = np.average(Ximageall, axis=1)
# df_oc["vstd"] = np.std(Ximageall, axis=1)
df_oc.reset_index(drop=True, inplace=True)
df_pred.reset_index(drop=True, inplace=True)
df_out = pd.concat([df_oc, df_pred], axis=1)
df_out.to_csv(output_file, sep=" ", index=False)
