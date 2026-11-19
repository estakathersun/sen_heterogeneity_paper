import sys
import os
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from scipy import sparse
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

import h5py

from SenCID.api import SenCID


output_dir = "../."
# os.makedirs(output_dir, exist_ok=True)
# path = ".././0.1_adata_downsampled.h5ad"

path = "../adata_all_harm_ribo.h5ad"
global_obs = pd.read_csv(
    ".././0_adata obs - reduced and filt.csv",
    index_col=0,
)

# Для obs и var есть важный нюанс: это HDF5-группы, которые нужно корректно преобразовать в pandas.DataFrame.
# Самый простой вариант — использовать внутренние функции старого anndata только для этих групп, не вызывая чтение всего объекта:
#
with h5py.File(path, "r") as f:
    # X
    X_group = f["X"]
    X = sparse.csr_matrix(
        (X_group["data"][:], X_group["indices"][:], X_group["indptr"][:]),
        shape=tuple(X_group.attrs["shape"]),
    )

    # OBS
    obs_data = {}
    for col in f["obs"].keys():
        obj = f["obs"][col]
        if isinstance(obj, h5py.Dataset):
            values = obj[:]
            # bytes -> str
            if values.dtype.kind == "S":
                values = values.astype(str)
            obs_data[col] = values
        elif isinstance(obj, h5py.Group):
            # categorical: codes + categories
            codes = obj["codes"][:]
            categories = obj["categories"][:]

            if categories.dtype.kind == "S":
                categories = categories.astype(str)

            # просто восстанавливаем значения
            values = categories[codes]
            obs_data[col] = values

    obs = pd.DataFrame(obs_data)

    # индекс клеток
    index_name = f["obs"].attrs["_index"]
    if isinstance(index_name, bytes):
        index_name = index_name.decode()
    obs.index = obs[index_name].astype(str)
    obs = obs.drop(columns=index_name)

    # VAR
    var_data = {}
    for col in f["var"].keys():
        obj = f["var"][col]
        if isinstance(obj, h5py.Dataset):
            values = obj[:]
            if values.dtype.kind == "S":
                values = values.astype(str)
            var_data[col] = values
        elif isinstance(obj, h5py.Group):
            codes = obj["codes"][:]
            categories = obj["categories"][:]
            if categories.dtype.kind == "S":
                categories = categories.astype(str)

            values = categories[codes]
            var_data[col] = values
    var = pd.DataFrame(var_data)

    # индекс генов
    index_name = f["var"].attrs["_index"]
    if isinstance(index_name, bytes):
        index_name = index_name.decode()
    var.index = var[index_name].astype(str)
    var = var.drop(columns=index_name)

# СОЗДАЁМ НОВЫЙ nnData
adata = ad.AnnData(X=X, obs=obs, var=var)

# для данных без даунсемплинга:
adata = adata[global_obs.index].copy()
adata.obs = global_obs.loc[adata.obs_names].copy()


all_dfs = []

for t in tqdm(set(adata.obs.Tissue_global)):
    adata_t = adata[adata.obs.Tissue_global == t]
    for ct in set(adata_t.obs.cell_type):
        adata_ct = adata_t[adata_t.obs.cell_type == ct]
        sc.pp.filter_genes(adata_ct, min_cells=1)

        # Run SenCID
        pred_dict_ind, recSID_ind, _ = SenCID(
            adata=adata_ct,
            # sidnums=[1, 2, 3, 4, 5, 6],
            denoising=False,  # DCA
            binarize=False,
            threads=16,
        )

        # Combine output
        for key in pred_dict_ind.keys():
            df = pred_dict_ind[key]
            df.columns = [f"{key}_{col}" for col in df.columns]
            recSID_ind = pd.concat([recSID_ind, df], axis=1)
            # recSID_ind['Tissue'] = t
            # recSID_ind['cell_type'] = ct


        all_dfs.append(recSID_ind)

merged_df = pd.concat(all_dfs)
merged_df.to_csv(f"{output_dir}/3_full_SenCID_all_res.csv")
