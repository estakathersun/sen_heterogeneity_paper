import sys, os, logging
import scanpy as sc
import numpy as np
import pandas as pd
from tqdm import tqdm

from DeepScence.api import DeepScence


# logging config
# ---------------------------
logger = logging.getLogger("pipeline")
logger.setLevel(logging.INFO)
# ❗ Отключаем проброс в root (это убирает дубли)
logger.propagate = False
# вывод в консоль:
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)

out_dir = '../.'
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    logger.info(f'создаю папку {out_dir}')
else:
    logger.info(f'{out_dir} уже существует')

# вывод в файл:
file_handler = logging.FileHandler(os.path.join(out_dir, '2_log_run_ds_only.txt'))
file_handler.setLevel(logging.INFO)
# 4. Формат логов
formatter = logging.Formatter(
    "%(asctime)s | %(levelname)s | %(message)s")
# 5. Подключить формат к handlers
console_handler.setFormatter(formatter)
file_handler.setFormatter(formatter)
# 6. “Прикрутить” handlers к logger
logger.addHandler(console_handler)
logger.addHandler(file_handler)


# ловит все ошибки:
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    logger.error("Uncaught exception",
                 exc_info=(exc_type, exc_value, exc_traceback))


sys.excepthook = handle_exception


# adata = sc.read_h5ad("../adata_all_harm_ribo.h5ad")
# obs = pd.read_csv(
#     ".././0_adata obs - reduced and filt.csv",
#     index_col=0,
# )
# adata = adata[obs.index].copy()
# adata.obs = obs.loc[adata.obs_names].copy()

adata = sc.read_h5ad('.././0.1_adata_downsampled.h5ad')

logger.info("AnnData loaded: %s cells, %s genes", adata.n_obs, adata.n_vars)

# adata = adata[adata.obs.Tissue_global.isin(['Heart_Cell', 'Heart_Nuclei'])]
sc.pp.filter_genes(adata, min_cells=1)



res_obs =[]
grouped = adata.obs.groupby(["Tissue_global", "cell_type"], observed=True)

for (t, ct), idx in tqdm(grouped.groups.items(), total=len(grouped)):
    adata_sub = adata[idx].copy()
    sc.pp.filter_genes(adata_sub, min_counts=1)

    n_batches = adata_sub.obs["GSE_id"].nunique()
    if n_batches >= 2:
        adata_sub.obs["b"] = adata_sub.obs["GSE_id"]
    adata_sub = DeepScence(adata_sub, binarize=False, dropout=True, verbose=True)  # предполагается: returns AnnData. binarize=True будет падать с ошибкой из-за внутр ошибки библиотеки

    res_obs.append(adata_sub.obs)

df = pd.concat(res_obs, axis=0, ignore_index=False)
df.to_csv(os.path.join(out_dir, '2_dw_obs_ds_only.csv'))

logger.info("Pipeline finished successfully")
