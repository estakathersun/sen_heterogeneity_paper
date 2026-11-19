import pandas as pd
import scanpy as sc
import decoupler as dc


#  ------------ load & prep adata ------------ 
# adata = sc.read_h5ad("../adata_all_harm_ribo.h5ad")
# obs = pd.read_csv(".././0_adata obs - reduced and filt.csv",index_col=0)
# adata = adata[obs.index].copy()
# adata.obs = obs.loc[adata.obs_names].copy()

# dw
adata = sc.read_h5ad(".././0.1_adata_downsampled.h5ad")

sc.pp.filter_genes(adata, min_cells=1)
sc.pp.normalize_total(adata,target_sum=1e4)
sc.pp.log1p(adata)

# ------------- gene sets ----------------
# prep cellage sigs
cellage = pd.read_csv(".././_gene_sets/cellage3.tsv",sep="\t",)
cellage = cellage[(cellage["Senescence Effect"] == "Induces") & (cellage["Cancer Cell"] == "No")]
cellage_genes = list(set(cellage["Gene symbol"]) & set(adata.var_names))
# Convert gene sets to the network format expected by decoupler (source-target weight)
ca = pd.DataFrame({
    "source": ["cellage"] * len(cellage_genes),
    "target": cellage_genes
})

# fridman_senescence_up sig
fr = dc.pp.read_gmt(
    ".././_gene_sets/FRIDMAN_SENESCENCE_UP.v2026.1.Hs.gmt"
)
# тут фильтрация не нужна, все гены из сигнатуры есть в данном adata

# prep senmayo sig
filepath = '.././_gene_sets/SAUL_SEN_MAYO.v2026.1.Hs.txt'
with open(filepath, 'r') as f:
    senmayo_genes = [line.strip() for line in f]
sm_genes = list(set(senmayo_genes) & set(adata.var_names))
sm = pd.DataFrame({
    "source": ["senmayo"] * len(sm_genes),
    "target": sm_genes
})

net = pd.concat([ca, fr, sm]) # n_genes = 222, 77, 124 resp.

# ---------run -----------------
dc.mt.aucell(
    data=adata,
    net=net,
    # tmin= фильтр на слишком маленькие сигнатуры из net
    # layer=
    raw=False, # if True, игнорирует layers
    empty=True, # удаляет пустые клетки и гены
    verbose=True,
    n_up=1000
)

adata.obsm['score_aucell'].to_csv('.././1_dw_all_aucell_scores_1000.csv')
