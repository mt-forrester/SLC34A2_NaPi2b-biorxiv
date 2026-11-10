"""Surfaceome screen for AT2-cell surface markers in public lung scRNA-seq data.

For every gene in the SURFY in silico human surfaceome (Bausch-Fluck et al., PNAS 2018),
and in each of 4 human and 5 mouse lung scRNA-seq datasets, compare AT2 cells against all
other annotated lung cells and compute:

  auroc      area under the ROC curve, gene expression as the sole classifier of AT2 vs rest
  tau        tissue-specificity index (Yanai et al. 2005) over the per-cell-type mean expression
  log2fc_ct  log2( AT2 mean / median of the other cell types' means ); ranked only among
             genes detected in >= 20% of AT2 cells

Outputs (in --out):
  ranking_<dataset>.csv               every detected surface gene with its metrics and ranks
  fig_tau_vs_auroc_<dataset>.png/.svg tau (x) vs AUROC (y, log(1-AUROC) scale)
  fig_tau_vs_log2fc_<dataset>.png/.svg tau (x) vs log2FC (y)
  slc34a2_summary_{human,mouse}.csv   SLC34A2 / Slc34a2 metrics and ranks per dataset
  fig_slc34a2_across_{human,mouse}_{auroc,log2fc}.png/.svg
  fig_dotplot_{hlca,hurskainen_normoxia}.png/.svg  marker dot plots (+ dotplot_values.csv)

Usage:
  bash download_data.sh                 # SURFY, HLCA annotations, GEO files (Census is pulled here)
  python surfaceome_screen.py           # all 9 datasets
  python surfaceome_screen.py --datasets hlca hurskainen_normoxia   # Figure 1 datasets
"""
import argparse
import gzip
import os
import shutil
import subprocess

import anndata as ad
import cellxgene_census as cc
import matplotlib
import numpy as np
import pandas as pd
import scipy.io
from scipy.sparse import csr_matrix
from sklearn.metrics import roc_auc_score

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CENSUS_VERSION = "2025-11-08"     # pinned CELLxGENE Census release
DETECT_FLOOR = 20                 # % of AT2 cells; log2FC is ranked only above this
CENSUS_AT2 = "pulmonary alveolar type 2 cell"   # Cell Ontology label used by Census

LUNG = "tissue_general=='lung'"
PRIMARY_LUNG = "is_primary_data==True and tissue_general=='lung'"

# One row per dataset. AT2 label and all other cell-type labels are used verbatim from
# each source. HLCA core is an integrated object (is_primary_data is False), so it takes
# the lung filter only; Tabula Sapiens is multi-organ, hence the tissue filter. In HLCA
# core the lung filter keeps the lung-parenchyma samples only (nose and airway samples are
# tissue_general 'nose' / 'respiratory system'), and cells are labelled with HLCA's own
# ann_finest_level, with AT2 = HLCA ann_level_3 'AT2' (= 'AT2' + 'AT2 proliferating').
# `keep` subsets cells after loading; `summary=False` keeps a dataset out of Fig S1.
DATASETS = [
    # ---- human (CELLxGENE Census) ----
    dict(key="hlca", species="human", source="census", at2="AT2", hlca_labels=True,
         dataset_id="066943a2-fdac-4b29-b348-40cede398e4e", obs_filter=LUNG,
         title="HLCA core (Sikkema 2023), lung parenchyma — healthy"),
    dict(key="tsapiens", species="human", source="census", at2=CENSUS_AT2,
         dataset_id="53d208b0-2cfd-4366-9866-c3c6114081bc", obs_filter=PRIMARY_LUNG,
         title="Tabula Sapiens lung (2022) — healthy"),
    dict(key="krasnow", species="human", source="census", at2=CENSUS_AT2,
         dataset_id="8c42cfd0-0b0a-46d5-910c-fc833d83c45e", obs_filter=PRIMARY_LUNG,
         title="Travaglini/Krasnow (2020) — healthy"),
    dict(key="bpd", species="human", source="census", at2=CENSUS_AT2,
         dataset_id="d68a8b48-abf0-4bd5-8834-155d34ec448a", obs_filter=PRIMARY_LUNG,
         title="Bronchopulmonary dysplasia + control"),
    # ---- mouse ----
    dict(key="lungmap", species="mouse", source="census", at2=CENSUS_AT2,
         dataset_id="047d57f2-4d14-45de-aa98-336c6f583750",
         obs_filter="development_stage in "
                    "['prime adult stage', '2-week-old stage', '7-day-old stage']",
         title="LungMAP postnatal (7d, 2wk, adult) — healthy"),
    dict(key="strunz_pbs", species="mouse", source="strunz", at2="AT2 cells",
         keep=("grouping", "PBS"), title="Strunz 2020 (GSE141259) — PBS control"),
    dict(key="strunz_bleo", species="mouse", source="strunz", at2="AT2 cells",
         title="Strunz 2020 (GSE141259) — bleomycin time course"),
    dict(key="angelidis", species="mouse", source="angelidis", at2="Type_2_pneumocytes",
         title="Angelidis 2019 (GSE124872) — aging"),
    dict(key="hurskainen", species="mouse", source="hurskainen", at2="AT2",
         title="Hurskainen 2021 (GSE151974) — normoxia + hyperoxia"),
    dict(key="hurskainen_normoxia", species="mouse", source="hurskainen", at2="AT2",
         keep=("oxygen", "Normoxia"), summary=False,
         title="Hurskainen 2021 (GSE151974) — normoxia only"),
]

# Dot plots (Figure 1 datasets): SLC34A2, AT2 identity markers, then current sort handles.
DOTPLOT = {
    "hlca": ["SLC34A2", "SFTPC", "ABCA3", "LAMP3", "EPCAM", "MUC1"],
    "hurskainen_normoxia": ["Slc34a2", "Sftpc", "Abca3", "Lamp3", "Epcam", "H2-Aa", "H2-Ab1"],
}
N_AT2_MARKERS = 4            # first 4 dot-plot genes are AT2 markers, the rest sort handles
DOT_MIN_N = 100              # cell types with fewer cells are left out of the dot plots
# Not surface proteins, so not screened; pulled only for the dot plots.
EXTRA_GENES = {"human": ["SFTPC"], "mouse": ["Sftpc"]}

FOCUS = {"human": "SLC34A2", "mouse": "Slc34a2"}
# Established AT2 isolation handles, overlaid on every scatter for reference.
# MUC1 is the likely HTII-280 antigen. Mouse MHC-II = I-A chains only (C57BL/6 are
# H2-b haplotype with a null H2-Ea, so I-E is not surface-expressed).
CANONICAL = {"human": {"EPCAM": "EpCAM", "MUC1": "MUC1 (HTII-280)"},
             "mouse": {"Epcam": "EpCAM", "H2-Aa": "MHC-II I-A (H2-Aa)",
                       "H2-Ab1": "MHC-II I-A (H2-Ab1)"}}


# ----------------------------------------------------------------------------
# Gene lists
# ----------------------------------------------------------------------------

def human_surfaceome(data_dir):
    """SURFY 'in silico surfaceome only' sheet: 2,886 proteins -> 2,799 unique gene symbols."""
    xlsx = os.path.join(data_dir, "surfy", "table_S3_surfaceome.xlsx")
    surf = pd.read_excel(xlsx, sheet_name="in silico surfaceome only", header=1)
    return sorted(surf["UniProt gene"].dropna().astype(str).unique())


def mouse_surfaceome(human_genes):
    """Map human surface genes to mouse: MGI homology first, else title-case symbol."""
    ortho = pd.read_csv(os.path.join(HERE, "reference", "human_mouse_orthologs_mgi.csv"))
    mgi = ortho.groupby("human")["mouse"].apply(list).to_dict()
    rows = []
    for hg in human_genes:
        for mg in mgi.get(hg, [hg.capitalize()]):
            rows.append((mg, hg))
    m = pd.DataFrame(rows, columns=["gene", "human_gene"]).drop_duplicates("gene")
    return sorted(m["gene"]), dict(zip(m["gene"], m["human_gene"]))


# ----------------------------------------------------------------------------
# Loaders. All return AnnData (cells x surface genes) with obs['cell_type'] and X =
# per-cell library-size-normalised expression (counts / total counts, no log), which is
# what the Census 'normalized' layer contains. Each is cached to <data>/cache/.
# ----------------------------------------------------------------------------

def libnorm(M):
    lib = np.asarray(M.sum(1)).ravel()
    lib[lib == 0] = 1.0
    return csr_matrix(M.multiply(1.0 / lib[:, None]))


def load_census(spec, genes, data_dir):
    organism = "homo_sapiens" if spec["species"] == "human" else "mus_musculus"
    with cc.open_soma(census_version=CENSUS_VERSION) as census:
        a = cc.get_anndata(
            census, organism=organism, measurement_name="RNA", X_name="normalized",
            obs_value_filter=f"dataset_id=='{spec['dataset_id']}' and {spec['obs_filter']}",
            var_value_filter=f"feature_name in {genes}",
            column_names={"obs": ["cell_type", "observation_joinid"], "var": ["feature_name"]})
    a.var_names = a.var["feature_name"].astype(str).values
    if spec.get("hlca_labels"):
        # HLCA's own annotations live in the source h5ad, matched on observation_joinid.
        src = ad.read_h5ad(os.path.join(data_dir, "hlca", "hlca_core.h5ad"), backed="r").obs
        src = src.set_index("observation_joinid")
        j = a.obs["observation_joinid"].astype(str).values
        lab = src.loc[j, "ann_finest_level"].astype(str).values
        a.obs["census_cell_type"] = a.obs["cell_type"].astype(str).values
        a.obs["cell_type"] = np.where(src.loc[j, "ann_level_3"].astype(str).values == "AT2",
                                      "AT2", lab)
    return a


def load_strunz(spec, genes, data_dir):
    g = os.path.join(data_dir, "geo", "GSE141259", "GSE141259_WholeLung_")
    gene_ids = pd.read_csv(g + "genes.txt.gz", header=None)[0].astype(str).values
    barcodes = pd.read_csv(g + "barcodes.txt.gz", header=None)[0].astype(str).values
    info = pd.read_csv(g + "cellinfo.csv.gz", index_col=0).reindex(barcodes)
    M = scipy.io.mmread(g + "rawcounts.mtx.gz").tocsr()
    if M.shape == (len(gene_ids), len(barcodes)):
        M = M.T.tocsr()
    a = ad.AnnData(X=libnorm(M), obs=pd.DataFrame(index=barcodes),
                   var=pd.DataFrame(index=gene_ids))
    a.obs["cell_type"] = info["cell.type"].astype(str).values
    a.obs["grouping"] = info["grouping"].astype(str).values
    a = a[a.obs["cell_type"] != "nan"]
    return a[:, [x for x in genes if x in set(gene_ids)]].copy()


def load_angelidis(spec, genes, data_dir):
    """The GEO .mtx has no row/column names; they are read from the companion .RData."""
    d = os.path.join(data_dir, "geo", "GSE124872")
    genes_f, bc_f = os.path.join(d, "genes.txt"), os.path.join(d, "barcodes.txt")
    if not os.path.exists(genes_f):
        rdata = os.path.join(d, "GSE124872_raw_counts_single_cell.RData")
        with gzip.open(rdata + ".gz") as fi, open(rdata, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        subprocess.run(["Rscript", "-e", f'library(Matrix); load("{rdata}"); o <- get(ls()[1]); '
                        f'writeLines(rownames(o), "{genes_f}"); '
                        f'writeLines(colnames(o), "{bc_f}")'], check=True)
    gene_ids = pd.read_csv(genes_f, header=None)[0].astype(str).values
    meta = pd.read_csv(os.path.join(d, "GSE124872_Angelidis_2018_metadata.csv.gz"), index_col=0)
    cols = pd.read_csv(bc_f, header=None)[0].astype(str).values
    assert [c.rsplit(":", 1)[-1] for c in (cols[0], cols[-1])] == \
        [c.rsplit(":", 1)[-1] for c in (meta.index[0], meta.index[-1])], "cell order mismatch"
    M = scipy.io.mmread(os.path.join(d, "GSE124872_raw_counts_single_cell.mtx.gz")).tocsr()
    if M.shape == (len(gene_ids), len(meta)):
        M = M.T.tocsr()                       # matrix columns follow metadata row order
    a = ad.AnnData(X=libnorm(M), obs=pd.DataFrame(index=meta.index),
                   var=pd.DataFrame(index=gene_ids))
    a.obs["cell_type"] = meta["celltype"].astype(str).values
    a = a[a.obs["cell_type"] != "low_quality_cells"]
    return a[:, [x for x in genes if x in set(gene_ids)]].copy()


def load_hurskainen(spec, genes, data_dir):
    """Dense genes x cells CSV, streamed keeping surface-gene rows only. Normalised by the
    authors' nCount_RNA. The authors' 'AT2 1' and 'AT2 2' clusters are merged as AT2."""
    d = os.path.join(data_dir, "geo", "GSE151974")
    meta = pd.read_csv(os.path.join(d, "GSE151974_cell_metadata_postfilter.csv.gz"), index_col=0)
    want, kept, rows = set(genes), [], []
    with gzip.open(os.path.join(d, "GSE151974_raw_umi_matrix_postfilter.csv.gz"), "rt") as fh:
        barcodes = np.array(fh.readline().rstrip("\n").split(",")[1:])
        for line in fh:
            i = line.find(",")
            if line[:i] in want:
                kept.append(line[:i])
                rows.append(np.array(line[i + 1:].split(","), dtype=np.float32))
    meta = meta.reindex(barcodes)
    lib = meta["nCount_RNA"].to_numpy(float)
    lib[lib == 0] = 1.0
    a = ad.AnnData(X=csr_matrix((np.vstack(rows) / lib[None, :]).T),
                   obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=kept))
    a.obs["cell_type"] = (meta["CellType"].astype(str)
                          .replace({"AT2 1": "AT2", "AT2 2": "AT2"}).to_numpy())
    a.obs["oxygen"] = meta["Oxygen"].astype(str).to_numpy()
    return a[a.obs["cell_type"] != "nan"].copy()


LOADERS = dict(census=load_census, strunz=load_strunz,
               angelidis=load_angelidis, hurskainen=load_hurskainen)


def load(spec, genes, data_dir):
    cache_key = spec["source"] if spec["source"] in ("strunz", "hurskainen") else spec["key"]
    path = os.path.join(data_dir, "cache", f"{cache_key}_surfaceome.h5ad")
    if os.path.exists(path):
        a = ad.read_h5ad(path)
    else:
        print(f"  building {path}")
        a = LOADERS[spec["source"]](spec, genes, data_dir)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        a.write_h5ad(path)
    a = a[:, a.var_names.isin(genes)]
    if spec.get("keep"):
        col, val = spec["keep"]
        a = a[a.obs[col].astype(str) == val].copy()
    if a.var_names.has_duplicates:     # duplicate symbols: keep the highest-total column
        tot = np.asarray(a.X.sum(0)).ravel()
        first = ~pd.Index(a.var_names[np.argsort(-tot)]).duplicated()
        a = a[:, np.sort(np.argsort(-tot)[first])].copy()
    return a


# ----------------------------------------------------------------------------
# Metrics
# ----------------------------------------------------------------------------

def tau(v):
    """Yanai tau: 0 = uniform across cell types, 1 = expressed in exactly one."""
    mx = v.max()
    return float((1 - v / mx).sum() / (len(v) - 1)) if mx > 0 else np.nan


def screen(a, at2_label):
    tot = np.asarray(a.X.sum(0)).ravel()
    a = a[:, tot > 0].copy()                  # genes never detected in this dataset
    genes = np.array(a.var_names)
    ct = a.obs["cell_type"].astype(str).values
    cats = pd.Categorical(ct)
    k = list(cats.categories).index(at2_label)
    n = len(ct)

    X = csr_matrix(a.X).tocsc()
    onehot = csr_matrix((np.ones(n), (np.arange(n), cats.codes)), shape=(n, len(cats.categories)))
    counts = np.bincount(cats.codes).astype(float)
    mean_ct = np.asarray((onehot.T @ X).todense()) / counts[:, None]
    det_ct = np.asarray((onehot.T @ (X > 0).astype(float)).todense())

    is_at2 = ct == at2_label
    eps = 1e-9
    others = np.delete(mean_ct, k, axis=0)
    top = np.argsort(-mean_ct, axis=0)
    df = pd.DataFrame({
        "gene": genes,
        "auroc": [roc_auc_score(is_at2, X[:, j].toarray().ravel()) for j in range(len(genes))],
        "tau": [tau(mean_ct[:, j]) for j in range(len(genes))],
        "log2fc_ct": np.log2((mean_ct[k] + eps) / (np.median(others, axis=0) + eps)),
        "pct_at2": det_ct[k] / counts[k] * 100,
        "pct_rest": (det_ct.sum(0) - det_ct[k]) / (n - counts[k]) * 100,
        "at2_is_top": top[0] == k,
        "top_offtarget": np.array(cats.categories)[np.where(top[0] == k, top[1], top[0])],
    })
    df = df.sort_values("auroc", ascending=False).reset_index(drop=True)
    df["rank_auroc"] = df["auroc"].rank(method="first", ascending=False).astype(int)
    pool = df["pct_at2"] >= DETECT_FLOOR
    df["rank_log2fc"] = np.nan
    df.loc[pool, "rank_log2fc"] = df.loc[pool, "log2fc_ct"].rank(method="first", ascending=False)
    meta = dict(n_cells=n, n_cell_types=len(cats.categories), n_at2=int(counts[k]),
                n_surface=len(df), n_log2_pool=int(pool.sum()))
    return df, meta


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------

def _fwd(a):   # AUROC -> -log10(1 - AUROC); ticks are labelled in plain AUROC
    return -np.log10(np.clip(1.0 - np.asarray(a, float), 1e-4, None))


def _inv(t):
    return 1.0 - 10.0 ** (-np.asarray(t, float))


def scatter(df, spec, metric, out_stem):
    focus, canon = FOCUS[spec["species"]], CANONICAL[spec["species"]]
    col = "auroc" if metric == "auroc" else "log2fc_ct"
    pool = df if metric == "auroc" else df[df["pct_at2"] >= DETECT_FLOOR]
    below = df[df["gene"].isin(canon) & ~df["gene"].isin(pool["gene"])]  # canonical, under floor
    top = pool.nlargest(15, col)
    foc = df[df["gene"] == focus].iloc[0]
    rank = int(foc["rank_auroc" if metric == "auroc" else "rank_log2fc"])

    fig, ax = plt.subplots(figsize=(11, 8))
    if metric == "auroc":
        ax.set_yscale("function", functions=(_fwd, _inv))
        lo, hi = 0.30, min(0.9997, 1 - (1 - pool[col].max()) / 3)
        ticks = [t for t in (0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99, 0.995, 0.998) if lo <= t <= hi]
        ax.set_yticks(ticks)
        ax.set_yticklabels([str(t) for t in ticks])
        ax.minorticks_off()
        ax.axhline(0.5, ls="--", color="red", alpha=0.5)
        ylabel = "AUROC, AT2 vs rest  (log(1 − AUROC) scale)"
    else:
        lo, hi = min(0.0, pool[col].min()), pool[col].max() * 1.08
        ax.axhline(0.0, ls="--", color="red", alpha=0.4)
        ylabel = f"log2 FC, AT2 mean / median cell type  (genes in ≥{DETECT_FLOOR}% of AT2)"

    ax.scatter(pool["tau"], pool[col].clip(lo, hi), s=14, c="#bbbbbb", edgecolors="none",
               label=f"surface genes ({len(pool)})")
    ax.scatter(top["tau"], top[col], s=70, c="#1f77b4", edgecolors="black", linewidths=0.5,
               label="top 15")
    for sub, color, lab in [(pool[pool["gene"].isin(canon)], "#eb6834", "canonical AT2 marker"),
                            (below, "#f6bfa4", f"canonical, < {DETECT_FLOOR}% of AT2")]:
        if len(sub):
            ax.scatter(sub["tau"], sub[col], s=70, c=color, edgecolors="white", label=lab)
            for _, r in sub.iterrows():
                ax.annotate(canon[r["gene"]], (r["tau"], r[col]), xytext=(7, -11),
                            textcoords="offset points", fontsize=8.5, fontweight="bold",
                            color="#eb6834")
    ax.scatter([foc["tau"]], [foc[col]], s=130, c="#e34948", edgecolors="white", label=focus)
    ax.annotate(focus, (foc["tau"], foc[col]), xytext=(-60, 10), textcoords="offset points",
                fontsize=11, fontweight="bold", color="red",
                arrowprops=dict(arrowstyle="-", color="red", lw=0.9))

    # top-15 labels stacked in a column at the left, with leader lines
    others = top[~top["gene"].isin([focus, *canon])]
    if metric == "auroc":
        ys = _inv(np.linspace(_fwd(min(0.99, hi)), _fwd(0.65), len(others)))
    else:
        ys = np.linspace(hi * 0.95, lo + (hi - lo) * 0.05, len(others))
    for y, (_, r) in zip(ys, others.iterrows()):
        ax.annotate(r["gene"], (r["tau"], r[col]), xytext=(0.63, y), textcoords="data",
                    fontsize=8.5, ha="right", va="center",
                    arrowprops=dict(arrowstyle="-", color="gray", lw=0.5, alpha=0.7))

    ax.set_xlim(0.5, 1.005)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("tau  (specificity across cell types)")
    ax.set_ylabel(ylabel)
    ax.set_title(f"{spec['title']}\n{focus} ranks #{rank} of {len(pool)} surface genes "
                 f"by {'AUROC' if metric == 'auroc' else 'log2 FC'}")
    ax.legend(loc="lower left")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(f"{out_stem}.{ext}", dpi=150, bbox_inches="tight")
    plt.close(fig)


def across_datasets(summary, species, metric, out_stem):
    col, rk, denom = (("auroc", "rank_auroc", "n_surface") if metric == "auroc"
                      else ("log2fc_ct", "rank_log2fc", "n_log2_pool"))
    s = summary.iloc[::-1].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(10, 0.8 * len(s) + 2))
    ax.barh(range(len(s)), s[col], color="#d62728", edgecolor="black", linewidth=0.5)
    ax.set_yticks(range(len(s)))
    ax.set_yticklabels(s["title"], fontsize=9)
    lo = 0.5 if metric == "auroc" else 0.0
    ax.set_xlim(lo, 1.0 if metric == "auroc" else s[col].max() * 1.35)
    for i, r in s.iterrows():
        ax.text(r[col], i, f"  #{int(r[rk])}/{int(r[denom])}, n AT2 = {r['n_at2']}",
                va="center", fontsize=8)
    ax.set_xlabel(f"{FOCUS[species]} {'AUROC' if metric == 'auroc' else 'log2 FC'} (AT2 vs rest)")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(f"{out_stem}.{ext}", dpi=150, bbox_inches="tight")
    plt.close(fig)


def dotplot(a, spec, out_stem):
    """Rows = cell types (n >= DOT_MIN_N), columns = DOTPLOT genes.
    colour = z-score of the cell type's mean expression, per gene across cell types;
    size = % of cells with the gene detected. z is per gene within this dataset only."""
    genes, at2 = DOTPLOT[spec["key"]], spec["at2"]
    names = list(a.var_names)
    X = a.X[:, [names.index(g) for g in genes]]
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    df = pd.DataFrame(X, columns=genes)
    df["cell_type"] = a.obs["cell_type"].astype(str).values
    g = df.groupby("cell_type")
    n = g.size()
    keep = n[n >= DOT_MIN_N].index
    mean = g.mean().loc[keep]
    pct = (g.agg(lambda v: (v > 0).mean()) * 100).loc[keep]
    z = (mean - mean.mean()) / mean.std(ddof=0)
    order = list(mean[genes[0]].sort_values(ascending=False).index)
    long = pd.concat({"mean_expr": mean.stack(), "z": z.stack(), "pct_detected": pct.stack()},
                     axis=1).reset_index()
    long.columns = ["cell_type", "gene", "mean_expr", "z", "pct_detected"]
    long["n_cells"] = long["cell_type"].map(n)

    cmap = LinearSegmentedColormap.from_list(
        "div", ["#184f95", "#3987e5", "#f0efec", "#e34948", "#9e2a2a"])  # blue-grey-red
    zmax = float(np.ceil(long["z"].abs().max()))
    norm = TwoSlopeNorm(vmin=-zmax, vcenter=0, vmax=zmax)
    s_max = 260
    row = {c: i for i, c in enumerate(order)}
    col = {gg: j for j, gg in enumerate(genes)}
    fig, ax = plt.subplots(figsize=(0.75 * len(genes) + 4.6, 0.24 * len(order) + 2.0))
    ax.scatter(long["gene"].map(col), long["cell_type"].map(row), s=long["pct_detected"] / 100 * s_max,
               c=long["z"], cmap=cmap, norm=norm, edgecolors="#7a7a7a", linewidths=0.4)
    tick = {"MUC1": "MUC1\n(HTII-280)", "H2-Aa": "H2-Aa\n(MHC-II)", "H2-Ab1": "H2-Ab1\n(MHC-II)"}
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels([tick.get(gg, gg) for gg in genes], fontsize=8.5)
    ax.get_xticklabels()[0].set_fontweight("bold")
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f"{c}  (n={n[c]:,})" for c in order], fontsize=7.5)
    for lab, c in zip(ax.get_yticklabels(), order):
        if c == at2:
            lab.set_fontweight("bold")
    ax.set_xlim(-0.6, len(genes) - 0.4)
    ax.set_ylim(len(order) - 0.5, -0.5)
    k = N_AT2_MARKERS - 0.5
    ax.axvline(k, color="#bbbbbb", lw=0.8)
    for x, txt in (((k - 0.5) / 2, "AT2 markers"), ((k + len(genes) - 0.5) / 2, "sort handles")):
        ax.annotate(txt, xy=(x, 1.0), xycoords=("data", "axes fraction"), xytext=(0, 34),
                    textcoords="offset points", ha="center", va="bottom", fontsize=8,
                    color="#555555", annotation_clip=False)
    ax.tick_params(length=0)
    ax.grid(color="#eeeeee", lw=0.6)
    ax.set_axisbelow(True)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(spec["title"], fontsize=11, pad=56)
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                      fraction=0.04, pad=0.03, shrink=0.35, anchor=(0, 1))
    cb.set_label("mean expression\n(z-score across cell types)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    cb.outline.set_visible(False)
    handles = [plt.scatter([], [], s=p / 100 * s_max, color="#999999", edgecolors="#7a7a7a",
                           linewidths=0.4) for p in (25, 50, 75, 100)]
    ax.legend(handles, ["25%", "50%", "75%", "100%"], title="% of cells\ndetected",
              loc="lower left", bbox_to_anchor=(1.02, 0.0), frameon=False, fontsize=7.5,
              title_fontsize=8, labelspacing=1.1, borderaxespad=0)
    for ext in ("png", "svg"):
        fig.savefig(f"{out_stem}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    long.insert(0, "dataset", spec["key"])
    return long


# ----------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--data", default=os.path.join(HERE, "data"))
    p.add_argument("--out", default=os.path.join(HERE, "results"))
    p.add_argument("--datasets", nargs="*", default=[d["key"] for d in DATASETS])
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)

    human = human_surfaceome(args.data)
    mouse, mouse_to_human = mouse_surfaceome(human)
    print(f"surfaceome: {len(human)} human genes -> {len(mouse)} mouse genes")

    summary, dots = [], []
    for spec in [d for d in DATASETS if d["key"] in args.datasets]:
        print(f"[{spec['key']}] {spec['title']}")
        genes = human if spec["species"] == "human" else mouse
        a = load(spec, genes + EXTRA_GENES[spec["species"]], args.data)
        df, meta = screen(a[:, a.var_names.isin(genes)], spec["at2"])
        if spec["key"] in DOTPLOT:
            dots.append(dotplot(a, spec, os.path.join(args.out, f"fig_dotplot_{spec['key']}")))
        if spec["species"] == "mouse":
            df.insert(1, "human_gene", df["gene"].map(mouse_to_human))
        df.round(4).to_csv(os.path.join(args.out, f"ranking_{spec['key']}.csv"), index=False)
        for m, tag in (("auroc", "auroc"), ("log2fc", "log2fc")):
            scatter(df, spec, m, os.path.join(args.out, f"fig_tau_vs_{tag}_{spec['key']}"))
        foc = df[df["gene"] == FOCUS[spec["species"]]].iloc[0]
        if spec.get("summary", True):
            summary.append(dict(key=spec["key"], species=spec["species"], title=spec["title"],
                                **foc[["auroc", "tau", "log2fc_ct", "pct_at2", "pct_rest",
                                       "rank_auroc", "rank_log2fc"]].to_dict(), **meta))
        print(f"  {meta['n_cells']} cells, {meta['n_at2']} AT2; {FOCUS[spec['species']]}: "
              f"AUROC {foc['auroc']:.3f} (#{foc['rank_auroc']}/{meta['n_surface']}), "
              f"log2FC {foc['log2fc_ct']:.2f} (#{int(foc['rank_log2fc'])}/{meta['n_log2_pool']})")

    if dots:
        pd.concat(dots).round(4).to_csv(os.path.join(args.out, "dotplot_values.csv"), index=False)
    summary = pd.DataFrame(summary)
    for sp, s in summary.groupby("species"):
        s.round(4).to_csv(os.path.join(args.out, f"slc34a2_summary_{sp}.csv"), index=False)
        for m in ("auroc", "log2fc"):
            across_datasets(s, sp, m, os.path.join(args.out, f"fig_slc34a2_across_{sp}_{m}"))


if __name__ == "__main__":
    main()
