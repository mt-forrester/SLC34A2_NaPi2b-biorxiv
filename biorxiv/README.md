# Surfaceome screen for AT2-cell surface markers

Code for the transcriptomic surfaceome analysis that nominated NaPi2b (*SLC34A2*) as an
alveolar type 2 (AT2) cell surface marker (Figure 1 and Figure S1). For every gene in the
SURFY human surfaceome, in each of 4 human and 5 mouse public lung scRNA-seq datasets, the
script scores how well that gene's expression separates AT2 cells from all other annotated
lung cells. It also draws marker dot plots for the two Figure 1 datasets.

## Quick start

```bash
pip install -r requirements.txt      # Python 3.12; R/Rscript also needed (see below)
bash download_data.sh                # SURFY, HLCA annotations, GEO files -> ./data (~6.1 GB)
python surfaceome_screen.py          # pulls Census data on first run; writes ./results
python surfaceome_screen.py --datasets hlca hurskainen_normoxia   # just the Figure 1 datasets
```

The first run caches each dataset to `data/cache/` (cells × surface genes only), so later
runs are fast. Pulling the HLCA core parenchyma cells from Census (333k cells) is the slow
step, about 30 minutes on a laptop. The 5.9 GB HLCA source file is read only for HLCA's
own cell annotations.

**R:** the GSE124872 matrix has no gene or cell names in the `.mtx`. They are read from the
companion `.RData` with `Rscript` (base R plus the `Matrix` package, which ships with R). Only the Angelidis dataset needs this.

## Outputs (`results/`)

| File | Contents |
|---|---|
| `fig_tau_vs_auroc_<dataset>.png/.svg` | τ (x) vs AUROC (y, log(1 − AUROC) scale; ticks read as AUROC). Figure 1: `hurskainen_normoxia` (mouse), `hlca` (human) |
| `fig_tau_vs_log2fc_<dataset>.png/.svg` | τ (x) vs log2 FC (AT2 mean / median cell-type mean), genes detected in ≥20% of AT2 |
| `fig_slc34a2_across_{human,mouse}_{auroc,log2fc}.png/.svg` | SLC34A2 value and rank in each dataset (Figure S1) |
| `fig_dotplot_{hlca,hurskainen_normoxia}.png/.svg` | Dot plot of SLC34A2, AT2 markers (SFTPC, ABCA3, LAMP3) and sort handles across cell types. Colour = z-score of mean expression per gene across cell types; size = % of cells detected. Values in `dotplot_values.csv` |
| `ranking_<dataset>.csv` | Every detected surface gene: `auroc`, `tau`, `log2fc_ct`, `pct_at2`, `pct_rest`, `at2_is_top`, `top_offtarget`, `rank_auroc`, `rank_log2fc` |
| `slc34a2_summary_{human,mouse}.csv` | SLC34A2 / Slc34a2 metrics, ranks and dataset sizes |

In every scatter, grey = all detected surface genes, blue = top 15 by that metric, red =
SLC34A2, orange = established AT2 isolation handles (EpCAM and MUC1, the likely
HTII-280 antigen, in human; EpCAM and MHC-II I-A in mouse).

## Data sources

Everything is public. Census datasets are pulled by `surfaceome_screen.py` from the
pinned **CELLxGENE Census release `2025-11-08`**. The GEO files and SURFY table are
fetched by `download_data.sh`.

### Surface-protein list

| Resource | Source | Used as |
|---|---|---|
| SURFY *in silico* human surfaceome (Bausch-Fluck et al., *PNAS* 2018; PMID 30373828) | `table_S3_surfaceome.xlsx` from <https://wollscheidlab.org/SURFY/>, sheet "in silico surfaceome only" | 2,886 predicted surface proteins → **2,799 unique gene symbols** (4 entries have no gene symbol, and 83 share a symbol with another entry) |
| MGI human–mouse homology | `HOM_MouseHumanSequence.rpt` from <https://www.informatics.jax.org/downloads/reports/> | `reference/human_mouse_orthologs_mgi.csv` (human, mouse symbol pairs). **Frozen copy** from the time of analysis: MGI updates this file, and the October 2026 release adds and drops 19 of 24,584 pairs (2 involve surface genes) |

Human surface genes are mapped to mouse through MGI homology (2,554 of 2,799 genes). The
other 245 fall back to the title-case symbol (`SLC34A2` → `Slc34a2`). A gene is screened in
a dataset only if it is present there and has a non-zero total count.

### Human datasets (CELLxGENE Census)

| Key | Dataset | Census `dataset_id` | Cell filter | Publication |
|---|---|---|---|---|
| `hlca` | An integrated cell atlas of the human lung in health and disease (core). Cell labels from HLCA's own `ann_finest_level`, read from the source h5ad (`download_data.sh`) | `066943a2-fdac-4b29-b348-40cede398e4e` | lung (= the lung-parenchyma samples; nose and airway samples excluded) | Sikkema et al., *Nat Med* 2023, doi:10.1038/s41591-023-02327-2 |
| `tsapiens` | Tabula Sapiens — All Cells | `53d208b0-2cfd-4366-9866-c3c6114081bc` | lung, primary cells | Tabula Sapiens Consortium, *Science* 2022, doi:10.1126/science.abl4896 |
| `krasnow` | Krasnow Lab Human Lung Cell Atlas, 10X | `8c42cfd0-0b0a-46d5-910c-fc833d83c45e` | lung, primary cells | Travaglini et al., *Nature* 2020, doi:10.1038/s41586-020-2922-4 |
| `bpd` | Bronchopulmonary dysplasia: single-nucleus RNA-seq (10x 3′ v3), BPD and age-matched control infants aged 4 months to 3 years (LungMAP **LMEX0000004400**; Xin Sun, Gloria Pryhuber, Allen Wang) | `d68a8b48-abf0-4bd5-8834-155d34ec448a` | lung, primary cells | No publication listed. LungMAP <https://www.lungmap.net/dataset?dataset_id=LMEX0000004400>, which links to CELLxGENE collection <https://cellxgene.cziscience.com/collections/3a5dbf8a-9b3e-4309-b4c5-d8a024f83734> (dataset version `759340c3-4596-40f6-88b9-b2b984686e06`); LungMAP project <https://data-browser.lungmap.net/projects/1977dc47-8414-4263-a870-6b0f207d8ab3> |

### Mouse datasets

| Key | Dataset | Source | Files | Publication |
|---|---|---|---|---|
| `lungmap` | Developmental single-cell atlas of the murine lung (LungMAP), postnatal stages only (P7, 2 wk, adult) | Census `047d57f2-4d14-45de-aa98-336c6f583750` | pulled by script | Zepp et al., *Science* 2021, doi:10.1126/science.abc3172 |
| `strunz_pbs`, `strunz_bleo` | Bleomycin injury time course, whole lung. PBS controls only, or all timepoints | GEO [GSE141259](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141259) | `GSE141259_WholeLung_{rawcounts.mtx,genes.txt,barcodes.txt,cellinfo.csv}.gz` | Strunz et al., *Nat Commun* 2020, doi:10.1038/s41467-020-17358-3 |
| `angelidis` | Aging mouse lung (Drop-seq) | GEO [GSE124872](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124872) | `GSE124872_raw_counts_single_cell.{mtx,RData}.gz`, `GSE124872_Angelidis_2018_metadata.csv.gz` | Angelidis et al., *Nat Commun* 2019, doi:10.1038/s41467-019-08831-9 |
| `hurskainen`, `hurskainen_normoxia` | Neonatal lung (P3–P14), normoxia + hyperoxia, or normoxia only (Figure 1) | GEO [GSE151974](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151974) | `GSE151974_raw_umi_matrix_postfilter.csv.gz`, `GSE151974_cell_metadata_postfilter.csv.gz` | Hurskainen et al., *Nat Commun* 2021, doi:10.1038/s41467-021-21865-2 |

Cell-type labels are used exactly as given by each source. AT2 is: HLCA `ann_level_3` =
`AT2` (its `AT2` + `AT2 proliferating`); the Census Cell Ontology label `pulmonary alveolar
type 2 cell` for the other Census datasets; GSE141259 `cell.type` = `AT2 cells`; GSE124872
`celltype` = `Type_2_pneumocytes`; GSE151974 `CellType` `AT2 1` + `AT2 2`. No cell types
were re-annotated or grouped beyond these source-defined AT2 classes.

HLCA uses its own annotations because the Census Cell Ontology mapping obscures them. For
example, HLCA's `AT0` becomes "pulmonary alveolar epithelial cell", and the nasal,
bronchial and subsegmental goblet subsets become three differently-worded goblet labels.

## Method in brief

See `METHODS.md` for the full methods. Expression is library-size normalised (counts
divided by each cell's total counts, not log-transformed). This is the Census `normalized`
layer, and the GEO data are normalised the same way. For each gene:

- **AUROC**: `sklearn.metrics.roc_auc_score(is_AT2, expression)` over all cells (AT2 vs all
  other cells). 0.5 = no separation, 1 = perfect.
- **τ** (Yanai et al. 2005): over the per-cell-type mean expression *x*,
  τ = Σ(1 − *x*ᵢ / max *x*) / (*n* − 1). 0 = uniform, 1 = restricted to one cell type. τ is
  not AT2-specific: it is high for any gene concentrated in *some* cell type.
- **log2 FC**: log2(AT2 mean / median of the other cell types' means). Each cell type
  counts once, so abundant populations do not dominate. It is ranked only among genes
  detected in ≥20% of AT2 cells, because the ratio is unstable for sparsely detected
  genes.

The log2 FC pool size depends on sequencing depth (36–554 genes across datasets), so its
rank is comparable within a dataset only. AUROC ranks are over all detected surface genes.

## Reproducibility notes

`results/` was produced by this exact code on 2026-10-29 from freshly downloaded GEO,
SURFY and HLCA files and Census release `2025-11-08`. Differences from the earlier working
analysis, all corrected here:

- **HLCA and Krasnow cell sets.** The earlier HLCA and Krasnow caches had been pulled
  before the `tissue_general == 'lung'` filter was added, so they also held non-lung cells:
  HLCA nose and airway samples (584,944 cells vs 333,468 parenchyma cells), and 4,669
  Krasnow cells outside lung. With the filter, HLCA SLC34A2 is AUROC 0.971, #1 of 2,217
  (was 0.975, #1 of 2,221). Krasnow is 0.971, #1 of 2,208 (was 0.972, #1 of 2,216).
- **Rounding.** The earlier analysis ranked log2 FC on values rounded to 2 decimals and
  applied the 20% detection floor to values rounded to 1 decimal. Here unrounded values
  are used, which changes a few log2 FC ranks and pool sizes by one.
- **Duplicate symbols.** The earlier HLCA pull held SIGLEC5 twice (two Ensembl IDs).
  Duplicates are now collapsed.
