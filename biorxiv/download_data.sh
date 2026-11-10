#!/usr/bin/env bash
# Download every input not fetched by surfaceome_screen.py itself.
# (The 5 CELLxGENE Census datasets are pulled on first run by the Python script.)
set -euo pipefail
DATA="${1:-$(dirname "$0")/data}"
GEO=https://ftp.ncbi.nlm.nih.gov/geo/series

get() {  # get <dest dir> <url>
  mkdir -p "$1"
  local f="$1/$(basename "$2")"
  [ -s "$f" ] || curl -fL --retry 3 -o "$f" "$2"
}

# SURFY in silico human surfaceome (Bausch-Fluck et al. 2018, PNAS), table S3
get "$DATA/surfy" https://wollscheidlab.org/SURFY/table_S3_surfaceome.xlsx

# HLCA core source file (5.9 GB), read only for HLCA's own cell annotations
# (ann_finest_level / ann_level_3); expression comes from Census. CELLxGENE dataset
# 066943a2-fdac-4b29-b348-40cede398e4e, collection 6f6d381a-7701-4781-935c-db10d30de293.
mkdir -p "$DATA/hlca"
[ -s "$DATA/hlca/hlca_core.h5ad" ] || curl -fL --retry 3 -o "$DATA/hlca/hlca_core.h5ad" \
  https://datasets.cellxgene.cziscience.com/688185ad-11c2-4172-a53a-f4f1f4076860.h5ad

# Strunz et al. 2020, mouse bleomycin whole lung
for f in barcodes.txt.gz genes.txt.gz cellinfo.csv.gz rawcounts.mtx.gz; do
  get "$DATA/geo/GSE141259" "$GEO/GSE141nnn/GSE141259/suppl/GSE141259_WholeLung_$f"
done

# Angelidis et al. 2019, mouse aging lung (Drop-seq)
for f in raw_counts_single_cell.mtx.gz raw_counts_single_cell.RData.gz Angelidis_2018_metadata.csv.gz; do
  get "$DATA/geo/GSE124872" "$GEO/GSE124nnn/GSE124872/suppl/GSE124872_$f"
done

# Hurskainen et al. 2021, neonatal mouse lung, normoxia vs hyperoxia
for f in raw_umi_matrix_postfilter.csv.gz cell_metadata_postfilter.csv.gz; do
  get "$DATA/geo/GSE151974" "$GEO/GSE151nnn/GSE151974/suppl/GSE151974_$f"
done
echo "done -> $DATA"
