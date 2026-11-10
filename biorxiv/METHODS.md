# Methods: transcriptomic surfaceome analysis

### Surfaceome analysis of public single-cell RNA-seq data

**Candidate surface proteins.** Candidates were restricted to proteins with potential for
cell-surface expression, defined by the SURFY *in silico* human surfaceome [1]. This list
has 2,886 predicted surface proteins (Table S3, sheet "in silico surfaceome only",
downloaded from https://wollscheidlab.org/SURFY/), which correspond to 2,799 unique gene
symbols. For mouse datasets, human genes were mapped to mouse orthologs with the Mouse
Genome Informatics (MGI) human–mouse homology table (HOM_MouseHumanSequence.rpt) [2]. This
covered 2,554 genes. The remaining 245 genes, which have no MGI entry, were matched by
symbol (e.g., *SLC34A2* → *Slc34a2*).

**Datasets.** We analysed four human and five mouse lung scRNA-seq datasets (table below).
Human data were the lung-parenchyma samples of the Human Lung Cell Atlas (HLCA) core [3]
(nasal and airway samples excluded), the lung portion of Tabula
Sapiens [4], the Travaglini et al. lung atlas [5], and a LungMAP single-nucleus RNA-seq
dataset of bronchopulmonary dysplasia (BPD) and age-matched control infants (4 months to
3 years; LungMAP LMEX0000004400, deposited by X. Sun, G. Pryhuber and A. Wang; accessed
through CELLxGENE collection 3a5dbf8a-9b3e-4309-b4c5-d8a024f83734). One mouse dataset
was the postnatal (P7, 2-week and adult) cells of the LungMAP developmental atlas [6]. The
others were a bleomycin injury time course [7] (GSE141259; analysed both as PBS controls
alone and as the full time course), an aging lung atlas [8] (GSE124872), and a neonatal
normoxia/hyperoxia dataset [9] (GSE151974; analysed both in full and as the normoxia arm
alone). The HLCA, Tabula Sapiens, Travaglini, BPD
and LungMAP data were obtained from CZ CELLxGENE Census [10] (release 2025-11-08),
restricted to lung tissue and, except for the integrated HLCA core, to primary cells. The
GEO datasets were built from the deposited raw count matrices and author cell annotations.
Datasets were included only if the non-AT2 compartment contained substantial alveolar
epithelium, including AT1 cells. In datasets dominated by immune cells, an AT2-versus-rest
comparison mainly measures epithelial-versus-immune separation.

**Cell annotations and normalisation.** Cell-type labels were used as provided by each
source, with no re-clustering or re-annotation. For HLCA, the atlas's own finest-level
annotations (ann_finest_level) were used, obtained from the CELLxGENE source file, and AT2
cells were those in its level-3 "AT2" class ("AT2" and "AT2 proliferating"). For the
other Census datasets, AT2 cells were those labelled "pulmonary alveolar type 2 cell"
(Cell Ontology). In the GEO datasets they were those labelled "AT2 cells" (GSE141259),
"Type_2_pneumocytes" (GSE124872), or "AT2 1" and "AT2 2" combined (GSE151974). All other annotated cells formed the comparison group. Expression was
library-size normalised by dividing each cell's counts by its total counts, without log
transformation. This matches the Census "normalized" layer, and the GEO data were
normalised the same way. Surface genes with zero total counts in a dataset were excluded.

**Specificity metrics.** Three metrics were computed for each surface gene in each
dataset. (i) *AUROC*: the area under the receiver operating characteristic curve, using
the gene's expression as the only classifier of AT2 versus all other cells
(scikit-learn `roc_auc_score` [11]). This measures how well the gene separates AT2 cells
from everything else. 0.5 = no separation and 1 = perfect separation. (ii) *Tau (τ)*:
the tissue-specificity index of Yanai et al. [12], computed over the mean expression of
each annotated cell type as τ = Σᵢ(1 − *x*ᵢ/*x*ₘₐₓ)/(*n* − 1). τ ranges from 0 (uniform
expression) to 1 (expression in a single cell type). (iii) *log2 fold change*: log2 of
the mean AT2 expression divided by the median of the other cell types' mean expression.
Each cell type is weighted equally regardless of its abundance. The log2 fold change was
ranked only among genes detected in ≥20% of AT2 cells, because the ratio is unstable for
sparsely detected genes. Genes were ranked by AUROC among all detected surface genes and
by log2 fold change among the ≥20%-detected subset. Results are shown as scatter plots of
τ against AUROC (plotted on a log(1 − AUROC) scale to resolve values near 1) or against
log2 fold change. Figure 1 shows the HLCA parenchyma and the Hurskainen normoxia arm, the
mouse dataset in which every lung compartment (AT1, airway, endothelial, mesenchymal and
immune cells) is well represented. SLC34A2's AUROC and rank in each dataset are summarised
in Figure S1.

**Marker dot plots.** For the two Figure 1 datasets, SLC34A2 was compared with the AT2
markers SFTPC, ABCA3 and LAMP3 and with current sorting handles (EPCAM and MUC1, the likely
HTII-280 antigen, in human; EPCAM and the MHC class II I-A chains H2-Aa and H2-Ab1 in
mouse). For each cell type with ≥100 cells, dot colour shows the z-score of the cell
type's mean expression, computed for each gene across cell types, and dot size shows the
percentage of cells in which the gene was detected. All analyses used Python 3.12 with cellxgene-census 1.18.0, anndata 0.12.18,
NumPy, pandas, SciPy, scikit-learn 1.9.0 and Matplotlib. Code is available at https://github.com/mt-forrester/SLC34A2_NaPi2b-biorxiv/tree/main/biorxiv.

---

### Datasets used in the surfaceome analysis

| Species | Dataset | Accession | Condition | Cells | Cell types | AT2 cells |
|---|---|---|---|---|---|---|
| Human | HLCA core, lung parenchyma [3] | Census 066943a2-fdac-4b29-b348-40cede398e4e | healthy | 333,468 | 58 | 60,424 |
| Human | Tabula Sapiens, lung [4] | Census 53d208b0-2cfd-4366-9866-c3c6114081bc | healthy | 65,847 | 34 | 11,594 |
| Human | Travaglini et al. [5] | Census 8c42cfd0-0b0a-46d5-910c-fc833d83c45e | healthy | 60,993 | 46 | 3,814 |
| Human | LungMAP BPD (snRNA-seq) | LungMAP LMEX0000004400; Census d68a8b48-abf0-4bd5-8834-155d34ec448a | BPD + age-matched control infants | 255,204 | 39 | 88,298 |
| Mouse | LungMAP, postnatal [6] | Census 047d57f2-4d14-45de-aa98-336c6f583750 | healthy | 32,929 | 17 | 1,760 |
| Mouse | Strunz et al., PBS [7] | GEO GSE141259 | healthy control | 6,798 | 39 | 2,342 |
| Mouse | Strunz et al., time course [7] | GEO GSE141259 | bleomycin injury | 29,268 | 39 | 2,711 |
| Mouse | Angelidis et al. [8] | GEO GSE124872 | aging | 14,318 | 33 | 4,199 |
| Mouse | Hurskainen et al. [9] | GEO GSE151974 | normoxia + hyperoxia | 61,839 | 33 | 6,155 |
| Mouse | Hurskainen et al., normoxia only [9] | GEO GSE151974 | healthy neonatal | 33,561 | 33 | 3,073 |

### References

1. Bausch-Fluck D, *et al.* The in silico human surfaceome. *Proc Natl Acad Sci USA* 2018;115:E10988–E10997. PMID 30373828.
2. Baldarelli RM, *et al.* Mouse Genome Informatics: an integrated knowledgebase system for the laboratory mouse. *Genetics* 2024;227:iyae031.
3. Sikkema L, *et al.* An integrated cell atlas of the lung in health and disease. *Nat Med* 2023;29:1563–1577.
4. The Tabula Sapiens Consortium. The Tabula Sapiens: a multiple-organ, single-cell transcriptomic atlas of humans. *Science* 2022;376:eabl4896.
5. Travaglini KJ, *et al.* A molecular cell atlas of the human lung from single-cell RNA sequencing. *Nature* 2020;587:619–625.
6. Zepp JA, *et al.* Genomic, epigenomic, and biophysical cues controlling the emergence of the lung alveolus. *Science* 2021;371:eabc3172.
7. Strunz M, *et al.* Alveolar regeneration through a Krt8+ transitional stem cell state that persists in human lung fibrosis. *Nat Commun* 2020;11:3559.
8. Angelidis I, *et al.* An atlas of the aging lung mapped by single cell transcriptomics and deep tissue proteomics. *Nat Commun* 2019;10:963.
9. Hurskainen M, *et al.* Single cell transcriptomic analysis of murine lung development on hyperoxia-induced damage. *Nat Commun* 2021;12:1565.
10. CZI Cell Science Program, *et al.* CZ CELLxGENE Discover: a single-cell data platform for scalable exploration, analysis and modeling of aggregated data. *Nucleic Acids Res* 2025;53:D886–D900.
11. Pedregosa F, *et al.* Scikit-learn: machine learning in Python. *J Mach Learn Res* 2011;12:2825–2830.
12. Yanai I, *et al.* Genome-wide midrange transcription profiles reveal expression level relationships in human tissue specification. *Bioinformatics* 2005;21:650–659.
