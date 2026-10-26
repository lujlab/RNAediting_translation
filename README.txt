RNA-editing study: code and experimental inputs
===============================================
The numbered Bash files describe the analyses; code/ contains their scripts,
data/ contains the bundled experimental inputs, and config/ contains input
configuration examples. Each experimental dataset is supplied once.
Generated figures, result tables and historical figure PDFs are not bundled.
The plotting commands generate outputs locally; no full-page composition step
is required. Table1 filter scatter panels and S10 species panels are saved separately.

Requirements
Use Python >=3.10 and install requirements.txt. R package requirements are listed below; code/08/setup/install_packages.R
installs the workflow 08 dependencies. Set RSCRIPT and PYTHON_BIN when needed.
Cairo graphics and the fonts used by a figure must be available for rendering.
Ribo-seq/RNA-seq/MS analyses additionally require the tools named in their Bash
entrypoints (including aligners, UMI tools, MaxQuant/DIA-NN and references).

R dependencies used across workflows (install those needed for your analysis):
  CRAN: data.table, dplyr, tidyr, tibble, readr, readxl, stringr, plyr,
  ggplot2, ggprism, ggpubr, cowplot, patchwork, gridExtra, gplots, ggseqlogo,
  pheatmap, RColorBrewer, scales, UpSetR, lme4, lmerTest, glmmTMB, VGAM,
  pspearman, psych and statmod. MASS is an R recommended package.
  Bioconductor: AnnotationDbi, Biobase, ChIPseeker, DEGreport, DESeq2, edgeR,
  clusterProfiler, GO.db, GOSemSim, GenomicFeatures, GenomicRanges, rtracklayer,
  org.Dm.eg.db, org.Hs.eg.db, TxDb.Hsapiens.UCSC.hg38.knownGene and txdbmaker.
  Also install riboWaltz for workflow 01; knitr/rmarkdown for the Figure 5 notebook.
The dependency lists are not a version-locked, clean-install-tested environment.

All-analysis entrypoints
01  Ribo-seq processing and Figure S3.
02  RNA-seq processing, codon counts and variants.
03  Editing-site confirmation and Figures 1/S1.
04  RNA/Ribo comparisons and Tables 1/S4/S5.
05  Stalling, gmRC, Figures 2/3 and associated supplementary panels.
06  MS processing, protein evidence and associated figures/tables.
07  Figures 4C-F, S15D-G, S16, S17, S18 and S27.
08  Genotype-dependent editing/expression analyses and Figure 5.
09  Human/mouse analyses, Figure S28 and Tables S12/S13.
10  Figure 6.

Configure section 00 and run the required numbered sections in order. These
are analysis recipes with alternative branches, not a single command that runs
all branches together. Intermediate counts, references and model inputs are
required between processing stages; they are not extra raw-data copies.
The 07 entrypoint calls S16 in section 02 and retains only final figures in
results/07/. S17/S18/S27 require the external inputs described in their scripts.
With existing human/mouse counts, run 09 section 00 followed by sections 05-07.
Table S10 has a standalone entrypoint: bash code/06/run_table_s10.sh.

Table1 filter scatter analysis
code/04/table1_filter/ contains the scatter analysis using the Table 1
minrep=2 eligibility rules. Section 07b of 04_RNA_Ribo_alledit_Table1_analysis.sh
calls its run.sh. TABLE1_FILTER_OUT sets the output directory; individual
panels use the Table1_filter_ filename prefix, with metric, species and sex.

Workflow 05 code directories
Figure numbers follow the manuscript and supplementary figure legends.
figure2_3_s6_edger/           Figure 2B, 3A and S6A: edgeR peak enrichment.
  peak_only/                 Peak models and panels called by the main Bash.
figure_s6_s8_pausepred/       Figure S6B/C and S8: PausePred and enrichment.
  enrichment/                Editing-site mapping and hairpin enrichment.
figure2_s7_metagene/          Figure 2C and S7A: metagene profiles.
figure3_s8_s11_rna_structure/ RNA hairpin inputs for Figure 3, S8 and S11.
figure_s10_codon_preference/ Figure S10: codon-preference models and plots.
figure2_3_s7_s10_s11_gmrc/    gmRC workflow and final models for Figure 2/3 and S7/S10/S11.
config/                     Shared sample tables and reference configuration.
These directories are under code/05/. Use 05_Stalling_analysis.sh for the
analysis order; shared inputs and scripts are kept in a single location.

Workflow 06 code directories
figure4_s14_diann/   Figure 4B/S14 DIA-NN proteoform fractions.
figure_s9_assets/   Figure S9 layout and font resources used by figure_s9_equalrep.py.
table_s10_multiedit/ Table S10 multi-editing-site RNA/RPF/MS evidence.
config/            Shared mass-spectrometry configuration examples.
These directories are under code/06/. Sections 08, 08B and 08C of
06_MS_analysis.sh call Table S10, Figure S9 and Figure 4B/S14, respectively.

Workflow 08 code directories
setup/               Input paths, sample metadata, package installation and input checks.
analysis/editing/    Editing-level models, clustering, editing groups and GO analysis.
analysis/expression/ RNA counts, DESeq2 models, expression clusters and enrichment.
analysis/common.R    Shared paths and analysis settings.
figures/             Standalone plots of auto-editing, protein output and expression changes.
figures/figure_s26/   Figure S26 fitting, grouping and plotting workflow.
All paths above are relative to code/08/. The main 08 Bash file calls the
numbered analysis scripts and figures/figure_s26/run.sh. The Figure 5 notebook
and the descriptive-named analysis/plotting scripts are run separately after
configuring their input paths; the pairwise UpSet script also requires the
RNA objects prepared by the Figure 5 notebook.
