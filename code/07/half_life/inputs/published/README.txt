The bundled Duan2017_Graveley2011_embryo_CDS_editing_sites_long.tsv is the
exact embryo CDS editing-site input selected for Figure S18, not a sequencing
count matrix. It was copied unchanged from the study half-life analysis.
Source measurements: Graveley et al. 2011 Supplementary Table 29; preparation
code: ../../scripts/extract_duan2017_embryo_sites.py (see --help for inputs).
The half-life launcher reads this filename directly. The three public expression/half-life tables used by config.R are also bundled:
Thomsenetal.normalised.classification_orig.txt, 2023_S1_Table.xlsx and
2023_S2_Table.xlsx. They are copied unchanged from the actual analysis inputs;
source URLs and SHA-256 values are in published_tables_manifest.tsv.
00_prepare_inputs.R uses these local files only; a missing file causes an error,
not a download. Reference genome GTF/FASTA files remain external inputs.
SHA-256: 73893d72b23dd76c66ec6362f611369cab0bdd0e45e77880ea3b61a970a6ec08
