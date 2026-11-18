Supplementary Figure S18

Run from the repository root:
  bash code/07/half_life/run_half_life.sh

Dependencies
  R: data.table and jsonlite
  Python 3.6 or later: PyMuPDF
  Set RSCRIPT or PYTHON to select an interpreter. PDF_PYTHON can select
  a separate Python interpreter for plotting.

Reference files
  HALF_LIFE_GTF: dmel-all-r6.12.gtf
  HALF_LIFE_GENOME_FASTA: dmel-all-chromosome-r6.12.fasta
  HALF_LIFE_TRANSCRIPT_FASTA: dmel-all-transcript-r6.12.fasta
  Set these variables to absolute paths. The script can also locate the
  files in a ref_genome directory in the repository or an ancestor.

Inputs and method
  The public embryonic editing table and the normalized expression time
  series from Thomsen et al. (2010) are in inputs/published/. The editing
  table uses dm6 coordinates. Half-lives are estimated from the expression
  time series. Genes with conflicting decay classes among mapped probe sets
  are excluded.

  For an edited gene, choose the longest CDS from transcripts whose CDS
  contains a detected editing site. For an unedited gene, choose the
  longest CDS. Break ties by mature transcript length and transcript ID.
  Map the CG9243 measurement once to CG43346 (FBgn0263051).

  For each site in the selected CDS, edited-read counts are reconstructed by
  rounding the product of reported depth and editing fraction. Counts and
  depths are pooled across embryonic intervals with positive reported editing.
  The pooled count divided by pooled depth gives the site fraction. The sum
  of site fractions divided by the CDS adenosine count gives the gene-level mean.

Outputs
  HALF_LIFE_OUT_DIR sets the analysis table directory, which defaults to
  code/07/half_life/output/duan2017_embryo_cds_min1/.
  HALF_LIFE_PDF sets the figure path, which defaults to
  output/Supplementary_Figure_18_r6_12_transcripts_editable.pdf.
  Panel E includes 145 gene-level measurements and 165 editing sites
  (Spearman rho = -0.117379; P = 0.159705).
