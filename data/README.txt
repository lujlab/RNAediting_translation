Experimental input inventory
============================
Each dataset below is stored once. Reference results, combined presentation
workbooks and complete historical archives are not included as extra copies.
Final figures are generated under results/. Required intermediate tables are
created while an analysis runs; standalone validation reports are not distributed.

07/20260322data-ITP.xlsx
  Figure 4C RNA-digest measurements and the standards/controls used to derive
  concentrations; read by code/07/figure4_CEF.R. Existing calibration and
  outlier-handling rules remain in the analysis code.

07/dual_plot.txt
  Figure 4D: 11 concentrations x 9 included technical wells. Read by
  code/07/figure4_CEF.R and code/07/standalone_panel_data.R. Both readers use
  this same file. Its 90-minute translation assay is separate from S16.

07/WB.txt
  Figure 4E: 6 concentrations x 3 measurements. Read by code/07/figure4_CEF.R.

07/Westernblot_and_luciferase.rds
  Figure S15D-G: four experimental tables. Both reporter scripts read the same
  RDS. Replicate measurements are used to recalculate plotted means/SDs;
  stored mean/SE columns are not treated as additional observations.

07/20260413-DECAY_data_R_3.xls.xlsx
07/20260410-DECAY_data_R_5.xls.xlsx
  Figure S16: 72 Ct measurements per workbook, 144 total, only 5/30/60 min.
  code/07/run_S16.sh reads these through extract_s16_qpcr.py, then runs
  figure_S16_qPCR.R. The extracted Ct table is temporary; the command exports only the final PDF/PNG.

10/phenotype_data.tsv
  Figure 6: 1239 observations from the eight plotted assays; read directly by
  code/10/figure6.R through 10_Phenotype_analysis.sh. Original assay/source_row,
  sample_id, batch and raw_value are retained. The value field records the
  existing plotting-unit conversion for that same observation, not another
  sample. Historical Egg Production and Starvation datasets are not included.
  Sample IDs were recovered from the laboratory Rdata; legacy labels were
  decoded from GB18030 to UTF-8 when they were not already UTF-8.

Figure 4F's 18 raw luciferase measurements are defined once in the data_F block
of code/07/figure4_CEF.R; there is no second source-data file for those values.

Shared configuration and public reference inputs
The RNA correction table is stored once at config/rna_count_scaling.tsv.
The published half-life/editing inputs under code/07/half_life/inputs/published
are used by the S18 analysis; their sources/checksums are documented there.
Required figure-layout/font assets are code resources, not extra observations.
