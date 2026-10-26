Final Figure S9 assets

layout.json: approved 84-codon row ordering, text coordinates and font sizes,
independent heatmap cell coordinates and colorbar rectangles.
fonts.json: embedded Arial regular/italic/bold and Symbol PDF font resources.
These are style resources, not analytical counts or old heatmap results.
The calculation starts from the complete upstream two-editable-site catalog
and allele counts; it does not filter by the historical 84-codon table.
A changed retained set raises an explicit error because the approved publication
layout contains these 84 rows. Update layout explicitly for a different dataset.

Run from the package root: bash code/06/run_figure_s9.sh
Override FS9_COUNTS_DIR, FS9_OUT_DIR, FS9_PYTHON as needed.
Runtime: Python >=3.10, standard library only. No /tmp dependencies, R, fitz,
reportlab, old FigureS9 directories, or source reference PDFs are needed.
Outputs: an editable PDF, source-data TSVs and cohort summary counts. Font editing to new glyphs may require
locally installed Arial. No embedded raster images or Illustrator private data.
