#!/usr/bin/env bash
# Example overrides; replace input paths with the locations of the study data.
export PYTHON_BIN=python3
export PLOT_PYTHON=python3
export RENDER_PYTHON=python3
export RSCRIPT_BIN=Rscript
export RSCRIPT=Rscript
export RNA_EDITING_DATA_ROOT=/path/to/RNA_editing_DATAprocessing_new
export GMRC_ANNOTATION_ROOT=/path/to/peak_find
export RESCUE_ROOT=/path/to/codon_editing_project
export MAMMALS_DATA_ROOT=/path/to/11.mammals
export MAMMALS_REFERENCE_ROOT=/path/to/mammalian_references
# Configure the per-workflow sample/configuration TSVs in code/ before execution.

export TABLE_S10_INPUTS=/path/to/upstream_MS_codon_evidence
export DIA_NN_ANALYSIS_DIR=/path/to/DIA_NN_analysis
export MS_PROJECT_ROOT=/path/to/MS_analysis_project
export EDITING_PROJECT_ROOT=/path/to/editing_project

# Workflow 03 literature-reference inputs (not downloaded automatically).
export PUBLISHED_EDITING_DIR=/path/to/editing_sites_PUBlist
export DUAN2017_FILE=/path/to/2017edit_all5925.txt
