#!/usr/bin/env Rscript

# Memory-efficient wrapper for All_metagene_mean_existingMean_any1libLm02.R.  It
# preserves the original statistic, but streams each mean-gmRC table and retains only
# transcripts containing a new-rule editing center.  Those are the only transcripts
# used by both edited centers and their matched controls.

source_script <- Sys.getenv("METAGENE_SOURCE_SCRIPT", unset = "/data/dongxn/RNAediting/codex/rescue_3R_18907675_1v1_pipeline_clean/09_Ipos/All_metagene_mean_existingMean_any1libLm02.R")
new_alledit_dir <- Sys.getenv("METAGENE_ALLEDIT_DIR", unset = "/data/dongxn/RNAediting/codex/rescue_3R_18907675_1v1_pipeline_clean/04_alledit_new/results")

stopifnot(file.exists(source_script), dir.exists(new_alledit_dir))
code <- readLines(source_script, warn = FALSE)

# The sourced analysis keeps all of its ancillary inputs relative to 09_Ipos, not
# relative to this wrapper in the MS workspace.
script_dir_line <- grep("^script_dir <- normalizePath", code)
stopifnot(length(script_dir_line) == 1L)
code <- append(code, sprintf("script_dir <- %s", deparse(normalizePath(dirname(source_script)))),
    after = script_dir_line)

# Route candidate editing calls to the counter-rule alledit results.
first_alledit <- grep("^alledit_dir <- file.path", code)
last_alledit <- grep("^alledit_dir <- normalizePath", code)
stopifnot(length(first_alledit) == 1L, length(last_alledit) == 1L, last_alledit > first_alledit)
code[first_alledit:(last_alledit - 1L)] <- sprintf("alledit_dir <- %s", deparse(normalizePath(new_alledit_dir)))
code[last_alledit] <- "alledit_dir <- normalizePath(alledit_dir, mustWork = TRUE)"

# Replace the full-table reader.  `current_edit_species` is set immediately before
# every call, so all retained rows belong to a transcript with an edited center.  The
# original code draws controls only from those same transcripts, hence this does not
# alter the eligible data universe.
reader_start <- grep("^read_mean_table <- function", code)
reader_end <- which(seq_along(code) > reader_start & code == "}")[1L]
stopifnot(length(reader_start) == 1L, !is.na(reader_end))
stream_reader <- c("read_mean_table <- function(file) {", "  allowed <- unique(sub(\"_.*\", \"\", any_edit_by_species[[current_edit_species]]))",
    "  if (!length(allowed)) stop(\"No editing-center transcripts for \", current_edit_species)",
    "  pattern_file <- tempfile(pattern = \"metagene_fbtr_\", fileext = \".txt\")", "  on.exit(unlink(pattern_file), add = TRUE)",
    "  writeLines(allowed, pattern_file)", "  awk_program <- \"FILENAME == ARGV[1] {keep[$0]=1; next} FNR == 1 {for(i=1;i<=NF;i++){h=$i; gsub(/^\\\"|\\\"$/, \\\"\\\", h); if(h==\\\"FBtr\\\")col=i} if(!col)exit 2; print; next} {k=$col; gsub(/^\\\"|\\\"$/, \\\"\\\", k); if(k in keep)print}\"",
    "  cmd <- sprintf(\"gzip -cd -- %s | awk -F '\\t' %s %s -\", shQuote(file), shQuote(awk_program), shQuote(pattern_file))",
    "  x <- data.table::fread(cmd = cmd, select = c(\"FBtr_codonpos\", \"FBtr\", \"codon\", \"AA\", \"gmRC\", \"chr\", \"pos\"), na.strings = c(\"\", \"NA\", \"na\", \"NaN\", \"nan\"))",
    "  x[, `:=`(FBtr_codonpos = as.character(FBtr_codonpos), FBtr = as.character(FBtr), codon = toupper(as.character(codon)), AA = as.character(AA), gmRC = suppressWarnings(as.numeric(gmRC)), chr = as.character(chr), pos = suppressWarnings(as.numeric(pos)))]",
    "  x <- x[!is.na(FBtr_codonpos) & FBtr_codonpos != \"\"]", "  x[, codon_pos := suppressWarnings(as.integer(sub(\".*_\", \"\", FBtr_codonpos)))]",
    "  x[, FBtr := data.table::fifelse(!is.na(FBtr) & FBtr != \"\", FBtr, sub(\"_.*\", \"\", FBtr_codonpos))]",
    "  x[, handle := data.table::fifelse(!is.na(chr) & chr != \"\" & is.finite(pos), paste(chr, format(pos, scientific = FALSE, trim = TRUE), sep = \"_\"), NA_character_)]",
    "  unique(x, by = \"FBtr_codonpos\")", "}")
code <- c(code[seq_len(reader_start - 1L)], stream_reader, code[seq.int(reader_end + 1L,
    length(code))])

# Let the streamed reader know the species of the unit currently processed.
species_line <- grep("^[[:space:]]+edit_species <- info\\$edit_species", code)
stopifnot(length(species_line) == 1L)
code <- append(code, "  current_edit_species <<- edit_species", after = species_line)

eval(parse(text = code), envir = globalenv())
