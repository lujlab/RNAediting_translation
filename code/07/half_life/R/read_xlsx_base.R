# Minimal dependency-free XLSX reader for ordinary rectangular worksheets.
# It is sufficient for the public PLOS S1/S2 workbooks used by this project.
xml_decode <- function(x) {
    x <- gsub("&lt;", "<", x, fixed = TRUE)
    x <- gsub("&gt;", ">", x, fixed = TRUE)
    x <- gsub("&quot;", "\"", x, fixed = TRUE)
    x <- gsub("&apos;", "'", x, fixed = TRUE)
    gsub("&amp;", "&", x, fixed = TRUE)
}
xml_blocks <- function(x, tag) {
    m <- gregexpr(paste0("(?s)<", tag, "(?:[[:space:]][^>]*)?>(.*?)</", tag, ">"), x, perl = TRUE)
    z <- regmatches(x, m)[[1]]
    if (!length(z) || identical(z, character(0)) || identical(z, ""))
        return(character())
    sub(paste0("(?s)^<", tag, "(?:[[:space:]][^>]*)?>"), "", sub(paste0("</", tag, ">$"),
        "", z, perl = TRUE), perl = TRUE)
}
xml_full_blocks <- function(x, tag) {
    m <- gregexpr(paste0("(?s)<", tag, "(?:[[:space:]][^>]*)?>.*?</", tag, ">"), x, perl = TRUE)
    z <- regmatches(x, m)[[1]]
    if (!length(z) || identical(z, character(0)) || identical(z, ""))
        character() else z
}
read_zip_text <- function(zipfile, member) {
    con <- unz(zipfile, member, open = "rb")
    on.exit(close(con), add = TRUE)
    paste(readLines(con, warn = FALSE, encoding = "UTF-8"), collapse = "")
}
excel_col_number <- function(ref) {
    s <- sub("[0-9]+$", "", ref)
    vapply(strsplit(s, ""), function(z) Reduce(function(a, b) a * 26 + match(b, LETTERS),
        z, init = 0), numeric(1))
}
read_xlsx_sheet_base <- function(file, sheet_index = 1L) {
    members <- unzip(file, list = TRUE)$Name
    shared <- character()
    if ("xl/sharedStrings.xml" %in% members) {
        si <- xml_blocks(read_zip_text(file, "xl/sharedStrings.xml"), "si")
        shared <- vapply(si, function(z) paste0(xml_decode(xml_blocks(z, "t")), collapse = ""),
            character(1))
    }
    member <- paste0("xl/worksheets/sheet", sheet_index, ".xml")
    if (!member %in% members)
        stop("Worksheet not found: ", member)
    rows <- xml_blocks(read_zip_text(file, member), "row")
    parsed <- lapply(rows, function(row) {
        cells <- xml_full_blocks(row, "c")
        if (!length(cells))
            return(setNames(character(), character()))
        refs <- sub(".*? r=\"([A-Z]+[0-9]+)\".*", "\\1", cells, perl = TRUE)
        types <- ifelse(grepl(" t=\"s\"", cells, fixed = TRUE), "s", ifelse(grepl(" t=\"inlineStr\"",
            cells, fixed = TRUE), "inline", "n"))
        vals <- vapply(seq_along(cells), function(i) {
            v <- xml_blocks(cells[i], "v")
            if (types[i] == "s" && length(v)) {
                j <- suppressWarnings(as.integer(v[1])) + 1L
                if (!is.na(j) && j <= length(shared))
                  return(shared[j])
            }
            if (types[i] == "inline")
                return(paste0(xml_decode(xml_blocks(cells[i], "t")), collapse = ""))
            if (length(v))
                xml_decode(v[1]) else ""
        }, character(1))
        setNames(vals, as.character(excel_col_number(refs)))
    })
    max_col <- max(vapply(parsed, function(z) if (length(z)) max(as.integer(names(z))) else 0L,
        integer(1)))
    mat <- matrix(NA_character_, nrow = length(parsed), ncol = max_col)
    for (i in seq_along(parsed)) if (length(parsed[[i]]))
        mat[i, as.integer(names(parsed[[i]]))] <- parsed[[i]]
    nonempty <- rowSums(!is.na(mat) & mat != "") > 0
    mat <- mat[nonempty, , drop = FALSE]
    # Choose the first row containing at least two non-empty cells as header.
    h <- which(rowSums(!is.na(mat) & mat != "") >= 2)[1]
    headers <- mat[h, ]
    headers[is.na(headers) | headers == ""] <- paste0("V", which(is.na(headers) | headers ==
        ""))
    headers <- make.unique(headers)
    out <- as.data.table(mat[(h + 1L):nrow(mat), , drop = FALSE])
    setnames(out, headers)
    out[rowSums(!is.na(out) & out != "") > 0]
}
