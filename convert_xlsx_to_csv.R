# ==============================================================================
# convert_xlsx_to_csv.R — stay-point .xlsx (ArcGIS "Table To Excel") -> CSV
#
#   Usage:  source("convert_xlsx_to_csv.R", encoding = "UTF-8")
#
#   The pipeline reads CSV directly and accepts small .xlsx files, but a stay-point
#   table exported from ArcGIS as .xlsx is often far too large for readxl (which
#   expands the whole sheet in memory). This script streams the XML inside the
#   .xlsx, keeps only the columns the pipeline needs, renames them to the standard
#   names, and appends the rows to a CSV. Memory stays around 1 GB.
#
#   * Several sheets are read in order and concatenated (a header row repeated on
#     later sheets is skipped).
#   * Each file in CFG$files is converted to <name>.csv in the same folder; files
#     that already have a .csv are skipped unless CFG$overwrite = TRUE.
#   * Output is written to <name>.csv.part and renamed on completion, so an
#     interrupted conversion never leaves a truncated .csv behind.
#   * A summary of all converted files (rows, users, period, days without rows,
#     rows without a town-block) is written to CFG$summary_file.
#   * Only .xlsx files with inline strings are supported (the format written by
#     ArcGIS); files with a sharedStrings table are refused.
#
#   Columns kept (output name <- source column, matched on the end of the name):
#     id <- id / datetime <- datetime / leaving_datetime <- leaving_datetime
#     lng <- lng / lat <- lat / accuracy <- accuracy / uid <- uid
#     PREF_NAME <- AddSpatialJoin.PREF_NAME / CITY_NAME <- AddSpatialJoin.CITY_NAME
#     NAME <- AddSpatialJoin.S_NAME (town-block name) / src_file <- source-file column (optional)
# ==============================================================================

library(data.table)

# ------------------------------------------------------------------------------
# 0. Settings
# ------------------------------------------------------------------------------

CFG <- list(
  dir       = "data",
  files     = c("staypoints_150m15min",       # names without extension; .xlsx is looked up in dir
                "staypoints_150m20min",
                "staypoints_200m15min",
                "staypoints_200m20min"),
  overwrite = FALSE,                          # TRUE: rebuild the .csv even if it exists
  check     = TRUE,                           # re-read each CSV and print a summary
  chunk_mb  = 32,                             # amount of uncompressed XML read at a time
  summary_file = "staypoints_conversion_summary.csv"   # written to dir
)

## output name = end of the source column name. join = TRUE marks the spatial-join columns
FIELDS <- list(
  id               = list(src = "id",               join = FALSE, required = TRUE),
  datetime         = list(src = "datetime",         join = FALSE, required = TRUE),
  leaving_datetime = list(src = "leaving_datetime", join = FALSE, required = TRUE),
  lng              = list(src = "lng",              join = FALSE, required = TRUE),
  lat              = list(src = "lat",              join = FALSE, required = TRUE),
  accuracy         = list(src = "accuracy",         join = FALSE, required = TRUE),
  uid              = list(src = "uid",              join = FALSE, required = TRUE),
  PREF_NAME        = list(src = "PREF_NAME",        join = TRUE,  required = TRUE),
  CITY_NAME        = list(src = "CITY_NAME",        join = TRUE,  required = TRUE),
  NAME             = list(src = "S_NAME",           join = TRUE,  required = TRUE),
  src_file         = list(src = c("src_file", "ファイル名"), join = FALSE, required = FALSE))

# ------------------------------------------------------------------------------
# 1. Helpers
# ------------------------------------------------------------------------------

msg <- function(...) cat(sprintf(...), "\n")

## decode XML character references (&#12501; &amp; ...) once per distinct value
unescape_xml <- function(x) {
  hit <- which(!is.na(x) & grepl("&", x, fixed = TRUE))
  if (!length(hit)) return(x)
  u <- unique(x[hit]); d <- u
  m <- gregexpr("&#[xX]?[0-9A-Fa-f]+;", d)
  regmatches(d, m) <- lapply(regmatches(d, m), function(v) vapply(v, function(e) {
    body <- substr(e, 3, nchar(e) - 1L)
    code <- if (substr(body, 1, 1) %in% c("x", "X")) strtoi(substr(body, 2, nchar(body)), 16L) else as.integer(body)
    intToUtf8(code) }, character(1), USE.NAMES = FALSE))
  for (p in list(c("&lt;", "<"), c("&gt;", ">"), c("&quot;", "\""), c("&apos;", "'"), c("&amp;", "&")))
    d <- gsub(p[1], p[2], d, fixed = TRUE)
  x[hit] <- d[match(x[hit], u)]
  x
}

## value of the cell in column `letter` for each row (rows = character vector of row XML)
cell_value <- function(rows, letter, ub) {
  pat <- paste0('<c r="', letter, '[0-9]+"[^>/]*(?:/>|>(.*?)</c>)')
  m <- regexpr(pat, rows, perl = TRUE, useBytes = ub)
  cs <- attr(m, "capture.start")[, 1]; cl <- attr(m, "capture.length")[, 1]
  v <- rep(NA_character_, length(rows)); ok <- m > 0 & cl > 0
  if (any(ok)) {
    inner <- substring(rows[ok], cs[ok], cs[ok] + cl[ok] - 1L)
    v[ok] <- gsub("<[^>]+>", "", inner, perl = TRUE, useBytes = ub)     # strip <is><t>..</t></is> / <v>..</v>
  }
  v
}

## header row (XML of one row) -> column letters and names
parse_header <- function(row, ub) {
  m <- gregexpr('<c r="([A-Z]+)[0-9]+"[^>/]*(?:/>|>(.*?)</c>)', row, perl = TRUE, useBytes = ub)[[1]]
  cs <- attr(m, "capture.start"); cl <- attr(m, "capture.length")
  letters <- substring(row, cs[, 1], cs[, 1] + cl[, 1] - 1L)
  names_  <- gsub("<[^>]+>", "", substring(row, cs[, 2], cs[, 2] + cl[, 2] - 1L), perl = TRUE)
  Encoding(names_) <- "UTF-8"
  data.table(letter = letters, name = unescape_xml(names_))
}

## decide which column letter feeds each output field
map_fields <- function(hdr) {
  out <- character()
  for (nm in names(FIELDS)) {
    f <- FIELDS[[nm]]; cand <- hdr
    hit <- cand[name %in% f$src]
    if (!nrow(hit)) hit <- cand[grepl(paste0("(^|[._])(", paste(f$src, collapse = "|"), ")$"), name)]
    if (nrow(hit) > 1) {                                           # several columns share the same ending
      pick <- hit[grepl("AddSpatialJoin", name) == f$join]
      if (nrow(pick)) hit <- pick
    }
    if (nm == "datetime") hit <- hit[!grepl("leaving_datetime$", name)]
    if (!nrow(hit)) {
      if (f$required) stop("column not found in the header: ", paste(f$src, collapse = " / "), "\n  header: ", paste(head(hdr$name, 70), collapse = ", "))
      next
    }
    out[nm] <- hit$letter[1]
    msg("  %-17s <- %s  (column %s)", nm, hit$name[1], hit$letter[1])
  }
  out
}

# ------------------------------------------------------------------------------
# 2. Convert one file:  .xlsx -> .csv.part -> .csv
# ------------------------------------------------------------------------------

convert_one <- function(xlsx, csv, chunk_mb = CFG$chunk_mb) {
  t0 <- Sys.time()
  part <- paste0(csv, ".part"); if (file.exists(part)) file.remove(part)
  zl <- utils::unzip(xlsx, list = TRUE)
  sheets <- zl[grepl("^xl/worksheets/sheet[0-9]+\\.xml$", zl$Name), ]
  sheets <- sheets[order(as.integer(gsub("\\D", "", basename(sheets$Name)))), ]
  if (!nrow(sheets)) stop("no worksheet found (not an .xlsx?): ", xlsx)
  msg("%s : %.0f MB, %d sheet(s)", basename(xlsx), file.size(xlsx) / 1024^2, nrow(sheets))
  if (any(grepl("sharedStrings", zl$Name)))
    stop("this .xlsx stores strings in a shared-string table (sharedStrings.xml), which this script does not support. Export the table from ArcGIS with Table To Excel.")

  END_ROW <- charToRaw("</row>")
  col_map <- NULL; header_names <- NULL; n_total <- 0L; first_write <- TRUE
  sheet_log <- list()

  for (si in seq_len(nrow(sheets))) {
    con <- unz(xlsx, sheets$Name[si], open = "rb")
    carry <- raw(0); n_sheet <- 0L; first_chunk <- TRUE; bytes_read <- 0
    repeat {
      new <- readBin(con, what = "raw", n = chunk_mb * 1024^2)
      bytes_read <- bytes_read + length(new)
      buf <- c(carry, new); eof <- length(new) == 0
      if (!length(buf)) break
      ends <- grepRaw(END_ROW, buf, fixed = TRUE, all = TRUE)
      if (!length(ends)) { if (eof) break; carry <- buf; next }
      cut <- ends[length(ends)] + length(END_ROW) - 1L
      chunk <- buf[seq_len(cut)]
      carry <- if (cut < length(buf)) buf[(cut + 1L):length(buf)] else raw(0)
      ub <- !any(chunk > as.raw(0x7f))                                 # pure ASCII: byte-wise processing is faster
      txt <- rawToChar(chunk); rm(chunk, buf, new)
      if (!ub) Encoding(txt) <- "UTF-8"
      rows <- strsplit(txt, "</row>", fixed = TRUE, useBytes = ub)[[1]]; rm(txt)
      if (!ub) Encoding(rows) <- "UTF-8"

      if (first_chunk) {                                               # first row of the sheet = header (check)
        first_chunk <- FALSE
        h <- parse_header(rows[1], ub)
        if (is.null(col_map)) {
          msg("sheet %d: %d header columns", si, nrow(h)); header_names <- h$name
          col_map <- map_fields(h)
          rows <- rows[-1]
        } else if (identical(h$name, header_names)) {
          rows <- rows[-1]                                             # repeated header on a later sheet
        } else msg("sheet %d: first row is not a header; read as data", si)
      }
      if (length(rows)) {
        out <- as.data.table(lapply(col_map, function(L) cell_value(rows, L, ub)))
        for (cc in names(out)) { v <- out[[cc]]; Encoding(v) <- "UTF-8"; set(out, j = cc, value = unescape_xml(v)) }
        fwrite(out, part, append = !first_write, col.names = first_write, bom = first_write)
        first_write <- FALSE
        n_sheet <- n_sheet + nrow(out); n_total <- n_total + nrow(out)
        msg("  sheet %d: %s rows (total %s)  %.0f min", si, format(n_sheet, big.mark = ","), format(n_total, big.mark = ","),
            as.numeric(difftime(Sys.time(), t0, units = "mins")))
        rm(out, rows)
      }
      if (eof) break
    }
    close(con)
    if (isTRUE(sheets$Length[si] > 0 && sheets$Length[si] < 4294967295 && bytes_read != sheets$Length[si]))
      stop(sprintf("sheet %d was not read to the end (%.0f / %.0f bytes). The file may be corrupt or still synchronizing.",
                   si, bytes_read, sheets$Length[si]))
    sheet_log[[si]] <- data.table(sheet = si, rows = n_sheet, xml_gb = round(bytes_read / 1024^3, 2))
  }
  sheet_log <- rbindlist(sheet_log)
  if (!n_total) stop("no data rows: ", xlsx)
  if (sheet_log[.N, rows] >= 1048575L)
    warning(basename(xlsx), ": the last sheet reaches the Excel row limit (1,048,576 rows); check that the table is not truncated.")
  if (file.exists(csv)) file.remove(csv)
  if (!file.rename(part, csv)) stop("rename failed: ", part, " -> ", csv)
  msg("written: %s  (%s rows, %.0f min)", csv, format(n_total, big.mark = ","), as.numeric(difftime(Sys.time(), t0, units = "mins")))
  print(sheet_log)
  invisible(sheet_log)
}

# ------------------------------------------------------------------------------
# 3. Check a CSV (re-read it and return a one-row summary)
# ------------------------------------------------------------------------------

check_csv <- function(csv, stem) {
  chk <- fread(csv, encoding = "UTF-8", select = c("id", "uid", "datetime", "PREF_NAME", "CITY_NAME", "NAME"))
  day <- as.Date(substr(chk$datetime, 1, 10))
  all_days <- seq(min(day, na.rm = TRUE), max(day, na.rm = TRUE), by = "day")
  miss <- setdiff(as.character(all_days), as.character(unique(day)))
  pref <- chk[, .(rows = .N), by = PREF_NAME][order(-rows)]
  msg("[%s] rows %s | users %s | duplicated id %s | missing id %s", stem, format(nrow(chk), big.mark = ","),
      format(uniqueN(chk$uid), big.mark = ","), format(sum(duplicated(chk$id)), big.mark = ","), sum(is.na(chk$id) | chk$id == ""))
  msg("[%s] period %s .. %s | days with data %d of %d", stem, min(day, na.rm = TRUE), max(day, na.rm = TRUE), uniqueN(day), length(all_days))
  if (length(miss)) msg("[%s] days without any row (%d): %s%s", stem, length(miss), paste(head(miss, 20), collapse = ", "), if (length(miss) > 20) ", ..." else "")
  msg("[%s] rows without a town-block (NAME empty): %s", stem, format(sum(is.na(chk$NAME) | chk$NAME == ""), big.mark = ","))
  print(pref); print(head(fread(csv, encoding = "UTF-8", nrows = 3)))
  data.table(file = stem, rows = nrow(chk), users = uniqueN(chk$uid),
             duplicated_id = sum(duplicated(chk$id)), missing_id = sum(is.na(chk$id) | chk$id == ""),
             first_day = min(day, na.rm = TRUE), last_day = max(day, na.rm = TRUE),
             days_with_data = uniqueN(day), days_in_period = length(all_days),
             days_without_rows = paste(miss, collapse = " "),
             rows_no_town_block = sum(is.na(chk$NAME) | chk$NAME == ""),
             top_pref = paste(sprintf("%s %s", head(pref$PREF_NAME, 3), format(head(pref$rows, 3), big.mark = ",", trim = TRUE)), collapse = "; "))
}

# ------------------------------------------------------------------------------
# 4. Main: process CFG$files in order
# ------------------------------------------------------------------------------

T0 <- Sys.time(); summary_rows <- list()
for (stem in CFG$files) {
  xlsx <- file.path(CFG$dir, paste0(stem, ".xlsx")); csv <- file.path(CFG$dir, paste0(stem, ".csv"))
  msg("\n==== %s ====", stem)
  status <- tryCatch({
    if (file.exists(csv) && !CFG$overwrite) {
      msg("skip: %s already exists (set CFG$overwrite = TRUE to rebuild)", basename(csv)); "skipped"
    } else if (!file.exists(xlsx)) {
      msg("missing: %s not found", basename(xlsx)); "missing"
    } else { convert_one(xlsx, csv); "converted" }
  }, error = function(e) { msg("ERROR: %s", conditionMessage(e)); "error" })
  if (CFG$check && status %in% c("converted", "skipped") && file.exists(csv)) {
    s <- tryCatch(check_csv(csv, stem), error = function(e) { msg("check ERROR: %s", conditionMessage(e)); NULL })
    if (!is.null(s)) { s[, status := status]; summary_rows[[stem]] <- s }
  } else summary_rows[[stem]] <- data.table(file = stem, status = status)
}
summary_tbl <- rbindlist(summary_rows, fill = TRUE)
fwrite(summary_tbl, file.path(CFG$dir, CFG$summary_file), bom = TRUE)
msg("\n==== summary -> %s  (%.0f min) ====", CFG$summary_file, as.numeric(difftime(Sys.time(), T0, units = "mins")))
print(summary_tbl[, intersect(c("file", "status", "rows", "users", "first_day", "last_day", "days_with_data", "days_in_period",
                                "days_without_rows", "rows_no_town_block", "duplicated_id"), names(summary_tbl)), with = FALSE])
