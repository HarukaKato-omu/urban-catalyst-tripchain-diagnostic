# ==============================================================================
# Urban catalyst or enclave? A GPS trip-chain diagnostic of public-facility
# connectivity in a Japanese suburban city center — analysis pipeline (v2.0.0)
#
# Replication code for:
#   Kato, H. "Urban catalyst or enclave? A GPS trip-chain diagnostic of
#   public-facility connectivity in a Japanese suburban city center"
#   (manuscript under review).
#
# USAGE
#   source("onikuru_tripchain_pipeline.R", encoding = "UTF-8")
#   The script runs from top to bottom for ONE stay-point specification, chosen
#   in CFG below (spec, stay_file, facility_file, out_dir). To run another
#   specification, define CFG_OVERRIDE before sourcing, e.g.
#     CFG_OVERRIDE <- list(spec = "200m20min", stay_file = "staypoints_200m20min",
#                          facility_file = "facility_staypoints_200m20min", out_dir = "output/200m20min")
#   Table S4 of the article is assembled from the tables/TableS4.csv column of
#   each run.
#
# INPUT (all in CFG$data_dir; see data/schema.csv; .csv preferred, .xlsx accepted)
#   staypoints_<spec>            stay-point table, one row per stay episode, already
#                                joined to town-block polygons (PREF_NAME, CITY_NAME, NAME)
#   facility_staypoints_<spec>   the subset of stays whose centroid lies inside the
#                                facility footprint (matched to the stay table by id)
#   residence                    residence mesh centroid per user (uid, home_mesh_lon,
#                                home_mesh_lat); used only if the stay table lacks home_mesh_*
#   The licensed GPS data cannot be redistributed; make_synthetic_inputs.R writes
#   schema-compatible synthetic files for an end-to-end test.
#
# OUTPUT (CFG$out_dir)
#   S*.csv                       intermediate results of every step (unrounded)
#   tables/Table1-3.csv          Tables 1-3 of the article
#   tables/TableS1-S6.csv        Supplementary Tables S1-S6 (S4 = this specification's
#                                column + home stays excluded + post-opening phases)
#   tables/TableNotes.csv        the numbers quoted in table notes and in the text
#   figs/Figure5-8, FigureS2     article figures (PDF + PNG); Figure4_nodes/edges.csv are
#                                the data for the GIS map (Figure 4); figure_map.csv maps
#                                article figures to files
#   diagnostics/                 internal tables and figures (CFG$write_diagnostics = TRUE)
#
# STRUCTURE
#   PART A   0 Config  1 Packages  2 Helpers  3 Load and preprocess  4 Zones
#            5 User attributes  6 Trip chains  7 Transitions  8 H1/H2  9 H3  10 H4
#            11 H5  12 H4 cluster bootstrap  13 Permutation test  14 Edge stability
#            15 Sequence analysis  16A Within-person stage 2  16B Multi-anchor benchmark
#   PART A2  16C-16M additional analyses (benchmark variants and rank bootstrap,
#            anchor-wise contrasts, distance model, difference-in-differences,
#            destinations, town-block benchmark, H2 user bootstrap, phase split)
#   PART B   B-1..B-3 diagnostic tables and figures; B-4 article tables and figures
#
# DESIGN NOTES (Section 3 of the article)
#   * Observation window 2023-10-11 to 2024-02-11; opening 2023-11-26; the New
#     Year holiday 2023-12-29 to 2024-01-03 is excluded.
#   * Two spatial definitions: the facility node (stays inside the building
#     footprint, post-opening only) and the ONIKURU-area (facility + Ekimae-3),
#     the unit of the pre/post comparisons.
#   * Daily trip chains are run-length-encoded zone sequences; the central-city
#     reach indicator is 1(unique central-city nodes >= 2) and is undirected.
#   * Town-block Fisher p-values are screening statistics on non-independent
#     transitions; user-level cluster bootstraps and the person-day permutation
#     test provide the inferential statements.
#   * Zone names are Japanese (UTF-8); the ROMA table maps them to the labels
#     used in the article. Source the script with encoding = "UTF-8".
# ==============================================================================

Sys.setenv(TZ = "Asia/Tokyo")
set.seed(20260611)

# ------------------------------------------------------------------------------
# 0. Config
# ------------------------------------------------------------------------------

CFG <- list(
  data_dir     = "data",
  spec         = "150m15min",
  stay_file    = "staypoints_150m15min",       # stay-point table (extension optional; .csv preferred over .xlsx)
  facility_file = "facility_staypoints_150m15min",  # stay points whose centroid falls inside the facility footprint
  home_file    = "residence",                   # residence mesh per user; used only if the stay-point table lacks home_mesh_*
  out_dir      = "output/150m15min",            # relative to the working directory (absolute path also accepted)

  onikuru_label      = "ONIKURU",
  onikuru_area_label = "ONIKURU_AREA",
  jr_area_label      = "JR_AREA",

  pref_name = "大阪府",
  city_name = "茨木市",
  onikuru_area_name = "駅前三丁目",
  jr_area_name      = "西駅前町",

  center_names = c(
    "駅前一丁目", "駅前二丁目", "駅前三丁目", "駅前四丁目",
    "春日一丁目", "西駅前町", "西中条町", "岩倉町", "片桐町",
    "元町", "大手町", "本町", "宮元町", "別院町", "永代町", "双葉町"
  ),
  station_cores = c("西駅前町", "永代町"),                  # the two station cores of the reading rule (JR side / Hankyu side)
  dest_groups = list(mall   = c("大手町", "別院町", "元町"),  # 16G: grouping of onward destinations
                     hankyu = c("永代町", "双葉町"),
                     jr     = c("西駅前町", "駅前一丁目")),

  open_date     = as.Date("2023-11-26"),
  baseline_from = as.Date("2023-10-11"),
  baseline_to   = as.Date("2023-11-25"),
  follow_from   = as.Date("2023-11-26"),
  follow_to     = as.Date("2024-02-11"),
  ny_exclude_from = as.Date("2023-12-29"),
  ny_exclude_to   = as.Date("2024-01-03"),

  accuracy_max_m = 50,
  secondary_windows_hours = c(2, 4),

  min_pre_observed_days = 5,
  low_propensity_rate_cut = 0,
  high_propensity_rate_cut = 0.25,

  B_boot = 500,
  B_perm = 1000,
  top_edges_for_stability = 30,

  min_anchor_traffic = 20,    # multi-anchor: minimum in+out strength
  min_anchor_days    = 10,    # multi-anchor: minimum number of anchor-days
  gis_window_km = 2.0,        # GIS: half-width (km) of the check plot centred on ONIKURU

  h2_min_inflow_edges  = 10,  # H2 figure: minimum outgoing transitions of a block to be drawn (block -> ONIKURU)
  h2_min_outflow_edges = 2,   # H2 figure: minimum transitions to a block to be drawn (ONIKURU -> block)

  run_sequence_analysis = TRUE,
  sequence_k = 4,

  # ---- residence, home stays, and additional analyses ----
  restrict_to_residents = TRUE,  # keep only users with a residence mesh (city residents)
  home_mesh_rule      = "first", # residence file has several rows per user: "first" row or "mode" (most frequent mesh)
  home_radius_m       = 200,     # a stay within this distance of the residence mesh centroid is a home stay
  home_block_radius_m = 150,     # home block = modal town-block of the stays within this distance of the centroid
  near_km             = 1.0,     # 16F: threshold for "living near" the facility (home to ONIKURU)
  min_days_ref_week   = 5,       # 16F: observed days required in the reference week of the event study
  B_rank = 2000,                 # 16C: replications of the rank bootstrap
  B_dest = 20000,                # 16G: replications of the multinomial destination simulation
  run_additional = TRUE,         # run PART A2 (additional analyses)
  sequence_max_n = 6000,         # 15: above this many sequences, identical sequences are collapsed before the distance matrix
  xlsx_max_mb = 80,              # an .xlsx larger than this is refused (supply a CSV instead)
  fig8_medians = TRUE,           # Figure 8: draw the anchor medians (dashed) and the ONIKURU bootstrap interval
  write_diagnostics = FALSE      # write the internal tables and figures to diagnostics/
)
if (exists("CFG_OVERRIDE") && is.list(CFG_OVERRIDE)) CFG[names(CFG_OVERRIDE)] <- CFG_OVERRIDE

dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Packages (fixest for the fixed-effect regressions; readxl only for .xlsx input)
# ------------------------------------------------------------------------------

REQ <- c("data.table", "lubridate", "igraph", "fixest", "TraMineR",
         "ggplot2", "patchwork")

install_if_missing <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) install.packages(miss, repos = "https://cloud.r-project.org")
  invisible(lapply(pkgs, library, character.only = TRUE))
}
install_if_missing(REQ)

# ------------------------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------------------------

msg <- function(...) cat(sprintf(...), "\n")
safe_fwrite <- function(x, file, ...) { data.table::fwrite(x, file = file, bom = TRUE, ...); msg("write: %s", file) }

## Read date-times as Japan wall-clock time. Accepts character, POSIXct (fread / readxl read as UTC)
## and Excel serial numbers; rows in a different format are re-parsed individually.
parse_dt_jst <- function(x) {
  if (inherits(x, "POSIXct")) return(lubridate::force_tz(lubridate::with_tz(x, "UTC"), "Asia/Tokyo"))
  if (inherits(x, "Date")) return(lubridate::force_tz(as.POSIXct(x), "Asia/Tokyo"))
  x <- trimws(as.character(x)); x[x %in% c("", "NA", "NaN")] <- NA_character_
  y <- suppressWarnings(lubridate::ymd_hms(x, tz = "Asia/Tokyo", quiet = TRUE))
  bad <- is.na(y) & !is.na(x)
  if (any(bad)) {
    num <- suppressWarnings(as.numeric(x[bad]))
    ser <- !is.na(num) & num > 20000 & num < 80000            # Excel serial number (days)
    z <- suppressWarnings(lubridate::parse_date_time(x[bad],
      orders = c("ymd HMS", "ymd HM", "ymd", "mdy HMS", "mdy HM"), tz = "Asia/Tokyo", quiet = TRUE))
    if (any(ser)) z[ser] <- lubridate::force_tz(
      as.POSIXct(round((num[ser] - 25569) * 86400, 3), origin = "1970-01-01", tz = "UTC"), "Asia/Tokyo")
    y[bad] <- z
  }
  y
}
clean_chr <- function(x, missing = "UNKNOWN") { y <- trimws(as.character(x)); y[is.na(y)|y==""|y=="NA"|y=="NaN"] <- missing; y }
mode_first <- function(x) {
  x_nonmiss <- x[!is.na(x)]; if (is.character(x_nonmiss)) x_nonmiss <- x_nonmiss[nzchar(trimws(x_nonmiss))]
  if (length(x_nonmiss) == 0L) return(x[NA_integer_][1L])
  ux <- unique(x_nonmiss); ux[which.max(tabulate(match(x_nonmiss, ux), nbins = length(ux)))]
}
make_zone_id <- function(pref, city, name) paste(clean_chr(pref), clean_chr(city), clean_chr(name), sep = "_")
row_normalize <- function(M) { P <- M; rs <- rowSums(P); idx <- rs > 0; if (any(idx)) P[idx,] <- P[idx,,drop=FALSE]/rs[idx]; P }

# ---- Input handling ----------------------------------------------------------

hav_m <- function(lon1, lat1, lon2, lat2) {            # distance between two points (m, haversine)
  p1 <- lat1 * pi / 180; p2 <- lat2 * pi / 180
  a <- sin((p2 - p1) / 2)^2 + cos(p1) * cos(p2) * sin((lon2 - lon1) * pi / 360)^2
  2 * 6371000 * asin(pmin(1, sqrt(a)))
}

## standard name = list(exact = accepted column names, suffix = accepted endings of prefixed column names)
## Longer names are matched first (so that lat does not capture home_mesh_lat).
COL_SPEC <- list(
  leaving_datetime = list(exact = "leaving_datetime", suffix = "leaving_datetime"),
  home_mesh_lon = list(exact = "home_mesh_lon", suffix = "home_mesh_lon"),
  home_mesh_lat = list(exact = "home_mesh_lat", suffix = "home_mesh_lat"),
  work_mesh_lon = list(exact = "work_mesh_lon", suffix = "work_mesh_lon"),
  work_mesh_lat = list(exact = "work_mesh_lat", suffix = "work_mesh_lat"),
  home_mesh     = list(exact = "home_mesh", suffix = "home_mesh"),
  work_mesh     = list(exact = "work_mesh", suffix = "work_mesh"),
  PREF_NAME     = list(exact = "PREF_NAME", suffix = "PREF_NAME"),
  CITY_NAME     = list(exact = "CITY_NAME", suffix = "CITY_NAME"),
  NAME          = list(exact = c("NAME", "S_NAME"), suffix = "S_NAME"),
  ONIKURUstay   = list(exact = "ONIKURUstay", suffix = "ONIKURUstay"),
  datetime      = list(exact = "datetime", suffix = "datetime"),
  accuracy      = list(exact = "accuracy", suffix = "accuracy"),
  uid           = list(exact = c("uid", "userid"), suffix = c("uid", "userid")),
  id            = list(exact = "id", suffix = "id"),
  lng           = list(exact = c("lng", "lon", "longitude"), suffix = c("lng", "lon")),
  lat           = list(exact = c("lat", "latitude"), suffix = "lat"),
  src_file      = list(exact = c("src_file", "ファイル名"), suffix = c("src_file", "ファイル名")))

match_cols <- function(nms) {                           # standard name -> actual column name
  out <- character(); free <- nms
  for (canon in names(COL_SPEC)) {
    sp <- COL_SPEC[[canon]]
    hit <- free[free %in% sp$exact]
    if (!length(hit)) hit <- free[tolower(free) %in% tolower(sp$exact)]
    if (!length(hit)) hit <- free[grepl(paste0("(_|\\.|の)(", paste(sp$suffix, collapse = "|"), ")$"), free)]
    if (!length(hit)) next
    out[canon] <- hit[1]; free <- setdiff(free, hit[1])
  }
  out
}
pick_cols <- function(nms, want) {                      # actual column names matching want (standard or actual names)
  m <- match_cols(nms)
  unique(c(unname(m[intersect(want, names(m))]), intersect(want, nms)))
}
std_names <- function(x, label) {
  m <- match_cols(names(x))
  for (canon in names(m)) {
    if (m[[canon]] == canon || canon %in% names(x)) next
    msg("[%s] standardized column name: %s -> %s", label, m[[canon]], canon)
    data.table::setnames(x, m[[canon]], canon)
  }
  x
}

## look for stem / stem.csv / stem.xlsx under dir (.csv preferred)
resolve_file <- function(dir, stem, required = TRUE) {
  if (is.null(stem) || is.na(stem) || !nzchar(stem)) return(NA_character_)
  cand <- c(paste0(stem, c(".csv", ".xlsx", ".txt", ".xls")), stem)
  cand <- c(file.path(dir, cand), cand)
  hit <- cand[file.exists(cand) & !dir.exists(cand)]
  if (length(hit)) return(hit[1])
  if (required) stop("file not found: ", stem, " (searched in: ", dir, ")")
  NA_character_
}

## read .csv / .xlsx and standardize column names; want keeps only the needed columns (large files).
read_table_any <- function(path, label, want = NULL) {
  ext <- tolower(tools::file_ext(path)); mb <- file.size(path) / 1024^2
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl", repos = "https://cloud.r-project.org")
    ## readxl expands the whole sheet in memory (a 62-column, 300,000-row, ~100 MB .xlsx needs > 6 GB).
    if (mb > CFG$xlsx_max_mb) stop("[", label, "] ", basename(path), " is a ", round(mb), " MB Excel file, too large to read in R",
      " (limit CFG$xlsx_max_mb = ", CFG$xlsx_max_mb, " MB).\n  Save the same table as CSV (UTF-8) in the same folder under the same name (extension .csv);",
      " the CSV is then read automatically.")
    sheets <- readxl::excel_sheets(path); xl_max <- if (ext == "xlsx") 1048575L else 65535L
    parts <- lapply(sheets, function(s) {
      x <- data.table::as.data.table(readxl::read_excel(path, sheet = s, col_types = "text"))
      msg("[%s] sheet '%s': %s rows", label, s, format(nrow(x), big.mark = ","))
      if (!is.null(want) && nrow(x)) x <- x[, pick_cols(names(x), want), with = FALSE]
      x
    })
    n <- vapply(parts, nrow, integer(1)); parts <- parts[n > 0]; n <- n[n > 0]
    if (!length(parts)) stop("[", label, "] no sheet with data rows: ", path)
    if (n[length(n)] >= xl_max) stop("[", label, "] the last sheet reaches the Excel row limit (", xl_max + 1L,
      " rows); the table is probably truncated, so export it again as CSV: ", path)
    same <- vapply(parts, function(x) identical(names(x), names(parts[[1]])), logical(1))
    if (!all(same)) msg("[%s] skipping %d sheet(s) whose columns differ from the first sheet", label, sum(!same))
    x <- data.table::rbindlist(parts[same], use.names = TRUE)
  } else {
    hdr <- names(data.table::fread(path, nrows = 0, encoding = "UTF-8"))
    sel <- if (is.null(want)) NULL else pick_cols(hdr, want)
    x <- data.table::fread(path, encoding = "UTF-8", showProgress = TRUE, select = sel, integer64 = "double")
  }
  msg("[%s] %s : %s rows, %d cols", label, basename(path), format(nrow(x), big.mark = ","), ncol(x))
  std_names(x, label)
}

calc_centrality <- function(edge_counts, target_node) {
  out <- data.table(target = target_node, n_nodes = 0L, n_edges = 0L, betweenness = NA_real_,
    out_degree = NA_integer_, in_degree = NA_integer_, out_strength = NA_real_, in_strength = NA_real_, pagerank = NA_real_)
  if (nrow(edge_counts) == 0) return(out)
  g <- igraph::graph_from_data_frame(edge_counts[, .(from, to, weight = N)], directed = TRUE)
  out[, `:=`(n_nodes = igraph::vcount(g), n_edges = igraph::ecount(g))]
  if (!(target_node %in% igraph::V(g)$name)) return(out)
  igraph::E(g)$dist <- 1 / pmax(igraph::E(g)$weight, 1e-9)
  pr <- igraph::page_rank(g, directed = TRUE, weights = igraph::E(g)$weight)$vector
  out[, `:=`(
    betweenness  = as.numeric(igraph::betweenness(g, v = target_node, directed = TRUE, weights = igraph::E(g)$dist)),
    out_degree   = as.integer(igraph::degree(g, v = target_node, mode = "out")),
    in_degree    = as.integer(igraph::degree(g, v = target_node, mode = "in")),
    out_strength = as.numeric(igraph::strength(g, vids = target_node, mode = "out", weights = igraph::E(g)$weight)),
    in_strength  = as.numeric(igraph::strength(g, vids = target_node, mode = "in", weights = igraph::E(g)$weight)),
    pagerank     = as.numeric(pr[target_node]))]
  out
}

# chain type for any anchor (the ONIKURU version is a thin wrapper)
classify_anchor_chain <- function(zone, in_center, is_anchor) {
  pos <- which(is_anchor); if (length(pos) == 0L) return(NA_character_)
  oc <- in_center & !is_anchor
  if (!any(oc, na.rm = TRUE)) return("onikuru_only_center")   # = anchor_only_center (legacy label)
  before <- any(vapply(pos, function(p) p > 1L && any(oc[seq_len(p-1L)], na.rm=TRUE), logical(1)))
  after  <- any(vapply(pos, function(p) p < length(zone) && any(oc[(p+1L):length(zone)], na.rm=TRUE), logical(1)))
  if (before && after) return("middle_bridge")
  if (after) return("post_bridge")
  if (before) return("final_after_center")
  "onikuru_only_center"
}
classify_onikuru_chain <- function(zone, in_center, is_onikuru) classify_anchor_chain(zone, in_center, is_onikuru)

decompose_ego <- function(tc_period, target, center_map, period_label, level_label) {
  out_e <- tc_period[from == target, .(node = to,   N)]
  in_e  <- tc_period[to   == target, .(node = from, N)]
  out_e <- merge(out_e, center_map, by.x = "node", by.y = "zone", all.x = TRUE)
  in_e  <- merge(in_e,  center_map, by.x = "node", by.y = "zone", all.x = TRUE)
  out_e[is.na(zone_center), zone_center := FALSE]; in_e[is.na(zone_center), zone_center := FALSE]
  ot <- sum(out_e$N); it <- sum(in_e$N)
  oc <- sum(out_e[zone_center == TRUE, N]); ic <- sum(in_e[zone_center == TRUE, N])
  data.table(level = level_label, period = period_label, target = target,
    n_out_edges = nrow(out_e), n_in_edges = nrow(in_e), out_strength = ot, in_strength = it,
    out_to_center = oc, out_to_noncenter = ot - oc, in_from_center = ic, in_from_noncenter = it - ic,
    out_center_share = if (ot>0) oc/ot else NA_real_, in_center_share = if (it>0) ic/it else NA_real_,
    flow_through = if (it>0) ot/it else NA_real_, onward_to_center = if (it>0) oc/it else NA_real_)
}
ego_edgelist <- function(tc_period, target, center_map, period_label, level_label) {
  out_e <- tc_period[from == target, .(direction = "out", node = to,   N)]
  in_e  <- tc_period[to   == target, .(direction = "in",  node = from, N)]
  e <- rbind(out_e, in_e); if (nrow(e) == 0) return(data.table())
  e <- merge(e, center_map, by.x = "node", by.y = "zone", all.x = TRUE)
  e[is.na(zone_center), zone_center := FALSE]
  e[, `:=`(level = level_label, period = period_label, target = target,
           node_type = fifelse(zone_center, "center", "noncenter"))]
  setorder(e, direction, -N); e[, .(level, period, target, direction, node, node_type, N)]
}

# ------------------------------------------------------------------------------
# 3. Load and preprocess  (+ sample-flow counts)
#    3a stay points  3b facility flag (facility file matched by id)  3c residence (joined by uid)
#    3d dates and periods  3e accuracy  3f drop pre-opening stays inside the footprint  3g observed-day diagnostics
# ------------------------------------------------------------------------------

ATTR_COLS <- c("sex","birthyear","job","position","address","driver","cartype","transportation",
  "answerdate2","maritalstatus","familylivingtogether","childlivingtogether",
  "youngestbirthyear","residencestatus","recentacademichistory","personalincome","householdincome")
HOME_COLS <- c("home_mesh","home_mesh_lon","home_mesh_lat","work_mesh","work_mesh_lon","work_mesh_lat")
STAY_COLS <- c("id","uid","datetime","leaving_datetime","lng","lat","accuracy",
               "PREF_NAME","CITY_NAME","NAME","ONIKURUstay","src_file", HOME_COLS, ATTR_COLS)

# ---- 3a. Stay points ----
stay_path <- resolve_file(CFG$data_dir, CFG$stay_file)
dt <- read_table_any(stay_path, "stay", want = STAY_COLS)

need_cols <- c("uid", "datetime", "leaving_datetime", "PREF_NAME", "CITY_NAME", "NAME")
missing_cols <- setdiff(need_cols, names(dt))
if (length(missing_cols) > 0) stop(
  "required columns missing from the stay-point file: ", paste(missing_cols, collapse = ", "),
  "\n  PREF_NAME / CITY_NAME / NAME (S_NAME) come from the spatial join with the town-block polygons.",
  "\n  A table of raw stay points (id, datetime, lng, lat, ...) must be spatially joined before it is passed to the pipeline.")

flow <- list()
flow_add <- function(lbl) flow[[length(flow)+1]] <<- data.table(step = lbl, rows = nrow(dt), uids = uniqueN(dt$uid))
flow_add("Raw stay points (file)")

dt[, uid := clean_chr(uid, missing = NA_character_)]; dt <- dt[!is.na(uid)]; flow_add("Valid uid")
for (cc in intersect(c("lng", "lat", "accuracy"), names(dt))) dt[, (cc) := suppressWarnings(as.numeric(get(cc)))]

# ---- 3b. Facility flag: match the facility stay-point file ----
oni_path  <- resolve_file(CFG$data_dir, CFG$facility_file, required = FALSE)
oni_check <- data.table(item = "source", value = "ONIKURUstay column in the stay-point file")
if (!is.na(oni_path)) {
  oni <- read_table_any(oni_path, "onikuru", want = c("id", "uid", "datetime"))
  if ("id" %in% names(dt) && "id" %in% names(oni)) {
    key_used <- "id"
    k17 <- tolower(trimws(as.character(dt$id))); k18 <- tolower(trimws(as.character(oni$id)))
  } else {
    key_used <- "uid + datetime"
    kf <- function(u, t) paste(clean_chr(u), format(parse_dt_jst(t), "%Y-%m-%d %H:%M:%S"))
    k17 <- kf(dt$uid, dt$datetime); k18 <- kf(oni$uid, oni$datetime)
  }
  flag18 <- k17 %in% k18
  if (!any(flag18)) stop("no row of the facility stay-point file was found in the stay-point table (key: ", key_used,
    "). Check that both files come from the same stay-point extraction (same specification and version).")
  flag17 <- if ("ONIKURUstay" %in% names(dt)) clean_chr(dt$ONIKURUstay, missing = "") == CFG$onikuru_label else NULL
  oni_check <- data.table(
    item  = c("source", "key", "rows_in_onikuru_file", "stay_rows_flagged", "onikuru_rows_absent_from_stay_file",
              "existing_ONIKURUstay_rows", "existing_flag_disagreements"),
    value = c(basename(oni_path), key_used, nrow(oni), sum(flag18), sum(!(k18 %in% k17)),
              if (is.null(flag17)) NA else sum(flag17), if (is.null(flag17)) NA else sum(flag17 != flag18)))
  dt[, ONIKURUstay := fifelse(flag18, CFG$onikuru_label, NA_character_)]
  msg("ONIKURU flag: %d stay rows flagged from %s (key = %s)", sum(flag18), basename(oni_path), key_used)
  rm(oni, k17, k18, flag18, flag17)
} else if (!("ONIKURUstay" %in% names(dt))) {
  stop("cannot identify facility stays: ", CFG$facility_file, " was not found and the stay-point table has no ONIKURUstay column.")
}
safe_fwrite(oni_check, file.path(CFG$out_dir, "S2d_onikuru_flag_check.csv"))

# ---- 3c. Residence mesh: join the residence file if absent, and restrict to residents ----
home_multi_users <- NA_integer_
if (!all(c("home_mesh_lon", "home_mesh_lat") %in% names(dt))) {
  home_path <- resolve_file(CFG$data_dir, CFG$home_file, required = FALSE)
  if (!is.na(home_path)) {
    hm <- read_table_any(home_path, "home", want = c("uid", HOME_COLS))
    hm[, uid := clean_chr(uid, missing = NA_character_)]; hm <- hm[!is.na(uid)]
    ## The residence file may hold several rows per user (one per estimation date). "first" = first row,
    ## "mode" = the most frequent residence mesh.
    hm[, .mesh := paste(home_mesh_lon, home_mesh_lat)]
    home_multi_users <- hm[, .(k = uniqueN(.mesh)), by = uid][k > 1, .N]
    if (identical(CFG$home_mesh_rule, "mode")) {
      hm[, .n := .N, by = .(uid, .mesh)]; hm[, .row := seq_len(.N)]
      setorder(hm, uid, -.n, .row)
    }
    hm <- unique(hm, by = "uid")
    msg("home mesh: %d of %d users have more than one home mesh in the residence file (rule = %s)",
        home_multi_users, nrow(hm), CFG$home_mesh_rule)
    add <- intersect(HOME_COLS, names(hm))
    dt[hm, on = "uid", (add) := mget(paste0("i.", add))]
    msg("home mesh: joined %s (%d users in the residence file)", basename(home_path), nrow(hm)); rm(hm)
  }
}
HAS_HOME <- all(c("home_mesh_lon", "home_mesh_lat") %in% names(dt))
if (HAS_HOME) {
  for (cc in intersect(c("home_mesh_lon","home_mesh_lat","work_mesh_lon","work_mesh_lat"), names(dt)))
    dt[, (cc) := suppressWarnings(as.numeric(get(cc)))]
  if (isTRUE(CFG$restrict_to_residents)) {
    dt <- dt[is.finite(home_mesh_lon) & is.finite(home_mesh_lat)]
    flow_add("Residents (home mesh in the residence data)")
  }
} else if (isTRUE(CFG$restrict_to_residents)) {
  stop("no residence mesh (home_mesh_lon / home_mesh_lat): the stay-point table has no such columns and ", CFG$home_file, " was not found in ", CFG$data_dir,
       ".\n  Set CFG$restrict_to_residents = FALSE to run without the resident restriction (the sample then differs from the article).")
} else warning("no residence mesh: the resident restriction and the analyses that use home location and distance are skipped.")

# ---- 3d. Dates and periods ----
dt[, arrival := parse_dt_jst(datetime)]; dt[, departure := parse_dt_jst(leaving_datetime)]
dt <- dt[!is.na(arrival)]; flow_add("Parsed datetime")
dt[is.na(departure) | departure < arrival, departure := arrival]
dt[, duration_min := as.numeric(difftime(departure, arrival, units = "mins"))]
dt[, day := as.IDate(as.Date(arrival, tz = "Asia/Tokyo"))]
dt <- dt[day >= CFG$baseline_from & day <= CFG$follow_to]; flow_add("In study window")
dt <- dt[!(day >= CFG$ny_exclude_from & day <= CFG$ny_exclude_to)]; flow_add("Exclude New Year")
dt[, period := fifelse(day <= CFG$baseline_to, "pre", "post")]
dt[, post := as.integer(period == "post")]
dt[, week := as.integer(floor(as.integer(as.Date(day) - CFG$baseline_from) / 7L))]
dt[, dow := lubridate::wday(as.Date(day), week_start = 1)]

# ---- 3e. Accuracy ----
if ("accuracy" %in% names(dt)) {
  dt <- dt[!is.na(accuracy) & accuracy < CFG$accuracy_max_m]
} else warning("accuracy column is absent; accuracy filter skipped.")
flow_add(sprintf("Accuracy < %d m", CFG$accuracy_max_m))

for (cc in c("PREF_NAME", "CITY_NAME", "NAME")) dt[, (cc) := clean_chr(get(cc))]

# ---- 3f. Drop pre-opening stays inside the footprint (all later counts refer to the remaining rows) ----
dt[, is_onikuru_raw := clean_chr(ONIKURUstay, missing = "") == CFG$onikuru_label]
dt[is.na(is_onikuru_raw), is_onikuru_raw := FALSE]
n_pre_onikuru <- dt[is_onikuru_raw == TRUE & period == "pre", .N]
dt <- dt[!(is_onikuru_raw == TRUE & period == "pre")]; flow_add("Remove pre-opening ONIKURU pts")
dt[, is_onikuru := is_onikuru_raw & period == "post"]
msg("removed pre-opening ONIKURU facility points: %s", n_pre_onikuru)

# ---- 3g. Observed-day diagnostics ----
cal <- data.table(day = as.IDate(seq(CFG$baseline_from, CFG$follow_to, by = "day")))
cal[, period := fifelse(day <= CFG$baseline_to, "pre", "post")]
cal[, ny_excluded := day >= CFG$ny_exclude_from & day <= CFG$ny_exclude_to]
obs_days <- dt[, .(rows = .N, uids = uniqueN(uid), onikuru_rows = sum(is_onikuru)), by = day][cal, on = "day"]
for (cc in c("rows", "uids", "onikuru_rows")) obs_days[is.na(get(cc)), (cc) := 0L]
obs_days[, observed := rows > 0]
setorder(obs_days, day)
safe_fwrite(obs_days[, .(day, period, ny_excluded, observed, rows, uids, onikuru_rows)],
            file.path(CFG$out_dir, "S2b_observed_days.csv"))
gap_days <- obs_days[ny_excluded == FALSE & observed == FALSE]
obs_gaps <- if (nrow(gap_days)) {
  gap_days[, run := cumsum(c(1L, diff(as.integer(day)) != 1L))]
  gap_days[, .(from = min(day), to = max(day), n_days = .N, period = paste(unique(period), collapse = "/")), by = run][, run := NULL][]
} else data.table(from = as.IDate(character()), to = as.IDate(character()), n_days = integer(), period = character())
safe_fwrite(obs_gaps, file.path(CFG$out_dir, "S2c_observation_gaps.csv"))
obs_sum <- obs_days[ny_excluded == FALSE, .(calendar_days = .N, observed_days = sum(observed)), by = period]
n_obs <- function(pp, col) { v <- obs_sum[period == pp][[col]]; if (length(v)) v else 0L }
msg("observed days: pre %d of %d, post %d of %d (New Year excluded)",
    n_obs("pre","observed_days"), n_obs("pre","calendar_days"), n_obs("post","observed_days"), n_obs("post","calendar_days"))
if (nrow(obs_gaps)) {
  msg("!! %d day(s) in the observation window have no stay points at all (S2c_observation_gaps.csv):", sum(obs_gaps$n_days))
  print(obs_gaps)
  warning("days without data inside the observation window; check S2c_observation_gaps.csv for missing daily files.")
}

pref_tab <- dt[, .(rows = .N), by = PREF_NAME][order(-rows)][, share := round(rows / sum(rows), 4)][]
safe_fwrite(pref_tab, file.path(CFG$out_dir, "S2e_prefecture_composition.csv"))
if (pref_tab[PREF_NAME != CFG$pref_name & PREF_NAME != "UNKNOWN", sum(rows)] == 0)
  msg("!! no stay points outside %s. If the input was clipped to the prefecture, trips beyond it do not appear in the chains.", CFG$pref_name)

diag_s2 <- data.table(
  item = c("rows","uids","date_min","date_max","pre_rows","post_rows","duration_lt_15min_rows",
           "pre_calendar_days","pre_observed_days","post_calendar_days","post_observed_days",
           "missing_days_in_window","removed_pre_opening_onikuru_rows","rows_outside_prefecture"),
  value = c(nrow(dt), uniqueN(dt$uid), as.character(min(dt$day)), as.character(max(dt$day)),
            dt[period=="pre", .N], dt[period=="post", .N], dt[duration_min < 15, .N],
            n_obs("pre","calendar_days"), n_obs("pre","observed_days"),
            n_obs("post","calendar_days"), n_obs("post","observed_days"),
            sum(obs_gaps$n_days), n_pre_onikuru, dt[PREF_NAME != CFG$pref_name, .N]))
safe_fwrite(diag_s2, file.path(CFG$out_dir, "S2_data_diagnostics.csv"))

# ------------------------------------------------------------------------------
# 4. Zone assignment  (+ zone centroids for GIS, + home / work stays)
# ------------------------------------------------------------------------------

dt[, is_ibaraki := PREF_NAME == CFG$pref_name & CITY_NAME == CFG$city_name]
dt[, base_zone_id := make_zone_id(PREF_NAME, CITY_NAME, NAME)]
dt[, base_zone_label := paste(CITY_NAME, NAME, sep = "_")]
dt[, zone_id := fifelse(is_onikuru, CFG$onikuru_label, base_zone_id)]
dt[, zone_label := fifelse(is_onikuru, CFG$onikuru_label, base_zone_label)]
dt[, is_onikuru_area := is_ibaraki & NAME == CFG$onikuru_area_name]
dt[, is_jr_area := is_ibaraki & NAME == CFG$jr_area_name]
dt[, zone_area_id := fifelse(is_onikuru | is_onikuru_area, CFG$onikuru_area_label,
                             fifelse(is_jr_area, CFG$jr_area_label, base_zone_id))]
dt[, zone_area_label := fifelse(is_onikuru | is_onikuru_area, CFG$onikuru_area_label,
                                fifelse(is_jr_area, CFG$jr_area_label, base_zone_label))]
dt[, in_center := is_onikuru | (is_ibaraki & NAME %in% CFG$center_names)]
if (dt[in_center == TRUE & is_onikuru == FALSE, .N] == 0)
  stop("no stay point falls in a central town-block. Check the values and encoding of PREF_NAME / CITY_NAME / NAME.")

zone_diag <- data.table(
  item = c("ONIKURU_facility_post","ONIKURU_AREA_pre","ONIKURU_AREA_post",
           "JR_AREA_pre","JR_AREA_post","center_pre","center_post"),
  value = c(dt[is_onikuru == TRUE & period == "post", .N],
            dt[zone_area_id == CFG$onikuru_area_label & period == "pre", .N],
            dt[zone_area_id == CFG$onikuru_area_label & period == "post", .N],
            dt[zone_area_id == CFG$jr_area_label & period == "pre", .N],
            dt[zone_area_id == CFG$jr_area_label & period == "post", .N],
            dt[in_center == TRUE & period == "pre", .N],
            dt[in_center == TRUE & period == "post", .N]))
safe_fwrite(zone_diag, file.path(CFG$out_dir, "S3_zone_diagnostics.csv"))
if (dt[is_onikuru == TRUE & NAME != CFG$onikuru_area_name, .N] > 0)
  msg("!! %d facility-flagged stay points lie in a town-block other than %s.", dt[is_onikuru == TRUE & NAME != CFG$onikuru_area_name, .N], CFG$onikuru_area_name)

# zone centroids (if stay-point coordinates are available)
HAS_XY <- all(c("lng", "lat") %in% names(dt))
zone_coords <- NULL
if (HAS_XY) {
  dt[, .lon := lng]; dt[, .lat := lat]
  zone_coords <- dt[is.finite(.lon) & is.finite(.lat),
    .(lon = median(.lon), lat = median(.lat), is_center = any(in_center),
      is_onikuru = any(is_onikuru), n = .N), by = zone_id]
  safe_fwrite(zone_coords, file.path(CFG$out_dir, "S_zone_coords.csv"))
} else msg("GIS: coordinate columns (lng / lat) not found; zone centroids skipped")

# home and work stays (for the additional analyses; not used in the PART A definitions)
HAS_HOME <- HAS_HOME && HAS_XY
dt[, `:=`(is_home = FALSE, is_work = FALSE)]
if (HAS_HOME) {
  dt[, d_home := hav_m(.lon, .lat, home_mesh_lon, home_mesh_lat)]
  dt[, is_home := !is.na(d_home) & d_home <= CFG$home_radius_m & !is_onikuru]
  if (all(c("work_mesh_lon", "work_mesh_lat") %in% names(dt))) {
    dt[, d_work := hav_m(.lon, .lat, work_mesh_lon, work_mesh_lat)]
    dt[, is_work := !is.na(d_work) & d_work <= CFG$home_radius_m & !is_onikuru & !is_home]
  }
}

# ------------------------------------------------------------------------------
# 5. UID attributes
# ------------------------------------------------------------------------------

attr_cols <- intersect(c(
  "sex","birthyear","job","position","address","driver","cartype","transportation",
  "answerdate2","maritalstatus","familylivingtogether","childlivingtogether",
  "youngestbirthyear","residencestatus","recentacademichistory","personalincome",
  "householdincome","home_mesh","home_mesh_lon","home_mesh_lat",
  "work_mesh","work_mesh_lon","work_mesh_lat"), names(dt))
uid_attr <- dt[, lapply(.SD, mode_first), by = uid, .SDcols = attr_cols]
if ("birthyear" %in% names(uid_attr)) {
  uid_attr[, birthyear_num := suppressWarnings(as.integer(birthyear))]
  uid_attr[, age := lubridate::year(CFG$open_date) - birthyear_num]
  uid_attr[age < 10 | age > 100, age := NA_integer_]
} else uid_attr[, age := NA_integer_]
safe_fwrite(uid_attr, file.path(CFG$out_dir, "S3_uid_attributes.csv"))

# ------------------------------------------------------------------------------
# 6. Build trip chains  (+ finalize Sample Flow)
# ------------------------------------------------------------------------------

build_chain <- function(d, zone_col, label_col) {
  dd <- copy(d); dd[, z := get(zone_col)]; dd[, z_label := get(label_col)]
  setorder(dd, uid, day, arrival, departure)
  dd[, run_id := data.table::rleid(z), by = .(uid, day)]
  cp <- dd[, .(period = period[1], post = post[1], week = week[1], dow = dow[1],
               zone = z[1], zone_label = z_label[1],
               arrival = min(arrival, na.rm = TRUE), departure = max(departure, na.rm = TRUE),
               duration_min = sum(duration_min, na.rm = TRUE),
               in_center = any(in_center, na.rm = TRUE), is_onikuru = any(is_onikuru, na.rm = TRUE)),
           by = .(uid, day, run_id)]
  setorder(cp, uid, day, arrival, departure); cp[, ord := seq_len(.N), by = .(uid, day)]; cp[]
}

chain_facility <- build_chain(dt, "zone_id", "zone_label")
chain_area     <- build_chain(dt, "zone_area_id", "zone_area_label")

chain_day <- chain_facility[, .(
  chain_len = .N, unique_nodes = uniqueN(zone),
  has_onikuru = any(is_onikuru), has_center = any(in_center),
  center_stays = sum(in_center), unique_center_nodes = uniqueN(zone[in_center])),
  by = .(uid, day, period, post, week, dow)]
post_visitors <- chain_day[period == "post" & has_onikuru == TRUE, unique(uid)]
chain_day[, visitor := uid %in% post_visitors]

flow_tbl <- rbindlist(flow)
flow_tbl <- rbind(flow_tbl, data.table(
  step = c("Trip-chain person-days", "Post person-days", "ONIKURU post visitors (uid)"),
  rows = c(nrow(chain_day), chain_day[period=="post", .N], length(post_visitors)),
  uids = c(uniqueN(chain_day$uid), chain_day[period=="post", uniqueN(uid)], length(post_visitors))))
safe_fwrite(flow_tbl, file.path(CFG$out_dir, "S_sample_flow.csv"))

safe_fwrite(chain_facility, file.path(CFG$out_dir, "S4_chain_points_facility.csv"))
safe_fwrite(chain_area, file.path(CFG$out_dir, "S4_chain_points_area.csv"))
safe_fwrite(chain_day, file.path(CFG$out_dir, "S4_person_day_chain_summary.csv"))
chain_summary <- chain_day[, .(person_days = .N, uids = uniqueN(uid),
  mean_chain_len = mean(chain_len), mean_unique_nodes = mean(unique_nodes),
  onikuru_days = sum(has_onikuru), center_days = sum(has_center)), by = period]
safe_fwrite(chain_summary, file.path(CFG$out_dir, "S4_chain_summary_by_period.csv"))

# ------------------------------------------------------------------------------
# 7. Transitions
# ------------------------------------------------------------------------------

make_edges <- function(cp) {
  out <- cp[, { if (.N >= 2L) .(from = zone[-.N], to = zone[-1L]) else NULL }, by = .(uid, day, period, post, week)]
  if (nrow(out) == 0) return(data.table(uid=character(), day=as.IDate(character()), period=character(),
    post=integer(), week=integer(), from=character(), to=character()))
  out[]
}
edges_facility <- make_edges(chain_facility); edges_area <- make_edges(chain_area)
safe_fwrite(edges_facility, file.path(CFG$out_dir, "S5_edges_facility.csv"))
safe_fwrite(edges_area, file.path(CFG$out_dir, "S5_edges_area.csv"))
tc_facility <- edges_facility[, .(N = .N), by = .(period, from, to)]
tc_area     <- edges_area[, .(N = .N), by = .(period, from, to)]
safe_fwrite(tc_facility, file.path(CFG$out_dir, "S5_transition_counts_facility.csv"))
safe_fwrite(tc_area, file.path(CFG$out_dir, "S5_transition_counts_area.csv"))

# ------------------------------------------------------------------------------
# 8. H1 / H2
#    H2 = every central town-block <-> ONIKURU_AREA (both directions: inflow and outflow)
#         inflow : share of transitions leaving block z that go to ONIKURU_AREA
#         outflow: share of transitions leaving ONIKURU_AREA that go to block z
#         per block + aggregate (ALL_CENTER); pre/post Fisher tests with BH correction
# ------------------------------------------------------------------------------

all_periods <- data.table(period = c("pre", "post"))

# central town-blocks in the area zoning (ONIKURU_AREA itself is the destination and is excluded)
center_zones_area <- unique(chain_area[in_center == TRUE, zone])
center_zones_area <- setdiff(center_zones_area, CFG$onikuru_area_label)
# display label of an area zone (JR_AREA / ONIKURU_AREA as is; otherwise strip the prefecture and city)
area_zone_label <- function(z) ifelse(
  z %in% c(CFG$onikuru_area_label, CFG$jr_area_label), z, sub("^[^_]+_[^_]+_", "", z))

# build the pre/post 2x2 table from agg(period, out_edges, hit_edges) and run the Fisher test
.fisher_prepost <- function(agg) {
  agg <- agg[all_periods, on = "period"]
  agg[is.na(out_edges), out_edges := 0L]; agg[is.na(hit_edges), hit_edges := 0L]
  agg[, p_hit := fifelse(out_edges > 0, hit_edges / out_edges, NA_real_)]
  pa <- agg[period=="pre",  hit_edges]; pn <- agg[period=="pre",  out_edges - hit_edges]
  qa <- agg[period=="post", hit_edges]; qn <- agg[period=="post", out_edges - hit_edges]
  fp <- NA_real_; orr <- NA_real_
  if ((pa + pn) > 0 && (qa + qn) > 0 && (pa + qa) > 0) {
    ft <- tryCatch(fisher.test(matrix(c(pa, pn, qa, qn), nrow = 2, byrow = TRUE)), error = function(e) NULL)
    if (!is.null(ft)) { fp <- ft$p.value; orr <- unname(ft$estimate) }
  }
  data.table(
    pre_out_edges  = agg[period=="pre",  out_edges], pre_hit_edges  = agg[period=="pre",  hit_edges], pre_share  = agg[period=="pre",  p_hit],
    post_out_edges = agg[period=="post", out_edges], post_hit_edges = agg[period=="post", hit_edges], post_share = agg[period=="post", p_hit],
    share_change = agg[period=="post", p_hit] - agg[period=="pre", p_hit], fisher_p = fp, odds_ratio = orr)
}

# one-row summary for block z and direction
h2_zone_one <- function(z, direction) {
  if (direction == "inflow") {      # z -> ONIKURU : share of the outflow of z that goes to ONIKURU
    agg <- edges_area[from == z, .(out_edges = .N, hit_edges = sum(to == CFG$onikuru_area_label)), by = period]
  } else {                           # ONIKURU -> z : share of the outflow of ONIKURU that goes to z
    agg <- edges_area[from == CFG$onikuru_area_label, .(out_edges = .N, hit_edges = sum(to == z)), by = period]
  }
  cbind(data.table(zone = z, zone_label = area_zone_label(z), direction = direction), .fisher_prepost(agg))
}

h2_center <- rbindlist(lapply(center_zones_area, function(z)
  rbind(h2_zone_one(z, "inflow"), h2_zone_one(z, "outflow"))), fill = TRUE)
# Benjamini-Hochberg (within each direction; NA excluded)
h2_center[, fisher_p_BH := {
  v <- rep(NA_real_, .N); ok <- !is.na(fisher_p)
  if (any(ok)) v[ok] <- p.adjust(fisher_p[ok], method = "BH"); v }, by = direction]

# aggregate over all central town-blocks
agg_inflow  <- edges_area[from %in% center_zones_area,
  .(out_edges = .N, hit_edges = sum(to == CFG$onikuru_area_label)), by = period]
agg_outflow <- edges_area[from == CFG$onikuru_area_label,
  .(out_edges = .N, hit_edges = sum(to %in% center_zones_area)), by = period]
h2_agg <- rbind(
  cbind(data.table(zone = "ALL_CENTER", zone_label = "ALL_CENTER", direction = "inflow"),  .fisher_prepost(agg_inflow)),
  cbind(data.table(zone = "ALL_CENTER", zone_label = "ALL_CENTER", direction = "outflow"), .fisher_prepost(agg_outflow)))
h2_agg[, fisher_p_BH := NA_real_]

h2_all <- rbind(h2_agg, h2_center, fill = TRUE)
setorder(h2_all, direction, -post_share)
safe_fwrite(h2_all, file.path(CFG$out_dir, "S6_H2_center_to_ONIKURU_AREA_prepost.csv"))

# backward compatibility: also save the inflow row of the JR station core under the old file name
h2_jr <- h2_center[zone_label == CFG$jr_area_label & direction == "inflow"]
safe_fwrite(h2_jr, file.path(CFG$out_dir, "S6_H2_JR_to_ONIKURU_AREA_prepost.csv"))

central_area <- rbindlist(lapply(c("pre","post"), function(pp)
  cbind(period = pp, calc_centrality(tc_area[period == pp, .(from, to, N)], CFG$onikuru_area_label))), fill = TRUE)
central_facility_post <- calc_centrality(tc_facility[period == "post", .(from, to, N)], CFG$onikuru_label)
central_facility_post[, period := "post"]
safe_fwrite(central_area, file.path(CFG$out_dir, "S6_H1_centrality_ONIKURU_AREA.csv"))
safe_fwrite(central_facility_post, file.path(CFG$out_dir, "S6_H1_centrality_ONIKURU_facility_post.csv"))

# ------------------------------------------------------------------------------
# 8b. outflow decomposition
# ------------------------------------------------------------------------------

zone_center_map_fac  <- chain_facility[, .(zone_center = any(in_center, na.rm = TRUE)), by = .(zone)]
zone_center_map_area <- chain_area[,     .(zone_center = any(in_center, na.rm = TRUE)), by = .(zone)]
ego_fac_post  <- decompose_ego(tc_facility[period=="post", .(from,to,N)], CFG$onikuru_label, zone_center_map_fac, "post", "facility")
ego_fac_edges <- ego_edgelist(tc_facility[period=="post", .(from,to,N)], CFG$onikuru_label, zone_center_map_fac, "post", "facility")
ego_area <- rbindlist(lapply(c("pre","post"), function(pp)
  decompose_ego(tc_area[period==pp, .(from,to,N)], CFG$onikuru_area_label, zone_center_map_area, pp, "area")), fill = TRUE)
ego_area_edges <- rbindlist(lapply(c("pre","post"), function(pp)
  ego_edgelist(tc_area[period==pp, .(from,to,N)], CFG$onikuru_area_label, zone_center_map_area, pp, "area")), fill = TRUE)
ego_decomp <- rbindlist(list(ego_fac_post, ego_area), fill = TRUE)
ego_edges  <- rbindlist(list(ego_fac_edges, ego_area_edges), fill = TRUE)
safe_fwrite(ego_decomp, file.path(CFG$out_dir, "S6b_onikuru_outflow_decomposition.csv"))
safe_fwrite(ego_edges,  file.path(CFG$out_dir, "S6b_onikuru_ego_edgelist.csv"))

# ------------------------------------------------------------------------------
# 9. H3
# ------------------------------------------------------------------------------

pre_propensity <- chain_day[period == "pre", .(pre_observed_days = .N, pre_center_days = sum(has_center),
  pre_center_day_rate = mean(has_center), pre_mean_chain_len = mean(chain_len)), by = uid]
pre_propensity <- pre_propensity[pre_observed_days >= CFG$min_pre_observed_days]
pre_propensity[, pre_center_group := fifelse(pre_center_day_rate <= CFG$low_propensity_rate_cut, "low_pre_center",
  fifelse(pre_center_day_rate >= CFG$high_propensity_rate_cut, "high_pre_center", "middle_pre_center"))]
h3 <- chain_day[period == "post"][pre_propensity, on = "uid", nomatch = 0]
h3_summary <- h3[, .(uids = uniqueN(uid), person_days = .N, onikuru_days = sum(has_onikuru),
  onikuru_day_rate = mean(has_onikuru), visitor_uids = uniqueN(uid[has_onikuru == TRUE]),
  mean_chain_len = mean(chain_len), mean_unique_center_nodes = mean(unique_center_nodes)),
  by = pre_center_group][order(pre_center_group)]
safe_fwrite(pre_propensity, file.path(CFG$out_dir, "S7_H3_pre_center_propensity.csv"))
safe_fwrite(h3_summary, file.path(CFG$out_dir, "S7_H3_generative_summary.csv"))

## H3 logistic regression (crude / adjusted; OR per +0.10 of the pre-opening central-city use rate)
h3_frame <- merge(
  pre_propensity[, .(uid, pre_observed_days, pre_center_day_rate)],
  chain_day[period == "post", .(post_observed_days = .N, became_visitor = as.integer(any(has_onikuru))), by = uid], by = "uid")
tidy_logit <- function(m, label, scale = 0.10) {
  co <- summary(m)$coefficients; co <- co[rownames(co) != "(Intercept)", , drop = FALSE]
  s <- ifelse(rownames(co) == "pre_center_day_rate", scale, ifelse(rownames(co) == "log_km", log(2), 1))
  data.table(model = label, term = rownames(co),
    unit = ifelse(rownames(co) == "pre_center_day_rate", sprintf("per +%.2f", scale),
           ifelse(rownames(co) == "log_km", "per doubling of distance", "per day")),
    OR = exp(s * co[, 1]), OR_low = exp(s * (co[, 1] - 1.96 * co[, 2])), OR_high = exp(s * (co[, 1] + 1.96 * co[, 2])),
    p = co[, 4], n_users = stats::nobs(m), n_visitors = sum(m$y))
}
h3_logit <- data.table(); h3_m1 <- NULL
if (nrow(h3_frame) > 0 && uniqueN(h3_frame$became_visitor) == 2) {
  h3_m0 <- glm(became_visitor ~ pre_center_day_rate, data = h3_frame, family = binomial())
  h3_m1 <- glm(became_visitor ~ pre_center_day_rate + pre_observed_days + post_observed_days, data = h3_frame, family = binomial())
  h3_logit <- rbind(tidy_logit(h3_m0, "crude"), tidy_logit(h3_m1, "adjusted"))
  safe_fwrite(h3_logit, file.path(CFG$out_dir, "S7_H3_logistic.csv"))
}

# ------------------------------------------------------------------------------
# 10. H4
# ------------------------------------------------------------------------------

oni_first <- chain_facility[is_onikuru == TRUE, .(oni_ord = min(ord),
  oni_arrival = arrival[which.min(ord)], oni_departure = departure[which.min(ord)]), by = .(uid, day)]
post_sec_same_day <- chain_facility[oni_first, on = .(uid, day), allow.cartesian = TRUE][
  ord > oni_ord & in_center == TRUE & is_onikuru == FALSE,
  .(post_secondary_same_day = 1L, n_post_secondary_same_day = uniqueN(zone), first_secondary_arrival = min(arrival)), by = .(uid, day)]
chain_day[post_sec_same_day, on = .(uid, day), `:=`(post_secondary_same_day = i.post_secondary_same_day,
  n_post_secondary_same_day = i.n_post_secondary_same_day, first_secondary_arrival = i.first_secondary_arrival)]
chain_day[is.na(post_secondary_same_day), post_secondary_same_day := 0L]
chain_day[is.na(n_post_secondary_same_day), n_post_secondary_same_day := 0L]
for (w in CFG$secondary_windows_hours) {
  col <- paste0("post_secondary_", w, "h"); ncol <- paste0("n_post_secondary_", w, "h")
  tmp <- chain_facility[oni_first, on = .(uid, day), allow.cartesian = TRUE][
    ord > oni_ord & in_center == TRUE & is_onikuru == FALSE & arrival <= oni_departure + lubridate::hours(w),
    .(tmp_secondary = 1L, tmp_n_secondary = uniqueN(zone)), by = .(uid, day)]
  setnames(tmp, c("tmp_secondary","tmp_n_secondary"), c(col, ncol))
  chain_day[tmp, on = .(uid, day), `:=`(tmp_secondary_join = get(paste0("i.", col)), tmp_n_secondary_join = get(paste0("i.", ncol)))]
  setnames(chain_day, "tmp_secondary_join", col); setnames(chain_day, "tmp_n_secondary_join", ncol)
  chain_day[is.na(get(col)), (col) := 0L]; chain_day[is.na(get(ncol)), (ncol) := 0L]
}
h4_onikuru_days <- chain_day[period == "post" & has_onikuru == TRUE, .(
  onikuru_days = .N, same_day_secondary_rate = mean(post_secondary_same_day),
  same_day_secondary_n_mean = mean(n_post_secondary_same_day),
  secondary_2h_rate = if ("post_secondary_2h" %in% names(chain_day)) mean(post_secondary_2h) else NA_real_,
  secondary_4h_rate = if ("post_secondary_4h" %in% names(chain_day)) mean(post_secondary_4h) else NA_real_,
  mean_chain_len = mean(chain_len), mean_unique_center_nodes = mean(unique_center_nodes))]
h4_within_visitor <- chain_day[period == "post" & visitor == TRUE, .(
  person_days = .N, mean_chain_len = mean(chain_len), mean_unique_nodes = mean(unique_nodes),
  mean_unique_center_nodes = mean(unique_center_nodes), center_day_rate = mean(has_center)), by = has_onikuru]
safe_fwrite(h4_onikuru_days, file.path(CFG$out_dir, "S8_H4_ordered_secondary_summary.csv"))
safe_fwrite(h4_within_visitor, file.path(CFG$out_dir, "S8_H4_visitor_day_comparison.csv"))

# ------------------------------------------------------------------------------
# 11. H5
# ------------------------------------------------------------------------------

chain_type <- chain_facility[, .(chain_type = classify_onikuru_chain(zone, in_center, is_onikuru)), by = .(uid, day)]
chain_day[chain_type, on = .(uid, day), chain_type := i.chain_type]
chain_type_summary <- chain_day[period == "post" & has_onikuru == TRUE, .N, by = chain_type]
chain_type_summary[, share := N / sum(N)]; setorder(chain_type_summary, -share)
safe_fwrite(chain_type_summary, file.path(CFG$out_dir, "S9_H5_chain_type_summary.csv"))

h5_indicators <- data.table(
  indicator = c("post_onikuru_days","onikuru_only_center_share","final_after_center_share",
    "post_bridge_share","middle_bridge_share","bridge_total_share","same_day_secondary_rate",
    "secondary_2h_rate","secondary_4h_rate","ONIKURU_out_degree_post","ONIKURU_out_strength_post",
    "ONIKURU_betweenness_post","ONIKURU_fac_out_center_share_post","ONIKURU_fac_onward_to_center_post",
    "ONIKURU_fac_flow_through_post","ONIKURU_AREA_out_center_share_pre","ONIKURU_AREA_out_center_share_post",
    "ONIKURU_AREA_onward_to_center_pre","ONIKURU_AREA_onward_to_center_post"),
  value = c(chain_day[period=="post" & has_onikuru==TRUE, .N],
    chain_type_summary[chain_type=="onikuru_only_center", share],
    chain_type_summary[chain_type=="final_after_center", share],
    chain_type_summary[chain_type=="post_bridge", share],
    chain_type_summary[chain_type=="middle_bridge", share],
    chain_type_summary[chain_type %in% c("post_bridge","middle_bridge"), sum(share)],
    h4_onikuru_days$same_day_secondary_rate, h4_onikuru_days$secondary_2h_rate, h4_onikuru_days$secondary_4h_rate,
    central_facility_post$out_degree, central_facility_post$out_strength, central_facility_post$betweenness,
    ego_fac_post$out_center_share, ego_fac_post$onward_to_center, ego_fac_post$flow_through,
    ego_area[period=="pre", out_center_share], ego_area[period=="post", out_center_share],
    ego_area[period=="pre", onward_to_center], ego_area[period=="post", onward_to_center]))
h5_indicators[grepl("share|rate", indicator) & is.na(value), value := 0]
safe_fwrite(h5_indicators, file.path(CFG$out_dir, "S10_integrated_catalyst_enclave_indicators.csv"))

# ------------------------------------------------------------------------------
# 12. Cluster bootstrap for H4
# ------------------------------------------------------------------------------

boot_h4 <- function(dat, B = 500) {
  dat <- dat[period == "post" & has_onikuru == TRUE]; ids <- unique(dat$uid)
  if (length(ids) < 2 || nrow(dat) == 0) return(data.table())
  res <- rbindlist(lapply(seq_len(B), function(b) {
    smp <- sample(ids, length(ids), replace = TRUE)
    w <- data.table(uid = names(table(smp)), weight = as.numeric(table(smp)))
    d <- dat[w, on = "uid", nomatch = 0]
    data.table(b = b, same_day_secondary_rate = weighted.mean(d$post_secondary_same_day, d$weight),
      secondary_2h_rate = if ("post_secondary_2h" %in% names(d)) weighted.mean(d$post_secondary_2h, d$weight) else NA_real_,
      secondary_4h_rate = if ("post_secondary_4h" %in% names(d)) weighted.mean(d$post_secondary_4h, d$weight) else NA_real_)
  }))
  ## point estimates are the sample values; intervals are bootstrap percentiles.
  point <- c(same_day_secondary_rate = mean(dat$post_secondary_same_day),
    secondary_2h_rate = if ("post_secondary_2h" %in% names(dat)) mean(dat$post_secondary_2h) else NA_real_,
    secondary_4h_rate = if ("post_secondary_4h" %in% names(dat)) mean(dat$post_secondary_4h) else NA_real_)
  out <- melt(res, id.vars="b", variable.name="metric", value.name="value")[
    , .(boot_mean = mean(value, na.rm=TRUE), ci_low = quantile(value,0.025,na.rm=TRUE),
        ci_high = quantile(value,0.975,na.rm=TRUE)), by = metric]
  out[, estimate := unname(point[as.character(metric)])]
  out[, .(metric, estimate, ci_low, ci_high, boot_mean)]
}
h4_boot_ci <- boot_h4(chain_day, CFG$B_boot)
safe_fwrite(h4_boot_ci, file.path(CFG$out_dir, "S11_H4_cluster_bootstrap_CI.csv"))

# ------------------------------------------------------------------------------
# 13. Network distance permutation test
# ------------------------------------------------------------------------------

if (nrow(edges_area) > 0 && length(unique(edges_area$period)) == 2) {
  states <- sort(unique(c(edges_area$from, edges_area$to))); ns <- length(states)
  E <- copy(edges_area)
  E[, fi := match(from, states)]; E[, ti := match(to, states)]
  E[, k := fi + (ti - 1L) * ns]; E[, udid := .GRP, by = .(uid, day)]
  fi_vec <- E$fi; udid <- E$udid
  cellf <- factor(E$k); kc <- as.integer(cellf); nc <- nlevels(cellf)
  cell_fi <- fi_vec[match(seq_len(nc), kc)]
  ud <- unique(E[, .(udid, period)]); setorder(ud, udid); post_ud <- as.integer(ud$period == "post")
  dist_from_postedge <- function(pe) {
    cp <- tabulate(kc[pe], nbins = nc); cq <- tabulate(kc[!pe], nbins = nc)
    rp <- tabulate(fi_vec[pe], nbins = ns); rq <- tabulate(fi_vec[!pe], nbins = ns)
    Pp <- cp / rp[cell_fi]; Pq <- cq / rq[cell_fi]; Pp[!is.finite(Pp)] <- 0; Pq[!is.finite(Pq)] <- 0
    sqrt(sum((Pq - Pp)^2))
  }
  obs_dist <- dist_from_postedge(post_ud[udid] == 1L)
  perm_vals <- vapply(seq_len(CFG$B_perm), function(b) dist_from_postedge(sample(post_ud)[udid] == 1L), numeric(1))
  perm_dist <- data.table(b = seq_len(CFG$B_perm), distance = perm_vals)
  perm_p <- mean(perm_dist$distance >= obs_dist, na.rm = TRUE)
} else { obs_dist <- NA_real_; perm_dist <- data.table(); perm_p <- NA_real_ }
network_test <- data.table(observed_distance = obs_dist, permutation_p = perm_p, B_perm = CFG$B_perm)
safe_fwrite(network_test, file.path(CFG$out_dir, "S12_network_distance_permutation_test.csv"))
safe_fwrite(perm_dist, file.path(CFG$out_dir, "S12_network_distance_permutation_distribution.csv"))

# ------------------------------------------------------------------------------
# 14. Edge stability
# ------------------------------------------------------------------------------

edge_counts_uid <- edges_facility[, .(N = .N), by = .(uid, period, from, to)]
top_edges <- edge_counts_uid[period == "post", .(N = sum(N)), by = .(from, to)][order(-N)][seq_len(min(.N, CFG$top_edges_for_stability))]
boot_edge_stability <- function(edge_uid, top_edges, B = 500) {
  ids <- unique(edge_uid$uid); if (length(ids) < 2 || nrow(top_edges) == 0) return(data.table())
  res <- rbindlist(lapply(seq_len(B), function(b) {
    smp <- sample(ids, length(ids), replace = TRUE)
    w <- data.table(uid = names(table(smp)), weight = as.numeric(table(smp)))
    e <- edge_uid[w, on = "uid", nomatch = 0]; e[, wN := N * weight]
    e_sum <- e[, .(boot_N = sum(wN)), by = .(period, from, to)]
    out <- top_edges[e_sum[period == "post"], on = .(from, to)]; out[is.na(boot_N), boot_N := 0]; out[, b := b]; out[, .(b, from, to, boot_N)]
  }))
  res[, .(mean_N = mean(boot_N), ci_low = quantile(boot_N,0.025), ci_high = quantile(boot_N,0.975),
          positive_rate = mean(boot_N > 0)), by = .(from, to)][order(-mean_N)]
}
edge_stability <- boot_edge_stability(edge_counts_uid, top_edges, CFG$B_boot)
safe_fwrite(edge_stability, file.path(CFG$out_dir, "S13_edge_stability_uid_bootstrap.csv"))

# ------------------------------------------------------------------------------
# 15. Sequence analysis
# ------------------------------------------------------------------------------

chain_strings <- chain_facility[, .(sequence = paste(zone, collapse = " -> "),
  sequence_label = paste(zone_label, collapse = " -> "), chain_len = .N,
  has_onikuru = any(is_onikuru), has_center = any(in_center)), by = .(uid, day, period)]
safe_fwrite(chain_strings, file.path(CFG$out_dir, "S14_sequence_strings.csv"))
if (isTRUE(CFG$run_sequence_analysis)) {
  seq_target <- chain_facility[period == "post" & uid %in% post_visitors]
  if (nrow(seq_target) > 0) {
    seq_wide <- dcast(seq_target, uid + day ~ ord, value.var = "zone", fill = "NONE")
    seq_cols <- setdiff(names(seq_wide), c("uid","day"))
    if (length(seq_cols) >= 2 && nrow(seq_wide) >= CFG$sequence_k) {
      if (nrow(seq_wide) <= CFG$sequence_max_n) {
        seq_obj <- TraMineR::seqdef(seq_wide[, ..seq_cols], right = "NONE")
        seq_dist <- TraMineR::seqdist(seq_obj, method = "OM", indel = 1, sm = "CONSTANT")
        hc <- hclust(as.dist(seq_dist), method = "ward.D2")
        seq_wide[, seq_cluster := cutree(hc, k = CFG$sequence_k)]
      } else {
        ## many sequences: collapse identical sequences to one row and cluster with the counts as weights
        ## (no person-day x person-day distance matrix; with many ties the partition need not match the method above)
        uq <- seq_wide[, .(n_seq = .N), by = seq_cols]
        msg("[15] %d sequences -> %d distinct sequences", nrow(seq_wide), nrow(uq))
        if (nrow(uq) <= CFG$sequence_max_n && nrow(uq) >= CFG$sequence_k) {
          seq_obj <- TraMineR::seqdef(uq[, ..seq_cols], right = "NONE")
          seq_dist <- TraMineR::seqdist(seq_obj, method = "OM", indel = 1, sm = "CONSTANT")
          hc <- hclust(as.dist(seq_dist), method = "ward.D2", members = uq$n_seq)
          uq[, seq_cluster := cutree(hc, k = CFG$sequence_k)]
          seq_wide[uq, on = seq_cols, seq_cluster := i.seq_cluster]
        } else {
          msg("[15] too many distinct sequences; clustering skipped (CFG$sequence_max_n = %d)", CFG$sequence_max_n)
          seq_wide[, seq_cluster := NA_integer_]
        }
      }
      safe_fwrite(seq_wide, file.path(CFG$out_dir, "S14_sequence_clusters.csv"))
      safe_fwrite(seq_wide[, .N, by = seq_cluster][, share := N / sum(N)], file.path(CFG$out_dir, "S14_sequence_cluster_summary.csv"))
    }
  }
}

# ------------------------------------------------------------------------------
# 16A. Within-person stage 2 (central-city-day baseline) + uid x dow fixed-effect regressions
# ------------------------------------------------------------------------------

pday_wp <- chain_facility[, .(chain_len = .N, dow = dow[1],
  has_onikuru = any(is_onikuru), has_center = any(in_center),
  n_center_nodes = uniqueN(zone[in_center])), by = .(uid, day, period, post)]

d1 <- pday_wp[period == "post" & uid %in% post_visitors]
d1[, daytype := fifelse(has_onikuru, "ONIKURU_day", "non_ONIKURU_day_all")]
within_stage1 <- d1[, .(person_days = .N, mean_n_center_nodes = mean(n_center_nodes),
  onward_rate = mean(n_center_nodes >= 2L), mean_chain_len = mean(chain_len)), by = daytype]
safe_fwrite(within_stage1, file.path(CFG$out_dir, "S17_within_stage1_descriptive.csv"))

d2 <- pday_wp[period == "post" & uid %in% post_visitors & has_center == TRUE]
d2[, daytype := fifelse(has_onikuru, "ONIKURU_day", "other_center_day")]
d2[, onward := as.integer(n_center_nodes >= 2L)]
within_stage2 <- d2[, .(person_days = .N, mean_n_center_nodes = mean(n_center_nodes),
  onward_rate = mean(onward), mean_chain_len = mean(chain_len)), by = daytype][order(-daytype)]
safe_fwrite(within_stage2, file.path(CFG$out_dir, "S17_within_stage2_descriptive.csv"))

within_stage2_reg <- data.table()
if (requireNamespace("fixest", quietly = TRUE) && uniqueN(d2$has_onikuru) == 2) {
  d2[, is_oni := as.integer(has_onikuru)]
  ext <- function(m, lbl) { ct <- as.data.table(summary(m)$coeftable, keep.rownames = "term")
    setnames(ct, c("term","est","se","t","p")); ct <- ct[term == "is_oni"]
    data.table(outcome = lbl, est = ct$est, se = ct$se, ci_low = ct$est-1.96*ct$se, ci_high = ct$est+1.96*ct$se, p = ct$p) }
  m_nodes  <- fixest::feols(n_center_nodes ~ is_oni | uid + dow, data = d2, cluster = ~uid)
  m_onward <- fixest::feols(onward ~ is_oni | uid + dow, data = d2, cluster = ~uid)
  m_nodes2 <- fixest::feols(n_center_nodes ~ is_oni + chain_len | uid + dow, data = d2, cluster = ~uid)
  within_stage2_reg <- rbindlist(list(ext(m_nodes,"n_center_nodes (total)"),
    ext(m_onward,"onward (>=2 centers)"), ext(m_nodes2,"n_center_nodes | chain_len controlled")), fill = TRUE)
  safe_fwrite(within_stage2_reg, file.path(CFG$out_dir, "S17_within_stage2_regression.csv"))
  msg("[16A] within-person stage-2 regressions (central-city-day baseline, uid + dow FE)"); print(within_stage2_reg)
}

# ------------------------------------------------------------------------------
# 16B. Multi-anchor benchmark (ONIKURU vs the station cores vs every central zone)
# ------------------------------------------------------------------------------

anchors <- data.table(anchor_id = c(CFG$onikuru_label, make_zone_id(CFG$pref_name, CFG$city_name, CFG$center_names)),
                      label = c("ONIKURU", CFG$center_names))
post_e <- tc_facility[period == "post", .(from, to, N)]
multi_anchor <- rbindlist(lapply(anchors$anchor_id, function(a)
  decompose_ego(post_e, a, zone_center_map_fac, "post", "facility")), fill = TRUE)
multi_anchor <- merge(multi_anchor, anchors, by.x = "target", by.y = "anchor_id", all.x = TRUE)
multi_anchor <- multi_anchor[(out_strength + in_strength) >= CFG$min_anchor_traffic]
multi_anchor[, is_ONIKURU := target == CFG$onikuru_label]
setorder(multi_anchor, -out_center_share)
safe_fwrite(multi_anchor[, .(label, is_ONIKURU, out_strength, in_strength,
  out_center_share = round(out_center_share,3), onward_to_center = round(onward_to_center,3),
  flow_through = round(flow_through,3))], file.path(CFG$out_dir, "S18_multi_anchor_decomposition.csv"))

cf_post <- chain_facility[period == "post"]
multi_chain <- rbindlist(lapply(anchors$anchor_id, function(a) {
  keys <- cf_post[, .(has_a = any(zone == a)), by = .(uid, day)][has_a == TRUE, .(uid, day)]
  if (!nrow(keys)) return(NULL)
  cc <- cf_post[keys, on = .(uid, day)][, .(ctype = classify_anchor_chain(zone, in_center, zone == a)), by = .(uid, day)]
  data.table(anchor_id = a, n_days = nrow(cc),
    connecting_share = mean(cc$ctype %in% c("post_bridge","middle_bridge")),
    terminal_share   = mean(cc$ctype %in% c("onikuru_only_center","final_after_center")))
}), fill = TRUE)
multi_chain <- merge(multi_chain, anchors, by = "anchor_id", all.x = TRUE)[n_days >= CFG$min_anchor_days]
multi_chain[, is_ONIKURU := anchor_id == CFG$onikuru_label]
setorder(multi_chain, -connecting_share)
safe_fwrite(multi_chain[, .(label, is_ONIKURU, n_days,
  connecting_share = round(connecting_share,3), terminal_share = round(terminal_share,3))],
  file.path(CFG$out_dir, "S18_multi_anchor_chain_share.csv"))

# ==============================================================================
# PART A2. Additional analyses
#   Keep the PART A definitions (the definitions of the article) and, with the same chains and chain types, add
#   (1) home stays excluded, (2) sampling uncertainty, (3) pre/post comparisons,
#   (4) town-block-level benchmark (16J), (5) H2 intervals (16K), (6) post-opening phases (16L).
# ==============================================================================

anchor_ids <- make_zone_id(CFG$pref_name, CFG$city_name, CFG$center_names)
ANCH       <- c(CFG$onikuru_label, anchor_ids)                  # ONIKURU + the 16 central town-blocks
core_ids   <- make_zone_id(CFG$pref_name, CFG$city_name, CFG$station_cores)
anchor_lab <- function(z) ifelse(z == CFG$onikuru_label, z, sub("^[^_]+_[^_]+_", "", z))
CONNECT    <- c("post_bridge", "middle_bridge")
fe_row <- function(m, term) {                                    # extract one coefficient from a fixest model
  ct <- as.data.frame(fixest::coeftable(m))
  if (!(term %in% rownames(ct))) return(list(est = NA_real_, se = NA_real_, p = NA_real_))
  list(est = ct[term, 1], se = ct[term, 2], p = ct[term, 4])
}

# ------------------------------------------------------------------------------
# 16C. Multi-anchor benchmark (as a function): article definition / home stays excluded + rank bootstrap
# ------------------------------------------------------------------------------

## chain type of each anchor-day; a vectorized version of the rule in classify_anchor_chain():
##   after  = another central zone is visited after the first stay at the anchor
##   before = another central zone is visited before the last stay at the anchor
anchor_days <- function(cf, anch = ANCH) {
  cz <- cf[in_center == TRUE, .(first = min(ord), last = max(ord)), by = .(uid, day, zone)]
  a  <- cz[zone %in% anch, .(uid, day, anchor = zone, a_first = first, a_last = last)]
  j  <- cz[a, on = .(uid, day), allow.cartesian = TRUE, nomatch = 0][zone != anchor,
          .(after = any(last > a_first), before = any(first < a_last)), by = .(uid, day, anchor)]
  out <- j[a[, .(uid, day, anchor)], on = .(uid, day, anchor)]
  out[is.na(after), after := FALSE]; out[is.na(before), before := FALSE]
  out[, ctype := fifelse(before & after, "middle_bridge", fifelse(after, "post_bridge",
                 fifelse(before, "final_after_center", "onikuru_only_center")))]
  out[, conn := ctype %in% CONNECT]
  out[]
}

anchor_benchmark <- function(d, variant) {
  cf <- build_chain(d, "zone_id", "zone_label")[period == "post"]
  ad <- anchor_days(cf)
  tc <- make_edges(cf)[, .(N = .N), by = .(from, to)]
  cmap <- cf[, .(zone_center = any(in_center, na.rm = TRUE)), by = .(zone)]
  dec <- rbindlist(lapply(ANCH, function(a) decompose_ego(tc, a, cmap, "post", "facility")), fill = TRUE)
  ch  <- ad[, .(n_days = .N, n_connecting = sum(conn), connecting_share = mean(conn)), by = .(target = anchor)]
  tab <- merge(dec[, .(target, out_strength, in_strength, out_center_share, onward_to_center, flow_through)], ch, by = "target")
  tab <- tab[(out_strength + in_strength) >= CFG$min_anchor_traffic & n_days >= CFG$min_anchor_days]
  tab[, `:=`(variant = variant, label = anchor_lab(target), is_ONIKURU = target == CFG$onikuru_label)]
  tab[, rank_connecting := frank(-connecting_share, ties.method = "min")]
  tab[, rank_outflow    := frank(-out_center_share, ties.method = "min")]
  setorder(tab, -connecting_share)
  list(variant = variant, table = tab, ad = ad, cf = cf)
}

## resample users with replacement and recompute the connecting-chain-share ranks
boot_rank <- function(res, B = CFG$B_rank, seed = 20261013) {
  keep <- res$table$target; io <- match(CFG$onikuru_label, keep)
  if (is.na(io) || length(keep) < 3) return(NULL)
  per <- res$ad[anchor %in% keep, .(days = .N, conn = sum(conn)), by = .(uid, anchor)]
  D <- dcast(per, uid ~ anchor, value.var = "days", fill = 0)
  C <- dcast(per, uid ~ anchor, value.var = "conn", fill = 0)
  Dm <- as.matrix(D[, keep, with = FALSE]); Cm <- as.matrix(C[, keep, with = FALSE]); n <- nrow(Dm)
  set.seed(seed)
  W <- stats::rmultinom(B, n, rep(1 / n, n))                    # n x B resampling counts
  S <- crossprod(W, Cm) / pmax(crossprod(W, Dm), 1e-9)          # B x anchors connecting-chain shares
  rk <- rowSums(S > S[, io]) + 1L
  ci <- match(intersect(core_ids, keep), keep)
  data.table(variant = res$variant, B = B, n_users = n,
    share_ci_low = unname(quantile(S[, io], 0.025)), share_ci_high = unname(quantile(S[, io], 0.975)),
    rank_ci_low = unname(quantile(rk, 0.025, type = 1)), rank_ci_high = unname(quantile(rk, 0.975, type = 1)),
    p_above_median = mean(S[, io] > apply(S, 1, median)),
    p_above_both_cores = if (length(ci)) mean(rowSums(S[, ci, drop = FALSE] < S[, io]) == length(ci)) else NA_real_)
}

bench_ms <- anchor_benchmark(dt, "manuscript definition")
bench_list <- list(bench_ms)
if (HAS_HOME) {
  bench_list[[2]] <- anchor_benchmark(dt[is_home == FALSE], "home stays excluded")
  bench_list[[3]] <- anchor_benchmark(dt[is_home == FALSE & is_work == FALSE], "home and work stays excluded")
}
bench_all  <- rbindlist(lapply(bench_list, `[[`, "table"), fill = TRUE)
bench_boot <- rbindlist(lapply(bench_list, boot_rank), fill = TRUE)
safe_fwrite(bench_all[, .(variant, label, is_ONIKURU, n_days, n_connecting, connecting_share,
  rank_connecting, out_strength, in_strength, out_center_share, rank_outflow, onward_to_center, flow_through)],
  file.path(CFG$out_dir, "S20_anchor_benchmark_variants.csv"))
safe_fwrite(bench_boot, file.path(CFG$out_dir, "S20_anchor_benchmark_bootstrap.csv"))

home_diag <- data.table(
  item = c("home_radius_m", "stays_total", "stays_home", "stays_home_share", "stays_work", "users_total", "users_with_home_stay",
           "central_stays_total", "central_stays_home", "home_mesh_rule", "users_with_multiple_home_meshes_in_residence_file"),
  value = c(CFG$home_radius_m, nrow(dt), dt[is_home == TRUE, .N], round(dt[, mean(is_home)], 4), dt[is_work == TRUE, .N],
            uniqueN(dt$uid), dt[is_home == TRUE, uniqueN(uid)], dt[in_center == TRUE, .N], dt[in_center == TRUE & is_home == TRUE, .N],
            CFG$home_mesh_rule, home_multi_users))
safe_fwrite(home_diag, file.path(CFG$out_dir, "S20_home_stay_diagnostics.csv"))

# ------------------------------------------------------------------------------
# 16D. Within-person contrast for every anchor (same equation as 16A) and the reference values of the proviso
#      sample = post-opening central-city days of users with at least one anchor-day
#      reach = 1(central zones that day >= 2); D = anchor visited that day
# ------------------------------------------------------------------------------

h4_anchorwise <- function(res) {
  cf <- res$cf; ad <- res$ad
  pday <- cf[, .(dow = dow[1], has_center = any(in_center), nc = uniqueN(zone[in_center])), by = .(uid, day)]
  pday[, reach := as.numeric(nc >= 2L)]
  rbindlist(lapply(res$table$target, function(a) {
    da <- ad[anchor == a, .(uid, day, D = 1)]
    x  <- merge(pday[uid %in% da$uid], da, by = c("uid", "day"), all.x = TRUE); x[is.na(D), D := 0]
    xc <- x[has_center == TRUE]; b <- xc[D == 0]
    fe <- list(est = NA_real_, se = NA_real_, p = NA_real_)
    if (uniqueN(xc$D) == 2) {
      m <- tryCatch(fixest::feols(reach ~ D | uid + dow, data = xc, cluster = ~uid), error = function(e) NULL)
      if (!is.null(m)) fe <- fe_row(m, "D")
    }
    data.table(variant = res$variant, target = a, label = anchor_lab(a), is_ONIKURU = a == CFG$onikuru_label,
      n_users = uniqueN(da$uid), n_anchor_days = nrow(da),
      reach_anchor_day = xc[D == 1, mean(reach)],
      reach_other_central_day = if (nrow(b)) mean(b$reach) else NA_real_,
      ref_size_biased = if (nrow(b)) b[, sum(nc * reach) / sum(nc)] else NA_real_,   # weighted by the number of central zones
      ref_any_central = x[D == 0, mean(has_center)],                                   # share of non-anchor days with a central-city stay
      fe_est = fe$est, fe_se = fe$se, fe_ci_low = fe$est - 1.96 * fe$se, fe_ci_high = fe$est + 1.96 * fe$se, fe_p = fe$p)
  }), fill = TRUE)
}
h4_anchor <- rbindlist(lapply(bench_list[seq_len(min(2L, length(bench_list)))], h4_anchorwise), fill = TRUE)
h4_anchor[, rank_fe := frank(-fe_est, ties.method = "min", na.last = "keep"), by = variant]
h4_anchor[, excess_over_size_biased := reach_anchor_day - ref_size_biased]
safe_fwrite(h4_anchor, file.path(CFG$out_dir, "S21_H4_anchorwise.csv"))

## table checking the reading rule (Section 3.3.5 of the article) numerically
reading_rule <- rbindlist(lapply(bench_list[seq_len(min(2L, length(bench_list)))], function(res) {
  t <- res$table; o <- t[is_ONIKURU == TRUE]; if (!nrow(o)) return(NULL)
  cores <- t[target %in% core_ids]; h <- h4_anchor[variant == res$variant & is_ONIKURU == TRUE]
  bt <- bench_boot[variant == res$variant]
  lo_ref <- min(h$ref_size_biased, h$ref_any_central)
  data.table(variant = res$variant, n_anchors = nrow(t), onikuru_visit_days = o$n_days,
    connecting_share = o$connecting_share, connecting_median = median(t$connecting_share),
    connecting_rank = o$rank_connecting, above_median_connecting = o$connecting_share > median(t$connecting_share),
    outflow_share = o$out_center_share, outflow_median = median(t$out_center_share),
    outflow_rank = o$rank_outflow, above_median_outflow = o$out_center_share > median(t$out_center_share),
    core_connecting_max = if (nrow(cores)) max(cores$connecting_share) else NA_real_,
    above_both_station_cores = nrow(cores) == length(core_ids) && all(o$connecting_share > cores$connecting_share),
    visit_day_reach = h$reach_anchor_day, ref_size_biased = h$ref_size_biased, ref_any_central = h$ref_any_central,
    proviso_holds = h$reach_anchor_day >= lo_ref,
    h4_fe = h$fe_est, h4_fe_ci_low = h$fe_ci_low, h4_fe_ci_high = h$fe_ci_high, h4_fe_rank = h$rank_fe,
    boot_share_ci_low = if (nrow(bt)) bt$share_ci_low else NA_real_, boot_share_ci_high = if (nrow(bt)) bt$share_ci_high else NA_real_,
    boot_rank_ci_low = if (nrow(bt)) bt$rank_ci_low else NA_real_, boot_rank_ci_high = if (nrow(bt)) bt$rank_ci_high else NA_real_,
    boot_p_above_median = if (nrow(bt)) bt$p_above_median else NA_real_,
    boot_p_above_both_cores = if (nrow(bt)) bt$p_above_both_cores else NA_real_)
}), fill = TRUE)
safe_fwrite(reading_rule, file.path(CFG$out_dir, "S19_reading_rule.csv"))
msg("[16C/D] reading rule"); print(t(reading_rule))

## onward destinations from ONIKURU (transitions per block and by group)
grp_of <- function(z) { nm <- anchor_lab(z); g <- rep("other", length(nm))
  for (k in names(CFG$dest_groups)) g[nm %in% CFG$dest_groups[[k]]] <- k; g }
onward_dest <- tc_facility[period == "post" & from == CFG$onikuru_label & to %in% anchor_ids,
  .(to, label = anchor_lab(to), group = grp_of(to), N)][order(-N)]
safe_fwrite(onward_dest, file.path(CFG$out_dir, "S6c_onikuru_onward_destinations.csv"))

if (isTRUE(CFG$run_additional)) {

# ------------------------------------------------------------------------------
# 16E. User table (home-to-ONIKURU distance, home block) and H3 with the distance term
# ------------------------------------------------------------------------------

users <- data.table(uid = unique(dt$uid)); users[, visitor := as.integer(uid %in% post_visitors)]
h3_distance <- data.table()
if (HAS_HOME) {
  oni_xy <- dt[is_onikuru == TRUE & is.finite(.lon), .(lon = median(.lon), lat = median(.lat))]
  users <- merge(users, dt[, .(home_lon = home_mesh_lon[1], home_lat = home_mesh_lat[1]), by = uid], by = "uid")
  users[, dist_m := hav_m(home_lon, home_lat, oni_xy$lon, oni_xy$lat)]
  ## home block = modal town-block of the 7 nearest stay points (all users) within home_block_radius_m
  ##              of the residence mesh centroid
  r <- CFG$home_block_radius_m; dlat <- r / 110940; dlon <- r / (111320 * cos(oni_xy$lat * pi / 180))
  pts <- dt[is_onikuru == FALSE & is.finite(.lon) & is.finite(.lat), .(plon = .lon, plat = .lat, pzone = base_zone_id)]
  hq  <- unique(users[is.finite(home_lon), .(home_lon, home_lat)])
  hq[, `:=`(lon_lo = home_lon - dlon, lon_hi = home_lon + dlon, lat_lo = home_lat - dlat, lat_hi = home_lat + dlat)]
  hz <- rbindlist(lapply(split(seq_len(nrow(hq)), ceiling(seq_len(nrow(hq)) / 500)), function(ix) {   # 500 meshes at a time
    cand <- pts[hq[ix], on = .(plon >= lon_lo, plon <= lon_hi, plat >= lat_lo, plat <= lat_hi), nomatch = 0, allow.cartesian = TRUE,
                .(home_lon, home_lat, plon = x.plon, plat = x.plat, pzone)]
    cand[, d := hav_m(plon, plat, home_lon, home_lat)]
    cand[d <= r][order(d), .(home_zone = mode_first(head(pzone, 7L))), by = .(home_lon, home_lat)]
  }))
  users <- merge(users, hz, by = c("home_lon", "home_lat"), all.x = TRUE)
  users[is.na(home_zone), home_zone := "UNKNOWN"]
  users[, home_in_center := home_zone %in% anchor_ids]
  rm(pts, hq, hz)

  f <- merge(h3_frame, users[, .(uid, dist_m)], by = "uid")[is.finite(dist_m) & dist_m > 0]
  f[, log_km := log(dist_m / 1000)]
  if (nrow(f) > 0 && uniqueN(f$became_visitor) == 2) {
    m1 <- glm(became_visitor ~ pre_center_day_rate + pre_observed_days + post_observed_days, data = f, family = binomial())
    m2 <- glm(became_visitor ~ pre_center_day_rate + log_km + pre_observed_days + post_observed_days, data = f, family = binomial())
    h3_distance <- rbind(tidy_logit(m1, "adjusted (same sample)"), tidy_logit(m2, "adjusted + log distance"))
    nd <- data.table(pre_center_day_rate = c(0, 0.5), log_km = log(median(f$dist_m) / 1000),
                     pre_observed_days = mean(f$pre_observed_days), post_observed_days = mean(f$post_observed_days))
    h3_pred <- data.table(pre_center_day_rate = nd$pre_center_day_rate,
      p_adjusted = predict(m1, nd, type = "response"), p_with_distance = predict(m2, nd, type = "response"))
    h3_comp <- data.table(
      item = c("users_in_frame", "visitors_in_frame", "visitors_with_zero_pre_central_use", "share_zero_among_visitors",
               "share_zero_among_all", "median_distance_km_visitors", "median_distance_km_non_visitors"),
      value = c(nrow(f), sum(f$became_visitor), f[became_visitor == 1, sum(pre_center_day_rate == 0)],
                f[became_visitor == 1, mean(pre_center_day_rate == 0)], f[, mean(pre_center_day_rate == 0)],
                f[became_visitor == 1, median(dist_m) / 1000], f[became_visitor == 0, median(dist_m) / 1000]))
    safe_fwrite(h3_distance, file.path(CFG$out_dir, "S22_H3_logistic_with_distance.csv"))
    safe_fwrite(h3_pred, file.path(CFG$out_dir, "S22_H3_predicted.csv"))
    safe_fwrite(h3_comp, file.path(CFG$out_dir, "S22_H3_visitor_composition.csv"))
  }
}
safe_fwrite(users, file.path(CFG$out_dir, "S22_users.csv"))

# ------------------------------------------------------------------------------
# 16F. Pre/post difference-in-differences (visitors vs non-visitors)
#      outcome oc = a non-home stay in a central zone other than ONIKURU on that day (0/1)
#      oc ~ visitor x post | uid + week relative to opening + day of week; user-clustered SE
#      sample = users observed in both periods
# ------------------------------------------------------------------------------

d_out <- if (HAS_HOME) dt[is_home == FALSE] else dt
cf_out <- build_chain(d_out, "zone_id", "zone_label")
panel <- cf_out[, .(dow = dow[1], has_oni = any(is_onikuru),
                    n_oc = uniqueN(zone[in_center & !is_onikuru])), by = .(uid, day, period)]
panel[, `:=`(oc = as.numeric(n_oc > 0), post = as.numeric(period == "post"), visitor = as.numeric(uid %in% post_visitors),
             rel_week = as.integer(floor(as.numeric(as.Date(day) - CFG$open_date) / 7)))]
if (HAS_HOME) panel <- merge(panel, users[, .(uid, dist_m)], by = "uid", all.x = TRUE)
both <- panel[, .(np = uniqueN(post)), by = uid][np == 2, uid]
x <- panel[uid %in% both]; x[, D := visitor * post]

did_tab <- data.table(); es <- data.table(); did_means <- data.table()
if (uniqueN(x$D) == 2) {
  m_did <- fixest::feols(oc ~ D | uid + rel_week + dow, data = x, cluster = ~uid); r1 <- fe_row(m_did, "D")
  did_tab <- data.table(model = "visitor x post", term = "D", est = r1$est, se = r1$se,
    ci_low = r1$est - 1.96 * r1$se, ci_high = r1$est + 1.96 * r1$se, p = r1$p,
    n_person_days = nrow(x), n_visitors = x[visitor == 1, uniqueN(uid)], n_non_visitors = x[visitor == 0, uniqueN(uid)])
  if (HAS_HOME) {
    x[, near_post := as.numeric(!is.na(dist_m) & dist_m <= CFG$near_km * 1000) * post]
    if (uniqueN(x$near_post) == 2) {
      m_nr <- fixest::feols(oc ~ near_post | uid + rel_week + dow, data = x, cluster = ~uid); r2 <- fe_row(m_nr, "near_post")
      did_tab <- rbind(did_tab, data.table(model = sprintf("home within %.1f km x post", CFG$near_km), term = "near_post",
        est = r2$est, se = r2$se, ci_low = r2$est - 1.96 * r2$se, ci_high = r2$est + 1.96 * r2$se, p = r2$p,
        n_person_days = nrow(x), n_visitors = x[!is.na(dist_m) & dist_m <= CFG$near_km * 1000, uniqueN(uid)],
        n_non_visitors = x[is.na(dist_m) | dist_m > CFG$near_km * 1000, uniqueN(uid)]))
    }
  }
  ## event study: visitor x week; reference week = the last pre-opening week with at least min_days_ref_week observed days
  wk <- x[, .(n_days = uniqueN(day)), by = rel_week][order(rel_week)]
  cand_w <- wk[rel_week < 0 & n_days >= CFG$min_days_ref_week, rel_week]
  ref_week <- if (length(cand_w)) max(cand_w) else wk[rel_week < 0, max(rel_week)]
  m_es <- fixest::feols(oc ~ i(rel_week, visitor, ref = ref_week) | uid + rel_week + dow, data = x, cluster = ~uid)
  ct <- as.data.table(as.data.frame(fixest::coeftable(m_es)), keep.rownames = "term")
  setnames(ct, c("term", "est", "se", "t", "p"))
  ct <- ct[grepl("^rel_week::-?[0-9]+:visitor$", term)]
  ct[, rel_week := as.integer(sub("^rel_week::(-?[0-9]+):visitor$", "\\1", term))]
  es <- merge(ct[, .(rel_week, est, se, ci_low = est - 1.96 * se, ci_high = est + 1.96 * se, p)],
              x[, .(n_days = uniqueN(day), visitor_person_days = sum(visitor)), by = rel_week], by = "rel_week", all = TRUE)
  es[rel_week == ref_week, `:=`(est = 0, se = 0, ci_low = 0, ci_high = 0)]
  pre_terms <- ct[rel_week < 0, term]; pre_p <- NA_real_
  if (length(pre_terms)) {
    b <- stats::coef(m_es)[pre_terms]; V <- stats::vcov(m_es)[pre_terms, pre_terms, drop = FALSE]
    w_stat <- tryCatch(as.numeric(t(b) %*% solve(V) %*% b), error = function(e) NA_real_)
    pre_p <- stats::pchisq(w_stat, df = length(b), lower.tail = FALSE)
  }
  did_tab[, `:=`(pretrend_p = c(pre_p, rep(NA_real_, .N - 1L)), ref_week = ref_week)]
  vp <- x[visitor == 1 & post == 1]
  did_means <- rbind(
    x[, .(mean_oc = mean(oc), person_days = .N, users = uniqueN(uid)), by = .(group = fifelse(visitor == 1, "visitors", "non-visitors"),
                                                                               period = fifelse(post == 1, "post", "pre"))],
    data.table(group = "visitors", period = c("post: share of days with an ONIKURU visit",
                                               "post: share of days with ONIKURU and no other central stay",
                                               "post: share of days with any central stay incl. ONIKURU"),
               mean_oc = c(mean(vp$has_oni), mean(vp$has_oni & vp$oc == 0), mean(vp$oc > 0 | vp$has_oni)),
               person_days = nrow(vp), users = uniqueN(vp$uid)))
  safe_fwrite(did_tab, file.path(CFG$out_dir, "S23_DID_estimates.csv"))
  safe_fwrite(es[order(rel_week)], file.path(CFG$out_dir, "S23_DID_event_study.csv"))
  safe_fwrite(did_means, file.path(CFG$out_dir, "S23_DID_group_means.csv"))
  msg("[16F] DID"); print(did_tab)
}

# ------------------------------------------------------------------------------
# 16G. Size and distance reference for onward destinations (Poisson gravity model)
#      fit N ~ origin + destination + log(distance) to the transitions among anchors other than ONIKURU,
#      and compare the expected and observed destinations of ONIKURU (own block excluded)
# ------------------------------------------------------------------------------

dest_tab <- data.table(); dest_elasticity <- NA_real_
if (HAS_XY) {
  loc <- dt[zone_id %in% ANCH & is.finite(.lon), .(lon = median(.lon), lat = median(.lat)), by = .(zone = zone_id)]
  az  <- intersect(ANCH, loc$zone)
  e_post <- make_edges(cf_out[period == "post"])
  tcg <- e_post[from %in% az & to %in% az, .(N = .N), by = .(from, to)]
  g <- CJ(from = az, to = az)[from != to]
  g <- merge(g, tcg, by = c("from", "to"), all.x = TRUE); g[is.na(N), N := 0L]
  g[loc, on = .(from = zone), `:=`(lon1 = i.lon, lat1 = i.lat)]; g[loc, on = .(to = zone), `:=`(lon2 = i.lon, lat2 = i.lat)]
  g[, ld := log(hav_m(lon1, lat1, lon2, lat2))]
  fit <- g[from != CFG$onikuru_label & to != CFG$onikuru_label]
  own <- make_zone_id(CFG$pref_name, CFG$city_name, CFG$onikuru_area_name)
  og  <- g[from == CFG$onikuru_label & to != own]
  if (nrow(fit) > 0 && sum(fit$N) > 0 && sum(og$N) > 0) {
    gm <- suppressWarnings(glm(N ~ factor(from) + factor(to) + ld, data = fit, family = poisson()))
    cf_g <- stats::coef(gm); dest_elasticity <- unname(cf_g["ld"])
    de <- setNames(rep(0, length(az)), az); nm <- names(cf_g)[startsWith(names(cf_g), "factor(to)")]
    de[sub("^factor\\(to\\)", "", nm)] <- cf_g[nm]; de[is.na(de)] <- 0
    og[, w := exp(de[to] + dest_elasticity * ld)]; og[, p_exp := w / sum(w)]; og[, group := grp_of(to)]
    tg <- og[, .(observed = sum(N), p_exp = sum(p_exp)), by = group][order(group)]
    n_on <- sum(tg$observed); set.seed(20261014)
    sim <- stats::rmultinom(CFG$B_dest, n_on, tg$p_exp)
    dest_tab <- tg[, .(group, observed, expected = n_on * p_exp, ratio = observed / (n_on * p_exp),
      p_at_least_observed = rowMeans(sim >= observed), p_at_most_observed = rowMeans(sim <= observed))]
    dest_tab[, `:=`(n_onward = n_on, distance_elasticity = dest_elasticity)]
    safe_fwrite(dest_tab, file.path(CFG$out_dir, "S24_destination_gravity_reference.csv"))
    safe_fwrite(og[, .(to, label = anchor_lab(to), group, observed = N, expected = sum(N) * p_exp, distance_m = exp(ld))][order(-observed)],
                file.path(CFG$out_dir, "S24_destination_by_block.csv"))
  }
}

# ------------------------------------------------------------------------------
# 16H. Pre/post change by town-block (area level, home stays excluded) and returning home on visit-days
# ------------------------------------------------------------------------------

ca_out <- build_chain(d_out, "zone_area_id", "zone_area_label")
pd_n <- ca_out[, .(person_days = uniqueN(paste(uid, day))), by = period]
zs <- unique(ca_out[in_center == TRUE, zone])
vd <- ca_out[zone %in% zs, .(visit_days = uniqueN(paste(uid, day))), by = .(zone, period)]
block_change <- dcast(vd, zone ~ period, value.var = "visit_days", fill = 0)
for (pp in c("pre", "post")) if (!(pp %in% names(block_change))) block_change[, (pp) := 0]
n_pd <- function(pp) { v <- pd_n[period == pp, person_days]; if (length(v)) v else NA_real_ }
block_change[, `:=`(pre_per_1000 = 1000 * pre / n_pd("pre"), post_per_1000 = 1000 * post / n_pd("post"))]
block_change[, change := fifelse(pre_per_1000 > 0, post_per_1000 / pre_per_1000 - 1, NA_real_)]
block_change[, rank_change := frank(-change, ties.method = "min", na.last = "keep")]
block_change[, label := area_zone_label(zone)]
setorder(block_change, -change)
safe_fwrite(block_change[, .(zone, label, pre_visit_days = pre, post_visit_days = post, pre_per_1000, post_per_1000, change, rank_change)],
            file.path(CFG$out_dir, "S25_block_prepost_change.csv"))

## visit-days with a central-city stay after ONIKURU, and those on which all such stays are home stays
od <- dt[period == "post"][order(uid, day, arrival, departure)]
od[, k := seq_len(.N), by = .(uid, day)]
od[, first_oni := if (any(is_onikuru)) min(k[is_onikuru]) else NA_integer_, by = .(uid, day)]
aft <- od[!is.na(first_oni) & k > first_oni & in_center == TRUE & is_onikuru == FALSE, .(all_home = all(is_home)), by = .(uid, day)]
home_return <- data.table(item = c("onikuru_visit_days", "days_with_later_central_stay", "of_which_all_later_central_stays_at_home",
                                   "visitors", "visitors_living_in_central_blocks"),
  value = c(nrow(unique(od[!is.na(first_oni), .(uid, day)])), nrow(aft), sum(aft$all_home), length(post_visitors),
            if ("home_in_center" %in% names(users)) users[visitor == 1, sum(home_in_center)] else NA))
safe_fwrite(home_return, file.path(CFG$out_dir, "S25_home_return_on_visit_days.csv"))

# ------------------------------------------------------------------------------
# 16J. Town-block-level benchmark (area level): before and after opening
#      Instead of comparing the building with town-blocks, place the ONIKURU-area (facility + Ekimae-3)
#      among the other central town-blocks as one unit and compare it with the same block before opening (= Ekimae-3).
# ------------------------------------------------------------------------------

AREA_UNITS <- c(CFG$onikuru_area_label, center_zones_area)       # ONIKURU_AREA + the other central town-blocks (including JR_AREA)

area_benchmark <- function(ca, variant) {
  rbindlist(lapply(c("pre", "post"), function(pp) {
    cf <- ca[period == pp]; if (!nrow(cf)) return(NULL)
    ad <- anchor_days(cf, AREA_UNITS)
    tc <- make_edges(cf)[, .(N = .N), by = .(from, to)]
    cmap <- cf[, .(zone_center = any(in_center, na.rm = TRUE)), by = .(zone)]
    dec <- rbindlist(lapply(AREA_UNITS, function(a) decompose_ego(tc, a, cmap, pp, "area")), fill = TRUE)
    ch  <- ad[, .(n_days = .N, n_connecting = sum(conn), connecting_share = mean(conn)), by = .(target = anchor)]
    tab <- merge(dec[, .(target, out_strength, in_strength, out_center_share, onward_to_center, flow_through)], ch, by = "target")
    tab <- tab[(out_strength + in_strength) >= CFG$min_anchor_traffic & n_days >= CFG$min_anchor_days]
    tab[, `:=`(variant = variant, period = pp, label = area_zone_label(target), is_ONIKURU_AREA = target == CFG$onikuru_area_label)]
    tab[, rank_connecting := frank(-connecting_share, ties.method = "min")]
    tab[, rank_outflow    := frank(-out_center_share, ties.method = "min")]
    tab[, n_units := .N]
    tab[]
  }), fill = TRUE)
}

## resample users with replacement and recompute the pre-to-post change and rank of the ONIKURU-area
area_change_boot <- function(ca, tab, variant, B = CFG$B_rank, seed = 20261015) {
  keep <- intersect(tab[period == "pre", target], tab[period == "post", target])
  io <- match(CFG$onikuru_area_label, keep); if (is.na(io) || length(keep) < 3) return(NULL)
  ad <- rbindlist(lapply(c("pre", "post"), function(pp) anchor_days(ca[period == pp], keep)[, period := pp]))
  per <- ad[, .(days = .N, conn = sum(conn)), by = .(uid, anchor, period)]
  users_b <- unique(per$uid)
  mk <- function(pp, v) { m <- matrix(0, length(users_b), length(keep)); x <- per[period == pp]
    m[cbind(match(x$uid, users_b), match(x$anchor, keep))] <- x[[v]]; m }
  D0 <- mk("pre", "days"); C0 <- mk("pre", "conn"); D1 <- mk("post", "days"); C1 <- mk("post", "conn"); n <- nrow(D0)
  set.seed(seed); W <- stats::rmultinom(B, n, rep(1 / n, n))
  S0 <- crossprod(W, C0) / pmax(crossprod(W, D0), 1e-9); S1 <- crossprod(W, C1) / pmax(crossprod(W, D1), 1e-9)
  dA <- S1[, io] - S0[, io]; dO <- apply((S1 - S0)[, -io, drop = FALSE], 1, median)
  r0 <- rowSums(S0 > S0[, io]) + 1L; r1 <- rowSums(S1 > S1[, io]) + 1L
  o0 <- tab[period == "pre" & target == CFG$onikuru_area_label]; o1 <- tab[period == "post" & target == CFG$onikuru_area_label]
  oth <- merge(tab[period == "pre" & target %in% keep[-io], .(target, s0 = connecting_share)],
               tab[period == "post" & target %in% keep[-io], .(target, s1 = connecting_share)], by = "target")
  q <- function(x, p) unname(quantile(x, p))
  data.table(variant = variant, n_units = length(keep),
    pre_share = o0$connecting_share, post_share = o1$connecting_share, change = o1$connecting_share - o0$connecting_share,
    change_ci_low = q(dA, .025), change_ci_high = q(dA, .975), p_change_positive = mean(dA > 0),
    median_change_other_blocks = median(oth$s1 - oth$s0),
    change_minus_others = (o1$connecting_share - o0$connecting_share) - median(oth$s1 - oth$s0),
    change_minus_others_ci_low = q(dA - dO, .025), change_minus_others_ci_high = q(dA - dO, .975),
    pre_rank = o0$rank_connecting, post_rank = o1$rank_connecting,
    pre_rank_ci = sprintf("%d-%d", quantile(r0, .025, type = 1), quantile(r0, .975, type = 1)),
    post_rank_ci = sprintf("%d-%d", quantile(r1, .025, type = 1), quantile(r1, .975, type = 1)),
    p_rank_improved = mean(r1 < r0),
    pre_outflow_share = o0$out_center_share, post_outflow_share = o1$out_center_share,
    pre_outflow_rank = o0$rank_outflow, post_outflow_rank = o1$rank_outflow,
    pre_days = o0$n_days, post_days = o1$n_days)
}

area_variants <- list(list(ca = chain_area, v = "manuscript definition"))
if (HAS_HOME) area_variants[[2]] <- list(ca = ca_out, v = "home stays excluded")
area_bench  <- rbindlist(lapply(area_variants, function(a) area_benchmark(a$ca, a$v)), fill = TRUE)
area_change <- rbindlist(lapply(area_variants, function(a) area_change_boot(a$ca, area_bench[variant == a$v], a$v)), fill = TRUE)
setorder(area_bench, variant, -period, -connecting_share)
safe_fwrite(area_bench[, .(variant, period, label, is_ONIKURU_AREA, n_days, n_connecting, connecting_share, rank_connecting,
  out_strength, in_strength, out_center_share, rank_outflow, onward_to_center, flow_through, n_units)],
  file.path(CFG$out_dir, "S26_area_benchmark_prepost.csv"))
safe_fwrite(area_change, file.path(CFG$out_dir, "S26_area_benchmark_change.csv"))
msg("[16J] area-level benchmark, pre -> post"); print(t(area_change))

# ------------------------------------------------------------------------------
# 16K. Interval estimates for H2: user-level bootstrap
#      Fisher tests treat transitions as independent; transitions of the same user are not, so
#      user-cluster intervals are attached to the pre/post changes in the inflow and outflow shares.
# ------------------------------------------------------------------------------

h2_boot <- data.table()
if (nrow(edges_area) > 0 && uniqueN(edges_area$period) == 2) {
  A <- CFG$onikuru_area_label
  cnt <- edges_area[, .(in_out = sum(from %in% center_zones_area), in_hit = sum(from %in% center_zones_area & to == A),
                        out_out = sum(from == A), out_hit = sum(from == A & to %in% center_zones_area)), by = .(uid, period)]
  w <- dcast(cnt, uid ~ period, value.var = c("in_out", "in_hit", "out_out", "out_hit"), fill = 0)
  M <- as.matrix(w[, -1]); n <- nrow(M)
  set.seed(20261016); W <- stats::rmultinom(CFG$B_rank, n, rep(1 / n, n)); Bm <- crossprod(W, M)   # B x 8
  one <- function(lbl, hit, out) {
    obs <- function(pp) sum(M[, paste0(hit, "_", pp)]) / sum(M[, paste0(out, "_", pp)])
    bs  <- function(pp) Bm[, paste0(hit, "_", pp)] / pmax(Bm[, paste0(out, "_", pp)], 1e-9)
    d <- bs("post") - bs("pre")
    data.table(direction = lbl, pre_share = obs("pre"), post_share = obs("post"), change = obs("post") - obs("pre"),
      change_ci_low = unname(quantile(d, .025)), change_ci_high = unname(quantile(d, .975)),
      p_boot_two_sided = min(1, 2 * min(mean(d <= 0), mean(d >= 0))),
      pre_edges = sum(M[, paste0(out, "_pre")]), post_edges = sum(M[, paste0(out, "_post")]), n_users = n, B = CFG$B_rank)
  }
  h2_boot <- rbind(one("inflow: central blocks -> ONIKURU_AREA", "in_hit", "in_out"),
                   one("outflow: ONIKURU_AREA -> central blocks", "out_hit", "out_out"))
  safe_fwrite(h2_boot, file.path(CFG$out_dir, "S27_H2_user_bootstrap.csv"))
  msg("[16K] H2 user-cluster bootstrap"); print(h2_boot)
}

# ------------------------------------------------------------------------------
# 16L. Key indicators by post-opening phase (is the pattern an opening-period surge?)
#      split at the New Year exclusion: phase 1 = opening to the day before the exclusion, phase 2 = the day after it onwards
# ------------------------------------------------------------------------------

phase_tab <- rbindlist(lapply(list(
  list(lbl = sprintf("phase 1: %s to %s", CFG$follow_from, CFG$ny_exclude_from - 1), from = CFG$follow_from, to = CFG$ny_exclude_from - 1),
  list(lbl = sprintf("phase 2: %s to %s", CFG$ny_exclude_to + 1, CFG$follow_to), from = CFG$ny_exclude_to + 1, to = CFG$follow_to)),
  function(ph) {
    d <- dt[period == "post" & day >= ph$from & day <= ph$to]
    if (!nrow(d) || !any(d$is_onikuru)) return(NULL)
    res <- anchor_benchmark(d, ph$lbl); t_ <- res$table; o <- t_[is_ONIKURU == TRUE]; if (!nrow(o)) return(NULL)
    bt <- boot_rank(res, seed = 20261017); h <- h4_anchorwise(res)
    h[, rank_fe := frank(-fe_est, ties.method = "min", na.last = "keep")]; ho <- h[is_ONIKURU == TRUE]
    ct <- res$ad[anchor == CFG$onikuru_label, .N, by = ctype]; cores <- t_[target %in% core_ids]
    data.table(phase = ph$lbl, observed_days = uniqueN(d$day), person_days = uniqueN(d[, .(uid, day)]),
      onikuru_visit_days = o$n_days, visitors = res$ad[anchor == CFG$onikuru_label, uniqueN(uid)],
      visit_days_per_1000_person_days = 1000 * o$n_days / uniqueN(d[, .(uid, day)]),
      onikuru_only_share = sum(ct[ctype == "onikuru_only_center", N]) / o$n_days,
      connecting_share = o$connecting_share, share_ci_low = if (!is.null(bt)) bt$share_ci_low else NA_real_,
      share_ci_high = if (!is.null(bt)) bt$share_ci_high else NA_real_,
      connecting_median = median(t_$connecting_share), connecting_rank = o$rank_connecting, n_anchors = nrow(t_),
      p_above_median = if (!is.null(bt)) bt$p_above_median else NA_real_,
      above_both_station_cores = nrow(cores) == length(core_ids) && all(o$connecting_share > cores$connecting_share),
      outflow_share = o$out_center_share, outflow_median = median(t_$out_center_share), outflow_rank = o$rank_outflow,
      flow_through = o$flow_through,
      visit_day_reach = ho$reach_anchor_day, ref_size_biased = ho$ref_size_biased, ref_any_central = ho$ref_any_central,
      h4_fe = ho$fe_est, h4_fe_ci_low = ho$fe_ci_low, h4_fe_ci_high = ho$fe_ci_high, h4_fe_rank = ho$rank_fe)
  }), fill = TRUE)
safe_fwrite(phase_tab, file.path(CFG$out_dir, "S28_phase_split.csv"))
msg("[16L] phases of the post-opening period"); print(t(phase_tab))

# ------------------------------------------------------------------------------
# 16M. Anchor-wise "new uptake" gradient (reference for H3)
#      Every ONIKURU visitor is a new user of the facility, so for each town-block the outcome is
#      whether a user with no pre-opening stay in that block stayed there after opening, and
#      the H3 equation (pre-opening central-city use rate, observed days, + distance from home) is fitted.
#      Any central-city destination shows a positive gradient, so ONIKURU is read within that distribution.
# ------------------------------------------------------------------------------

h3_anchor <- data.table()
if (nrow(h3_frame) > 0) {
  pre_stay  <- unique(chain_facility[period == "pre"  & zone %in% ANCH, .(uid, zone)])
  post_stay <- unique(chain_facility[period == "post" & zone %in% ANCH, .(uid, zone)])
  a_loc <- if (HAS_XY) dt[zone_id %in% ANCH & is.finite(.lon), .(lon = median(.lon), lat = median(.lat)), by = .(zone = zone_id)] else NULL
  h3_anchor <- rbindlist(lapply(ANCH, function(a) {
    f <- h3_frame[!(uid %in% pre_stay[zone == a, uid])]                 # users with no pre-opening stay in a
    f[, y := as.integer(uid %in% post_stay[zone == a, uid])]
    if (sum(f$y) < 20 || uniqueN(f$y) < 2) return(NULL)
    m1 <- glm(y ~ pre_center_day_rate + pre_observed_days + post_observed_days, data = f, family = binomial())
    co <- summary(m1)$coefficients["pre_center_day_rate", ]
    nd <- data.table(pre_center_day_rate = c(0, 0.5), pre_observed_days = mean(f$pre_observed_days), post_observed_days = mean(f$post_observed_days))
    p1 <- predict(m1, nd, type = "response")
    out <- data.table(target = a, label = anchor_lab(a), is_ONIKURU = a == CFG$onikuru_label,
      n_users = nrow(f), n_new_users = sum(f$y), uptake_rate = mean(f$y),
      OR = exp(0.1 * co[1]), OR_low = exp(0.1 * (co[1] - 1.96 * co[2])), OR_high = exp(0.1 * (co[1] + 1.96 * co[2])),
      p_at_zero = p1[1], p_at_half = p1[2], ratio_zero_to_half = p1[1] / p1[2],
      share_zero_among_new_users = f[y == 1, mean(pre_center_day_rate == 0)])
    if (HAS_HOME && !is.null(a_loc) && a %in% a_loc$zone && "home_lon" %in% names(users)) {
      g <- merge(f, users[, .(uid, home_lon, home_lat)], by = "uid")
      g[, log_km := log(pmax(hav_m(home_lon, home_lat, a_loc[zone == a, lon], a_loc[zone == a, lat]), 1) / 1000)]
      g <- g[is.finite(log_km)]
      if (uniqueN(g$y) == 2) {
        m2 <- glm(y ~ pre_center_day_rate + log_km + pre_observed_days + post_observed_days, data = g, family = binomial())
        c2 <- summary(m2)$coefficients
        nd2 <- cbind(nd, log_km = log(median(exp(g$log_km))))
        p2 <- predict(m2, nd2, type = "response")
        out[, `:=`(OR_dist = exp(0.1 * c2["pre_center_day_rate", 1]),
          OR_dist_low = exp(0.1 * (c2["pre_center_day_rate", 1] - 1.96 * c2["pre_center_day_rate", 2])),
          OR_dist_high = exp(0.1 * (c2["pre_center_day_rate", 1] + 1.96 * c2["pre_center_day_rate", 2])),
          OR_distance_doubling = exp(log(2) * c2["log_km", 1]),
          ratio_zero_to_half_dist = p2[1] / p2[2])]
      }
    }
    out
  }), fill = TRUE)
  if (nrow(h3_anchor)) {
    h3_anchor[, rank_OR_flattest := frank(OR, ties.method = "min")]                     # 1 = flattest gradient (broadest uptake)
    if ("OR_dist" %in% names(h3_anchor)) h3_anchor[, rank_OR_dist_flattest := frank(OR_dist, ties.method = "min", na.last = "keep")]
    h3_anchor[, rank_ratio_highest := frank(-ratio_zero_to_half, ties.method = "min")]  # 1 = highest relative uptake among zero-use users
    setorder(h3_anchor, OR)
    safe_fwrite(h3_anchor, file.path(CFG$out_dir, "S29_H3_anchorwise_new_uptake.csv"))
    msg("[16M] new-user uptake gradient by anchor"); print(h3_anchor[, .(label, n_new_users, OR = round(OR, 3), ratio = round(ratio_zero_to_half, 3), rank_OR_flattest)])
  }
}

# ------------------------------------------------------------------------------
# 16I. Summary of the additional analyses
# ------------------------------------------------------------------------------

kv <- function(section, metric, value) data.table(section = section, metric = metric, value = as.character(value))
add_summary <- rbindlist(c(
  list(kv("spec", "spec", CFG$spec), kv("spec", "stay_file", basename(stay_path)),
       kv("data", diag_s2$item, diag_s2$value), kv("home", home_diag$item, home_diag$value),
       kv("home", home_return$item, home_return$value)),
  lapply(seq_len(nrow(reading_rule)), function(i) {
    rr <- reading_rule[i]; cols <- setdiff(names(rr), "variant")
    kv(paste0("reading rule: ", rr$variant), cols, vapply(cols, function(cc) format(rr[[cc]], digits = 4), character(1))) }),
  if (nrow(h3_logit)) list(kv("H3", paste(h3_logit$model, h3_logit$term, h3_logit$unit),
       sprintf("%.3f [%.3f, %.3f]", h3_logit$OR, h3_logit$OR_low, h3_logit$OR_high))),
  if (nrow(h3_distance)) list(kv("H3 + distance", paste(h3_distance$model, h3_distance$term, h3_distance$unit),
       sprintf("%.3f [%.3f, %.3f]", h3_distance$OR, h3_distance$OR_low, h3_distance$OR_high))),
  if (nrow(h3_anchor) && any(h3_anchor$is_ONIKURU)) list(kv("H3 by anchor (new users)",
       c("ONIKURU OR per +0.10 [95% CI]", "rank, flattest gradient first (anchors)", "median OR of the other anchors",
         "ONIKURU p(zero)/p(0.5)", "median p(zero)/p(0.5) of the other anchors"),
       c(sprintf("%.3f [%.3f, %.3f]", h3_anchor[is_ONIKURU == TRUE, OR], h3_anchor[is_ONIKURU == TRUE, OR_low], h3_anchor[is_ONIKURU == TRUE, OR_high]),
         sprintf("%d (%d)", as.integer(h3_anchor[is_ONIKURU == TRUE, rank_OR_flattest]), nrow(h3_anchor)),
         sprintf("%.3f", h3_anchor[is_ONIKURU == FALSE, median(OR)]),
         sprintf("%.3f", h3_anchor[is_ONIKURU == TRUE, ratio_zero_to_half]),
         sprintf("%.3f", h3_anchor[is_ONIKURU == FALSE, median(ratio_zero_to_half)])))),
  if (nrow(did_tab)) list(kv("DID", did_tab$model, sprintf("%.4f [%.4f, %.4f]", did_tab$est, did_tab$ci_low, did_tab$ci_high)),
       kv("DID", c("pretrend_p", "event_study_ref_week"), c(format(did_tab$pretrend_p[1], digits = 3), did_tab$ref_week[1]))),
  if (nrow(dest_tab)) list(kv("destinations", paste0(dest_tab$group, " observed / expected"),
       sprintf("%d / %.1f", as.integer(dest_tab$observed), dest_tab$expected))),
  if (nrow(area_change)) list(kv(paste0("area benchmark: ", area_change$variant), "connecting share pre -> post [change, 95% CI]",
       sprintf("%.3f -> %.3f [%+.3f, %+.3f to %+.3f]", area_change$pre_share, area_change$post_share, area_change$change,
               area_change$change_ci_low, area_change$change_ci_high)),
       kv(paste0("area benchmark: ", area_change$variant), "rank pre -> post (units)",
          sprintf("%d -> %d (%d)", as.integer(area_change$pre_rank), as.integer(area_change$post_rank), as.integer(area_change$n_units))),
       kv(paste0("area benchmark: ", area_change$variant), "change minus median change of other blocks [95% CI]",
          sprintf("%+.3f [%+.3f, %+.3f]", area_change$change_minus_others, area_change$change_minus_others_ci_low, area_change$change_minus_others_ci_high))),
  if (nrow(h2_boot)) list(kv("H2 user bootstrap", h2_boot$direction,
       sprintf("%.4f -> %.4f, change %+.4f [%+.4f, %+.4f], p = %.3f", h2_boot$pre_share, h2_boot$post_share, h2_boot$change,
               h2_boot$change_ci_low, h2_boot$change_ci_high, h2_boot$p_boot_two_sided))),
  if (nrow(phase_tab)) list(kv(phase_tab$phase, "visit-days / connecting share (rank) / H4 FE",
       sprintf("%d / %.3f (%d of %d) / %.3f", as.integer(phase_tab$onikuru_visit_days), phase_tab$connecting_share,
               as.integer(phase_tab$connecting_rank), as.integer(phase_tab$n_anchors), phase_tab$h4_fe))),
  list(kv("block change", c("ONIKURU_AREA change in visit-days per 1000 person-days", "rank among central blocks", "blocks compared"),
       c(format(block_change[zone == CFG$onikuru_area_label, change], digits = 4),
         block_change[zone == CFG$onikuru_area_label, rank_change], nrow(block_change))))), fill = TRUE)
safe_fwrite(add_summary, file.path(CFG$out_dir, "S30_additional_summary.csv"))

}  # run_additional

# ------------------------------------------------------------------------------
# 17. Final summary + save
# ------------------------------------------------------------------------------

final_summary <- rbindlist(list(
  data.table(section = "data", metric = diag_s2$item, value = as.character(diag_s2$value)),
  data.table(section = "zone", metric = zone_diag$item, value = as.character(zone_diag$value)),
  data.table(section = "H2",
    metric = c("ALL_CENTER_inflow_share_pre","ALL_CENTER_inflow_share_post","ALL_CENTER_inflow_fisher_p",
               "ALL_CENTER_outflow_share_pre","ALL_CENTER_outflow_share_post","ALL_CENTER_outflow_fisher_p"),
    value = as.character(c(round(h2_agg[direction=="inflow",  pre_share],4), round(h2_agg[direction=="inflow",  post_share],4), round(h2_agg[direction=="inflow",  fisher_p],4),
                           round(h2_agg[direction=="outflow", pre_share],4), round(h2_agg[direction=="outflow", post_share],4), round(h2_agg[direction=="outflow", fisher_p],4)))),
  data.table(section = "H4", metric = h4_boot_ci$metric,
             value = paste0(round(h4_boot_ci$estimate,4), " [", round(h4_boot_ci$ci_low,4), ", ", round(h4_boot_ci$ci_high,4), "]")),
  data.table(section = "H5", metric = h5_indicators$indicator, value = as.character(h5_indicators$value)),
  data.table(section = "network", metric = c("observed_distance","permutation_p"), value = as.character(c(obs_dist, perm_p))),
  data.table(section = "within2", metric = within_stage2$daytype, value = as.character(round(within_stage2$onward_rate,4))),
  data.table(section = "H3", metric = if (nrow(h3_logit)) paste(h3_logit$model, h3_logit$term) else character(),
             value = if (nrow(h3_logit)) sprintf("OR %.3f [%.3f, %.3f] %s", h3_logit$OR, h3_logit$OR_low, h3_logit$OR_high, h3_logit$unit) else character()),
  data.table(section = "reading_rule", metric = paste(rep(reading_rule$variant, each = 4),
               c("connecting_rank", "above_median_both", "above_both_station_cores", "proviso_holds")),
             value = as.character(c(rbind(paste0(reading_rule$connecting_rank, "/", reading_rule$n_anchors),
               reading_rule$above_median_connecting & reading_rule$above_median_outflow,
               reading_rule$above_both_station_cores, reading_rule$proviso_holds))))), fill = TRUE)
safe_fwrite(final_summary, file.path(CFG$out_dir, "S99_final_summary.csv"))

saveRDS(list(CFG=CFG, dt=dt, uid_attr=uid_attr, chain_facility=chain_facility, chain_area=chain_area,
  chain_day=chain_day, edges_facility=edges_facility, edges_area=edges_area, tc_facility=tc_facility, tc_area=tc_area,
  h2_all=h2_all, h2_center=h2_center, h2_agg=h2_agg, h2_jr=h2_jr,
  h3_summary=h3_summary, h4_onikuru_days=h4_onikuru_days, chain_type_summary=chain_type_summary,
  h5_indicators=h5_indicators, ego_decomp=ego_decomp, ego_edges=ego_edges, network_test=network_test,
  flow_tbl=flow_tbl, zone_coords=zone_coords, within_stage1=within_stage1, within_stage2=within_stage2,
  within_stage2_reg=within_stage2_reg, multi_anchor=multi_anchor, multi_chain=multi_chain,
  obs_days=obs_days, obs_gaps=obs_gaps, h3_logit=h3_logit, bench_all=bench_all, bench_boot=bench_boot,
  h4_anchor=h4_anchor, reading_rule=reading_rule,
  additional=mget(intersect(c("users","h3_distance","did_tab","es","did_means","dest_tab","block_change","home_return","area_bench","area_change","h2_boot","phase_tab","h3_anchor","add_summary"), ls()))),
  file = file.path(CFG$out_dir, "pipeline_results.rds"))
msg("PIPELINE DONE. Outputs saved in: %s", CFG$out_dir)

# ==============================================================================
# PART B. TABLES (CSV) & FIGURES (PDF + PNG)
#   tables : Table01..Table11 + TableS1; additional analyses TableA1..A11
#   figures: main F3..F12 / supplement FS1, FS2 / additional analyses FA1..FA3
# ==============================================================================

OUT <- normalizePath(CFG$out_dir)
FIG_M <- file.path(OUT, "figs");   dir.create(FIG_M, showWarnings = FALSE, recursive = TRUE)   # article figures (B-4)
TAB_M <- file.path(OUT, "tables"); dir.create(TAB_M, showWarnings = FALSE, recursive = TRUE)   # article tables (B-4)
FIG <- file.path(OUT, "diagnostics", "figs");   dir.create(FIG, showWarnings = FALSE, recursive = TRUE)   # internal (B-1 to B-3)
TAB <- file.path(OUT, "diagnostics", "tables"); dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
WRITE_DIAG <- isTRUE(CFG$write_diagnostics)
has_pw <- requireNamespace("patchwork", quietly = TRUE)
rd <- function(f) { p <- file.path(OUT, f); if (file.exists(p)) fread(p) else NULL }
wt <- function(x, f) { if (!WRITE_DIAG && TAB != TAB_M) return(invisible(NULL)); fwrite(x, file.path(TAB, f), bom = TRUE); cat("table:", f, "\n") }
pdev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
png_type <- if (capabilities("cairo")) "cairo" else NULL
open_png <- function(file, w, h) { if (!is.null(png_type)) grDevices::png(file, width=w, height=h, units="in", res=200, type=png_type)
  else grDevices::png(file, width=w, height=h, units="in", res=200) }
save_pdf_gg <- function(p, name, w=11, h=4.6) {
  if (!WRITE_DIAG && FIG != FIG_M) return(invisible(NULL))
  ggplot2::ggsave(file.path(FIG, paste0(name,".pdf")), p, width=w, height=h, device=pdev)
  ggplot2::ggsave(file.path(FIG, paste0(name,".png")), p, width=w, height=h, dpi=200)
  cat("figure:", paste0(name," (pdf+png)"), "\n") }
save_base <- function(drawfun, name, w=9, h=6.5) {
  if (!WRITE_DIAG && FIG != FIG_M) return(invisible(NULL))
  pdev(file.path(FIG, paste0(name,".pdf")), width=w, height=h); drawfun(); dev.off()
  open_png(file.path(FIG, paste0(name,".png")), w, h); drawfun(); dev.off()
  cat("figure:", paste0(name," (pdf+png)"), "\n") }
ggplot2::theme_set(ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), plot.title = ggplot2::element_text(face = "bold"), legend.position = "top"))
C_CENTER <- "#2c7fb8"; C_NON <- "#d95f0e"; C_A <- "#41ab5d"; C_B <- "#807dba"; C_GREY <- "#888888"

ROMA <- c("ONIKURU"="ONIKURU","ONIKURU_AREA"="ONIKURU-area","JR_AREA"="Nishi-Ekimae","ALL_CENTER"="Central-city zones",   # JR_AREA = Nishi-Ekimae
  "駅前一丁目"="Ekimae-1","駅前二丁目"="Ekimae-2","駅前三丁目"="Ekimae-3","駅前四丁目"="Ekimae-4",
  "双葉町"="Futaba","松ケ本町"="Matsugamoto","永代町"="Eitai","元町"="Motomachi","新庄町"="Shinjo",
  "東中条町"="Higashi-Chujo","大手町"="Otemachi","中津町"="Nakatsu","別院町"="Betsuin","本町"="Honmachi",
  "舟木町"="Funaki","西中条町"="Nishi-Chujo","新中条町"="Shin-Chujo","春日一丁目"="Kasuga-1","城の前町"="Shironomae",
  "中穂積一丁目"="Nakahozumi-1","稲葉町"="Inaba","大池一丁目"="Oike-1","上郡二丁目"="Kamigori-2",
  "千里万博公園"="Banpaku-Park","西駅前町"="Nishi-Ekimae","岩倉町"="Iwakura","片桐町"="Katagiri","宮元町"="Miyamoto")
roma <- function(x) { base <- sub("^[^_]+_[^_]+_", "", x); key <- ifelse(x %in% names(ROMA), x, base)
  ifelse(key %in% names(ROMA), ROMA[key], key) }

# ------------------------------------------------------------------------------
# B-1. Tables
# ------------------------------------------------------------------------------

s2 <- rd("S2_data_diagnostics.csv"); cs <- rd("S4_chain_summary_by_period.csv")
zd <- rd("S3_zone_diagnostics.csv")
h2t <- rd("S6_H2_center_to_ONIKURU_AREA_prepost.csv")
ca <- rd("S6_H1_centrality_ONIKURU_AREA.csv"); cfac <- rd("S6_H1_centrality_ONIKURU_facility_post.csv"); dec <- rd("S6b_onikuru_outflow_decomposition.csv")
ct <- rd("S9_H5_chain_type_summary.csv"); ind <- rd("S10_integrated_catalyst_enclave_indicators.csv"); bc <- rd("S11_H4_cluster_bootstrap_CI.csv")
w2 <- rd("S17_within_stage2_descriptive.csv"); w2r <- rd("S17_within_stage2_regression.csv")
ma <- rd("S18_multi_anchor_decomposition.csv"); mc <- rd("S18_multi_anchor_chain_share.csv")

# Table01: sample descriptives
if (!is.null(s2) && !is.null(cs)) {
  gv <- function(it) s2[item == it, value][1]
  wt(data.table(item = c("Calendar days in window (New Year excluded)","Observed days","Stay points (rows)","Person-days","Unique uids","Mean chain length","Mean unique nodes","Center stay person-days","ONIKURU stay person-days"),
    pre = c(gv("pre_calendar_days"), gv("pre_observed_days"), gv("pre_rows"), cs[period=="pre",person_days], cs[period=="pre",uids], round(cs[period=="pre",mean_chain_len],3),
      round(cs[period=="pre",mean_unique_nodes],3), cs[period=="pre",center_days], cs[period=="pre",onikuru_days]),
    post = c(gv("post_calendar_days"), gv("post_observed_days"), gv("post_rows"), cs[period=="post",person_days], cs[period=="post",uids], round(cs[period=="post",mean_chain_len],3),
      round(cs[period=="post",mean_unique_nodes],3), cs[period=="post",center_days], cs[period=="post",onikuru_days])),
    "Table01_sample_descriptives.csv")
}

# Table02: zone diagnostics
if (!is.null(zd)) wt(zd, "Table02_zone_diagnostics.csv")

# Table03: H2 center <-> ONIKURU (all blocks x both directions + ALL_CENTER aggregate)
if (!is.null(h2t)) {
  h2fmt <- h2t[, .(zone_label, direction,
    pre_share = round(pre_share,3), post_share = round(post_share,3), share_change = round(share_change,3),
    pre_out_edges, post_out_edges, post_hit_edges,
    odds_ratio = round(odds_ratio,3), fisher_p = round(fisher_p,4), fisher_p_BH = round(fisher_p_BH,4))]
  setorder(h2fmt, direction, -post_share)
  wt(h2fmt, "Table03_H2_center_to_ONIKURU.csv")
}

# Table04: centrality + outflow decomposition
if (!is.null(ca)) {
  cen <- rbindlist(list(ca, cfac), fill=TRUE)[, .(level = ifelse(target=="ONIKURU","facility","area"), period, target,
    n_nodes, n_edges, betweenness = round(betweenness,1), out_degree, in_degree, out_strength, in_strength, pagerank = round(pagerank,5))]
  if (!is.null(dec)) cen <- dec[, .(level, period, target, out_center_share = round(out_center_share,3),
    in_center_share = round(in_center_share,3), onward_to_center = round(onward_to_center,3), flow_through = round(flow_through,3))][cen, on = .(level, period, target)]
  wt(cen, "Table04_centrality_decomposition.csv")
}

# Table05: chain types (terminal vs connecting)
if (!is.null(ct)) wt(ct[, .(chain_type, N, share = round(share,3), kind = fifelse(chain_type %in% c("post_bridge","middle_bridge"),"connecting","terminal"))], "Table05_chain_types.csv")

# Table06: ordered secondary stay (cluster bootstrap CI)
if (!is.null(bc)) wt(bc[, .(metric, estimate=round(estimate,4), ci_low=round(ci_low,4), ci_high=round(ci_high,4))], "Table06_secondary_bootstrapCI.csv")

# Table07 / Table08: within-person stage 2 (descriptive / regression)
if (!is.null(w2))  wt(w2,  "Table07_within_stage2.csv")
if (!is.null(w2r)) wt(w2r, "Table08_within_stage2_regression.csv")

# Table09 / Table10: multi-anchor benchmark (decomposition / chain share)
if (!is.null(ma)) wt(ma, "Table09_multi_anchor_decomposition.csv")
if (!is.null(mc)) wt(mc, "Table10_multi_anchor_chain_share.csv")

# Table11: integrated catalyst/enclave indicators (summary table)
if (!is.null(ind)) wt(ind[, .(indicator, value = round(value,4))], "Table11_integrated_indicators.csv")

# TableS1: sample flow (supplement)
if (!is.null(rd("S_sample_flow.csv"))) wt(rd("S_sample_flow.csv"), "TableS1_sample_flow.csv")

# ------------------------------------------------------------------------------
# B-2. Figures
# ------------------------------------------------------------------------------

# ---- F3 / F4 : ego network (romanized labels) ----
draw_ego <- function(el, title, summary_txt, topn = 12) {
  ind2 <- el[direction=="in"][order(-N)][seq_len(min(topn,.N))]; outd <- el[direction=="out"][order(-N)][seq_len(min(topn,.N))]
  col <- function(t) ifelse(t=="center", C_CENTER, C_NON)
  yin <- if (nrow(ind2)) seq(1,-1,length.out=nrow(ind2)) else numeric(0); yout <- if (nrow(outd)) seq(1,-1,length.out=nrow(outd)) else numeric(0)
  plot(NA, xlim=c(-2.3,2.3), ylim=c(-1.9,1.45), axes=FALSE, xlab="", ylab="", main=title)
  if (nrow(ind2)) for (i in seq_len(nrow(ind2))) segments(-1+.06,yin[i],-.05,0,col=col(ind2$node_type[i]),lwd=.4+ind2$N[i]*.4)
  if (nrow(outd)) for (i in seq_len(nrow(outd))) segments(.05,0,1-.06,yout[i],col=col(outd$node_type[i]),lwd=.4+outd$N[i]*.4)
  if (nrow(ind2)) { points(rep(-1,nrow(ind2)),yin,pch=19,cex=.8+ind2$N/8,col=col(ind2$node_type)); text(-1-.08,yin,paste0(roma(ind2$node)," (",ind2$N,")"),pos=2,cex=.65) }
  if (nrow(outd)) { points(rep(1,nrow(outd)),yout,pch=19,cex=.8+outd$N/8,col=col(outd$node_type)); text(1+.08,yout,paste0(roma(outd$node)," (",outd$N,")"),pos=4,cex=.65) }
  points(0,0,pch=19,cex=4.2,col="black"); text(0,0,roma(el$target[1]),col="white",cex=.65,font=2)
  text(-1,1.3,"IN (origins)",font=2); text(1,1.3,"OUT (destinations)",font=2); text(0,-1.55,summary_txt,cex=.7)
  legend("bottom", legend=c("center","non-center",roma(el$target[1])), col=c(C_CENTER,C_NON,"black"), pch=19, horiz=TRUE, bty="n", cex=.8)
}
eg <- rd("S6b_onikuru_ego_edgelist.csv")
sumtxt <- function(lv, pp) { if (is.null(dec)) return(""); r <- dec[level==lv & period==pp]; if (!nrow(r)) return("")
  sprintf("out-str %g  in-str %g | out_center_share %.2f  onward %.2f  flow_through %.2f", r$out_strength, r$in_strength, r$out_center_share, r$onward_to_center, r$flow_through) }
if (!is.null(eg)) {
  save_base(function() draw_ego(eg[level=="facility" & period=="post"], "Fig. 3  ONIKURU ego transition network (facility, post)", sumtxt("facility","post")), "F3_ego_network_facility_post", w=9, h=6.5)
  save_base(function() { par(mfrow=c(1,2))
    draw_ego(eg[level=="area" & period=="pre"], sprintf("Fig. 4(a) PRE (%s..%s)", CFG$baseline_from, format(CFG$baseline_to, "%m-%d")), sumtxt("area","pre"))
    draw_ego(eg[level=="area" & period=="post"], sprintf("Fig. 4(b) POST (%s..%s)", CFG$follow_from, format(CFG$follow_to, "%m-%d")), sumtxt("area","post")); par(mfrow=c(1,1)) }, "F4_area_prepost_network", w=14, h=7)
}

# ---- F5 : structural verdict (outflow decomposition + chain types) ----
if (!is.null(dec) && !is.null(ct)) {
  dec[, grp := paste0(level,"\n(",period,")")]
  da <- melt(dec, id.vars="grp", measure.vars=c("out_center_share","onward_to_center"), variable.name="metric", value.name="v")
  p5a <- ggplot2::ggplot(da, ggplot2::aes(grp, v, fill=metric)) + ggplot2::geom_col(position="dodge", width=.7) +
    ggplot2::geom_text(ggplot2::aes(label=sprintf("%.2f",v)), position=ggplot2::position_dodge(.7), vjust=-.3, size=3) +
    ggplot2::scale_fill_manual(values=c(out_center_share=C_CENTER, onward_to_center=C_A)) +
    ggplot2::labs(title="Fig. 5(a) ONIKURU outflow decomposition", x=NULL, y="share", fill=NULL) + ggplot2::ylim(0,.6)
  ct[, lab := factor(chain_type, levels=c("onikuru_only_center","final_after_center","post_bridge","middle_bridge"))]
  ct[, kind := fifelse(chain_type %in% c("post_bridge","middle_bridge"),"connecting","terminal")]
  p5b <- ggplot2::ggplot(ct, ggplot2::aes(lab, share, fill=kind)) + ggplot2::geom_col(width=.7) +
    ggplot2::geom_text(ggplot2::aes(label=sprintf("%.1f%%",100*share)), vjust=-.3, size=3) +
    ggplot2::scale_fill_manual(values=c(terminal=C_NON, connecting=C_CENTER)) +
    ggplot2::labs(title="Fig. 5(b) Chain types (terminal vs connecting)", x=NULL, y="share", fill=NULL) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=20, hjust=1)) + ggplot2::ylim(0,.85)
  save_pdf_gg(if (has_pw) patchwork::wrap_plots(p5a,p5b,nrow=1) else p5a, "F5_structural_verdict")
}

# ---- F6 : daily chain types ----
if (!is.null(ct)) {
  ctf <- copy(ct)
  lvl <- c("onikuru_only_center","final_after_center","post_bridge","middle_bridge")
  labmap <- c(onikuru_only_center = "ONIKURU only (terminal)",
              final_after_center  = "after other center (terminal)",
              post_bridge         = "ONIKURU -> center (connecting)",
              middle_bridge       = "center -> ONIKURU -> center (connecting)")
  ctf <- ctf[chain_type %in% lvl]
  ctf[, kind := fifelse(chain_type %in% c("post_bridge","middle_bridge"),"connecting","terminal")]
  ctf[, lab := factor(labmap[chain_type], levels = rev(labmap[lvl]))]
  term_sh <- ctf[kind=="terminal", sum(share)]; conn_sh <- ctf[kind=="connecting", sum(share)]
  p6 <- ggplot2::ggplot(ctf, ggplot2::aes(lab, share, fill = kind)) +
    ggplot2::geom_col(width = .68) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%% (n=%d)", 100*share, N)), hjust = -.08, size = 3.2) +
    ggplot2::scale_fill_manual(values = c(terminal = C_NON, connecting = C_CENTER), name = NULL) +
    ggplot2::coord_flip() + ggplot2::ylim(0, max(ctf$share)*1.30) +
    ggplot2::labs(title = "Fig. 6  Daily trip-chain types on ONIKURU visit-days (post)",
                  subtitle = sprintf("terminal %.1f%%   vs   connecting %.1f%%", 100*term_sh, 100*conn_sh),
                  x = NULL, y = "share of ONIKURU visit-days")
  save_pdf_gg(p6, "F6_chain_types", w = 9, h = 4.6)
}

# ---- F7 : within-person visitor day comparison ----
vc <- rd("S8_H4_visitor_day_comparison.csv")
if (!is.null(vc)) {
  m <- melt(vc, id.vars="has_onikuru", measure.vars=c("mean_chain_len","mean_unique_center_nodes","center_day_rate"), variable.name="metric", value.name="v")
  m[, metric := factor(metric, levels=c("mean_chain_len","mean_unique_center_nodes","center_day_rate"), labels=c("Chain length","Unique center nodes","Center-day rate"))]
  m[, day := fifelse(has_onikuru,"ONIKURU day","non-ONIKURU day")]
  p7a <- ggplot2::ggplot(m, ggplot2::aes(metric, v, fill=day)) + ggplot2::geom_col(position="dodge", width=.7) +
    ggplot2::geom_text(ggplot2::aes(label=sprintf("%.2f",v)), position=ggplot2::position_dodge(.7), vjust=-.3, size=3) +
    ggplot2::scale_fill_manual(values=c("non-ONIKURU day"=C_GREY,"ONIKURU day"=C_B)) +
    ggplot2::labs(title="Fig. 7(a) Same visitors: ONIKURU vs non-ONIKURU days", x=NULL, y="mean per person-day", fill=NULL)
  if (!is.null(bc)) {
    bc[, metric := factor(metric, levels=c("same_day_secondary_rate","secondary_2h_rate","secondary_4h_rate"), labels=c("same-day","<=2h","<=4h"))]
    p7b <- ggplot2::ggplot(bc, ggplot2::aes(metric, estimate)) + ggplot2::geom_pointrange(ggplot2::aes(ymin=ci_low, ymax=ci_high), color=C_A, size=.8) +
      ggplot2::geom_text(ggplot2::aes(label=sprintf("%.3f",estimate)), hjust=-.3, size=3) +
      ggplot2::labs(title="Fig. 7(b) Ordered secondary stay (cluster bootstrap 95% CI)", x=NULL, y="secondary-stay rate") + ggplot2::ylim(0,.25)
    save_pdf_gg(if (has_pw) patchwork::wrap_plots(p7a,p7b,nrow=1) else p7a, "F7_within_person")
  } else save_pdf_gg(p7a, "F7_within_person")
}

# ---- F8 : within-person stage2 (onward connection) ----
ws1 <- rd("S17_within_stage1_descriptive.csv"); ws2 <- rd("S17_within_stage2_descriptive.csv")
if (!is.null(ws1) && !is.null(ws2)) {
  comp <- rbindlist(list(
    ws1[, .(stage="vs all non-ONIKURU days", grp = fifelse(daytype == "ONIKURU_day","ONIKURU day","benchmark"), onward_rate)],   # exact match on the day type
    ws2[, .(stage="vs non-ONIKURU central days", grp = fifelse(daytype == "ONIKURU_day","ONIKURU day","benchmark"), onward_rate)]))
  p8 <- ggplot2::ggplot(comp, ggplot2::aes(grp, onward_rate, fill=grp)) + ggplot2::geom_col(width=.6) +
    ggplot2::geom_text(ggplot2::aes(label=sprintf("%.2f",onward_rate)), vjust=-.3, size=3) +
    ggplot2::facet_wrap(~stage) + ggplot2::scale_fill_manual(values=c("ONIKURU day"=C_B,"benchmark"=C_GREY), guide="none") +
    ggplot2::labs(title="Fig. 8  Within-person onward-connection (>=2 center nodes/day)", x=NULL, y="onward rate")
  save_pdf_gg(p8, "F8_within_stage2", w=9, h=4.6)
}

# ---- F9 : generative (H3 pre-propensity) ----
g3 <- rd("S7_H3_generative_summary.csv")
if (!is.null(g3)) {
  g3[, grp := factor(pre_center_group, levels=c("low_pre_center","middle_pre_center","high_pre_center"), labels=c("low (=0 pre)","middle","high (>=0.25)"))]
  p9a <- ggplot2::ggplot(g3, ggplot2::aes(grp, onikuru_day_rate)) + ggplot2::geom_col(fill=C_B, width=.6) +
    ggplot2::geom_text(ggplot2::aes(label=sprintf("%.4f\n(n=%d)", onikuru_day_rate, visitor_uids)), vjust=-.2, size=3) +
    ggplot2::labs(title="Fig. 9(a) ONIKURU uptake by pre-opening center propensity", x=NULL, y="ONIKURU visit-day rate (post)")
  p9b <- ggplot2::ggplot(g3, ggplot2::aes(grp, mean_unique_center_nodes)) + ggplot2::geom_col(fill=C_CENTER, width=.6) +
    ggplot2::geom_text(ggplot2::aes(label=sprintf("%.3f", mean_unique_center_nodes)), vjust=-.3, size=3) +
    ggplot2::labs(title="Fig. 9(b) Central-city reach by pre propensity", x=NULL, y="mean unique center nodes (post)")
  save_pdf_gg(if (has_pw) patchwork::wrap_plots(p9a,p9b,nrow=1) else p9a, "F9_generative")
}

# ---- F10 : multi-anchor benchmark ----
mad <- rd("S18_multi_anchor_decomposition.csv"); mac <- rd("S18_multi_anchor_chain_share.csv")
if (!is.null(mad)) {
  mad[, lab := factor(roma(label), levels = roma(label)[order(out_center_share)])]
  pA <- ggplot2::ggplot(mad, ggplot2::aes(lab, out_center_share, fill = is_ONIKURU)) + ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values=c("FALSE"=C_GREY,"TRUE"=C_NON), guide="none") +
    ggplot2::labs(title="Fig. 10(a) Outflow share to center, by anchor", x=NULL, y="out_center_share")
  if (!is.null(mac)) {
    mac[, lab := factor(roma(label), levels = roma(label)[order(connecting_share)])]
    pB <- ggplot2::ggplot(mac, ggplot2::aes(lab, connecting_share, fill = is_ONIKURU)) + ggplot2::geom_col() + ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values=c("FALSE"=C_GREY,"TRUE"=C_NON), guide="none") +
      ggplot2::labs(title="Fig. 10(b) Connecting-chain share, by anchor", x=NULL, y="connecting_share")
    save_pdf_gg(if (has_pw) patchwork::wrap_plots(pA,pB,nrow=1) else pA, "F10_multi_anchor", w=12, h=6)
  } else save_pdf_gg(pA, "F10_multi_anchor", w=7, h=6)
}

# ---- F11 : H2 flows by town-block (pre vs post, dumbbells) ----
h2dat <- rd("S6_H2_center_to_ONIKURU_AREA_prepost.csv")
if (!is.null(h2dat)) {
  # inflow: blocks with enough outgoing transitions; outflow: blocks with enough incoming transitions
  d_in  <- h2dat[direction=="inflow"  & zone!="ALL_CENTER" & (pre_out_edges + post_out_edges) >= CFG$h2_min_inflow_edges]
  d_out <- h2dat[direction=="outflow" & zone!="ALL_CENTER" & (pre_hit_edges + post_hit_edges) >= CFG$h2_min_outflow_edges]
  mk2 <- function(d, ttl, xlab) {
    if (!nrow(d)) return(NULL)
    d <- copy(d); d[, sig := !is.na(fisher_p_BH) & fisher_p_BH < 0.05]
    d[, zr := paste0(roma(zone_label), ifelse(sig, " *", ""))]; setorder(d, post_share)
    d[, zr := factor(zr, levels=unique(zr))]
    dl <- melt(d[, .(zr, pre_share, post_share)], id.vars="zr", variable.name="period", value.name="share")
    ggplot2::ggplot() +
      ggplot2::geom_segment(data=d, ggplot2::aes(y=zr, yend=zr, x=pre_share, xend=post_share), color=C_GREY, linewidth=.8) +
      ggplot2::geom_point(data=dl, ggplot2::aes(share, zr, color=period), size=2.7) +
      ggplot2::scale_color_manual(values=c(pre_share=C_GREY, post_share=C_NON), labels=c("pre","post"), name=NULL) +
      ggplot2::labs(title=ttl, subtitle="* = BH-adjusted Fisher p < .05", x=xlab, y=NULL)
  }
  pIn  <- mk2(d_in,  "Fig. 11(a) Town-block -> ONIKURU (share of block outflow)", "share to ONIKURU")
  pOut <- mk2(d_out, "Fig. 11(b) ONIKURU -> town-block (share of ONIKURU outflow)", "share from ONIKURU")
  ps <- Filter(Negate(is.null), list(pIn, pOut))
  if (length(ps) == 2 && has_pw) save_pdf_gg(patchwork::wrap_plots(ps[[1]], ps[[2]], nrow=1), "F11_H2_center_flows", w=13, h=6.5)
  else if (length(ps) >= 1) save_pdf_gg(ps[[1]], "F11_H2_center_flows", w=8, h=6.5)
}

# ---- F12 : GIS geographic transition network ----
zc <- rd("S_zone_coords.csv"); tcf <- rd("S5_transition_counts_facility.csv")
if (!is.null(zc) && !is.null(tcf)) {
  oni <- zc[zone_id == CFG$onikuru_label]
  if (nrow(oni)) {
    dlat <- CFG$gis_window_km/111; dlon <- CFG$gis_window_km/(111*cos(oni$lat*pi/180))
    win <- zc[abs(lat-oni$lat) <= dlat & abs(lon-oni$lon) <= dlon]; kz <- win$zone_id
    e <- tcf[period=="post" & from %in% kz & to %in% kz & N >= 2]
    e <- merge(e, zc[, .(zone_id, lon, lat)], by.x="from", by.y="zone_id"); setnames(e, c("lon","lat"), c("x1","y1"))
    e <- merge(e, zc[, .(zone_id, lon, lat)], by.x="to", by.y="zone_id");   setnames(e, c("lon","lat"), c("x2","y2"))
    e[, oedge := from==CFG$onikuru_label | to==CFG$onikuru_label]
    draw_gis <- function() {
      asp <- 1/cos(oni$lat*pi/180)
      plot(win$lon, win$lat, type="n", asp=asp, xlab="lon", ylab="lat", main="Fig. 12  Geographic transition network (central area, post)")
      if (nrow(e)) segments(e$x1,e$y1,e$x2,e$y2, lwd=0.3+e$N*0.25, col=adjustcolor(ifelse(e$oedge,"#444444","#bbbbbb"), alpha.f=0.6))
      ncol <- ifelse(win$is_onikuru,"black", ifelse(win$is_center, C_CENTER, C_NON))
      ncex <- ifelse(win$is_onikuru, 2.4, 0.6 + pmin(win$n,200)/120)
      points(win$lon, win$lat, pch=19, col=ncol, cex=ncex)
      points(oni$lon, oni$lat, pch=19, col="black", cex=2.6); text(oni$lon, oni$lat, "ONIKURU", pos=3, font=2, cex=0.8)
      legend("topleft", legend=c("ONIKURU","center zone","non-center zone","ONIKURU edge","other edge"),
        pch=c(19,19,19,NA,NA), lty=c(NA,NA,NA,1,1), lwd=c(NA,NA,NA,2,2), col=c("black",C_CENTER,C_NON,"#444444","#bbbbbb"), bty="n", cex=0.75)
    }
    save_base(draw_gis, "F12_gis_transition_network", w=8, h=8)
  }
} else cat("F12_gis skipped: no coordinates (S_zone_coords.csv). Provide lng / lat columns to draw it.\n")

# ---- FS1 : sample flow diagram (CONSORT style, supplement) ----
ft <- rd("S_sample_flow.csv")
if (!is.null(ft)) {
  draw_flow <- function() {
    n <- nrow(ft); plot(NA, xlim=c(0,10), ylim=c(-0.5, n*1.4+0.5), axes=FALSE, xlab="", ylab="", main="Fig. S1  Sample flow diagram")
    for (i in seq_len(n)) { y <- (n-i)*1.4 + 0.4
      rect(1.2, y-0.55, 8.8, y+0.55, col="#eef3f8", border="#7f95b0", lwd=1.2)
      text(5, y+0.22, ft$step[i], font=2, cex=0.82)
      text(5, y-0.25, sprintf("rows = %s   |   uids = %s", format(ft$rows[i], big.mark=","), format(ft$uids[i], big.mark=",")), cex=0.75)
      if (i < n) arrows(5, y-0.55, 5, y-0.85, length=0.08, lwd=1.2) }
  }
  save_base(draw_flow, "FS1_sample_flow", w=7.5, h=max(7, nrow(ft)*1.05))
}

# ---- FS2 : permutation null distribution ----
pdd <- rd("S12_network_distance_permutation_distribution.csv"); ntt <- rd("S12_network_distance_permutation_test.csv")
if (!is.null(pdd) && !is.null(ntt)) {
  obs <- ntt$observed_distance[1]; pv <- ntt$permutation_p[1]
  pS2 <- ggplot2::ggplot(pdd, ggplot2::aes(distance)) + ggplot2::geom_histogram(bins=40, fill=C_GREY, alpha=.7) +
    ggplot2::geom_vline(xintercept=obs, color=C_NON, linewidth=1) +
    ggplot2::annotate("text", x=obs, y=Inf, vjust=1.5, hjust=-.05, label=sprintf("observed = %.1f\np = %.3f", obs, pv), color=C_NON) +
    ggplot2::labs(title="Fig. S2  Permutation null distribution of pre/post network distance", x="Frobenius distance (permuted)", y="count")
  save_pdf_gg(pS2, "FS2_permutation_null", w=8)
}

# ------------------------------------------------------------------------------
# B-3. Tables and figures of the additional analyses (TableA1..A11 / FA1..FA3)
# ------------------------------------------------------------------------------

rnd <- function(x, d = 3) { x <- copy(x); for (cc in names(x)) if (is.double(x[[cc]])) set(x, j = cc, value = round(x[[cc]], d)); x }
a_rr <- rd("S19_reading_rule.csv"); a_bv <- rd("S20_anchor_benchmark_variants.csv"); a_bb <- rd("S20_anchor_benchmark_bootstrap.csv")
a_h4 <- rd("S21_H4_anchorwise.csv"); a_h3 <- rd("S7_H3_logistic.csv"); a_h3d <- rd("S22_H3_logistic_with_distance.csv")
a_did <- rd("S23_DID_estimates.csv"); a_es <- rd("S23_DID_event_study.csv"); a_dm <- rd("S23_DID_group_means.csv")
a_dest <- rd("S24_destination_gravity_reference.csv"); a_blk <- rd("S25_block_prepost_change.csv")

if (!is.null(a_rr)) wt(rnd(a_rr), "TableA1_reading_rule.csv")
if (!is.null(a_bv)) {
  ta2 <- rnd(a_bv)
  if (!is.null(a_bb) && nrow(a_bb)) ta2 <- merge(ta2, rnd(a_bb)[, .(variant, is_ONIKURU = TRUE, share_ci_low, share_ci_high,
    rank_ci_low, rank_ci_high, p_above_median, p_above_both_cores)], by = c("variant", "is_ONIKURU"), all.x = TRUE, sort = FALSE)
  wt(ta2, "TableA2_benchmark_variants.csv")
}
if (!is.null(a_h4)) wt(rnd(a_h4)[, .(variant, label, is_ONIKURU, n_users, n_anchor_days, reach_anchor_day, reach_other_central_day,
  ref_size_biased, ref_any_central, fe_est, fe_ci_low, fe_ci_high, rank_fe)], "TableA3_H4_anchorwise.csv")
if (!is.null(a_h3)) wt(rnd(rbindlist(list(a_h3, a_h3d), fill = TRUE)), "TableA4_H3_logistic.csv")
if (!is.null(a_did)) { wt(rnd(a_did, 4), "TableA5_DID.csv"); if (!is.null(a_dm)) wt(rnd(a_dm, 4), "TableA5b_DID_group_means.csv") }
if (!is.null(a_dest)) wt(rnd(a_dest), "TableA6_destinations.csv")
if (!is.null(a_blk)) wt(rnd(a_blk), "TableA7_block_change.csv")
a_ab <- rd("S26_area_benchmark_prepost.csv"); a_ac <- rd("S26_area_benchmark_change.csv")
a_h2b <- rd("S27_H2_user_bootstrap.csv"); a_ph <- rd("S28_phase_split.csv")
if (!is.null(a_ab)) wt(rnd(a_ab), "TableA8_area_benchmark_prepost.csv")
if (!is.null(a_ac)) wt(rnd(a_ac), "TableA8b_area_benchmark_change.csv")
if (!is.null(a_h2b)) wt(rnd(a_h2b, 4), "TableA9_H2_user_bootstrap.csv")
if (!is.null(a_ph)) wt(rnd(a_ph), "TableA10_phase_split.csv")
a_h3a <- rd("S29_H3_anchorwise_new_uptake.csv")
if (!is.null(a_h3a)) wt(rnd(a_h3a), "TableA11_H3_anchorwise_new_uptake.csv")

# ---- FA1 : benchmark (article definition / home stays excluded) ----
if (!is.null(a_bv)) {
  bv <- a_bv[variant %in% c("manuscript definition", "home stays excluded")]
  bv[, variant := factor(variant, levels = c("manuscript definition", "home stays excluded"))]
  ord <- bv[variant == levels(variant)[1]][order(connecting_share), label]
  bv[, lab := factor(roma(label), levels = unique(roma(c(ord, setdiff(label, ord)))))]
  med <- bv[, .(m = median(connecting_share)), by = variant]
  pA1 <- ggplot2::ggplot(bv, ggplot2::aes(lab, connecting_share, fill = is_ONIKURU)) + ggplot2::geom_col(width = .72) +
    ggplot2::geom_hline(data = med, ggplot2::aes(yintercept = m), linetype = 2, color = "#444444") +
    ggplot2::coord_flip() + ggplot2::facet_wrap(~variant) +
    ggplot2::scale_fill_manual(values = c("FALSE" = C_GREY, "TRUE" = C_NON), guide = "none") +
    ggplot2::labs(title = "Fig. A1  Connecting-chain share by anchor", subtitle = "dashed line = anchor median",
                  x = NULL, y = "connecting-chain share")
  if (!is.null(a_bb) && nrow(a_bb)) {
    bb <- merge(bv[is_ONIKURU == TRUE, .(variant, lab)], a_bb[, .(variant, share_ci_low, share_ci_high)], by = "variant")
    bb[, variant := factor(variant, levels = levels(bv$variant))]
    pA1 <- pA1 + ggplot2::geom_errorbar(data = bb, ggplot2::aes(x = lab, ymin = share_ci_low, ymax = share_ci_high),
                                        width = .3, inherit.aes = FALSE)
  }
  save_pdf_gg(pA1, "FA1_benchmark_variants", w = 11, h = 6)
}

# ---- FA2 : event study ----
if (!is.null(a_es) && nrow(a_es[!is.na(est)])) {
  pA2 <- ggplot2::ggplot(a_es[!is.na(est)], ggplot2::aes(rel_week, est)) +
    ggplot2::geom_hline(yintercept = 0, color = "#444444") + ggplot2::geom_vline(xintercept = -0.5, linetype = 2, color = "#444444") +
    ggplot2::geom_linerange(ggplot2::aes(ymin = ci_low, ymax = ci_high), color = C_CENTER, linewidth = .7) +
    ggplot2::geom_point(color = C_CENTER, size = 2.2) +
    ggplot2::labs(title = "Fig. A2  Visitors vs non-visitors: other central-city stays by week",
                  subtitle = "difference relative to the reference week; person, week and day-of-week fixed effects; 95% CI",
                  x = "week relative to opening (0 = opening week)", y = "difference in share of days")
  save_pdf_gg(pA2, "FA2_event_study", w = 9, h = 4.8)
}

# ---- FA3 : onward destinations, observed vs expected ----
if (!is.null(a_dest)) {
  dl <- melt(a_dest[, .(group, observed = as.numeric(observed), expected)], id.vars = "group", variable.name = "kind", value.name = "n")
  pA3 <- ggplot2::ggplot(dl, ggplot2::aes(group, n, fill = kind)) + ggplot2::geom_col(position = "dodge", width = .7) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", n)), position = ggplot2::position_dodge(.7), vjust = -.3, size = 3) +
    ggplot2::scale_fill_manual(values = c(observed = C_NON, expected = C_GREY), name = NULL) +
    ggplot2::labs(title = "Fig. A3  Onward destinations from ONIKURU: observed vs size-and-distance reference",
                  x = NULL, y = "onward transitions")
  save_pdf_gg(pA3, "FA3_destinations", w = 8, h = 4.8)
}

if (WRITE_DIAG) cat("\nDiagnostic tables -> ", TAB, "\n  Diagnostic figures (pdf+png) -> ", FIG, "\n", sep="")
# ------------------------------------------------------------------------------
# B-4. Tables and figures in the layout of the article (written to tables/ and figs/)
#      CSV tables and figures with the rows, columns, order and notation of the article (Tables 1-3, S1-S6).
#        tables: tables/Table1.csv .. Table3.csv, TableS1.csv .. TableS6.csv
#            tables/TableS4.csv           ... this specification + home stays excluded + post-opening phases (one column per run)
#            tables/TableNotes.csv        ... numbers quoted in table notes and in the text (medians, intervals, ranks)
#        figures: figs/Figure4_nodes.csv, Figure4_edges.csv, Figure4_reference, Figure5, Figure6, Figure7, Figure8, FigureS2
#            figs/figure_map.csv          ... article figure numbers and files
#      Values are computed from the unrounded S files and rounded in the tables; minus sign U+2212, thousands separated.
# ------------------------------------------------------------------------------

TAB <- TAB_M; FIG <- FIG_M   # from here on: article tables and figures
MINUS <- "−"
fnum <- function(x, d = 3) { s <- sprintf(paste0("%.", d, "f"), x); sub("^-", MINUS, s) }
fsgn <- function(x, d = 3) { s <- sprintf(paste0("%+.", d, "f"), x); sub("^-", MINUS, s) }
fint <- function(x) formatC(round(x), format = "d", big.mark = ",")
fci  <- function(est, lo, hi, d = 3) sprintf("%s [%s, %s]", fnum(est, d), fnum(lo, d), fnum(hi, d))
fp   <- function(p) ifelse(is.na(p), "–", ifelse(p < 0.0005, "<0.001", ifelse(p < 0.10, sprintf("%.3f", p), sprintf("%.2f", p))))
fp3  <- function(p) ifelse(is.na(p), "–", ifelse(p < 0.0005, "<0.001", sprintf("%.3f", p)))
DASH <- "—"; NOTES <- list()
note <- function(key, value, where = "") NOTES[[length(NOTES) + 1L]] <<- data.table(item = key, value = as.character(value), where = where)
mrd  <- function(f) { x <- rd(f); if (is.null(x)) msg("[B-4] missing: %s", f); x }
ord_num <- function(k) c("first","second","third","fourth","fifth","sixth","seventh","eighth","ninth","tenth",
                         "eleventh","twelfth","thirteenth","fourteenth","fifteenth","sixteenth","seventeenth")[k]

m_s2 <- mrd("S2_data_diagnostics.csv"); m_cs <- mrd("S4_chain_summary_by_period.csv"); m_zd <- mrd("S3_zone_diagnostics.csv")
m_fl <- mrd("S_sample_flow.csv"); m_ca <- mrd("S6_H1_centrality_ONIKURU_AREA.csv"); m_cf <- mrd("S6_H1_centrality_ONIKURU_facility_post.csv")
m_dec <- mrd("S6b_onikuru_outflow_decomposition.csv"); m_h2 <- mrd("S6_H2_center_to_ONIKURU_AREA_prepost.csv"); m_h2b <- mrd("S27_H2_user_bootstrap.csv")
m_w1 <- mrd("S17_within_stage1_descriptive.csv"); m_w2 <- mrd("S17_within_stage2_descriptive.csv"); m_wr <- mrd("S17_within_stage2_regression.csv")
m_bc <- mrd("S11_H4_cluster_bootstrap_CI.csv"); m_ct <- mrd("S9_H5_chain_type_summary.csv"); m_perm <- mrd("S12_network_distance_permutation_test.csv")
m_bv <- mrd("S20_anchor_benchmark_variants.csv"); m_bb <- mrd("S20_anchor_benchmark_bootstrap.csv"); m_hd <- mrd("S20_home_stay_diagnostics.csv")
m_rr <- mrd("S19_reading_rule.csv"); m_h4a <- mrd("S21_H4_anchorwise.csv"); m_h3a <- mrd("S29_H3_anchorwise_new_uptake.csv")
m_h3 <- mrd("S7_H3_logistic.csv"); m_h3d <- mrd("S22_H3_logistic_with_distance.csv"); m_h3p <- mrd("S22_H3_predicted.csv"); m_h3c <- mrd("S22_H3_visitor_composition.csv")
m_did <- mrd("S23_DID_estimates.csv"); m_dm <- mrd("S23_DID_group_means.csv"); m_es <- mrd("S23_DID_event_study.csv")
m_ab <- mrd("S26_area_benchmark_prepost.csv"); m_ac <- mrd("S26_area_benchmark_change.csv"); m_ph <- mrd("S28_phase_split.csv")
gv <- function(it) as.numeric(m_s2[item == it, value][1]); zv <- function(it) as.numeric(m_zd[item == it, value][1])
obs_pre <- gv("pre_observed_days"); obs_post <- gv("post_observed_days")

# ---- Table 1 : sample characteristics ----
if (!is.null(m_s2) && !is.null(m_cs) && !is.null(m_zd) && !is.null(m_fl)) {
  pr <- m_cs[period == "pre"]; po <- m_cs[period == "post"]
  users_all <- m_fl[step == "Remove pre-opening ONIKURU pts", uids]; vis <- m_fl[step == "ONIKURU post visitors (uid)", uids]
  r <- function(...) as.list(c(...))
  M1 <- rbindlist(list(
    r("Analytic sample", "", "", "", ""),
    r("Observed days:", "(Days)", obs_pre, obs_post, obs_pre + obs_post),
    r("Analytic stay points (after all filters):", "(Stay points)", fint(gv("pre_rows")), fint(gv("post_rows")), fint(gv("pre_rows") + gv("post_rows"))),
    r("Daily trip chains:", "(Person-days)", fint(pr$person_days), fint(po$person_days), fint(pr$person_days + po$person_days)),
    r("Unique individuals observed:", "(Users)", fint(pr$uids), fint(po$uids), fint(users_all)),
    r("Trip-chain characteristics", "", "", "", ""),
    r("Mean chain length:", "(Nodes/person-day)", fnum(pr$mean_chain_len), fnum(po$mean_chain_len), DASH),
    r("Mean unique nodes per day:", "(Nodes/person-day)", fnum(pr$mean_unique_nodes), fnum(po$mean_unique_nodes), DASH),
    r("Central-city exposure", "", "", "", ""),
    r("Person-days with at least one central-city stay:", "(Person-days)", fint(pr$center_days), fint(po$center_days), fint(pr$center_days + po$center_days)),
    r("ONIKURU exposure", "", "", "", ""),
    r("ONIKURU visit-days:", "(Person-days)", "0", fint(po$onikuru_days), fint(po$onikuru_days)),
    r("ONIKURU visitors:", "(Unique users)", "0", fint(vis), fint(vis)),
    r("Zone-level diagnostics", "", "", "", ""),
    r("ONIKURU-facility:", "(Stay points)", DASH, fint(zv("ONIKURU_facility_post")), fint(zv("ONIKURU_facility_post"))),
    r("ONIKURU-area:", "(Stay points)", fint(zv("ONIKURU_AREA_pre")), fint(zv("ONIKURU_AREA_post")), fint(zv("ONIKURU_AREA_pre") + zv("ONIKURU_AREA_post"))),
    r("Central-city zones:", "(Stay points)", fint(zv("center_pre")), fint(zv("center_post")), fint(zv("center_pre") + zv("center_post")))))
  setnames(M1, c("Measure", "Unit", "Pre-opening", "Post-opening", "Overall")); wt(M1, "Table1.csv")
  note("t1_missing_post_days", gv("post_calendar_days") - obs_post, "Table 1 note")
  note("t1_person_days_per_calendar_day_pre_post", sprintf("%.0f / %.0f", pr$person_days / obs_pre, po$person_days / obs_post), "4.1")
  note("t1_mean_observed_days_per_user", fnum((pr$person_days + po$person_days) / users_all, 1), "4.1")
  note("t1_area_stays_per_1000_person_days_pre_post", sprintf("%.1f / %.1f", 1000 * zv("ONIKURU_AREA_pre") / pr$person_days, 1000 * zv("ONIKURU_AREA_post") / po$person_days), "4.1")
}

# ---- Table 2 : centrality, outflow decomposition, connecting-chain share ----
if (!is.null(m_ca) && !is.null(m_cf) && !is.null(m_dec) && !is.null(m_ab) && !is.null(m_rr)) {
  cpre <- m_ca[period == "pre"]; cpost <- m_ca[period == "post"]
  dpre <- m_dec[level == "area" & period == "pre"]; dpost <- m_dec[level == "area" & period == "post"]; dfac <- m_dec[level == "facility" & period == "post"]
  bn <- function(x) x$betweenness / ((x$n_nodes - 1) * (x$n_nodes - 2))
  abm <- m_ab[variant == "manuscript definition" & is_ONIKURU_AREA == TRUE]
  cs_pre <- abm[period == "pre", n_connecting / n_days]; cs_post <- abm[period == "post", n_connecting / n_days]
  cs_fac <- m_rr[variant == "manuscript definition", connecting_share]
  two <- function(a, b, d = 0) sprintf("%s / %s", fnum(a, d), fnum(b, d)); dtwo <- function(a, b, d = 0) sprintf("%s / %s", fsgn(a, d), fsgn(b, d))
  M2 <- rbindlist(list(
    list("Betweenness (normalized):", fnum(bn(cpre), 4), fnum(bn(cpost), 4), fsgn(bn(cpost) - bn(cpre), 4), fnum(bn(m_cf), 4)),
    list("Out-degree / In-degree:", two(cpre$out_degree, cpre$in_degree), two(cpost$out_degree, cpost$in_degree), dtwo(cpost$out_degree - cpre$out_degree, cpost$in_degree - cpre$in_degree), two(m_cf$out_degree, m_cf$in_degree)),
    list("Out-strength / In-strength:", two(cpre$out_strength, cpre$in_strength), two(cpost$out_strength, cpost$in_strength), dtwo(cpost$out_strength - cpre$out_strength, cpost$in_strength - cpre$in_strength), two(m_cf$out_strength, m_cf$in_strength)),
    list("per observed day:", two(cpre$out_strength / obs_pre, cpre$in_strength / obs_pre, 1), two(cpost$out_strength / obs_post, cpost$in_strength / obs_post, 1),
         dtwo(cpost$out_strength / obs_post - cpre$out_strength / obs_pre, cpost$in_strength / obs_post - cpre$in_strength / obs_pre, 1), two(m_cf$out_strength / obs_post, m_cf$in_strength / obs_post, 1)),
    list("PageRank:", fnum(cpre$pagerank), fnum(cpost$pagerank), fsgn(cpost$pagerank - cpre$pagerank), fnum(m_cf$pagerank)),
    list("Central-city outflow share:", fnum(dpre$out_center_share), fnum(dpost$out_center_share), fsgn(dpost$out_center_share - dpre$out_center_share), fnum(dfac$out_center_share)),
    list("Central-city onward-flow ratio:", fnum(dpre$onward_to_center), fnum(dpost$onward_to_center), fsgn(dpost$onward_to_center - dpre$onward_to_center), fnum(dfac$onward_to_center)),
    list("Flow-through ratio (outflow/inflow):", fnum(dpre$flow_through), fnum(dpost$flow_through), fsgn(dpost$flow_through - dpre$flow_through), fnum(dfac$flow_through)),
    list("Connecting-chain share:", fnum(cs_pre), fnum(cs_post), fsgn(cs_post - cs_pre), fnum(cs_fac))))
  setnames(M2, c("", "ONIKURU-area: pre-opening", "ONIKURU-area: post-opening", "Change", "ONIKURU-facility: post-opening")); wt(M2, "Table2.csv")
  note("t2_observed_days_pre_post", sprintf("%d and %d", obs_pre, obs_post), "Table 2 note")
  note("t2_betweenness_raw_pre_post_fac", sprintf("%s / %s / %s", fint(cpre$betweenness), fint(cpost$betweenness), fint(m_cf$betweenness)), "4.2")
  note("t2_n_nodes_pre_post", sprintf("%s / %s", fint(cpre$n_nodes), fint(cpost$n_nodes)), "4.2")
  note("t2_facility_out_to_center_of_out", sprintf("%d of %d", dfac$out_to_center, dfac$out_strength), "4.2")
  if (!is.null(m_ac)) { acm <- m_ac[variant == "manuscript definition"]
    note("t2_area_connecting_change_ci", sprintf("%s [%s, %s]", fsgn(acm$change), fsgn(acm$change_ci_low), fsgn(acm$change_ci_high)), "Table 2 note / Table S6")
    note("t2_area_connecting_rank_pre_post", sprintf("%d of %d to %d of %d", acm$pre_rank, acm$n_units, acm$post_rank, acm$n_units), "4.5") }
}

# ---- Table 3 : within-person comparisons ----
if (!is.null(m_w1) && !is.null(m_w2) && !is.null(m_wr)) {
  o <- m_w2[daytype == "ONIKURU_day"]; al <- m_w1[daytype == "non_ONIKURU_day_all"]; ce <- m_w2[daytype == "other_center_day"]
  g <- function(pat) m_wr[grepl(pat, outcome, fixed = TRUE)][1]
  r_reach <- g("onward"); r_nodes <- g("n_center_nodes (total)"); r_ctrl <- g("chain_len")
  ref_size <- (ce$mean_n_center_nodes - (1 - ce$onward_rate)) / ce$mean_n_center_nodes; ref_any <- ce$person_days / al$person_days
  M3 <- rbindlist(list(
    list("Central-city reach rate (≥2 central-city nodes)", fnum(o$onward_rate), fnum(al$onward_rate), fnum(ce$onward_rate), fci(r_reach$est, r_reach$ci_low, r_reach$ci_high), fp3(r_reach$p)),
    list("Central-city nodes/day", fnum(o$mean_n_center_nodes), fnum(al$mean_n_center_nodes), fnum(ce$mean_n_center_nodes), fci(r_nodes$est, r_nodes$ci_low, r_nodes$ci_high), fp3(r_nodes$p)),
    list("(+ chain-length control)", DASH, DASH, DASH, fci(r_ctrl$est, r_ctrl$ci_low, r_ctrl$ci_high), fp3(r_ctrl$p)),
    list("Daily chain length", fnum(o$mean_chain_len), fnum(al$mean_chain_len), fnum(ce$mean_chain_len), DASH, DASH),
    list("n (person-days)", fint(o$person_days), fint(al$person_days), fint(ce$person_days), DASH, DASH)))
  setnames(M3, c("", "ONIKURU", "All non-ONIKURU", "Non-ONIKURU central", "Adjusted [95% CI]", "p-value")); wt(M3, "Table3.csv")
  note("t3_ref_size_biased_formula", sprintf("(%s − %s)/%s = %s", fnum(ce$mean_n_center_nodes), fnum(1 - ce$onward_rate), fnum(ce$mean_n_center_nodes), fnum(ref_size)), "Table 3 note")
  note("t3_ref_any_central_formula", sprintf("%s/%s = %s", fint(ce$person_days), fint(al$person_days), fnum(ref_any)), "Table 3 note")
  if (!is.null(m_ct)) { cs <- function(k) m_ct[chain_type == k, share]
    note("t3_reach_decomposition_onward_vs_preceding", sprintf("%s = %s + %s", fnum(o$onward_rate), fnum(cs("post_bridge") + cs("middle_bridge")), fnum(cs("final_after_center"))), "4.4") }
  if (!is.null(m_bc)) for (k in m_bc$metric) note(paste0("t3_", k), fci(m_bc[metric == k, estimate], m_bc[metric == k, ci_low], m_bc[metric == k, ci_high]), "4.4 / Figure 7c")
}

# ---- Table S1 : bidirectional transitions, all central town-blocks ----
if (!is.null(m_h2)) {
  h <- copy(m_h2); h[, lab := ifelse(zone_label == "ALL_CENTER", "Central-city zones", roma(zone_label))]
  one <- function(d) { d <- copy(d); zero <- (d$pre_hit_edges + d$post_hit_edges) == 0
    d[, `:=`(pre = fnum(pre_share), post = fnum(post_share), delta = ifelse(zero, "0.000", fsgn(post_share - pre_share)),
             p = ifelse(zero, "–", fp3(fisher_p)), bh = ifelse(zero | zone_label == "ALL_CENTER", "–", fp3(fisher_p_BH)))]
    d[, .(lab, pre, post, delta, p, bh)] }
  inn <- one(h[direction == "inflow"]); out <- one(h[direction == "outflow"])
  ord <- h[direction == "inflow"][order(-(zone_label == "ALL_CENTER"), -post_share), lab]
  S1 <- merge(inn, out, by = "lab", suffixes = c("_in", "_out"))[match(ord, lab)]
  setnames(S1, c("Zone", "Inflow: Pre share", "Inflow: Post share", "Inflow: Delta", "Inflow: Fisher p", "Inflow: BH p",
                 "Outflow: Pre share", "Outflow: Post share", "Outflow: Delta", "Outflow: Fisher p", "Outflow: BH p")); wt(S1, "TableS1.csv")
  ai <- h[direction == "inflow" & zone_label == "ALL_CENTER"]; ao <- h[direction == "outflow" & zone_label == "ALL_CENTER"]
  note("s1_inflow_denominators_pre_post", sprintf("%s and %s", fint(ai$pre_out_edges), fint(ai$post_out_edges)), "Table S1 note")
  note("s1_outflow_denominators_pre_post", sprintf("%s and %s", fint(ao$pre_out_edges), fint(ao$post_out_edges)), "Table S1 note")
  note("s1_aggregate_inflow_counts", sprintf("%d/%s -> %d/%s, Fisher p = %s", ai$pre_hit_edges, fint(ai$pre_out_edges), ai$post_hit_edges, fint(ai$post_out_edges), fp3(ai$fisher_p)), "4.2")
  note("s1_aggregate_outflow_counts", sprintf("%d/%s -> %d/%s, Fisher p = %s", ao$pre_hit_edges, fint(ao$pre_out_edges), ao$post_hit_edges, fint(ao$post_out_edges), fp3(ao$fisher_p)), "4.2")
  sig_in <- h[direction == "inflow" & zone_label != "ALL_CENTER" & !is.na(fisher_p_BH) & fisher_p_BH < 0.05]
  sig_out <- h[direction == "outflow" & zone_label != "ALL_CENTER" & !is.na(fisher_p_BH) & fisher_p_BH < 0.05]
  note("s1_blocks_BH_significant_inflow", if (nrow(sig_in)) paste(sprintf("%s (%d/%d -> %d/%d, BH p = %s)", sig_in$lab, sig_in$pre_hit_edges, sig_in$pre_out_edges, sig_in$post_hit_edges, sig_in$post_out_edges, fp3(sig_in$fisher_p_BH)), collapse = "; ") else "none", "Table S1 note / 4.2")
  note("s1_blocks_BH_significant_outflow", if (nrow(sig_out)) paste(sig_out$lab, collapse = "; ") else "none", "Table S1 note / 4.2")
  if (!is.null(m_h2b)) { bi <- m_h2b[grepl("^inflow", direction)]; bo <- m_h2b[grepl("^outflow", direction)]
    note("s1_user_bootstrap_inflow_change_ci", sprintf("%s [%s, %s]", fsgn(bi$change), fsgn(bi$change_ci_low), fsgn(bi$change_ci_high)), "Table S1 note / 4.2")
    note("s1_user_bootstrap_outflow_change_ci", sprintf("%s [%s, %s]", fsgn(bo$change), fsgn(bo$change_ci_low), fsgn(bo$change_ci_high)), "Table S1 note / 4.2")
    note("s1_user_bootstrap_B_users", sprintf("B = %s, %s users", fint(bi$B), fint(bi$n_users)), "Table S1 note") }
}

# ---- Table S3 : multi-anchor benchmark (H5 + H4 + H3) ----
if (!is.null(m_bv) && !is.null(m_h4a) && !is.null(m_h3a)) {
  bm <- m_bv[variant == "manuscript definition"][, .(label, is_ONIKURU, n_days, conn = n_connecting / n_days, ocs = round(out_center_share * out_strength) / out_strength,
                                                     onw = round(out_center_share * out_strength) / in_strength, ft = out_strength / in_strength)]   # recover the counts from the rounded share and recompute
  h4 <- m_h4a[variant == "manuscript definition"][, .(label, reach = reach_anchor_day, fe_est, fe_ci_low, fe_ci_high, rank_fe)]
  h3 <- m_h3a[, .(label, n_new_users, n_users, OR, OR_low, OR_high, rank_OR_flattest)]
  S2 <- merge(merge(bm, h4, by = "label"), h3, by = "label")[order(-conn)]
  S2[, en := roma(label)]
  out <- S2[, .(Anchor = en, `Anchor-days, n` = fint(n_days), `Central-city outflow share` = fnum(ocs), `Connecting-chain share` = fnum(conn),
                `Reach rate` = fnum(reach), `FE contrast [95% CI]` = fci(fe_est, fe_ci_low, fe_ci_high),
                `New users / at risk` = sprintf("%s / %s", fint(n_new_users), fint(n_users)), `OR [95% CI]` = fci(OR, OR_low, OR_high))]
  wt(out, "TableS3.csv")
  on <- S2[is_ONIKURU == TRUE]; oth <- S2[is_ONIKURU == FALSE]
  note("s2_median_outflow_share", fnum(median(S2$ocs)), "Table S3 note / 4.5"); note("s2_median_connecting_share", fnum(median(S2$conn)), "Table S3 note / 4.5")
  note("s2_median_onward_flow_ratio_not_tabulated", fnum(median(S2$onw)), "(reference only)")
  note("s2_onikuru_rank_connecting", sprintf("%d of %d", rank(-S2$conn, ties.method = "min")[S2$is_ONIKURU], nrow(S2)), "4.5")
  note("s2_onikuru_rank_outflow", sprintf("%d of %d", rank(-S2$ocs, ties.method = "min")[S2$is_ONIKURU], nrow(S2)), "4.5")
  note("s2_station_cores_connecting", paste(sprintf("%s %s", roma(CFG$station_cores), fnum(S2[label %in% CFG$station_cores, conn][match(CFG$station_cores, S2[label %in% CFG$station_cores, label])])), collapse = "; "), "4.5")
  big <- oth[n_days >= 1000 & conn < on$conn & ocs < on$ocs, en]
  note("s2_large_anchors_exceeded_on_both", if (length(big)) paste(big, collapse = ", ") else "none", "Table S3 note")
  note("s2_spearman_connecting_vs_days", fnum(cor(S2$conn, S2$n_days, method = "spearman"), 2), "4.5")
  note("s2_spearman_connecting_vs_outflow", fnum(cor(S2$conn, S2$ocs, method = "spearman"), 2), "4.5")
  note("s2_anchors_self_contained_above_80pct", sprintf("%d of %d", sum((1 - S2$conn) > 0.80), nrow(S2)), "4.5")
  cmp <- S2[n_days >= on$n_days / 2 & n_days <= on$n_days * 2]   # anchors with 1/2 to 2 times the anchor-days of ONIKURU
  note("s2_comparable_size_band_days", sprintf("%s\u2013%s anchor-days (%d anchors)", fint(on$n_days / 2), fint(on$n_days * 2), nrow(cmp)), "4.5")
  note("s2_rank_among_comparable_size", sprintf("%d of %d", rank(-cmp$conn, ties.method = "min")[cmp$is_ONIKURU], nrow(cmp)), "4.5")
  if (!is.null(m_bb)) { b <- m_bb[variant == "manuscript definition"]
    note("s2_onikuru_share_boot_ci", sprintf("%s–%s", fnum(b$share_ci_low), fnum(b$share_ci_high)), "Table S3 note")
    note("s2_onikuru_rank_boot_interval", sprintf("%d–%d", b$rank_ci_low, b$rank_ci_high), "Table S3 note")
    note("s2_p_above_median_and_cores", sprintf("%.4f / %.4f", b$p_above_median, b$p_above_both_cores), "Table S3 note") }
  note("s2_onikuru_rank_fe", sprintf("%s of %d", ord_num(on$rank_fe), nrow(S2)), "Table S3 note / 4.4")
  note("s2_onikuru_rank_OR_flattest", sprintf("%s flattest of %d", ord_num(on$rank_OR_flattest), nrow(S2)), "Table S3 note / 4.3")
  note("s2_median_OR_other_anchors", fnum(median(oth$OR), 2), "Table S3 note / 4.3")
  note("s2_all_gradients_positive", all(S2$OR_low > 1 | S2$OR > 1), "Table S3 note")
}

# ---- Table S4 : this specification + home stays excluded + phases (one column per run) ----
if (!is.null(m_rr) && !is.null(m_dec) && !is.null(m_w2) && !is.null(m_wr) && !is.null(m_h2)) {
  rm_ <- m_rr[variant == "manuscript definition"]; rh_ <- m_rr[variant == "home stays excluded"]
  ph1 <- if (!is.null(m_ph) && nrow(m_ph) >= 1) m_ph[1] else NULL; ph2 <- if (!is.null(m_ph) && nrow(m_ph) >= 2) m_ph[2] else NULL
  dfac <- m_dec[level == "facility" & period == "post"]; dpost <- m_dec[level == "area" & period == "post"]
  ai <- m_h2[direction == "inflow" & zone_label == "ALL_CENTER"]; ao <- m_h2[direction == "outflow" & zone_label == "ALL_CENTER"]
  o <- m_w2[daytype == "ONIKURU_day"]; ce <- m_w2[daytype == "other_center_day"]; r_reach <- m_wr[grepl("onward", outcome, fixed = TRUE)][1]
  hv <- if (!is.null(m_bv)) m_bv[variant == "home stays excluded" & label == CFG$onikuru_label] else NULL
  hh4 <- if (!is.null(m_h4a)) m_h4a[variant == "home stays excluded" & label == CFG$onikuru_label] else NULL
  acm <- if (!is.null(m_ac)) m_ac[variant == "manuscript definition"] else NULL; ach <- if (!is.null(m_ac)) m_ac[variant == "home stays excluded"] else NULL
  po <- m_cs[period == "post"]
  rule <- function(x) if (is.null(x) || !nrow(x)) NA else (isTRUE(x$above_median_connecting) && isTRUE(x$above_median_outflow) && isTRUE(x$above_both_station_cores) && isTRUE(x$proviso_holds))
  rule_ph <- function(q) if (is.null(q)) NA else (q$connecting_share > q$connecting_median && q$outflow_share > q$outflow_median && isTRUE(q$above_both_station_cores) &&
                                                   q$visit_day_reach >= min(q$ref_size_biased, q$ref_any_central))
  yn <- function(v) ifelse(is.na(v), DASH, ifelse(v, "Yes", "No"))
  pq <- function(q, f) if (is.null(q)) DASH else f(q)
  nz <- function(x) if (is.null(x) || !nrow(x)) DASH else x
  row <- function(label, main, home = DASH, p1 = DASH, p2 = DASH) list(label, main, home, p1, p2)
  S3 <- rbindlist(list(
    row("Post-opening observed days:", as.character(obs_post), as.character(obs_post), pq(ph1, function(q) as.character(q$observed_days)), pq(ph2, function(q) as.character(q$observed_days))),
    row("Stay points (rows):", fint(gv("pre_rows") + gv("post_rows")), if (!is.null(m_hd)) fint(gv("pre_rows") + gv("post_rows") - as.numeric(m_hd[item == "stays_home", value])) else DASH),
    row("ONIKURU visit person-days:", fint(po$onikuru_days), fint(po$onikuru_days), pq(ph1, function(q) fint(q$onikuru_visit_days)), pq(ph2, function(q) fint(q$onikuru_visit_days))),
    row("Visit-days per 1,000 person-days:", fnum(1000 * po$onikuru_days / po$person_days, 1), DASH, pq(ph1, function(q) fnum(q$visit_days_per_1000_person_days, 1)), pq(ph2, function(q) fnum(q$visit_days_per_1000_person_days, 1))),
    row("H1 outflow share", "", "", "", ""),
    row("    ONIKURU-area central-city:", fnum(dpost$out_center_share), if (!is.null(ach)) fnum(ach$post_outflow_share) else DASH),
    row("    ONIKURU-facility central-city:", fnum(dfac$out_center_share), if (!is.null(hv) && nrow(hv)) fnum(round(hv$out_center_share * hv$out_strength) / hv$out_strength) else DASH, pq(ph1, function(q) fnum(q$outflow_share)), pq(ph2, function(q) fnum(q$outflow_share))),
    row("    ONIKURU-facility flow-through ratio:", fnum(dfac$flow_through), if (!is.null(hv) && nrow(hv)) fnum(hv$out_strength / hv$in_strength) else DASH, pq(ph1, function(q) fnum(q$flow_through)), pq(ph2, function(q) fnum(q$flow_through))),
    row("H2 inflow share", "", "", "", ""),
    row("    Pre-period:", fnum(ai$pre_share)), row("    Post-period:", fnum(ai$post_share)), row("    Fisher p:", fp3(ai$fisher_p)),
    row("H2 outflow share", "", "", "", ""),
    row("    Pre-period:", fnum(ao$pre_share)), row("    Post-period:", fnum(ao$post_share)), row("    Fisher p:", fp3(ao$fisher_p)),
    row("H3 uptake odds ratio per +0.10 [95% CI]:", if (!is.null(m_h3)) { a <- m_h3[model == "adjusted" & term == "pre_center_day_rate"]; fci(a$OR, a$OR_low, a$OR_high) } else DASH),
    row("H4 within-person reach rate", "", "", "", ""),
    row("    ONIKURU central day:", fnum(o$onward_rate), if (!is.null(hh4) && nrow(hh4)) fnum(hh4$reach_anchor_day) else DASH, pq(ph1, function(q) fnum(q$visit_day_reach)), pq(ph2, function(q) fnum(q$visit_day_reach))),
    row("    non-ONIKURU central day:", fnum(ce$onward_rate), if (!is.null(hh4) && nrow(hh4)) fnum(hh4$reach_other_central_day) else DASH),
    row("    Reference, size-biased:", fnum(rm_$ref_size_biased), fnum(rh_$ref_size_biased), pq(ph1, function(q) fnum(q$ref_size_biased)), pq(ph2, function(q) fnum(q$ref_size_biased))),
    row("    Reference, any central day:", fnum(rm_$ref_any_central), fnum(rh_$ref_any_central), pq(ph1, function(q) fnum(q$ref_any_central)), pq(ph2, function(q) fnum(q$ref_any_central))),
    row("    Fixed-effect contrast [95% CI]:", fci(r_reach$est, r_reach$ci_low, r_reach$ci_high), fci(rh_$h4_fe, rh_$h4_fe_ci_low, rh_$h4_fe_ci_high),
        pq(ph1, function(q) fci(q$h4_fe, q$h4_fe_ci_low, q$h4_fe_ci_high)), pq(ph2, function(q) fci(q$h4_fe, q$h4_fe_ci_low, q$h4_fe_ci_high))),
    row("H4 same-day secondary-stay rate [95% CI]:", if (!is.null(m_bc)) { b <- m_bc[metric == "same_day_secondary_rate"]; fci(b$estimate, b$ci_low, b$ci_high) } else DASH),
    row("H5 ONIKURU-only-chain share:", if (!is.null(m_ct)) fnum(m_ct[chain_type == "onikuru_only_center", share]) else DASH, DASH, pq(ph1, function(q) fnum(q$onikuru_only_share)), pq(ph2, function(q) fnum(q$onikuru_only_share))),
    row("H5 connecting-chain share [95% CI]:", fci(rm_$connecting_share, rm_$boot_share_ci_low, rm_$boot_share_ci_high), fci(rh_$connecting_share, rh_$boot_share_ci_low, rh_$boot_share_ci_high),
        pq(ph1, function(q) fci(q$connecting_share, q$share_ci_low, q$share_ci_high)), pq(ph2, function(q) fci(q$connecting_share, q$share_ci_low, q$share_ci_high))),
    row("    Anchor median:", fnum(rm_$connecting_median), fnum(rh_$connecting_median), pq(ph1, function(q) fnum(q$connecting_median)), pq(ph2, function(q) fnum(q$connecting_median))),
    row(sprintf("    Rank among %d anchors:", rm_$n_anchors), as.character(rm_$connecting_rank), as.character(rh_$connecting_rank), pq(ph1, function(q) as.character(q$connecting_rank)), pq(ph2, function(q) as.character(q$connecting_rank))),
    row("Reading rule met:", yn(rule(rm_)), yn(rule(rh_)), yn(rule_ph(ph1)), yn(rule_ph(ph2))),
    row("Difference-in-differences [95% CI]:", if (!is.null(m_did)) { d <- m_did[model == "visitor x post"]; fci(d$est, d$ci_low, d$ci_high) } else DASH),
    row("ONIKURU-area connecting-chain share, pre / post:", if (!is.null(acm)) sprintf("%s / %s", fnum(acm$pre_share), fnum(acm$post_share)) else DASH,
        if (!is.null(ach)) sprintf("%s / %s", fnum(ach$pre_share), fnum(ach$post_share)) else DASH),
    row("Permutation p:", if (!is.null(m_perm)) fnum(m_perm$permutation_p) else DASH)))
  setnames(S3, c("row", CFG$spec, "home_stays_excluded", "phase_1", "phase_2"))
  if (!is.null(ph1)) note("s3_phase_labels", paste(m_ph$phase, collapse = " | "), "Table S4 header")
  wt(S3, "TableS4.csv")
}

# ---- Table S2 : H3 logistic, Model 1 / Model 2 ----
if (!is.null(m_h3) && !is.null(m_h3d)) {
  g1 <- function(t) m_h3[model == "adjusted" & term == t]; g2 <- function(t) m_h3d[model == "adjusted + log distance" & term == t]
  orc <- function(x) if (nrow(x)) fci(x$OR, x$OR_low, x$OR_high) else DASH; pc <- function(x) if (nrow(x)) fp(x$p) else DASH
  S4 <- rbindlist(list(
    list("Pre-opening central-city use rate (per +0.10)", orc(g1("pre_center_day_rate")), pc(g1("pre_center_day_rate")), orc(g2("pre_center_day_rate")), pc(g2("pre_center_day_rate"))),
    list("Distance from home to ONIKURU (per doubling)", DASH, DASH, orc(g2("log_km")), pc(g2("log_km"))),
    list("Pre-opening observed days (per day)", orc(g1("pre_observed_days")), pc(g1("pre_observed_days")), orc(g2("pre_observed_days")), pc(g2("pre_observed_days"))),
    list("Post-opening observed days (per day)", orc(g1("post_observed_days")), pc(g1("post_observed_days")), orc(g2("post_observed_days")), pc(g2("post_observed_days")))))
  setnames(S4, c("Predictor", "Model 1: OR [95% CI]", "Model 1: p", "Model 2 (with distance): OR [95% CI]", "Model 2: p")); wt(S4, "TableS2.csv")
  a <- g1("pre_center_day_rate"); note("s4_n_users_visitors", sprintf("n = %s users, of whom %s became visitors", fint(a$n_users), fint(a$n_visitors)), "Table S2 note / 4.3")
  if (!is.null(m_h3p)) { p0 <- m_h3p[pre_center_day_rate == 0]; p5 <- m_h3p[pre_center_day_rate == 0.5]
    note("s4_predicted_model1_at_0_and_0.5", sprintf("%s and %s", fnum(p0$p_adjusted), fnum(p5$p_adjusted)), "Table S2 note / 4.3")
    note("s4_predicted_model2_at_0_and_0.5", sprintf("%s and %s", fnum(p0$p_with_distance), fnum(p5$p_with_distance)), "Table S2 note / 4.3") }
  if (!is.null(m_h3c)) { cv <- function(k) as.numeric(m_h3c[item == k, value])
    note("s4_visitors_with_zero_pre_central_use", sprintf("%d of %d (%.1f%%) vs %.1f%% of all users", cv("visitors_with_zero_pre_central_use"), cv("visitors_in_frame"), 100 * cv("share_zero_among_visitors"), 100 * cv("share_zero_among_all")), "Table S2 note / 4.3")
    note("s4_median_distance_km_visitors_nonvisitors", sprintf("%.2f / %.2f", cv("median_distance_km_visitors"), cv("median_distance_km_non_visitors")), "4.3") }
}

# ---- Table S5 : difference-in-differences ----
if (!is.null(m_did) && !is.null(m_dm)) {
  gm <- function(g, p) m_dm[group == g & period == p]
  d1 <- m_did[model == "visitor x post"]; d2 <- m_did[grepl("^home within", model)]
  vz <- function(p) m_dm[group == "visitors" & grepl(p, period, fixed = TRUE)]
  S5 <- rbindlist(list(
    list("Panel A. Share of person-days with a stay in another central-city zone", "Pre-opening", "Post-opening", "Person-days (pre / post)"),
    list(sprintf("ONIKURU visitors (n = %s)", fint(d1$n_visitors)), fnum(gm("visitors", "pre")$mean_oc, 4), fnum(gm("visitors", "post")$mean_oc, 4), sprintf("%s / %s", fint(gm("visitors", "pre")$person_days), fint(gm("visitors", "post")$person_days))),
    list(sprintf("Non-visitors (n = %s)", fint(d1$n_non_visitors)), fnum(gm("non-visitors", "pre")$mean_oc, 4), fnum(gm("non-visitors", "post")$mean_oc, 4), sprintf("%s / %s", fint(gm("non-visitors", "pre")$person_days), fint(gm("non-visitors", "post")$person_days))),
    list("Panel B. Estimates", "Estimate", "95% CI", "p"),
    list("Visitor × post-opening", fnum(d1$est, 4), sprintf("[%s, %s]", fnum(d1$ci_low, 4), fnum(d1$ci_high, 4)), fp(d1$p)),
    if (nrow(d2)) list(sprintf("Home within %.1f km of ONIKURU × post-opening (n = %s)", CFG$near_km, fint(d2$n_visitors)), fnum(d2$est, 4), sprintf("[%s, %s]", fnum(d2$ci_low, 4), fnum(d2$ci_high, 4)), fp(d2$p)) else NULL,
    list("Joint test of pre-opening weekly coefficients", DASH, DASH, fp(d1$pretrend_p)),
    list("Panel C. Visitors' post-opening person-days", "Share", "", ""),
    list("with an ONIKURU visit", fnum(vz("share of days with an ONIKURU visit")$mean_oc), "", ""),
    list("with an ONIKURU visit and no other central-city stay", fnum(vz("no other central stay")$mean_oc), "", ""),
    list("with any central-city stay, ONIKURU included", fnum(vz("any central stay")$mean_oc), "", "")))
  setnames(S5, c("", "col_1", "col_2", "col_3")); wt(S5, "TableS5.csv")
  note("s5_did_sentence", sprintf("%s [95%% CI: %s, %s]", fnum(d1$est), fnum(d1$ci_low), fnum(d1$ci_high)), "4.4 (one sentence) / Table S4")
  note("s5_sample_person_days", fint(d1$n_person_days), "Table S5 note"); note("s5_reference_week", d1$ref_week, "Figure S2 caption")
  if (!is.null(m_es)) note("s5_last_week_observed_days", sprintf("week %d: %d day(s)", m_es[.N, rel_week], m_es[.N, n_days]), "Figure S2 caption")
}

# ---- Table S6 : town-block-level benchmark, pre and post ----
if (!is.null(m_ab)) {
  a <- m_ab[variant == "manuscript definition"]
  w <- dcast(a[, .(label, period, conn = n_connecting / n_days, rank = rank_connecting, ocs = round(out_center_share * out_strength) / out_strength)], label ~ period, value.var = c("conn", "rank", "ocs"))
  w <- w[order(rank_post)]; w[, Unit := ifelse(label == CFG$onikuru_area_label, "ONIKURU-area", roma(label))]
  S6 <- w[, .(Unit, `Connecting-chain share: Pre` = fnum(conn_pre), `Connecting-chain share: Post` = fnum(conn_post),
              `Rank: Pre` = as.character(rank_pre), `Rank: Post` = as.character(rank_post),
              `Central-city outflow share: Pre` = fnum(ocs_pre), `Central-city outflow share: Post` = fnum(ocs_post))]
  wt(S6, "TableS6.csv")
  if (!is.null(m_ac)) for (v in m_ac$variant) { x <- m_ac[variant == v]; k <- if (v == "manuscript definition") "s6" else "s6_home_excluded"
    note(paste0(k, "_area_days_pre_post"), sprintf("%s and %s", fint(x$pre_days), fint(x$post_days)), "Table S6 note")
    note(paste0(k, "_change_ci"), sprintf("%s (%s, %s)", fsgn(x$change), fsgn(x$change_ci_low), fsgn(x$change_ci_high)), "Table S6 note")
    note(paste0(k, "_change_minus_other_blocks_ci"), sprintf("%s (%s, %s)", fsgn(x$change_minus_others), fsgn(x$change_minus_others_ci_low), fsgn(x$change_minus_others_ci_high)), "Table S6 note")
    note(paste0(k, "_median_change_other_blocks"), fsgn(x$median_change_other_blocks), "4.5")
    note(paste0(k, "_rank_pre_to_post"), sprintf("%d (%s) to %d (%s)", x$pre_rank, gsub("-", "–", x$pre_rank_ci), x$post_rank, gsub("-", "–", x$post_rank_ci)), "Table S6 note")
    note(paste0(k, "_outflow_rank_pre_to_post"), sprintf("%d to %d", x$pre_outflow_rank, x$post_outflow_rank), "4.5") }
}

# ---- TableNotes ----
if (length(NOTES)) wt(rbindlist(NOTES), "TableNotes.csv")

# ---- article figures (no titles, panel letters, romanized labels) ----
th_m <- ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom",
          plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0), plot.title.position = "plot")

## Figure 6 : H2 dumbbells
if (!is.null(m_h2)) {
  d_in  <- m_h2[direction == "inflow"  & zone_label != "ALL_CENTER" & (pre_out_edges + post_out_edges) >= CFG$h2_min_inflow_edges]
  d_out <- m_h2[direction == "outflow" & zone_label != "ALL_CENTER" & (pre_hit_edges + post_hit_edges) >= CFG$h2_min_outflow_edges]
  mk6 <- function(d, ttl) { if (!nrow(d)) return(NULL); d <- copy(d); d[, sig := !is.na(fisher_p_BH) & fisher_p_BH < 0.05]
    d[, zr := paste0(roma(zone_label), ifelse(sig, " *", ""))]; setorder(d, post_share); d[, zr := factor(zr, levels = unique(zr))]
    dl <- melt(d[, .(zr, pre_share, post_share)], id.vars = "zr", variable.name = "period", value.name = "share")
    ggplot2::ggplot() + ggplot2::geom_segment(data = d, ggplot2::aes(y = zr, yend = zr, x = pre_share, xend = post_share), color = C_GREY, linewidth = .8) +
      ggplot2::geom_point(data = dl, ggplot2::aes(share, zr, color = period), size = 2.7) +
      ggplot2::scale_color_manual(values = c(pre_share = C_GREY, post_share = C_NON), labels = c("Pre-period", "Post-period"), name = NULL) +
      ggplot2::labs(title = ttl, x = NULL, y = NULL) + th_m }
  p6a <- mk6(d_in, "(a) Inflow"); p6b <- mk6(d_out, "(b) Outflow")
  if (!is.null(p6a) && !is.null(p6b) && has_pw) save_pdf_gg(patchwork::wrap_plots(p6a, p6b, nrow = 1) + patchwork::plot_layout(guides = "collect") & ggplot2::theme(legend.position = "bottom"), "Figure6", w = 10, h = 4.6)
}

## Figure 7 : within-person connectivity. (a) reach rate by day type with the two H4 reference lines; (b) ordered secondary stays
if (!is.null(m_w1) && !is.null(m_w2) && !is.null(m_bc)) {
  o <- m_w2[daytype == "ONIKURU_day"]; al <- m_w1[daytype == "non_ONIKURU_day_all"]; ce <- m_w2[daytype == "other_center_day"]
  refs2 <- c((ce$mean_n_center_nodes - (1 - ce$onward_rate)) / ce$mean_n_center_nodes, ce$person_days / al$person_days)
  d7 <- data.table(grp = factor(c("All non-ONIKURU\ndays", "Non-ONIKURU\ncentral days", "ONIKURU\nvisit-days"),
                                levels = c("All non-ONIKURU\ndays", "Non-ONIKURU\ncentral days", "ONIKURU\nvisit-days")),
                   v = c(al$onward_rate, ce$onward_rate, o$onward_rate), n = c(al$person_days, ce$person_days, o$person_days), oni = c(FALSE, FALSE, TRUE))
  p7a <- ggplot2::ggplot(d7, ggplot2::aes(grp, v, fill = oni)) + ggplot2::geom_col(width = .62) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f\n(n = %s)", v, format(n, big.mark = ",", trim = TRUE))), vjust = -.25, size = 3, lineheight = .9) +
    ggplot2::geom_hline(yintercept = refs2, linetype = c("22", "42"), color = "#444444") +
    ggplot2::annotate("text", x = 0.45, y = refs2[1], label = sprintf("Reference: node-count-weighted baseline (%.3f)", refs2[1]), hjust = 0, vjust = 1.4, size = 2.8, color = "#444444") +
    ggplot2::annotate("text", x = 0.45, y = refs2[2], label = sprintf("Reference: any central-city stay (%.3f)", refs2[2]), hjust = 0, vjust = -.5, size = 2.8, color = "#444444") +
    ggplot2::scale_fill_manual(values = c("FALSE" = C_GREY, "TRUE" = C_B), guide = "none") + ggplot2::scale_y_continuous(limits = c(0, max(0.34, 1.25 * max(c(d7$v, refs2)))), expand = c(0, 0)) +
    ggplot2::labs(title = "(a) Central-city reach rate by day type", x = NULL, y = "Central-city reach rate") + th_m
  bc <- copy(m_bc); bc[, metric := factor(metric, levels = c("same_day_secondary_rate", "secondary_2h_rate", "secondary_4h_rate"), labels = c("same-day", "≤2 h", "≤4 h"))]
  p7b <- ggplot2::ggplot(bc, ggplot2::aes(metric, estimate)) + ggplot2::geom_pointrange(ggplot2::aes(ymin = ci_low, ymax = ci_high), color = C_A, size = .7) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", estimate)), hjust = -.35, size = 3.2) + ggplot2::scale_y_continuous(limits = c(0, max(0.25, 1.3 * max(bc$ci_high)))) +
    ggplot2::labs(title = "(b) Ordered secondary stays after ONIKURU", x = NULL, y = "Secondary-stay rate") + th_m
  if (has_pw) save_pdf_gg(patchwork::wrap_plots(p7a, p7b, nrow = 1, widths = c(1.25, 1)), "Figure7", w = 9.5, h = 4.2)
}

## Figure 8 : multi-anchor benchmark (a) outflow share (b) connecting-chain share, medians + ONIKURU CI
if (!is.null(m_bv)) {
  bm <- m_bv[variant == "manuscript definition"][, .(label, is_ONIKURU, conn = n_connecting / n_days, ocs = round(out_center_share * out_strength) / out_strength)]; bm[, en := roma(label)]
  mk8 <- function(v, ttl, xlab, ci = NULL) { d <- copy(bm); d[, lab := factor(en, levels = en[order(get(v))])]
    p <- ggplot2::ggplot(d, ggplot2::aes(lab, get(v), fill = is_ONIKURU)) + ggplot2::geom_col(width = .72) + ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = c("FALSE" = C_GREY, "TRUE" = C_NON), guide = "none") + ggplot2::labs(title = ttl, x = NULL, y = xlab) + th_m
    if (isTRUE(CFG$fig8_medians)) p <- p + ggplot2::geom_hline(yintercept = median(d[[v]]), linetype = 2, color = "#444444")
    if (!is.null(ci) && isTRUE(CFG$fig8_medians)) p <- p + ggplot2::geom_errorbar(data = d[is_ONIKURU == TRUE], ggplot2::aes(x = lab, ymin = ci[1], ymax = ci[2]), width = .3, inherit.aes = FALSE)
    p }
  ci <- if (!is.null(m_bb)) { b <- m_bb[variant == "manuscript definition"]; c(b$share_ci_low, b$share_ci_high) } else NULL
  p8a <- mk8("ocs", "(a) Outflow share to center, by anchor", "Central-city outflow share"); p8b <- mk8("conn", "(b) Connecting-chain share, by anchor", "Connecting-chain share", ci)
  if (has_pw) save_pdf_gg(patchwork::wrap_plots(p8a, p8b, nrow = 1), "Figure8", w = 10, h = 5)
}

## Figure S2 : event study
if (!is.null(m_es) && nrow(m_es[!is.na(est)])) {
  pS2 <- ggplot2::ggplot(m_es[!is.na(est)], ggplot2::aes(rel_week, est)) + ggplot2::geom_hline(yintercept = 0, color = "#444444") +
    ggplot2::geom_vline(xintercept = -0.5, linetype = 2, color = "#444444") +
    ggplot2::geom_linerange(ggplot2::aes(ymin = ci_low, ymax = ci_high), color = C_CENTER, linewidth = .7) + ggplot2::geom_point(color = C_CENTER, size = 2.2) +
    ggplot2::scale_y_continuous(labels = function(x) sub("^-", MINUS, format(x))) +
    ggplot2::labs(x = "Week relative to opening (0 = opening week)", y = "Difference in share of days") + th_m
  save_pdf_gg(pS2, "FigureS2", w = 8, h = 4.2)
}

## write out the zone names that have no romanized label (Japanese would remain in figure labels)
eg_all <- rd("S6b_onikuru_ego_edgelist.csv")
nm_all <- unique(c(if (!is.null(eg_all) && "node" %in% names(eg_all)) eg_all$node, if (!is.null(m_h2)) m_h2$zone_label, if (!is.null(m_bv)) m_bv$label))
nm_base <- sub("^[^_]+_[^_]+_", "", nm_all); un <- unique(nm_base[!(nm_all %in% names(ROMA)) & !(nm_base %in% names(ROMA))])
if (length(un)) { fwrite(data.table(name = un), file.path(TAB, "unmapped_names.csv"), bom = TRUE); msg("[B-4] %d place names have no romanization (tables/unmapped_names.csv): %s", length(un), paste(head(un, 12), collapse = ", ")) }


## Figure 5 : ego transition networks (a) area pre (b) area post (c) facility post (article layout)
eg_m <- rd("S6b_onikuru_ego_edgelist.csv")
if (!is.null(eg_m) && nrow(eg_m)) {
  draw_ego_m <- function(el, title, topn = 12, nmax = 60) {
    ind2 <- el[direction == "in"][order(-N)][seq_len(min(topn, .N))]; outd <- el[direction == "out"][order(-N)][seq_len(min(topn, .N))]
    col <- function(t) ifelse(t == "center", C_CENTER, C_NON); cx <- function(n) 0.7 + 2.2 * sqrt(pmin(n, nmax) / nmax)
    yin <- if (nrow(ind2)) seq(1, -1, length.out = max(nrow(ind2), 2))[seq_len(nrow(ind2))] else numeric(0)
    yout <- if (nrow(outd)) seq(1, -1, length.out = max(nrow(outd), 2))[seq_len(nrow(outd))] else numeric(0)
    par(mar = c(0.5, 0.5, 2.2, 0.5)); plot(NA, xlim = c(-2.6, 2.6), ylim = c(-1.35, 1.45), axes = FALSE, xlab = "", ylab = "")
    title(main = title, adj = 0, font.main = 2, cex.main = 1.05)
    if (nrow(ind2)) for (i in seq_len(nrow(ind2))) segments(-1, yin[i], 0, 0, col = col(ind2$node_type[i]), lwd = 0.5 + 6 * ind2$N[i] / nmax)
    if (nrow(outd)) for (i in seq_len(nrow(outd))) segments(0, 0, 1, yout[i], col = col(outd$node_type[i]), lwd = 0.5 + 6 * outd$N[i] / nmax)
    if (nrow(ind2)) { points(rep(-1, nrow(ind2)), yin, pch = 19, cex = cx(ind2$N), col = col(ind2$node_type)); text(-1.1, yin, sprintf("%s (%d)", roma(ind2$node), ind2$N), pos = 2, cex = 0.72) }
    if (nrow(outd)) { points(rep(1, nrow(outd)), yout, pch = 19, cex = cx(outd$N), col = col(outd$node_type)); text(1.1, yout, sprintf("%s (%d)", roma(outd$node), outd$N), pos = 4, cex = 0.72) }
    points(0, 0, pch = 21, cex = 3.2, bg = "#e41a1c", col = "black"); segments(0, 0.12, 0, 0.62); text(0, 0.68, roma(el$target[1]), cex = 0.75, font = 2)
    text(-1, 1.32, "IN (origins)", font = 2, cex = 0.85); text(1, 1.32, "OUT (destinations)", font = 2, cex = 0.85)
  }
  nmax5 <- max(eg_m$N); sz5 <- pretty(c(1, nmax5), 3)[-1]
  draw_fig5 <- function() {
    layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE))
    draw_ego_m(eg_m[level == "area" & period == "pre"], "(a) ONIKURU-area in pre-period", nmax = nmax5)
    draw_ego_m(eg_m[level == "area" & period == "post"], "(b) ONIKURU-area in post-period", nmax = nmax5)
    draw_ego_m(eg_m[level == "facility" & period == "post"], "(c) ONIKURU-facility in post-period", nmax = nmax5)
    par(mar = c(0.5, 0.5, 0.5, 0.5)); plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
    legend("center", legend = c("Town-block of the central city", "Town-block outside the central city", sprintf("Transition volume: %d", sz5)),
           pch = 19, col = c(C_CENTER, C_NON, rep("#777777", length(sz5))), pt.cex = c(1.4, 1.4, 0.7 + 2.2 * sqrt(pmin(sz5, nmax5) / nmax5)),
           bty = "o", box.col = "#444444", cex = 0.9, y.intersp = 1.5)
  }
  save_base(draw_fig5, "Figure5", w = 12, h = 9)
}

## Figure 4 : data for the GIS map (post-opening transition network), plus a check plot
zc4 <- rd("S_zone_coords.csv"); tc4 <- rd("S5_transition_counts_facility.csv")
if (!is.null(zc4) && !is.null(tc4)) {
  nodes4 <- zc4[, .(zone_id, label = roma(zone_id), lon, lat, is_center, is_onikuru, n_stays = n)]
  e4 <- tc4[period == "post" & from != to & N >= 2]
  e4 <- merge(e4, zc4[, .(zone_id, lon, lat)], by.x = "from", by.y = "zone_id"); setnames(e4, c("lon", "lat"), c("lon_from", "lat_from"))
  e4 <- merge(e4, zc4[, .(zone_id, lon, lat)], by.x = "to", by.y = "zone_id");   setnames(e4, c("lon", "lat"), c("lon_to", "lat_to"))
  e4[, `:=`(from_label = roma(from), to_label = roma(to), from_onikuru = from == CFG$onikuru_label, to_onikuru = to == CFG$onikuru_label,
            class = cut(N, c(1, 2, 5, 10, 20, Inf), labels = c("2", "3-5", "6-10", "11-20", "20+"), right = TRUE))]
  setorder(e4, -N)
  fwrite(nodes4, file.path(FIG, "Figure4_nodes.csv"), bom = TRUE)
  fwrite(e4[, .(from, to, from_label, to_label, N, class, from_onikuru, to_onikuru, lon_from, lat_from, lon_to, lat_to)], file.path(FIG, "Figure4_edges.csv"), bom = TRUE)
  cat("figure data: Figure4_nodes.csv, Figure4_edges.csv (post-period transitions, N >= 2; class = legend bins)\n")
  oni4 <- zc4[zone_id == CFG$onikuru_label]
  if (nrow(oni4)) {
    dlat <- CFG$gis_window_km / 111; dlon <- CFG$gis_window_km / (111 * cos(oni4$lat * pi / 180))
    win <- zc4[abs(lat - oni4$lat) <= dlat & abs(lon - oni4$lon) <= dlon]; ew <- e4[N >= 3 & from %in% win$zone_id & to %in% win$zone_id]
    draw_f4 <- function() {
      par(mar = c(1, 1, 1, 1)); plot(win$lon, win$lat, type = "n", asp = 1 / cos(oni4$lat * pi / 180), xlab = "", ylab = "", axes = FALSE); box()
      if (nrow(ew)) segments(ew$lon_from, ew$lat_from, ew$lon_to, ew$lat_to, lwd = 0.4 + as.integer(ew$class) * 0.6,
                             col = adjustcolor(ifelse(ew$from_onikuru, "#d6336c", "#e6c200"), alpha.f = 0.55))
      points(win$lon, win$lat, pch = 19, col = ifelse(win$is_center, C_CENTER, C_NON), cex = 1.1)
      points(oni4$lon, oni4$lat, pch = 19, col = "#e41a1c", cex = 2.4); text(oni4$lon, oni4$lat, "ONIKURU", pos = 3, font = 2, cex = 0.8)
      legend("topleft", legend = c("Town block of center zone", "Town block of non-center zone", "3-5", "6-10", "11-20", "20+", "from ONIKURU"),
             pch = c(19, 19, NA, NA, NA, NA, NA), lty = c(NA, NA, 1, 1, 1, 1, 1), lwd = c(NA, NA, 1.6, 2.2, 2.8, 3.4, 2.2),
             col = c(C_CENTER, C_NON, rep("#e6c200", 4), "#d6336c"), bty = "n", cex = 0.75, title = "Number of transitions")
    }
    save_base(draw_f4, "Figure4_reference", w = 8, h = 7)
  }
} else cat("Figure 4 data skipped: S_zone_coords.csv or S5_transition_counts_facility.csv not found\n")

fwrite(data.table(manuscript = c("Figure 1", "Figure 2", "Figure 3", "Figure 4", "Figure 5", "Figure 6", "Figure 7", "Figure 8", "Figure S1", "Figure S2"),
  file = c("(photo; not from the pipeline)", "(GIS map; not from the pipeline)", "(framework diagram; not from the pipeline)",
           "Figure4_nodes.csv + Figure4_edges.csv, drawn in GIS; Figure4_reference is a check plot", "Figure5", "Figure6",
           "Figure7", "Figure8 (CFG$fig8_medians)",
           "(photos; not from the pipeline)", "FigureS2")), file.path(FIG, "figure_map.csv"), bom = TRUE)
msg("[B-4] manuscript tables (tables/Table1-3, TableS1-S6, TableNotes) and figures (figs/Figure4-8, FigureS2) written")

cat("\nDONE.\n",
    "  Intermediate results (S*.csv): ", OUT, "\n",
    "  Article tables:  ", TAB_M, "  (Table1-3, TableS1-S6, TableNotes)\n",
    "  Article figures: ", FIG_M, "  (Figure5-8, FigureS2; Figure4_nodes/edges.csv for the GIS map; figure_map.csv)\n",
    if (WRITE_DIAG) paste0("  Diagnostics:     ", file.path(OUT, "diagnostics"), "\n") else "",
    "  Table S4 of the article: place the tables/TableS4.csv column of each specification side by side.\n", sep = "")
