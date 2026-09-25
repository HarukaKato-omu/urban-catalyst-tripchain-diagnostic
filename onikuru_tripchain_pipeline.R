# ==============================================================================
# Urban catalyst or enclave? A sequence-based diagnostic of public-facility
# connectivity using GPS trip chains and ego networks - full pipeline (v1.0.0)
#
# Replication code for:
#   Kato, H. "Urban catalyst or enclave? A sequence-based diagnostic of
#   public-facility connectivity using GPS trip chains and ego networks"
#   (manuscript under review).
#
# NOTE ON THIS PUBLIC VERSION
#   v1.0.0 is a cleaned public release that merges two internal scripts that
#   produced the manuscript results: the trip-chain pipeline (internal v3.7,
#   2026-06-18) and the H3 logistic-regression script (H3_section9_logistic.R,
#   2026-06). Relative to the internal versions: comments and console messages
#   were translated into English; the two scripts were merged so that H3 uses
#   the same person-day frame as H4/H5; the four-specification sensitivity
#   analysis (Table 7) was turned into a loop over stay-point files; two
#   exploratory steps that are not reported in the article (optimal-matching
#   sequence clustering and a top-edge stability bootstrap) were removed; the
#   fixed-effect regressions fall back to lm() + cluster-robust standard errors
#   when fixest is not installed; and execution was reorganized into
#   checkpoint-style step functions. The statistical logic of every reported
#   quantity is unchanged.
#
# WHAT THIS SCRIPT DOES (article elements in brackets)
#   S1  Sample characteristics ......................... Table 1
#   S2  H1: area-level centrality and facility-level
#       ego-network decomposition ...................... Table 2; Figure 5
#   S3  H2: bidirectional transition shares,
#       Fisher + Benjamini-Hochberg .................... Table 3; Table S1; Figure 6
#   S4  H3: generative uptake, logistic regression ..... Table 4
#   S5  H4: within-person comparison, fixed-effect
#       regressions, ordered secondary stays ........... Table 5; Figure 7
#   S6  H5: chain-type classification and
#       multi-anchor benchmark ......................... Table 6; Table S2; Figure 8
#   S7  Pre/post network-distance permutation test ..... Section 4.2; Figure S1
#   S8  Sensitivity to stay-point definitions .......... Table 7
#   Figure 4 (geographic transition network) is drawn in GIS software from the
#   edge list exported by S2 (report/figure4_transition_edges.csv). Figures 1-3
#   are not produced by this script.
#
# INPUT
#   data/staypoints_150m15min.csv   (primary specification; not redistributed)
#   data/staypoints_150m20min.csv   (sensitivity)
#   data/staypoints_200m15min.csv   (sensitivity)
#   data/staypoints_200m20min.csv   (sensitivity)
#   Each file is a stay-point table (one row per stay episode) derived from
#   the raw GPS records licensed from Geo-Technologies Inc.; see
#   data/schema.csv for the column definitions. The upstream detection step
#   is not part of this script: raw fixes with horizontal accuracy >= 50 m
#   were dropped, stay episodes were detected with scikit-mobility
#   stay_locations() (spatial_radius_km = 0.15 or 0.20, minutes_for_a_stop =
#   15 or 20, leaving_time = TRUE), and each stay centroid was assigned to
#   its town-block (cho/chome) polygon and tested against the ONIKURU
#   building footprint (Section 3.2.2 of the article). Neither the raw GPS
#   data nor the stay-point tables can be redistributed. A schema-compatible
#   synthetic example can be generated with make_synthetic_example().
#
# USAGE
#   source("onikuru_tripchain_pipeline.R", encoding = "UTF-8")   # loads helpers, prints status
#   run_all(force = TRUE)                     # full pipeline (all specs present)
#   ## End-to-end smoke test on synthetic data (does NOT reproduce the
#   ## published values):
#   make_synthetic_example(); options(ONIKURU_SYNTHETIC = TRUE,
#     ONIKURU_B_BOOT = 100, ONIKURU_B_PERM = 100); run_all(force = TRUE)
#
# DESIGN NOTES (aligned with Section 3 of the article)
#   * Observation window 2023-10-11 to 2024-02-11; opening 2023-11-26;
#     New Year holiday 2023-12-29 to 2024-01-03 excluded.
#   * Two spatial definitions: the facility node (ONIKURU building polygon,
#     post-opening only) and the ONIKURU-area (facility + Ekimae-3), which is
#     the unit of the pre/post comparisons (H1 area level, H2, permutation).
#   * Daily trip chains are run-length-encoded zone sequences (Eq. 1); the
#     central-city reach indicator is 1(unique central-city nodes >= 2) (Eq. 2)
#     and is undirected. Ordered secondary stays (H4) count stays in another
#     central-city zone that begin after the ONIKURU stay; the 2 h and 4 h
#     windows are measured from the departure time of the ONIKURU stay.
#   * H2 town-block p-values are Fisher screening statistics on non-independent
#     transitions; system-level inference is the person-day permutation test.
#   * Within-person estimates use individual and day-of-week fixed effects with
#     user-clustered standard errors (Eq. 3); the H4 secondary-stay CIs are
#     user-level cluster bootstraps.
#   * Every stochastic step is seeded deterministically from its label
#     (seed_for()), so steps can be re-run in isolation.
# ==============================================================================

# ==============================================================================
# S0 Configuration
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(igraph)
  library(ggplot2)
})
Sys.setenv(TZ = "Asia/Tokyo")
## Zone names are Japanese (UTF-8). Source this file with encoding = "UTF-8"
## on non-UTF-8 locales; all string constants are normalized with enc2utf8().
if (!isTRUE(l10n_info()[["UTF-8"]]))
  message("Note: non-UTF-8 locale detected; source() this file with encoding = \"UTF-8\".")

DIR_DATA <- "data"
DIR_OUT  <- "output"; dir.create(DIR_OUT, showWarnings = FALSE)
DIR_REP  <- "report"; dir.create(DIR_REP, showWarnings = FALSE)

BASE_SEED <- 20260611

## Stay-point specifications (radius / minimum duration used in
## the upstream scikit-mobility step). The first entry is the primary specification.
SPECS <- list(
  `150m15min` = list(file = "staypoints_150m15min.csv", radius_m = 150, min_minutes = 15),
  `150m20min` = list(file = "staypoints_150m20min.csv", radius_m = 150, min_minutes = 20),
  `200m15min` = list(file = "staypoints_200m15min.csv", radius_m = 200, min_minutes = 15),
  `200m20min` = list(file = "staypoints_200m20min.csv", radius_m = 200, min_minutes = 20)
)
MAIN_SPEC <- "150m15min"

spec_path <- function(spec) {
  f <- SPECS[[spec]]$file
  if (isTRUE(getOption("ONIKURU_SYNTHETIC", FALSE))) f <- sub("^staypoints_", "synthetic_", f)
  file.path(DIR_DATA, f)
}

## Study period
OPEN_DATE     <- as.Date("2023-11-26")
BASELINE_FROM <- as.Date("2023-10-11"); BASELINE_TO <- as.Date("2023-11-25")
FOLLOW_FROM   <- as.Date("2023-11-26"); FOLLOW_TO   <- as.Date("2024-02-11")
NY_FROM       <- as.Date("2023-12-29"); NY_TO       <- as.Date("2024-01-03")

## Filters
ACCURACY_MAX_M <- 50          # horizontal accuracy filter (re-applied if column present)
RESIDENT_PREF  <- enc2utf8("大阪府")   # residency restriction, applied only if home_city exists
RESIDENT_CITY  <- enc2utf8("茨木市")

## Zones (Japanese names as they appear in the stay-point files)
ONIKURU_LABEL      <- "ONIKURU"
ONIKURU_AREA_LABEL <- "ONIKURU_AREA"
ONIKURU_AREA_NAME  <- enc2utf8("駅前三丁目")          # Ekimae-3, merged with the facility
JR_CORE_NAME       <- enc2utf8("西駅前町")            # Nishi-Ekimae (JR-Ibaraki Station area)
HANKYU_CORE_NAME   <- enc2utf8("永代町")              # Eitai (Hankyu-Ibaraki-City Station area)
CENTER_NAMES <- enc2utf8(c(
  "駅前一丁目", "駅前二丁目", "駅前三丁目", "駅前四丁目",
  "春日一丁目", "西駅前町", "西中条町", "岩倉町", "片桐町",
  "元町", "大手町", "本町", "宮元町", "別院町", "永代町", "双葉町"))

ROMA <- c("ONIKURU" = "ONIKURU", "ONIKURU_AREA" = "ONIKURU-area",
  "駅前一丁目" = "Ekimae-1", "駅前二丁目" = "Ekimae-2", "駅前三丁目" = "Ekimae-3",
  "駅前四丁目" = "Ekimae-4", "双葉町" = "Futaba", "松ケ本町" = "Matsugamoto",
  "永代町" = "Eitai", "元町" = "Motomachi", "新庄町" = "Shinjo",
  "東中条町" = "Higashi-Chujo", "大手町" = "Otemachi", "中津町" = "Nakatsu",
  "別院町" = "Betsuin", "本町" = "Honmachi", "舟木町" = "Funaki",
  "西中条町" = "Nishi-Chujo", "春日一丁目" = "Kasuga-1", "城の前町" = "Shironomae",
  "中穂積一丁目" = "Nakahozumi-1", "稲葉町" = "Inaba", "大池一丁目" = "Oike-1",
  "上郡二丁目" = "Kamigori-2", "千里万博公園" = "Banpaku-Park",
  "西駅前町" = "Nishi-Ekimae", "岩倉町" = "Iwakura", "片桐町" = "Katagiri",
  "宮元町" = "Miyamoto", "小川町" = "Ogawa", "新中条町" = "Shin-Chujo",
  "南春日丘一丁目" = "Minami-Kasugaoka-1", "中村町" = "Nakamura",
  "五日市一丁目" = "Itsukaichi-1", "沢良宜西一丁目" = "Sawaragi-Nishi-1",
  "紺屋町" = "Konya (Takatsuki)")
names(ROMA) <- enc2utf8(names(ROMA))
roma <- function(x) {
  x <- enc2utf8(x)
  base <- sub("^[^_]+_[^_]+_", "", x)
  key  <- ifelse(x %in% names(ROMA), x, base)
  unname(ifelse(key %in% names(ROMA), ROMA[key], key))
}

## Analysis parameters
SECONDARY_WINDOWS_H   <- c(2, 4)  # ordered secondary-stay windows (hours after leaving ONIKURU)
MIN_PRE_OBSERVED_DAYS <- 5        # H3: minimum pre-opening person-days
MIN_ANCHOR_TRAFFIC    <- 20       # H5 benchmark: minimum in+out transitions per anchor
MIN_ANCHOR_DAYS       <- 10       # H5 benchmark: minimum anchor-days for chain shares
b_boot <- function() getOption("ONIKURU_B_BOOT", 500)
b_perm <- function() getOption("ONIKURU_B_PERM", 1000)

## Deterministic per-step seeding
seed_for <- function(label) {
  v <- utf8ToInt(label)
  as.integer((BASE_SEED + sum(v * seq_along(v))) %% .Machine$integer.max)
}

# ==============================================================================
# Helpers
# ==============================================================================
msg <- function(...) message(sprintf(...))
write_out <- function(x, name, dir = DIR_OUT) {
  p <- file.path(dir, name); fwrite(x, p, bom = TRUE); msg("write: %s", p); invisible(p)
}
clean_chr <- function(x, missing = "UNKNOWN") {
  y <- trimws(as.character(x)); y[is.na(y) | y == "" | y == "NA" | y == "NaN"] <- missing; y
}
parse_dt_jst <- function(x) {
  y <- suppressWarnings(lubridate::ymd_hms(x, tz = "Asia/Tokyo", quiet = TRUE))
  if (all(is.na(y))) y <- suppressWarnings(lubridate::parse_date_time(
    x, orders = c("ymd HMS", "ymd HMS OS", "ymd HM", "ymd"), tz = "Asia/Tokyo", quiet = TRUE))
  y
}
make_zone_id <- function(pref, city, name) paste(clean_chr(pref), clean_chr(city), clean_chr(name), sep = "_")
zone_label <- function(z) roma(z)

save_gg <- function(p, stem, w = 11, h = 4.6) {
  ggsave(file.path(DIR_REP, paste0(stem, ".pdf")), p, width = w, height = h,
         device = if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf)
  ggsave(file.path(DIR_REP, paste0(stem, ".png")), p, width = w, height = h, dpi = 200)
  msg("figure: %s (pdf+png)", stem)
}
save_base <- function(draw, stem, w = 9, h = 6.5) {
  dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  dev(file.path(DIR_REP, paste0(stem, ".pdf")), width = w, height = h); draw(); dev.off()
  grDevices::png(file.path(DIR_REP, paste0(stem, ".png")), width = w, height = h, units = "in", res = 200,
                 type = if (capabilities("cairo")) "cairo" else NULL)
  draw(); dev.off()
  msg("figure: %s (pdf+png)", stem)
}
C_CENTER <- "#2c7fb8"; C_NON <- "#d95f0e"; C_A <- "#41ab5d"; C_B <- "#807dba"; C_GREY <- "#888888"
theme_set(theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        legend.position = "top"))

# ==============================================================================
# Preprocessing: stay points -> zones -> daily chains -> transitions
# ==============================================================================
prepare_staypoints <- function(spec = MAIN_SPEC) {
  path <- spec_path(spec)
  if (!file.exists(path)) stop("Stay-point file not found: ", path, " (see README.md)")
  dt <- fread(path, encoding = "UTF-8", showProgress = FALSE)
  need <- c("uid", "datetime", "leaving_datetime", "PREF_NAME", "CITY_NAME", "NAME", "ONIKURUstay")
  miss <- setdiff(need, names(dt))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))

  flow <- list(); add_flow <- function(lbl) flow[[length(flow) + 1]] <<-
    data.table(step = lbl, stay_points = nrow(dt), users = uniqueN(dt$uid))
  add_flow("Stay points in file")

  dt[, uid := clean_chr(uid, NA_character_)]; dt <- dt[!is.na(uid)]
  ## Residency restriction (the article's extract was pre-restricted; applied here if available)
  if ("home_city" %in% names(dt)) {
    dt[, home_city := enc2utf8(clean_chr(home_city))]
    dt <- dt[home_city == RESIDENT_CITY]; add_flow(sprintf("Residents of %s", RESIDENT_CITY))
  } else msg("home_city column absent: assuming the file is already restricted to %s residents", RESIDENT_CITY)

  dt[, arrival := parse_dt_jst(datetime)]; dt[, departure := parse_dt_jst(leaving_datetime)]
  dt <- dt[!is.na(arrival)]
  dt[is.na(departure) | departure < arrival, departure := arrival]
  dt[, duration_min := as.numeric(difftime(departure, arrival, units = "mins"))]
  dt[, day := as.IDate(as.Date(arrival, tz = "Asia/Tokyo"))]
  dt <- dt[day >= BASELINE_FROM & day <= FOLLOW_TO]; add_flow("In observation window")
  dt <- dt[!(day >= NY_FROM & day <= NY_TO)];        add_flow("New Year holiday excluded")
  dt[, period := fifelse(day <= BASELINE_TO, "pre", "post")]
  dt[, dow := lubridate::wday(as.Date(day), week_start = 1)]

  if ("accuracy" %in% names(dt)) {
    dt[, accuracy := suppressWarnings(as.numeric(accuracy))]
    dt <- dt[!is.na(accuracy) & accuracy < ACCURACY_MAX_M]
    add_flow(sprintf("Accuracy < %d m", ACCURACY_MAX_M))
  }
  for (cc in c("PREF_NAME", "CITY_NAME", "NAME")) dt[, (cc) := enc2utf8(clean_chr(get(cc)))]

  ## Zone assignment (two levels)
  dt[, is_onikuru_raw := clean_chr(ONIKURUstay, "") == ONIKURU_LABEL]
  n_pre_oni <- dt[is_onikuru_raw & period == "pre", .N]
  if (n_pre_oni > 0) msg("removed %d pre-opening points inside the ONIKURU polygon", n_pre_oni)
  dt <- dt[!(is_onikuru_raw & period == "pre")]
  dt[, is_onikuru := is_onikuru_raw & period == "post"]
  dt[, is_ibaraki := PREF_NAME == RESIDENT_PREF & CITY_NAME == RESIDENT_CITY]
  dt[, base_zone := make_zone_id(PREF_NAME, CITY_NAME, NAME)]
  dt[, zone_fac  := fifelse(is_onikuru, ONIKURU_LABEL, base_zone)]
  dt[, is_onikuru_area := is_ibaraki & NAME == ONIKURU_AREA_NAME]
  dt[, zone_area := fifelse(is_onikuru | is_onikuru_area, ONIKURU_AREA_LABEL, base_zone)]
  dt[, in_center := is_onikuru | (is_ibaraki & NAME %in% CENTER_NAMES)]
  add_flow("Analytic stay points")

  list(spec = spec, dt = dt, flow = rbindlist(flow))
}

build_chain <- function(d, zone_col) {
  dd <- copy(d); dd[, z := get(zone_col)]
  setorder(dd, uid, day, arrival, departure)
  dd[, run_id := rleid(z), by = .(uid, day)]
  cp <- dd[, .(period = period[1], dow = dow[1], zone = z[1],
               arrival = min(arrival), departure = max(departure),
               in_center = any(in_center), is_onikuru = any(is_onikuru)),
           by = .(uid, day, run_id)]
  setorder(cp, uid, day, arrival, departure); cp[, ord := seq_len(.N), by = .(uid, day)]; cp[]
}
make_edges <- function(cp) {
  e <- cp[, { if (.N >= 2L) .(from = zone[-.N], to = zone[-1L]) else NULL }, by = .(uid, day, period)]
  if (!nrow(e)) e <- data.table(uid = character(), day = as.IDate(character()),
                                period = character(), from = character(), to = character())
  e[]
}

prepare_context <- function(spec = MAIN_SPEC) {
  pp <- prepare_staypoints(spec); dt <- pp$dt
  chain_fac  <- build_chain(dt, "zone_fac")
  chain_area <- build_chain(dt, "zone_area")
  chain_day  <- chain_fac[, .(chain_len = .N, unique_nodes = uniqueN(zone),
                              has_onikuru = any(is_onikuru), has_center = any(in_center),
                              unique_center_nodes = uniqueN(zone[in_center]), dow = dow[1]),
                          by = .(uid, day, period)]
  post_visitors <- chain_day[period == "post" & has_onikuru, unique(uid)]
  chain_day[, visitor := uid %in% post_visitors]
  edges_fac  <- make_edges(chain_fac);  edges_area <- make_edges(chain_area)
  tc_fac  <- edges_fac[, .(N = .N), by = .(period, from, to)]
  tc_area <- edges_area[, .(N = .N), by = .(period, from, to)]
  center_map_fac  <- chain_fac[,  .(zone_center = any(in_center)), by = zone]
  center_map_area <- chain_area[, .(zone_center = any(in_center)), by = zone]
  center_zones_area <- setdiff(unique(chain_area[in_center == TRUE, zone]), ONIKURU_AREA_LABEL)
  list(spec = spec, dt = dt, flow = pp$flow, chain_fac = chain_fac, chain_area = chain_area,
       chain_day = chain_day, post_visitors = post_visitors, edges_fac = edges_fac,
       edges_area = edges_area, tc_fac = tc_fac, tc_area = tc_area,
       center_map_fac = center_map_fac, center_map_area = center_map_area,
       center_zones_area = center_zones_area)
}

# ==============================================================================
# Network primitives (H1, H5)
# ==============================================================================
calc_centrality <- function(tc, target) {
  out <- data.table(target = target, n_nodes = 0L, n_edges = 0L, betweenness = NA_real_,
                    out_degree = NA_integer_, in_degree = NA_integer_,
                    out_strength = NA_real_, in_strength = NA_real_, pagerank = NA_real_)
  if (!nrow(tc)) return(out)
  g <- graph_from_data_frame(tc[, .(from, to, weight = N)], directed = TRUE)
  out[, `:=`(n_nodes = vcount(g), n_edges = ecount(g))]
  if (!(target %in% V(g)$name)) return(out)
  E(g)$dist <- 1 / pmax(E(g)$weight, 1e-9)
  pr <- page_rank(g, directed = TRUE, weights = E(g)$weight)$vector
  out[, `:=`(
    betweenness  = as.numeric(betweenness(g, v = target, directed = TRUE, weights = E(g)$dist)),
    out_degree   = as.integer(degree(g, v = target, mode = "out")),
    in_degree    = as.integer(degree(g, v = target, mode = "in")),
    out_strength = as.numeric(strength(g, vids = target, mode = "out", weights = E(g)$weight)),
    in_strength  = as.numeric(strength(g, vids = target, mode = "in",  weights = E(g)$weight)),
    pagerank     = as.numeric(pr[target]))]
  out
}

## Ego-network decomposition: central-city outflow share, onward-flow ratio,
## flow-through ratio (Section 3.3.2)
decompose_ego <- function(tc, target, center_map) {
  out_e <- merge(tc[from == target, .(node = to,   N)], center_map, by.x = "node", by.y = "zone", all.x = TRUE)
  in_e  <- merge(tc[to   == target, .(node = from, N)], center_map, by.x = "node", by.y = "zone", all.x = TRUE)
  out_e[is.na(zone_center), zone_center := FALSE]; in_e[is.na(zone_center), zone_center := FALSE]
  ot <- sum(out_e$N); it <- sum(in_e$N)
  oc <- sum(out_e[zone_center == TRUE, N]); ic <- sum(in_e[zone_center == TRUE, N])
  data.table(target = target, out_strength = ot, in_strength = it,
             out_to_center = oc, in_from_center = ic,
             out_center_share = if (ot > 0) oc / ot else NA_real_,
             onward_to_center = if (it > 0) oc / it else NA_real_,
             flow_through     = if (it > 0) ot / it else NA_real_)
}
ego_edgelist <- function(tc, target, center_map) {
  e <- rbind(tc[from == target, .(direction = "out", node = to,   N)],
             tc[to   == target, .(direction = "in",  node = from, N)])
  if (!nrow(e)) return(data.table())
  e <- merge(e, center_map, by.x = "node", by.y = "zone", all.x = TRUE)
  e[is.na(zone_center), zone_center := FALSE]
  e[, node_type := fifelse(zone_center, "center", "noncenter")]
  setorder(e, direction, -N); e[, .(target = target, direction, node, node_label = zone_label(node), node_type, N)]
}

## Chain-type classification by the anchor's position within the central-city
## portion of the chain (Section 3.3.5). Labels follow Table 6.
classify_anchor_chain <- function(zone, in_center, is_anchor) {
  pos <- which(is_anchor); if (!length(pos)) return(NA_character_)
  oc  <- in_center & !is_anchor
  if (!any(oc)) return("anchor_only")
  before <- any(vapply(pos, function(p) p > 1L && any(oc[seq_len(p - 1L)]), logical(1)))
  after  <- any(vapply(pos, function(p) p < length(zone) && any(oc[(p + 1L):length(zone)]), logical(1)))
  if (before && after) return("bridge")
  if (after)  return("post_connection")
  if (before) return("terminal")
  "anchor_only"
}
CHAIN_TYPE_LABELS <- c(anchor_only = "ONIKURU-only chain", terminal = "ONIKURU-terminal chain",
                       post_connection = "Post-ONIKURU connection", bridge = "ONIKURU-bridge chain")
CHAIN_TYPE_CLASS  <- c(anchor_only = "Self-contained", terminal = "Self-contained",
                       post_connection = "Onward-linked", bridge = "Onward-linked")

# ==============================================================================
# S1 Sample characteristics (Table 1)
# ==============================================================================
run_sample <- function(ctx) {
  dt <- ctx$dt; cd <- ctx$chain_day
  per <- function(f, p) f(p)
  tab <- data.table(
    measure = c("Stay points retained after the accuracy filter", "Daily trip chains (person-days)",
                "Unique individuals observed", "Mean chain length", "Mean unique nodes per day",
                "Person-days with at least one central-city stay", "ONIKURU visit-days",
                "ONIKURU visitors", "ONIKURU-facility stay points", "ONIKURU-area stay points",
                "Central-city stay points"),
    pre  = c(dt[period == "pre", .N], cd[period == "pre", .N], cd[period == "pre", uniqueN(uid)],
             round(cd[period == "pre", mean(chain_len)], 3), round(cd[period == "pre", mean(unique_nodes)], 3),
             cd[period == "pre", sum(has_center)], 0, 0, NA,
             dt[period == "pre" & zone_area == ONIKURU_AREA_LABEL, .N], dt[period == "pre" & in_center, .N]),
    post = c(dt[period == "post", .N], cd[period == "post", .N], cd[period == "post", uniqueN(uid)],
             round(cd[period == "post", mean(chain_len)], 3), round(cd[period == "post", mean(unique_nodes)], 3),
             cd[period == "post", sum(has_center)], cd[period == "post", sum(has_onikuru)],
             length(ctx$post_visitors), dt[is_onikuru == TRUE, .N],
             dt[period == "post" & zone_area == ONIKURU_AREA_LABEL, .N], dt[period == "post" & in_center, .N]),
    overall = c(nrow(dt), nrow(cd), uniqueN(cd$uid), NA, NA, cd[, sum(has_center)],
                cd[, sum(has_onikuru)], length(ctx$post_visitors), dt[is_onikuru == TRUE, .N],
                dt[zone_area == ONIKURU_AREA_LABEL, .N], dt[in_center == TRUE, .N]))
  write_out(tab, "table1_sample_characteristics.csv")
  write_out(ctx$flow, "s1_sample_flow.csv")
  invisible(tab)
}

# ==============================================================================
# S2 H1: centrality and ego-network decomposition (Table 2; Figure 5; Figure 4 export)
# ==============================================================================
run_h1 <- function(ctx) {
  cen_area <- rbindlist(lapply(c("pre", "post"), function(p)
    cbind(level = "ONIKURU-area", period = p,
          calc_centrality(ctx$tc_area[period == p, .(from, to, N)], ONIKURU_AREA_LABEL))))
  cen_fac <- cbind(level = "ONIKURU-facility", period = "post",
                   calc_centrality(ctx$tc_fac[period == "post", .(from, to, N)], ONIKURU_LABEL))
  ego_area <- rbindlist(lapply(c("pre", "post"), function(p)
    cbind(level = "ONIKURU-area", period = p,
          decompose_ego(ctx$tc_area[period == p, .(from, to, N)], ONIKURU_AREA_LABEL, ctx$center_map_area))))
  ego_fac <- cbind(level = "ONIKURU-facility", period = "post",
                   decompose_ego(ctx$tc_fac[period == "post", .(from, to, N)], ONIKURU_LABEL, ctx$center_map_fac))
  cen <- rbind(cen_area, cen_fac); ego <- rbind(ego_area, ego_fac)
  tab2 <- merge(cen, ego[, .(level, period, out_center_share, onward_to_center, flow_through)],
                by = c("level", "period"))
  tab2[, `:=`(betweenness = round(betweenness, 0), pagerank = round(pagerank, 4),
              out_center_share = round(out_center_share, 3), onward_to_center = round(onward_to_center, 3),
              flow_through = round(flow_through, 3))]
  write_out(tab2, "table2_centrality_decomposition.csv")

  el <- rbind(
    cbind(level = "ONIKURU-area", period = "pre",
          ego_edgelist(ctx$tc_area[period == "pre",  .(from, to, N)], ONIKURU_AREA_LABEL, ctx$center_map_area)),
    cbind(level = "ONIKURU-area", period = "post",
          ego_edgelist(ctx$tc_area[period == "post", .(from, to, N)], ONIKURU_AREA_LABEL, ctx$center_map_area)),
    cbind(level = "ONIKURU-facility", period = "post",
          ego_edgelist(ctx$tc_fac[period == "post",  .(from, to, N)], ONIKURU_LABEL, ctx$center_map_fac)),
    fill = TRUE)
  write_out(el, "s2_ego_edgelists.csv")

  ## Figure 4 export: post-opening facility-level edges with zone centroids (if coordinates exist)
  if (all(c("lon", "lat") %in% names(ctx$dt))) {
    zc <- ctx$dt[is.finite(lon) & is.finite(lat),
                 .(lon = median(lon), lat = median(lat), is_center = any(in_center), n = .N), by = zone_fac]
    e4 <- ctx$tc_fac[period == "post"]
    e4 <- merge(e4, zc[, .(zone_fac, lon1 = lon, lat1 = lat)], by.x = "from", by.y = "zone_fac")
    e4 <- merge(e4, zc[, .(zone_fac, lon2 = lon, lat2 = lat)], by.x = "to",   by.y = "zone_fac")
    e4[, `:=`(from_label = zone_label(from), to_label = zone_label(to),
              onikuru_edge = from == ONIKURU_LABEL | to == ONIKURU_LABEL)]
    write_out(e4, "figure4_transition_edges.csv", dir = DIR_REP)
    write_out(zc[, .(zone = zone_fac, label = zone_label(zone_fac), lon, lat, is_center, n)],
              "figure4_zone_centroids.csv", dir = DIR_REP)
  }
  fig5_ego_networks(el, ego)
  invisible(tab2)
}

fig5_ego_networks <- function(el, ego, topn = 12) {
  draw_panel <- function(e, title, topn) {
    ind <- e[direction == "in"][order(-N)][seq_len(min(topn, .N))]
    out <- e[direction == "out"][order(-N)][seq_len(min(topn, .N))]
    col <- function(t) ifelse(t == "center", C_CENTER, C_NON)
    yin <- if (nrow(ind)) seq(1, -1, length.out = max(nrow(ind), 2))[seq_len(nrow(ind))] else numeric(0)
    yout <- if (nrow(out)) seq(1, -1, length.out = max(nrow(out), 2))[seq_len(nrow(out))] else numeric(0)
    plot(NA, xlim = c(-2.4, 2.4), ylim = c(-1.5, 1.45), axes = FALSE, xlab = "", ylab = "", main = title)
    if (nrow(ind)) for (i in seq_len(nrow(ind))) segments(-0.94, yin[i], -0.05, 0, col = col(ind$node_type[i]), lwd = 0.4 + ind$N[i] * 0.4)
    if (nrow(out)) for (i in seq_len(nrow(out))) segments(0.05, 0, 0.94, yout[i], col = col(out$node_type[i]), lwd = 0.4 + out$N[i] * 0.4)
    if (nrow(ind)) { points(rep(-1, nrow(ind)), yin, pch = 19, cex = 0.8 + ind$N / 8, col = col(ind$node_type))
      text(-1.08, yin, paste0(ind$node_label, " (", ind$N, ")"), pos = 2, cex = 0.65) }
    if (nrow(out)) { points(rep(1, nrow(out)), yout, pch = 19, cex = 0.8 + out$N / 8, col = col(out$node_type))
      text(1.08, yout, paste0(out$node_label, " (", out$N, ")"), pos = 4, cex = 0.65) }
    points(0, 0, pch = 19, cex = 4.2, col = "#c0392b"); text(0, 0.28, zone_label(e$target[1]), cex = 0.7, font = 2)
    text(-1, 1.3, "IN (origins)", font = 2); text(1, 1.3, "OUT (destinations)", font = 2)
  }
  draw <- function() {
    layout(matrix(c(1, 2, 3, 3), 2, 2, byrow = TRUE))
    draw_panel(el[level == "ONIKURU-area" & period == "pre"],  "(a) ONIKURU-area, pre-opening", topn)
    draw_panel(el[level == "ONIKURU-area" & period == "post"], "(b) ONIKURU-area, post-opening", topn)
    draw_panel(el[level == "ONIKURU-facility"], "(c) ONIKURU-facility, post-opening", topn)
    legend("bottom", legend = c("Town-block of city-center", "Town-block outside city-center"),
           col = c(C_CENTER, C_NON), pch = 19, horiz = TRUE, bty = "n", cex = 0.85)
    layout(1)
  }
  save_base(draw, "figure5_ego_networks", w = 13, h = 11)
}

# ==============================================================================
# S3 H2: bidirectional transition shares (Table 3; Table S1; Figure 6)
# ==============================================================================
fisher_prepost <- function(agg) {
  agg <- merge(data.table(period = c("pre", "post")), agg, by = "period", all.x = TRUE)
  agg[is.na(out_edges), out_edges := 0L]; agg[is.na(hit_edges), hit_edges := 0L]
  agg[, share := fifelse(out_edges > 0, hit_edges / out_edges, NA_real_)]
  pa <- agg[period == "pre", hit_edges]; pn <- agg[period == "pre", out_edges - hit_edges]
  qa <- agg[period == "post", hit_edges]; qn <- agg[period == "post", out_edges - hit_edges]
  fp <- NA_real_
  if ((pa + pn) > 0 && (qa + qn) > 0 && (pa + qa) > 0) {
    ft <- tryCatch(fisher.test(matrix(c(pa, pn, qa, qn), 2, byrow = TRUE)), error = function(e) NULL)
    if (!is.null(ft)) fp <- ft$p.value
  }
  data.table(pre_out_edges = agg[period == "pre", out_edges], pre_hit_edges = agg[period == "pre", hit_edges],
             pre_share = agg[period == "pre", share], post_out_edges = agg[period == "post", out_edges],
             post_hit_edges = agg[period == "post", hit_edges], post_share = agg[period == "post", share],
             delta = agg[period == "post", share] - agg[period == "pre", share], fisher_p = fp)
}
h2_table <- function(ctx) {
  ea <- ctx$edges_area; cz <- ctx$center_zones_area
  one <- function(z, direction) {
    agg <- if (direction == "inflow")
      ea[from == z, .(out_edges = .N, hit_edges = sum(to == ONIKURU_AREA_LABEL)), by = period]
    else ea[from == ONIKURU_AREA_LABEL, .(out_edges = .N, hit_edges = sum(to == z)), by = period]
    cbind(data.table(zone = z, zone_label = zone_label(z), direction = direction), fisher_prepost(agg))
  }
  blocks <- rbindlist(lapply(cz, function(z) rbind(one(z, "inflow"), one(z, "outflow"))))
  blocks[, fisher_p_BH := { v <- rep(NA_real_, .N); ok <- !is.na(fisher_p)
    if (any(ok)) v[ok] <- p.adjust(fisher_p[ok], method = "BH"); v }, by = direction]
  agg_in  <- ea[from %in% cz, .(out_edges = .N, hit_edges = sum(to == ONIKURU_AREA_LABEL)), by = period]
  agg_out <- ea[from == ONIKURU_AREA_LABEL, .(out_edges = .N, hit_edges = sum(to %in% cz)), by = period]
  aggr <- rbind(cbind(data.table(zone = "ALL_CENTER", zone_label = "Central-city zones", direction = "inflow"),  fisher_prepost(agg_in)),
                cbind(data.table(zone = "ALL_CENTER", zone_label = "Central-city zones", direction = "outflow"), fisher_prepost(agg_out)))
  aggr[, fisher_p_BH := NA_real_]
  out <- rbind(aggr, blocks, fill = TRUE); setorder(out, direction, -post_share); out[]
}
run_h2 <- function(ctx) {
  h2 <- h2_table(ctx)
  fmt <- copy(h2)[, `:=`(pre_share = round(pre_share, 3), post_share = round(post_share, 3),
                         delta = round(delta, 3), fisher_p = round(fisher_p, 3), fisher_p_BH = round(fisher_p_BH, 3))]
  write_out(fmt, "tableS1_h2_bidirectional_all_blocks.csv")
  cores <- c(make_zone_id(RESIDENT_PREF, RESIDENT_CITY, JR_CORE_NAME),
             make_zone_id(RESIDENT_PREF, RESIDENT_CITY, HANKYU_CORE_NAME))
  write_out(fmt[zone %in% c("ALL_CENTER", cores)], "table3_h2_bidirectional_main.csv")
  fig6_dumbbell(h2)
  invisible(h2)
}
fig6_dumbbell <- function(h2) {
  mk <- function(d, ttl, xlab) {
    d <- copy(d[zone != "ALL_CENTER" & !is.na(pre_share) & !is.na(post_share)])
    if (!nrow(d)) return(NULL)
    d[, sig := !is.na(fisher_p_BH) & fisher_p_BH < 0.05]
    d[, zr := paste0(zone_label, ifelse(sig, " *", ""))]; setorder(d, post_share); d[, zr := factor(zr, levels = unique(zr))]
    dl <- melt(d[, .(zr, pre_share, post_share)], id.vars = "zr", variable.name = "period", value.name = "share")
    ggplot() +
      geom_segment(data = d, aes(y = zr, yend = zr, x = pre_share, xend = post_share), color = C_GREY, linewidth = 0.8) +
      geom_point(data = dl, aes(share, zr, color = period), size = 2.7) +
      scale_color_manual(values = c(pre_share = C_GREY, post_share = C_NON), labels = c("Pre-period", "Post-period"), name = NULL) +
      labs(title = ttl, subtitle = "* BH-adjusted Fisher p < 0.05", x = xlab, y = NULL)
  }
  p_in  <- mk(h2[direction == "inflow"],  "(a) Inflow: town-block -> ONIKURU-area", "share of the town-block's outgoing transitions")
  p_out <- mk(h2[direction == "outflow"], "(b) Outflow: ONIKURU-area -> town-block", "share of the ONIKURU-area's outgoing transitions")
  ps <- Filter(Negate(is.null), list(p_in, p_out))
  if (length(ps) == 2 && requireNamespace("patchwork", quietly = TRUE))
    save_gg(patchwork::wrap_plots(ps[[1]], ps[[2]], nrow = 1), "figure6_h2_dumbbell", w = 13, h = 6.5)
  else if (length(ps)) save_gg(ps[[1]], "figure6_h2_dumbbell", w = 8, h = 6.5)
}

# ==============================================================================
# S4 H3: generative uptake, logistic regression (Table 4)
# ==============================================================================
h3_frame <- function(ctx) {
  cd <- ctx$chain_day
  pre <- cd[period == "pre", .(pre_observed_days = .N, pre_center_day_rate = mean(has_center)), by = uid]
  pre <- pre[pre_observed_days >= MIN_PRE_OBSERVED_DAYS]
  post <- cd[period == "post", .(post_observed_days = .N, became_visitor = as.integer(any(has_onikuru))), by = uid]
  merge(pre, post, by = "uid")
}
run_h3 <- function(ctx) {
  d <- h3_frame(ctx)
  msg("H3 frame: n = %d users, visitors = %d", nrow(d), sum(d$became_visitor))
  m0 <- glm(became_visitor ~ pre_center_day_rate, data = d, family = binomial())
  m1 <- glm(became_visitor ~ pre_center_day_rate + pre_observed_days + post_observed_days, data = d, family = binomial())
  tidy <- function(m, label, scale) {
    co <- summary(m)$coefficients; co <- co[rownames(co) != "(Intercept)", , drop = FALSE]
    s <- ifelse(rownames(co) == "pre_center_day_rate", scale, 1)
    data.table(model = label, term = rownames(co),
               unit = ifelse(rownames(co) == "pre_center_day_rate", sprintf("per +%.2f", scale), "per day"),
               OR = exp(s * co[, 1]), OR_low = exp(s * (co[, 1] - 1.96 * co[, 2])),
               OR_high = exp(s * (co[, 1] + 1.96 * co[, 2])), p = co[, 4])
  }
  tab4 <- rbind(tidy(m0, "crude", 0.10), tidy(m1, "adjusted", 0.10))
  tab4[, `:=`(OR = round(OR, 3), OR_low = round(OR_low, 3), OR_high = round(OR_high, 3), p = signif(p, 3))]
  write_out(tab4, "table4_h3_logistic.csv")
  nd <- data.table(pre_center_day_rate = c(0, 0.10, 0.25, 0.50),
                   pre_observed_days = mean(d$pre_observed_days), post_observed_days = mean(d$post_observed_days))
  lp <- predict(m1, newdata = nd, type = "link", se.fit = TRUE)
  nd[, `:=`(p_hat = plogis(lp$fit), p_low = plogis(lp$fit - 1.96 * lp$se.fit), p_high = plogis(lp$fit + 1.96 * lp$se.fit))]
  nd[, n_users := nrow(d)]
  write_out(nd[, .(pre_center_day_rate, p_hat = round(p_hat, 4), p_low = round(p_low, 4), p_high = round(p_high, 4), n_users)],
            "s4_h3_predicted_probabilities.csv")
  if (requireNamespace("mgcv", quietly = TRUE)) {   # dose-response check (not reported as a table)
    g <- mgcv::gam(became_visitor ~ s(pre_center_day_rate, k = 5) + pre_observed_days + post_observed_days,
                   data = d, family = binomial())
    st <- summary(g)$s.table
    write_out(data.table(term = "s(pre_center_day_rate)", edf = round(st[1, "edf"], 2), p = signif(st[1, "p-value"], 3)),
              "s4_h3_gam_spline_check.csv")
  }
  invisible(list(model = m1, table = tab4, predicted = nd))
}

# ==============================================================================
# S5 H4: within-person comparison and ordered secondary stays (Table 5; Figure 7)
# ==============================================================================
person_day_frame <- function(ctx) {
  cf <- ctx$chain_fac
  pd <- cf[, .(chain_len = .N, dow = dow[1], has_onikuru = any(is_onikuru), has_center = any(in_center),
               n_center_nodes = uniqueN(zone[in_center])), by = .(uid, day, period)]
  pd[, reach := as.integer(n_center_nodes >= 2L)]; pd[, is_oni := as.integer(has_onikuru)]
  pd
}
## Ordered secondary stays: central-city stay (not ONIKURU) that begins after the
## first ONIKURU stay of the day; windows are measured from the ONIKURU departure.
secondary_stays <- function(ctx) {
  cf <- ctx$chain_fac
  oni <- cf[is_onikuru == TRUE, .(oni_ord = min(ord), oni_departure = departure[which.min(ord)]), by = .(uid, day)]
  j <- cf[oni, on = .(uid, day), allow.cartesian = TRUE][ord > oni_ord & in_center & !is_onikuru]
  out <- oni[, .(uid, day)]
  out[, same_day := as.integer(paste(uid, day) %in% paste(j$uid, j$day))]
  for (w in SECONDARY_WINDOWS_H) {
    jw <- j[arrival <= oni_departure + lubridate::hours(w)]
    out[, (paste0("within_", w, "h")) := as.integer(paste(uid, day) %in% paste(jw$uid, jw$day))]
  }
  out[]
}
cluster_boot_rates <- function(sec, B = b_boot(), seed = seed_for("h4_boot")) {
  metrics <- setdiff(names(sec), c("uid", "day")); ids <- unique(sec$uid)
  set.seed(seed)
  res <- rbindlist(lapply(seq_len(B), function(b) {
    smp <- sample(ids, length(ids), replace = TRUE)
    w <- data.table(uid = names(table(smp)), weight = as.numeric(table(smp)))
    d <- sec[w, on = "uid", nomatch = 0]
    as.data.table(lapply(metrics, function(m) weighted.mean(d[[m]], d$weight)))[, b := b][]
  }))
  setnames(res, c(metrics, "b"))
  point <- sec[, lapply(.SD, mean), .SDcols = metrics]
  data.table(metric = metrics, estimate = as.numeric(point[1, ]),
             ci_low = sapply(metrics, function(m) quantile(res[[m]], 0.025)),
             ci_high = sapply(metrics, function(m) quantile(res[[m]], 0.975)), B = B)
}
fit_fe <- function(y, rhs, d, label) {
  if (requireNamespace("fixest", quietly = TRUE)) {
    m <- fixest::feols(as.formula(paste(y, "~", rhs, "| uid + dow")), data = d, cluster = ~uid)
    ct <- as.data.frame(fixest::coeftable(m))["is_oni", ]
    est <- ct[[1]]; se <- ct[[2]]; p <- ct[[4]]
  } else {
    m <- lm(as.formula(paste(y, "~", rhs, "+ factor(uid) + factor(dow)")), data = d)
    V <- sandwich::vcovCL(m, cluster = d$uid)
    est <- unname(coef(m)["is_oni"]); se <- sqrt(V["is_oni", "is_oni"]); p <- 2 * pnorm(-abs(est / se))
  }
  data.table(outcome = label, estimate = est, se = se, ci_low = est - 1.96 * se, ci_high = est + 1.96 * se, p = p,
             engine = if (requireNamespace("fixest", quietly = TRUE)) "fixest::feols" else "lm + sandwich::vcovCL")
}
run_h4 <- function(ctx) {
  pd <- person_day_frame(ctx)[period == "post" & uid %in% ctx$post_visitors]
  d_all <- copy(pd)[, daytype := fifelse(has_onikuru, "ONIKURU", "All non-ONIKURU")]
  d_cen <- pd[has_center == TRUE][, daytype := fifelse(has_onikuru, "ONIKURU", "Non-ONIKURU central")]
  desc <- rbind(d_all[, .(baseline = "all non-ONIKURU days", daytype, reach, n_center_nodes, chain_len)],
                d_cen[, .(baseline = "non-ONIKURU central days", daytype, reach, n_center_nodes, chain_len)])[
    , .(person_days = .N, reach_rate = mean(reach), mean_center_nodes = mean(n_center_nodes), mean_chain_len = mean(chain_len)),
    by = .(baseline, daytype)]
  fe <- rbind(fit_fe("reach", "is_oni", d_cen, "Central-city reach rate"),
              fit_fe("n_center_nodes", "is_oni", d_cen, "Central-city nodes/day"),
              fit_fe("n_center_nodes", "is_oni + chain_len", d_cen, "Central-city nodes/day (+ chain-length control)"))
  sec <- secondary_stays(ctx); boot <- cluster_boot_rates(sec)
  write_out(desc[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], "table5_within_person_means.csv")
  write_out(fe[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)], "table5_within_person_fe.csv")
  write_out(boot[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)], "s5_h4_secondary_stay_bootstrap.csv")
  fig7_within_person(desc, boot)
  invisible(list(desc = desc, fe = fe, boot = boot))
}
fig7_within_person <- function(desc, boot) {
  d <- copy(desc)[, panel := fifelse(baseline == "all non-ONIKURU days", "(a) vs all non-ONIKURU days", "(b) vs non-ONIKURU central days")]
  d[, grp := fifelse(daytype == "ONIKURU", "ONIKURU day", "benchmark")]
  p1 <- ggplot(d, aes(grp, reach_rate, fill = grp)) + geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%.3f", reach_rate)), vjust = -0.3, size = 3) + facet_wrap(~panel) +
    scale_fill_manual(values = c("ONIKURU day" = C_B, benchmark = C_GREY), guide = "none") +
    labs(x = NULL, y = "Central-city reach rate") + ylim(0, max(d$reach_rate) * 1.2)
  b <- copy(boot)[, metric := factor(metric, levels = c("same_day", paste0("within_", SECONDARY_WINDOWS_H, "h")),
                                     labels = c("same-day", paste0("<=", SECONDARY_WINDOWS_H, " h")))]
  p2 <- ggplot(b, aes(metric, estimate)) + geom_pointrange(aes(ymin = ci_low, ymax = ci_high), color = C_A, size = 0.8) +
    geom_text(aes(label = sprintf("%.3f", estimate)), hjust = -0.3, size = 3) +
    labs(title = "(c) Ordered secondary stays", x = NULL, y = "Secondary-stay rate") + ylim(0, max(b$ci_high) * 1.2)
  if (requireNamespace("patchwork", quietly = TRUE))
    save_gg(patchwork::wrap_plots(p1, p2, nrow = 1, widths = c(2, 1)), "figure7_within_person", w = 12, h = 4.6)
  else save_gg(p1, "figure7_within_person")
}

# ==============================================================================
# S6 H5: chain types and multi-anchor benchmark (Table 6; Table S2; Figure 8)
# ==============================================================================
chain_types <- function(ctx) {
  cf <- ctx$chain_fac
  ct <- cf[, .(chain_type = classify_anchor_chain(zone, in_center, is_onikuru)), by = .(uid, day, period)]
  ct <- ct[period == "post" & !is.na(chain_type), .N, by = chain_type]
  ct <- merge(data.table(chain_type = names(CHAIN_TYPE_LABELS)), ct, by = "chain_type", all.x = TRUE)
  ct[is.na(N), N := 0L]; ct[, share := N / sum(N)]
  ct[, `:=`(label = CHAIN_TYPE_LABELS[chain_type], classification = CHAIN_TYPE_CLASS[chain_type])]
  ct[match(names(CHAIN_TYPE_LABELS), chain_type)]
}
multi_anchor <- function(ctx) {
  anchors <- c(ONIKURU_LABEL, make_zone_id(RESIDENT_PREF, RESIDENT_CITY, CENTER_NAMES))
  post_e <- ctx$tc_fac[period == "post", .(from, to, N)]
  dec <- rbindlist(lapply(anchors, function(a) decompose_ego(post_e, a, ctx$center_map_fac)))
  cf <- ctx$chain_fac[period == "post"]
  ch <- rbindlist(lapply(anchors, function(a) {
    keys <- cf[, .(has_a = any(zone == a)), by = .(uid, day)][has_a == TRUE, .(uid, day)]
    if (!nrow(keys)) return(NULL)
    cc <- cf[keys, on = .(uid, day)][, .(ctype = classify_anchor_chain(zone, in_center, zone == a)), by = .(uid, day)]
    data.table(target = a, anchor_days = nrow(cc),
               connecting_share = mean(cc$ctype %in% c("post_connection", "bridge")))
  }))
  out <- merge(dec, ch, by = "target", all.x = TRUE)
  out <- out[(out_strength + in_strength) >= MIN_ANCHOR_TRAFFIC & !is.na(anchor_days) & anchor_days >= MIN_ANCHOR_DAYS]
  out[, `:=`(anchor = zone_label(target), is_ONIKURU = target == ONIKURU_LABEL)]
  setorder(out, -connecting_share)
  out[, .(anchor, is_ONIKURU, anchor_days, out_center_share, onward_to_center, flow_through, connecting_share)]
}
run_h5 <- function(ctx) {
  ct <- chain_types(ctx); ma <- multi_anchor(ctx)
  tab6 <- rbind(ct[, .(chain_type = label, n = N, share = round(share, 3), classification)],
                data.table(chain_type = c("All self-contained chains", "All onward-linked chains"),
                           n = c(ct[classification == "Self-contained", sum(N)], ct[classification == "Onward-linked", sum(N)]),
                           share = round(c(ct[classification == "Self-contained", sum(share)], ct[classification == "Onward-linked", sum(share)]), 3),
                           classification = ""))
  rank_oni <- which(ma$is_ONIKURU)
  msg("H5: ONIKURU connecting-chain share rank %s of %d anchors", ifelse(length(rank_oni), rank_oni, "NA"), nrow(ma))
  write_out(tab6, "table6_chain_types.csv")
  write_out(ma[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], "tableS2_multi_anchor.csv")
  fig8_multi_anchor(ma)
  invisible(list(chain_types = ct, multi_anchor = ma))
}
fig8_multi_anchor <- function(ma) {
  bar <- function(d, var, ttl, xlab) {
    d <- copy(d)[order(get(var))][, anchor := factor(anchor, levels = anchor)]
    ggplot(d, aes(anchor, .data[[var]], fill = is_ONIKURU)) + geom_col() + coord_flip() +
      scale_fill_manual(values = c(`FALSE` = C_GREY, `TRUE` = C_NON), guide = "none") + labs(title = ttl, x = NULL, y = xlab)
  }
  pa <- bar(ma, "connecting_share", "(a) Connecting-chain share, by anchor", "Connecting-chain share")
  pb <- bar(ma, "out_center_share", "(b) Central-city outflow share, by anchor", "Central-city outflow share")
  if (requireNamespace("patchwork", quietly = TRUE))
    save_gg(patchwork::wrap_plots(pa, pb, nrow = 1), "figure8_multi_anchor", w = 12, h = 6)
  else save_gg(pa, "figure8_multi_anchor", w = 7, h = 6)
}

# ==============================================================================
# S7 Pre/post network-distance permutation test (Section 3.3.6; Figure S1)
# ==============================================================================
permutation_test <- function(ctx, B = b_perm(), seed = seed_for("permutation")) {
  E <- copy(ctx$edges_area)
  if (!nrow(E) || uniqueN(E$period) < 2) return(list(observed = NA_real_, p = NA_real_, null = data.table()))
  states <- sort(unique(c(E$from, E$to))); ns <- length(states)
  E[, fi := match(from, states)]; E[, ti := match(to, states)]
  E[, k := fi + (ti - 1L) * ns]; E[, udid := .GRP, by = .(uid, day)]
  cellf <- factor(E$k); kc <- as.integer(cellf); nc <- nlevels(cellf)
  fi_vec <- E$fi; udid <- E$udid; cell_fi <- fi_vec[match(seq_len(nc), kc)]
  ud <- unique(E[, .(udid, period)]); setorder(ud, udid); post_ud <- as.integer(ud$period == "post")
  dist_fn <- function(pe) {
    cp <- tabulate(kc[pe], nbins = nc); cq <- tabulate(kc[!pe], nbins = nc)
    rp <- tabulate(fi_vec[pe], nbins = ns); rq <- tabulate(fi_vec[!pe], nbins = ns)
    Pp <- cp / rp[cell_fi]; Pq <- cq / rq[cell_fi]; Pp[!is.finite(Pp)] <- 0; Pq[!is.finite(Pq)] <- 0
    sqrt(sum((Pq - Pp)^2))
  }
  obs <- dist_fn(post_ud[udid] == 1L)
  set.seed(seed)
  perm <- vapply(seq_len(B), function(b) dist_fn(sample(post_ud)[udid] == 1L), numeric(1))
  list(observed = obs, p = mean(perm >= obs), null = data.table(b = seq_len(B), distance = perm), B = B)
}
run_permutation <- function(ctx) {
  pt <- permutation_test(ctx)
  write_out(data.table(observed_distance = pt$observed, permutation_p = pt$p, B = pt$B), "s7_network_permutation_test.csv")
  if (nrow(pt$null)) {
    p <- ggplot(pt$null, aes(distance)) + geom_histogram(bins = 40, fill = C_GREY, alpha = 0.7) +
      geom_vline(xintercept = pt$observed, color = C_NON, linewidth = 1) +
      annotate("text", x = pt$observed, y = Inf, vjust = 1.5, hjust = -0.05,
               label = sprintf("observed = %.2f\np = %.3f", pt$observed, pt$p), color = C_NON) +
      labs(title = "Permutation null distribution of the pre/post transition-matrix distance",
           x = "Frobenius distance (permuted period labels, person-day level)", y = "count")
    save_gg(p, "figureS1_permutation_null", w = 8)
  }
  invisible(pt)
}

# ==============================================================================
# S8 Sensitivity to stay-point definitions (Table 7)
# ==============================================================================
## One column of Table 7 for a given context; reuses the S2-S7 primitives.
core_metrics <- function(ctx) {
  ego_area_post <- decompose_ego(ctx$tc_area[period == "post", .(from, to, N)], ONIKURU_AREA_LABEL, ctx$center_map_area)
  ego_fac_post  <- decompose_ego(ctx$tc_fac[period == "post",  .(from, to, N)], ONIKURU_LABEL, ctx$center_map_fac)
  h2 <- h2_table(ctx)
  pd <- person_day_frame(ctx)[period == "post" & uid %in% ctx$post_visitors & has_center == TRUE]
  sec <- secondary_stays(ctx); boot <- cluster_boot_rates(sec, seed = seed_for(paste0("h4_boot_", ctx$spec)))
  ct <- chain_types(ctx)
  pt <- permutation_test(ctx, seed = seed_for(paste0("permutation_", ctx$spec)))
  g <- function(x, d = 3) round(x, d)
  data.table(metric = c(
    "Stay points (rows)", "ONIKURU visit person-days",
    "H1 ONIKURU-area central-city outflow share (post)", "H1 ONIKURU-facility central-city outflow share",
    "H1 ONIKURU-facility flow-through ratio",
    "H2 inflow share, pre", "H2 inflow share, post", "H2 inflow Fisher p",
    "H2 outflow share, pre", "H2 outflow share, post", "H2 outflow Fisher p",
    "H4 reach rate, ONIKURU central day", "H4 reach rate, non-ONIKURU central day",
    "H4 same-day secondary-stay rate", "H4 same-day secondary-stay 95% CI",
    "H5 ONIKURU-only-chain share", "H5 connecting-chain share", "Permutation p"),
    value = c(
      nrow(ctx$dt), ctx$chain_day[period == "post", sum(has_onikuru)],
      g(ego_area_post$out_center_share), g(ego_fac_post$out_center_share), g(ego_fac_post$flow_through),
      g(h2[zone == "ALL_CENTER" & direction == "inflow",  pre_share]), g(h2[zone == "ALL_CENTER" & direction == "inflow",  post_share]),
      signif(h2[zone == "ALL_CENTER" & direction == "inflow",  fisher_p], 3),
      g(h2[zone == "ALL_CENTER" & direction == "outflow", pre_share]), g(h2[zone == "ALL_CENTER" & direction == "outflow", post_share]),
      signif(h2[zone == "ALL_CENTER" & direction == "outflow", fisher_p], 3),
      g(pd[has_onikuru == TRUE, mean(reach)]), g(pd[has_onikuru == FALSE, mean(reach)]),
      g(boot[metric == "same_day", estimate]),
      sprintf("[%.3f, %.3f]", boot[metric == "same_day", ci_low], boot[metric == "same_day", ci_high]),
      g(ct[chain_type == "anchor_only", share]), g(ct[classification == "Onward-linked", sum(share)]),
      g(pt$p)))
}
run_sensitivity <- function(specs = names(SPECS), ctx_main = NULL) {
  cols <- list()
  for (s in specs) {
    if (!file.exists(spec_path(s))) { msg("SKIP sensitivity spec %s: %s not found", s, spec_path(s)); next }
    ctx <- if (!is.null(ctx_main) && s == ctx_main$spec) ctx_main else prepare_context(s)
    msg("S8 sensitivity: %s", s)
    cm <- core_metrics(ctx); setnames(cm, "value", s); cols[[s]] <- cm
  }
  if (!length(cols)) { msg("S8: no stay-point files found; Table 7 not written"); return(invisible(NULL)) }
  tab7 <- Reduce(function(a, b) merge(a, b, by = "metric", sort = FALSE), cols)
  write_out(tab7, "table7_sensitivity.csv")
  invisible(tab7)
}

# ==============================================================================
# Synthetic example (schema-compatible; does NOT reproduce the published values)
# ==============================================================================
make_synthetic_example <- function(n_users = 400, seed = 20260925, specs = names(SPECS), dir = DIR_DATA) {
  set.seed(seed)
  noncenter <- enc2utf8(c("松ケ本町", "新庄町", "東中条町", "中津町", "舟木町", "城の前町", "中穂積一丁目", "稲葉町",
                 "大池一丁目", "上郡二丁目", "南春日丘一丁目", "沢良宜西一丁目", "中村町", "五日市一丁目"))
  zones <- data.table(NAME = c(CENTER_NAMES, noncenter, enc2utf8(c("千里万博公園", "紺屋町"))),
                      CITY_NAME = c(rep(RESIDENT_CITY, length(CENTER_NAMES) + length(noncenter)), enc2utf8(c("吹田市", "高槻市"))),
                      PREF_NAME = RESIDENT_PREF)
  zones[, in_center := NAME %in% CENTER_NAMES]
  zones[, `:=`(lon = 135.5685 + rnorm(.N, 0, 0.006), lat = 34.8165 + rnorm(.N, 0, 0.005))]
  zones[NAME == ONIKURU_AREA_NAME, `:=`(lon = 135.5685, lat = 34.8165)]
  days <- seq(BASELINE_FROM, FOLLOW_TO, by = "day"); days <- days[!(days >= NY_FROM & days <= NY_TO)]
  users <- data.table(uid = sprintf("U%05d", seq_len(n_users)),
                      p_center = rbeta(n_users, 1.2, 4), p_obs = 0.15 + 0.45 * rbeta(n_users, 3, 3))
  zones_c <- zones[in_center == TRUE, NAME]; zones_n <- zones[in_center == FALSE, NAME]
  gen_one <- function(spec_seed, min_minutes) {
    set.seed(spec_seed)
    acc <- list(); k <- 0L
    for (i in seq_len(n_users)) {
      u <- users[i]
      obs_days <- days[runif(length(days)) < u$p_obs]
      for (d in as.character(obs_days)) {
        dd <- as.Date(d); post <- dd >= OPEN_DATE
        n_st <- 1L + rpois(1, 0.6)
        zs <- character(n_st); oni <- logical(n_st)
        in_c <- runif(1) < u$p_center
        for (s in seq_len(n_st)) {
          if (s > 1) in_c <- if (in_c) runif(1) < 0.55 else runif(1) < 0.10
          pool <- if (in_c) zones_c else zones_n; zs[s] <- pool[sample(length(pool), 1)]
        }
        if (post && runif(1) < (0.002 + 0.04 * u$p_center)) {
          j <- sample(n_st, 1); zs[j] <- ONIKURU_AREA_NAME; oni[j] <- TRUE
          if (j < n_st && runif(1) < 0.35) zs[j + 1] <- zones_c[sample(length(zones_c), 1)]
        }
        start <- as.POSIXct(paste(d, "08:00:00"), tz = "Asia/Tokyo") + runif(1, 0, 8 * 3600)
        dur <- pmax(min_minutes, round(rexp(n_st, 1 / 45) + min_minutes))
        arr <- start + cumsum(c(0, runif(n_st - 1, 40, 150) * 60 + dur[-n_st] * 60))
        dep <- arr + dur * 60
        k <- k + 1L
        acc[[k]] <- list(uid = rep(u$uid, n_st), arrival = arr, departure = dep, NAME = zs, oni = oni)
      }
    }
    out <- data.table(uid = unlist(lapply(acc, `[[`, "uid")),
                      arrival = as.POSIXct(unlist(lapply(acc, `[[`, "arrival")), origin = "1970-01-01", tz = "Asia/Tokyo"),
                      departure = as.POSIXct(unlist(lapply(acc, `[[`, "departure")), origin = "1970-01-01", tz = "Asia/Tokyo"),
                      NAME = unlist(lapply(acc, `[[`, "NAME")), oni = unlist(lapply(acc, `[[`, "oni")))
    out <- merge(out, zones[, .(NAME, PREF_NAME, CITY_NAME, lon, lat)], by = "NAME", sort = FALSE)
    setorder(out, uid, arrival)
    out[, `:=`(datetime = format(arrival, "%Y-%m-%d %H:%M:%S"), leaving_datetime = format(departure, "%Y-%m-%d %H:%M:%S"),
               lon = round(lon + rnorm(.N, 0, 0.0008), 6), lat = round(lat + rnorm(.N, 0, 0.0006), 6),
               accuracy = round(runif(.N, 5, 45), 1), ONIKURUstay = fifelse(oni, ONIKURU_LABEL, ""),
               home_city = RESIDENT_CITY)]
    out[, .(uid, datetime, leaving_datetime, lon, lat, accuracy, PREF_NAME, CITY_NAME, NAME, ONIKURUstay, home_city)]
  }
  for (s in specs) {
    out <- gen_one(seed + match(s, names(SPECS)), SPECS[[s]]$min_minutes)
    p <- file.path(dir, sub("^staypoints_", "synthetic_", SPECS[[s]]$file))
    fwrite(out, p, bom = TRUE); msg("synthetic example: %s (%d rows, %d users)", p, nrow(out), uniqueN(out$uid))
  }
  invisible(TRUE)
}

# ==============================================================================
# Execution (checkpoint style)
# ==============================================================================
pipeline_status <- function() {
  need_out <- c("table1_sample_characteristics.csv", "table2_centrality_decomposition.csv",
                "table3_h2_bidirectional_main.csv", "tableS1_h2_bidirectional_all_blocks.csv",
                "table4_h3_logistic.csv", "table5_within_person_means.csv", "table5_within_person_fe.csv",
                "table6_chain_types.csv", "tableS2_multi_anchor.csv", "s7_network_permutation_test.csv",
                "table7_sensitivity.csv")
  need_rep <- c("figure5_ego_networks.pdf", "figure6_h2_dumbbell.pdf", "figure7_within_person.pdf",
                "figure8_multi_anchor.pdf", "figureS1_permutation_null.pdf")
  f <- c(file.path(DIR_OUT, need_out), file.path(DIR_REP, need_rep))
  data.table(file = f, exists = file.exists(f))
}

run_all <- function(force = FALSE, specs = names(SPECS)) {
  st <- pipeline_status(); done <- function(x) all(st$exists[basename(st$file) %in% x])
  ctx <- NULL; get_ctx <- function() { if (is.null(ctx)) ctx <<- prepare_context(MAIN_SPEC); ctx }
  if (force || !done("table1_sample_characteristics.csv")) run_sample(get_ctx()) else message("SKIP S1 sample")
  if (force || !done(c("table2_centrality_decomposition.csv", "figure5_ego_networks.pdf"))) run_h1(get_ctx()) else message("SKIP S2 H1")
  if (force || !done(c("table3_h2_bidirectional_main.csv", "figure6_h2_dumbbell.pdf"))) run_h2(get_ctx()) else message("SKIP S3 H2")
  if (force || !done("table4_h3_logistic.csv")) run_h3(get_ctx()) else message("SKIP S4 H3")
  if (force || !done(c("table5_within_person_fe.csv", "figure7_within_person.pdf"))) run_h4(get_ctx()) else message("SKIP S5 H4")
  if (force || !done(c("table6_chain_types.csv", "figure8_multi_anchor.pdf"))) run_h5(get_ctx()) else message("SKIP S6 H5")
  if (force || !done("s7_network_permutation_test.csv")) run_permutation(get_ctx()) else message("SKIP S7 permutation")
  if (force || !done("table7_sensitivity.csv")) run_sensitivity(specs, ctx_main = get_ctx()) else message("SKIP S8 sensitivity")
  writeLines(capture.output(sessionInfo()), file.path(DIR_OUT, "sessionInfo.txt"))
  invisible(pipeline_status())
}

print(pipeline_status())
message("v1.0.0 loaded. First run: run_all(force = TRUE). ",
        "Synthetic end-to-end test (does NOT reproduce the published values): ",
        "make_synthetic_example(); options(ONIKURU_SYNTHETIC = TRUE, ONIKURU_B_BOOT = 100, ",
        "ONIKURU_B_PERM = 100); run_all(force = TRUE)")
