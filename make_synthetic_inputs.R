# ==============================================================================
# make_synthetic_inputs.R — schema-compatible synthetic input files
#
#   Writes the three input files that onikuru_tripchain_pipeline.R expects, so
#   that the pipeline can be executed end-to-end without the licensed GPS data:
#
#     data/staypoints_<spec>.csv           stay-point table (one row per stay episode)
#     data/facility_staypoints_<spec>.csv  the subset of stays inside the facility footprint
#     data/residence.csv                   residence mesh centroid per user
#
#   The simulated values follow the column definitions in data/schema.csv and
#   the structure the pipeline relies on (observation window, central town-block
#   names, a facility that exists only after the opening date), but they are
#   random: they bear NO relation to the published results.
#
#   Usage:
#     source("make_synthetic_inputs.R", encoding = "UTF-8")
#     make_synthetic_inputs()                       # 150m15min only
#     make_synthetic_inputs(specs = c("150m15min", "150m20min", "200m15min", "200m20min"))
#
#   The alternative specifications are derived from the same simulated episodes
#   (a 20-min minimum duration drops the shorter episodes; a 200 m radius merges
#   nothing in this simulation), which is enough to exercise the configuration
#   mechanism of the pipeline.
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

make_synthetic_inputs <- function(dir = "data", specs = "150m15min", n_users = 600, seed = 20261015) {
  set.seed(seed); dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  Sys.setenv(TZ = "Asia/Tokyo")

  ## ---- study constants (identical to CFG in the pipeline) -------------------
  pref <- "大阪府"; city <- "茨木市"
  center <- c("駅前一丁目", "駅前二丁目", "駅前三丁目", "駅前四丁目", "春日一丁目", "西駅前町", "西中条町", "岩倉町",
              "片桐町", "元町", "大手町", "本町", "宮元町", "別院町", "永代町", "双葉町")
  facility_block <- "駅前三丁目"
  open_date <- as.Date("2023-11-26"); win_from <- as.Date("2023-10-11"); win_to <- as.Date("2024-02-11")
  ny_from <- as.Date("2023-12-29"); ny_to <- as.Date("2024-01-03")

  ## ---- zones: 16 central blocks on a small grid between the two stations, other blocks around them ----
  grid <- expand.grid(x = seq(135.560, 135.577, length.out = 4), y = seq(34.812, 34.823, length.out = 4))
  zones <- data.table(NAME = center, lon = grid$x, lat = grid$y, center = TRUE, CITY_NAME = city, PREF_NAME = pref)
  other <- c("松ケ本町", "中津町", "東中条町", "新庄町", "上郡二丁目", "城の前町", "春日二丁目", "舟木町", "大池一丁目",
             "中穂積一丁目", "下穂積二丁目", "稲葉町", "若草町", "南春日丘一丁目", "総持寺一丁目", "沢良宜西一丁目",
             "美沢町", "三島丘一丁目", "五日市一丁目", "玉櫛一丁目")
  ang <- seq(0, 2 * pi, length.out = length(other) + 1)[-1]; rad <- runif(length(other), 0.012, 0.035)
  zones <- rbind(zones, data.table(NAME = other, lon = 135.5685 + rad * cos(ang) / cos(34.82 * pi / 180), lat = 34.8175 + rad * sin(ang),
                                   center = FALSE, CITY_NAME = city, PREF_NAME = pref))
  outside <- data.table(NAME = c("梅田", "芥川町", "千里中央"), lon = c(135.498, 135.618, 135.495), lat = c(34.705, 34.853, 34.808),
                        center = FALSE, CITY_NAME = c("大阪市北区", "高槻市", "豊中市"), PREF_NAME = pref)
  zones <- rbind(zones, outside); zones[, zid := .I]
  facility <- list(lon = 135.5695, lat = 34.8162)               # inside the facility block

  ## ---- users: residence mesh, central-city propensity, facility uptake --------
  users <- data.table(uid = sprintf("u%04d", seq_len(n_users)))
  home_zone <- sample(zones[CITY_NAME == city, zid], n_users, replace = TRUE, prob = ifelse(zones[CITY_NAME == city, center], 0.4, 1))
  users[, `:=`(home_zid = home_zone,
               home_mesh_lon = zones$lon[home_zone] + rnorm(n_users, 0, 0.0012),
               home_mesh_lat = zones$lat[home_zone] + rnorm(n_users, 0, 0.0010),
               propensity = rbeta(n_users, 1.2, 3),             # share of days with a central-city stay
               n_days = pmin(100L, 5L + rpois(n_users, 28)),
               has_work = runif(n_users) < 0.45)]
  users[, dist_fac_km := sqrt(((home_mesh_lon - facility$lon) * 91.3)^2 + ((home_mesh_lat - facility$lat) * 111)^2)]
  users[, p_visit := plogis(-1.6 + 2.5 * propensity - 0.35 * dist_fac_km)]
  users[, visitor := runif(.N) < p_visit]
  users[, work_zid := ifelse(has_work, sample(zones$zid, .N, replace = TRUE), NA_integer_)]
  no_residence <- sample(users$uid, round(0.08 * n_users))      # users absent from the residence file (dropped by the pipeline)

  ## ---- daily chains ---------------------------------------------------------
  all_days <- seq(win_from, win_to, by = "day"); all_days <- all_days[!(all_days >= ny_from & all_days <= ny_to)]
  fac_zid <- -1L
  acc <- list(uid = character(), start = numeric(), dur = numeric(), zid = integer(), lon = numeric(), lat = numeric(), day = character())
  push <- function(uid, start, dur, zid, lon, lat, day) {
    acc$uid[[length(acc$uid) + 1]] <<- uid; acc$start[[length(acc$start) + 1]] <<- as.numeric(start); acc$dur[[length(acc$dur) + 1]] <<- dur
    acc$zid[[length(acc$zid) + 1]] <<- zid; acc$lon[[length(acc$lon) + 1]] <<- lon; acc$lat[[length(acc$lat) + 1]] <<- lat; acc$day[[length(acc$day) + 1]] <<- day
  }
  for (i in seq_len(n_users)) {
    u <- users[i]; days <- sort(sample(all_days, u$n_days))
    for (d in as.list(days)) {
      d <- as.Date(d); post <- d >= open_date; dstr <- format(d, "%Y%m%d")
      t0 <- as.numeric(as.POSIXct(paste(d, "06:30:00"), tz = "Asia/Tokyo")) + runif(1, 0, 7200)
      ## morning at home (fixed start), then a sequence of stays with gaps of 15-90 min
      dur0 <- runif(1, 240, 420); push(u$uid, t0 - dur0 * 60 - runif(1, 900, 5400), dur0, u$home_zid, u$home_mesh_lon + rnorm(1, 0, 3e-4), u$home_mesh_lat + rnorm(1, 0, 3e-4), dstr)
      seq_z <- integer(); seq_d <- numeric()
      if (u$has_work && runif(1) < 0.75) { seq_z <- c(seq_z, u$work_zid); seq_d <- c(seq_d, runif(1, 240, 540)) }
      if (runif(1) < u$propensity) { k <- sample(1:2, 1, prob = c(0.8, 0.2)); zc <- sample(zones[center == TRUE, zid], k)
        seq_z <- c(seq_z, zc); seq_d <- c(seq_d, rlnorm(k, log(45), 0.6)) }
      if (post && u$visitor && runif(1) < 0.15) {
        seq_z <- c(seq_z, fac_zid); seq_d <- c(seq_d, rlnorm(1, log(75), 0.5))
        if (runif(1) < 0.25) { seq_z <- c(seq_z, sample(zones[center == TRUE & NAME != facility_block, zid], 1)); seq_d <- c(seq_d, rlnorm(1, log(40), 0.5)) }
      }
      if (runif(1) < 0.35) { seq_z <- c(seq_z, sample(zones$zid, 1)); seq_d <- c(seq_d, rlnorm(1, log(35), 0.6)) }
      home_back <- runif(1) < 0.85
      cur <- t0
      for (j in seq_along(seq_z)) {
        cur <- cur + runif(1, 15, 90) * 60; z <- seq_z[j]
        if (z == fac_zid) push(u$uid, cur, seq_d[j], z, facility$lon + rnorm(1, 0, 4e-5), facility$lat + rnorm(1, 0, 4e-5), dstr)
        else push(u$uid, cur, seq_d[j], z, zones$lon[z] + rnorm(1, 0, 4e-4), zones$lat[z] + rnorm(1, 0, 3.5e-4), dstr)
        cur <- cur + seq_d[j] * 60
      }
      if (home_back) push(u$uid, cur + runif(1, 15, 90) * 60, runif(1, 120, 400), u$home_zid, u$home_mesh_lon + rnorm(1, 0, 3e-4), u$home_mesh_lat + rnorm(1, 0, 3e-4), dstr)
    }
  }
  sp <- data.table(uid = unlist(acc$uid), start = unlist(acc$start), duration_min = unlist(acc$dur), zid = unlist(acc$zid),
                   lng = unlist(acc$lon), lat = unlist(acc$lat), day = unlist(acc$day))
  sp[, in_facility := zid == fac_zid]
  sp[, `:=`(PREF_NAME = ifelse(in_facility, pref, zones$PREF_NAME[pmax(zid, 1L)]), CITY_NAME = ifelse(in_facility, city, zones$CITY_NAME[pmax(zid, 1L)]),
            NAME = ifelse(in_facility, facility_block, zones$NAME[pmax(zid, 1L)]))]
  sp[, datetime := format(as.POSIXct(start, origin = "1970-01-01", tz = "Asia/Tokyo"), "%Y-%m-%d %H:%M:%S")]
  sp[, leaving_datetime := format(as.POSIXct(start + duration_min * 60, origin = "1970-01-01", tz = "Asia/Tokyo"), "%Y-%m-%d %H:%M:%S")]
  sp[, accuracy := round(rlnorm(.N, log(18), 0.55), 1)]
  sp[, src_file := sprintf("synthetic_%s.csv", day)]
  setorder(sp, uid, start); sp[, id := sprintf("sp%07d", .I)]

  ## ---- write files per specification ----------------------------------------
  for (s in specs) {
    min_dur <- if (grepl("20min", s)) 20 else 15
    x <- sp[duration_min >= min_dur]
    stay_file <- file.path(dir, sprintf("staypoints_%s.csv", s))
    fac_file  <- file.path(dir, sprintf("facility_staypoints_%s.csv", s))
    fwrite(x[, .(id, uid, datetime, leaving_datetime, lng, lat, accuracy, PREF_NAME, CITY_NAME, NAME, src_file)], stay_file, bom = TRUE)
    fwrite(x[in_facility == TRUE, .(id, uid, datetime, leaving_datetime, lng, lat)], fac_file, bom = TRUE)
    cat(sprintf("synthetic %s: %s rows, %d users -> %s (%d facility stays -> %s)\n", s, format(nrow(x), big.mark = ","), uniqueN(x$uid),
                stay_file, x[in_facility == TRUE, .N], fac_file))
  }
  res <- users[!uid %in% no_residence, .(uid, home_mesh = sprintf("mesh_%s", uid), home_mesh_lon = round(home_mesh_lon, 6), home_mesh_lat = round(home_mesh_lat, 6))]
  fwrite(res, file.path(dir, "residence.csv"), bom = TRUE)
  cat(sprintf("synthetic residence file: %d users -> %s\n", nrow(res), file.path(dir, "residence.csv")))
  invisible(list(staypoints = sp, users = users))
}

if (sys.nframe() == 0L) make_synthetic_inputs()
