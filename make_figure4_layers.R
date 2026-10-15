# ==============================================================================
# make_figure4_layers.R — layer tables for drawing Figure 4 in a GIS
#
#   Input : <FIG>/Figure4_nodes.csv and Figure4_edges.csv written by the
#           pipeline (post-opening transitions, N >= 2, with the legend classes)
#   Output: <FIG>/figure4_gis/
#     nodes_center.csv   panel (a): the 16 central town-blocks + the facility node
#     edges_center.csv   panel (a): transitions with both ends in the city center (N >= 2)
#     nodes_city.csv     panel (b): zones inside the study city with at least
#                        N_MIN_NODE_B stay points (central zones always kept)
#     edges_city.csv     panel (b): transitions with both ends inside the city (N >= N_MIN_EDGE_B)
#
#   Each edge table carries columns for drawing the two directions of a pair
#   side by side (use `side` as the sign of a line offset in the GIS):
#     pair_id    common key of the two directions between the same pair of zones
#     side       1 / -1 (direction within the pair)
#     both_ways  TRUE if the reverse transition is also present
#     is_center_from / is_center_to   both ends in the city center
#     width_class 1..5 = legend bins "2", "3-5", "6-10", "11-20", "20+"
#
#   Usage:  source("make_figure4_layers.R", encoding = "UTF-8")
#           (set FIG to the figs/ folder of the specification to draw)
# ==============================================================================
library(data.table)

FIG          <- "output/150m15min/figs"     # figs/ folder written by the pipeline
OUT          <- file.path(FIG, "figure4_gis")
ZONE_PREFIX  <- "大阪府_茨木市_"            # zone_id prefix of the study city (prefecture_city_, as built by the pipeline)
FACILITY_ID  <- "ONIKURU"                   # zone_id of the facility node (CFG$onikuru_label in the pipeline)
N_MIN_EDGE_B <- 3                           # panel (b): minimum number of transitions drawn (panel (a) keeps N >= 2)
N_MIN_NODE_B <- 20                          # panel (b): minimum number of stay points for a zone to be drawn

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
nodes <- fread(file.path(FIG, "Figure4_nodes.csv"), encoding = "UTF-8")
edges <- fread(file.path(FIG, "Figure4_edges.csv"), encoding = "UTF-8")
cat(sprintf("nodes %d (center %d) | edges %d\n", nrow(nodes), sum(nodes$is_center), nrow(edges)))

## pair the two directions of each zone pair (columns for offsetting them in the GIS)
add_pair <- function(e) {
  e[, pair_id := paste(pmin(from, to), pmax(from, to), sep = " | ")]
  e[, side := ifelse(from < to, 1L, -1L)]
  e[, both_ways := .N == 2L, by = pair_id]
  e[, is_center_from := from %in% nodes[is_center == TRUE, zone_id]]
  e[, is_center_to   := to   %in% nodes[is_center == TRUE, zone_id]]
  e[, width_class := as.integer(factor(class, levels = c("2", "3-5", "6-10", "11-20", "20+")))]   # 1..5 (line-width steps)
  e
}

in_city <- function(id) startsWith(id, ZONE_PREFIX) | id == FACILITY_ID

## (a) city center
nodes_a <- nodes[is_center == TRUE]
edges_a <- add_pair(edges[from %in% nodes_a$zone_id & to %in% nodes_a$zone_id])
## (b) whole city
nodes_b <- nodes[in_city(zone_id)][n_stays >= N_MIN_NODE_B | is_center == TRUE]
edges_b <- add_pair(edges[in_city(from) & in_city(to) & N >= N_MIN_EDGE_B])
edges_b <- edges_b[from %in% nodes_b$zone_id & to %in% nodes_b$zone_id]

files <- c(nodes_a = "nodes_center.csv", edges_a = "edges_center.csv", nodes_b = "nodes_city.csv", edges_b = "edges_city.csv")
for (nm in names(files)) {
  x <- get(nm); f <- file.path(OUT, files[[nm]])
  fwrite(x, f, bom = TRUE); cat(sprintf("%-8s %5d rows -> %s\n", nm, nrow(x), f))
}
cat("\n(a) edges by class:\n"); print(edges_a[, .N, by = class][order(class)])
cat("(b) edges by class:\n"); print(edges_b[, .N, by = class][order(class)])
cat("\nfacility out-edges in (a):", edges_a[from_onikuru == TRUE, .N], " in-edges:", edges_a[to_onikuru == TRUE, .N], "\n")
