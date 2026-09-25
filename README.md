# Urban catalyst or enclave? A sequence-based diagnostic of public-facility connectivity using GPS trip chains and ego networks

Replication code for:

> Kato, H. (2026). Urban catalyst or enclave? A sequence-based diagnostic of
> public-facility connectivity using GPS trip chains and ego networks.
> *Manuscript under review.*

This repository contains the R pipeline used to evaluate whether ONIKURU, a
public multifunctional facility opened in November 2023 in the suburban city
center of Ibaraki City (Osaka Metropolitan Area), operated as a connective or
generative urban catalyst or as a behavioral enclave. From consent-based
smartphone GPS stay points for Ibaraki City residents (190,323 stay points,
7,390 users, October 2023–February 2024), the pipeline reconstructs daily trip
chains and directed stay-transition networks before and after the opening and
reproduces every model-based quantity in the article: area-level centrality
and facility-level ego-network decomposition (H1), bidirectional
transition-share tests (H2), the logistic model of generative uptake (H3),
within-person fixed-effect comparisons and ordered secondary stays (H4),
chain-type classification and the 17-anchor benchmark (H5), the person-day
permutation test, and the four-specification sensitivity analysis.

## Repository structure

```
├── onikuru_tripchain_pipeline.R   # Full analysis pipeline (v1.0.0)
├── data/
│   └── schema.csv                 # Column definitions of the stay-point files (no raw data)
├── .gitignore                     # Prevents accidental data/output commits
├── CITATION.cff
├── LICENSE
└── README.md
```

Running the pipeline creates `output/` (tables as CSV) and `report/`
(figures as PDF + PNG, plus the Figure 4 edge export). These are not tracked
in the repository.

## Data availability

The stay-point files (`data/staypoints_<spec>.csv`) and the raw GPS records
from which they are derived are **not** redistributed. The GPS trajectory data
are licensed from Geo-Technologies Inc. under a data-use agreement that does
not permit redistribution and contain personal location information collected
from consenting smartphone users; access requires a formal application to the
provider (see the Data availability statement of the article).
`data/schema.csv` documents every column of the stay-point files (name, type,
unit, definition, source, and access conditions) so that researchers who
obtain an equivalent extract can run the pipeline without modification.

The upstream step that turns raw GPS fixes into stay points is not part of
this repository. As described in Section 3.2.2 of the article, fixes with
horizontal accuracy ≥ 50 m were dropped, stay episodes were detected with
`scikit-mobility` (`stay_locations()`, spatial radius 150 or 200 m, minimum
duration 15 or 20 min), and each stay centroid was assigned to its town-block
(chō/chōme) polygon from the public e-Stat small-area boundary data and tested
against the ONIKURU building footprint.

`make_synthetic_example()` writes schema-compatible synthetic stay-point files
(`data/synthetic_<spec>.csv`, 400 simulated users) so that the pipeline can be
executed end-to-end and its outputs inspected without the licensed data. The
synthetic values bear no relation to the published results.

## Requirements

- R ≥ 4.3 (syntax verified on R 4.3.3)
- Required packages: `data.table`, `lubridate`, `igraph`, `ggplot2`
- Optional packages: `fixest` (fixed-effect regressions; without it the
  pipeline falls back to `lm()` with `sandwich::vcovCL` cluster-robust
  standard errors, which requires `sandwich`), `patchwork` (multi-panel
  figures), `mgcv` (spline check of the H3 dose-response, not reported as a
  table)

```r
install.packages(c("data.table", "lubridate", "igraph", "ggplot2",
                   "fixest", "patchwork", "mgcv", "sandwich"))
```

## How to run

1. Obtain an equivalent GPS extract, detect stay points as described above,
   arrange the result as documented in `data/schema.csv`, and save it as
   `data/staypoints_150m15min.csv` (primary specification). The three
   sensitivity files (`150m20min`, `200m15min`, `200m20min`) are optional;
   Table 7 is built from whichever files are present.
2. `source("onikuru_tripchain_pipeline.R", encoding = "UTF-8")` — loads the
   configuration and helpers, then prints the checkpoint status.
3. `run_all(force = TRUE)` — full pipeline (cluster bootstrap B = 500,
   permutations B = 1,000). `run_all(force = FALSE)` skips steps whose outputs
   already exist, so an interrupted run can be resumed.

End-to-end test on synthetic data (does **not** reproduce the published
values):

```r
source("onikuru_tripchain_pipeline.R", encoding = "UTF-8")
make_synthetic_example()                       # writes data/synthetic_<spec>.csv
options(ONIKURU_SYNTHETIC = TRUE, ONIKURU_B_BOOT = 100, ONIKURU_B_PERM = 100)
run_all(force = TRUE)
```

Approximate runtime: about 10–20 minutes for the full run on a standard
laptop with the article's extract; the permutation test and the
four-specification sensitivity loop dominate. The synthetic test runs in
under a minute.

## Output ↔ article mapping

| Output | Article element |
| --- | --- |
| `output/table1_sample_characteristics.csv` | Table 1 |
| `output/table2_centrality_decomposition.csv` | Table 2 |
| `output/table3_h2_bidirectional_main.csv` (aggregate and the two station cores) | Table 3 |
| `output/tableS1_h2_bidirectional_all_blocks.csv` | Table S1 |
| `output/table4_h3_logistic.csv`, `output/s4_h3_predicted_probabilities.csv` | Table 4 and its note |
| `output/table5_within_person_means.csv`, `output/table5_within_person_fe.csv` | Table 5 |
| `output/s5_h4_secondary_stay_bootstrap.csv` | Section 4.4 (ordered secondary stays) |
| `output/table6_chain_types.csv` | Table 6 |
| `output/tableS2_multi_anchor.csv` | Table S2 |
| `output/s7_network_permutation_test.csv` | Section 4.2 (system-level inference) |
| `output/table7_sensitivity.csv` | Table 7 |
| `report/figure5_ego_networks` | Figure 5 |
| `report/figure6_h2_dumbbell` | Figure 6 |
| `report/figure7_within_person` | Figure 7 |
| `report/figure8_multi_anchor` | Figure 8 |
| `report/figureS1_permutation_null` | Supplementary figure |
| `report/figure4_transition_edges.csv`, `report/figure4_zone_centroids.csv` | Figure 4 (drawn in GIS from this export) |

Figures 1–3 (facility photograph, study area, analytical framework) are not
produced by the pipeline. Tables 3 and S1 report BH-adjusted Fisher p-values
as screening statistics; the person-day permutation test is the system-level
inference (Section 3.3.2).

## Reproducibility notes

- The pipeline merges the two internal scripts that produced the manuscript
  results (the trip-chain pipeline and the H3 logistic-regression script);
  H3 uses the same person-day frame as H4 and H5. Two exploratory steps not
  reported in the article (optimal-matching sequence clustering and a
  top-edge stability bootstrap) were removed. No estimation logic was changed.
- Every stochastic step (cluster bootstrap, permutation test, synthetic
  example) is seeded deterministically from its label (`seed_for()`), so any
  step can be re-run in isolation and results do not depend on execution
  order.
- The ordered secondary-stay windows (≤ 2 h, ≤ 4 h) are measured from the
  departure time of the ONIKURU stay; the same-day rate has no window.
- Fixed-effect estimates use `fixest::feols` with user-clustered standard
  errors when `fixest` is installed; the `lm()` fallback yields the same point
  estimates with CR1 cluster-robust standard errors, which may differ from
  the `fixest` defaults at the third decimal.
- Zone names are Japanese (UTF-8). Source the pipeline with
  `encoding = "UTF-8"` on non-UTF-8 locales; all constants are normalized with
  `enc2utf8()`.
- `output/sessionInfo.txt` records the exact R and package versions of each
  run.

## Citation

If you use this code, please cite the paper above and the archived version of
this repository:

> [Zenodo version DOI — inserted after the first GitHub release]

## License

MIT — see `LICENSE`. The license applies to the code only and does not extend
to the GPS data, which remain subject to the data-use agreement with
Geo-Technologies Inc.

## Contact

Haruka Kato, Associate Professor, Osaka Metropolitan University
[institutional email or researchmap URL]
