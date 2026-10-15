# Urban catalyst or enclave? A GPS trip-chain diagnostic of public-facility connectivity in a Japanese suburban city center

Replication code for:

> Kato, H. (2026). Urban catalyst or enclave? A GPS trip-chain diagnostic of
> public-facility connectivity in a Japanese suburban city center.
> *Manuscript under review.*

The R pipeline in this repository evaluates whether ONIKURU, a public
multifunctional facility opened in November 2023 in the suburban city center of
Ibaraki City (Osaka Metropolitan Area), operated as a connective or generative
urban catalyst or as a behavioral enclave. From consent-based smartphone GPS
stay points of city residents (October 2023–February 2024) it reconstructs
daily trip chains and directed stay-transition networks before and after the
opening and produces every model-based quantity in the article: area-level
centrality and the facility-level ego-network decomposition (H1), bidirectional
transition-share tests with user-level bootstrap intervals (H2), the logistic
model of generative uptake with and without a distance term (H3),
within-person fixed-effect comparisons against two reference values and the
anchor-wise gradient (H4), chain-type classification and the multi-anchor and
town-block benchmarks with the reading rule (H5), the person-day permutation
test, the difference-in-differences and phase-split checks, and the figures.

## Repository structure

```
├── onikuru_tripchain_pipeline.R   # Analysis pipeline (v2.0.0): one stay-point specification per run
├── make_synthetic_inputs.R        # Schema-compatible synthetic inputs for an end-to-end test
├── convert_xlsx_to_csv.R          # Optional: stream a large stay-point .xlsx (ArcGIS export) to CSV
├── make_figure4_layers.R          # Optional: GIS layer tables for Figure 4 from the pipeline output
├── data/
│   └── schema.csv                 # Column definitions of the three input files (no data)
├── .gitignore                     # Keeps data and outputs out of the repository
├── CITATION.cff
├── LICENSE
└── README.md
```

Running the pipeline creates `output/<spec>/` (intermediate results, article
tables and figures). Outputs and data are not tracked.

## Data availability

The input files (`data/staypoints_<spec>.csv`, `data/facility_staypoints_<spec>.csv`,
`data/residence.csv`) and the raw GPS records from which they are derived are
**not** redistributed. The GPS trajectory data are licensed from a commercial
provider under a data-use agreement that does not permit redistribution and
contain personal location information collected from consenting smartphone
users; access requires a formal application to the provider (see the Data
availability statement of the article). `data/schema.csv` documents every
column of the three files (name, type, unit, definition, source and access
conditions) so that researchers who obtain an equivalent extract can run the
pipeline without modification.

The upstream step that turns raw GPS fixes into stay points is not part of this
repository. As described in Section 3 of the article, fixes with horizontal
accuracy ≥ 50 m were dropped, stay episodes were detected with `scikit-mobility`
(`stay_locations()`; spatial radius 150 m and minimum duration 15 min in the
primary specification), each stay centroid was assigned to its town-block
(chō/chōme) polygon from the public e-Stat small-area boundary data, and the
stays whose centroid falls inside the facility building footprint were
extracted to the facility stay-point file.

`make_synthetic_inputs.R` writes schema-compatible synthetic versions of the
three files (600 simulated users) so that the pipeline can be executed
end-to-end and its outputs inspected without the licensed data. The synthetic
values bear no relation to the published results.

## Requirements

- R ≥ 4.3
- Packages: `data.table`, `lubridate`, `igraph`, `fixest`, `TraMineR`,
  `ggplot2`, `patchwork`; `readxl` only if an input file is supplied as
  `.xlsx`. Missing packages are installed from CRAN when the script is sourced.

```r
install.packages(c("data.table", "lubridate", "igraph", "fixest", "TraMineR",
                   "ggplot2", "patchwork"))
```

## How to run

1. Obtain an equivalent GPS extract, detect stay points as described above, and
   arrange the result as documented in `data/schema.csv`:
   `data/staypoints_150m15min.csv`, `data/facility_staypoints_150m15min.csv`
   and `data/residence.csv` (CSV preferred; `.xlsx` is accepted for small
   files, and `convert_xlsx_to_csv.R` streams a large ArcGIS `.xlsx` export to
   CSV).
2. Set the working directory to the repository root and run

   ```r
   source("onikuru_tripchain_pipeline.R", encoding = "UTF-8")
   ```

   The script runs from top to bottom for the specification set in `CFG`
   (`spec`, `stay_file`, `facility_file`, `out_dir`; the default is the
   150 m / 15 min primary specification) and writes to `output/150m15min/`.

3. To run another stay-point specification, either edit `CFG` or define
   `CFG_OVERRIDE` before sourcing; every element of `CFG` can be overridden
   this way:

   ```r
   CFG_OVERRIDE <- list(spec = "200m20min",
                        stay_file = "staypoints_200m20min",
                        facility_file = "facility_staypoints_200m20min",
                        out_dir = "output/200m20min")
   source("onikuru_tripchain_pipeline.R", encoding = "UTF-8")
   ```

   Each run writes its own column of the sensitivity table
   (`tables/TableS4.csv`); Table S4 of the article places the columns of the
   four specifications (150 m / 15 min, 150 m / 20 min, 200 m / 15 min,
   200 m / 20 min) side by side.

End-to-end test on synthetic data (does **not** reproduce the published values):

```r
source("make_synthetic_inputs.R", encoding = "UTF-8")
make_synthetic_inputs()                        # writes the three files to data/
source("onikuru_tripchain_pipeline.R", encoding = "UTF-8")
```

The synthetic test runs in about a minute. A run on the article's extract takes
considerably longer; the sequence analysis (Step 15), the permutation test
(Step 13) and the bootstraps dominate. `CFG$run_sequence_analysis = FALSE`
skips Step 15, which is not used by any article table.

## Output ↔ article mapping

All article elements are written to `output/<spec>/tables/` and
`output/<spec>/figs/`; the numbers quoted in table notes and in the text are
collected in `tables/TableNotes.csv` with the place where each is used.

| Output | Article element |
| --- | --- |
| `tables/Table1.csv` | Table 1 (sample characteristics) |
| `tables/Table2.csv` | Table 2 (H1: centrality and ego-network decomposition) |
| `tables/Table3.csv` | Table 3 (H4: within-person comparison) |
| `tables/TableS1.csv` | Table S1 (H2: bidirectional transition shares by town-block) |
| `tables/TableS2.csv` | Table S2 (H3: logistic model of uptake, with and without distance) |
| `tables/TableS3.csv` | Table S3 (H5: multi-anchor benchmark) |
| `tables/TableS4.csv` | Table S4 (this specification's column, home stays excluded, post-opening phases) |
| `tables/TableS5.csv` | Table S5 (difference-in-differences) |
| `tables/TableS6.csv` | Table S6 (town-block-level benchmark, pre and post) |
| `tables/TableNotes.csv` | Numbers in table notes and in Sections 4.1–4.5 |
| `figs/Figure4_nodes.csv`, `figs/Figure4_edges.csv` | Figure 4 (drawn in GIS from this export; `Figure4_reference` is a check plot, and `make_figure4_layers.R` prepares the GIS layers) |
| `figs/Figure5` | Figure 5 (ego networks) |
| `figs/Figure6` | Figure 6 (H2 dumbbell plot) |
| `figs/Figure7` | Figure 7 (H4 within-person contrast and ordered secondary stays) |
| `figs/Figure8` | Figure 8 (multi-anchor benchmark) |
| `figs/FigureS2` | Figure S2 (event study) |
| `figs/figure_map.csv` | Article figure → file map |

Figures 1–3 and S1 (facility photographs, study area, analytical framework)
are not produced by the pipeline. Town-block Fisher p-values in Table S1 are
screening statistics on non-independent transitions; the user-level cluster
bootstraps and the person-day permutation test provide the inferential
statements (Section 3 of the article).

The intermediate results of every step are also written as `S*.csv`
(unrounded), and `pipeline_results.rds` holds the main R objects of the run.
`CFG$write_diagnostics = TRUE` additionally writes the internal diagnostic
tables and figures to `output/<spec>/diagnostics/`.

## Reproducibility notes

- One run covers one stay-point specification; the sensitivity analysis is
  four runs with different `CFG` settings, not a loop inside the script.
- Every stochastic step (user-level cluster bootstraps, rank bootstrap,
  permutation test, synthetic inputs) is seeded, so results do not depend on
  the execution order.
- Fixed-effect estimates use `fixest::feols` with user-clustered standard
  errors.
- The ordered secondary-stay windows (≤ 2 h, ≤ 4 h) are measured from the
  departure time of the facility stay; the same-day rate has no window.
- Zone names are Japanese (UTF-8) and are mapped to the romanized labels of the
  article by the `ROMA` table in the script; zones without a label are listed
  in `tables/unmapped_names.csv`. Source the scripts with `encoding = "UTF-8"`
  on non-UTF-8 locales.
- The study constants (observation window, opening date, central town-blocks,
  station cores, facility labels) are set in `CFG` at the top of the script.

## Citation

If you use this code, please cite the article above and the archived version of
this repository (`CITATION.cff`).

## License

MIT — see `LICENSE`. The license applies to the code only and does not extend
to the GPS data, which remain subject to the data-use agreement with the
provider.

## Contact

Haruka Kato, Osaka Metropolitan University
