# Provisioning alternation in dunnock pairs and trios

Data, code and manuscript source for:

> Paschoaletto, I. P., Holtmann, B., Lara, C. E., Nakagawa, S. & Santos, E. S. A.
> *Who's gonna feed the kids? Provisioning alternation in dunnock pairs and
> trios is largely explained by each parent's own visiting pattern.* Preprint,
> EcoEvoRxiv (DOI to be added).

Archived on Zenodo: version 1.0.0,
[10.5281/zenodo.24043042](https://doi.org/10.5281/zenodo.24043042); all
versions, [10.5281/zenodo.24043041](https://doi.org/10.5281/zenodo.24043041).

We test whether dunnock (*Prunella modularis*) carers in socially monogamous
pairs and socially polyandrous trios alternate their provisioning visits more
than expected from each bird's own visiting pattern, using five null models
(after Ihle et al. 2019), same-day pseudo-pairs, and a power simulation. The
data are 4,376 provisioning visits in 297 nest watches at 163 nests, filmed in
the Dunedin Botanic Garden, New Zealand, over seven breeding seasons
(2009/10–2015/16).

## Contents

```text
data/raw/       visit-level data, one CSV per social mating system (see data/README.md)
R/              functions: alternation metrics, null models, models, power, pseudo-pairs
analysis/       numbered scripts, run in order by analysis/run_all.R
tests/          testthat tests of the alternation metrics, null models and power code
build/          results.json: every number reported in the manuscript
manuscript/     paper.qmd (Quarto source of the manuscript), references, figures
```

## Reproducing the analyses

Requirements: R (≥ 4.4; we used 4.6.0) with glmmTMB, DHARMa, ggplot2,
patchwork, scales, jsonlite, knitr and testthat; Quarto (≥ 1.6) to render the
manuscript, which uses Typst for the PDF.

From the repository root:

```bash
Rscript analysis/run_all.R
```

This runs, in order:

| Script | What it does | Time |
|---|---|---|
| `02_prepare_data.R` | builds the visit and watch tables from `data/raw/` | seconds |
| `03_null_models.R` | observed alternation and null datasets B–F (1,000 iterations) | ~7 min |
| `04_models.R` | observed-vs-null models and sensitivity analyses | ~30 s |
| `04a_power.R` | power simulation (parallel; resumes from `output/results/power_parts/`) | ~40 min on 8 cores |
| `04b_pseudo_pairs.R` | same-day pseudo-pairs | minutes |
| `05_figures.R` | figures in `manuscript/figures/` | seconds |
| `06_build_results.R` | writes `build/results.json` | seconds |

and then renders `manuscript/paper.qmd`. Every number in the manuscript is read
from `build/results.json`; none is typed by hand, and a `claims` chunk stops the
render if a re-run changes a qualitative statement in the text.

Run the tests with:

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

## Licence

Code: MIT. Data and manuscript text: CC BY 4.0. See `LICENSE.md`.
