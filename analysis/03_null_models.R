# Observed alternation and the five null datasets (B-F) for every scenario
#
# Main run: 1,000 iterations with the watch duration used in the paper.
# Sensitivity run: dataset B with the first-to-last-arrival span as the duration.
# Run from the repo root: Rscript analysis/03_null_models.R [n_iter]

source("R/alternation.R")
source("R/data.R")
source("R/nulls.R")

args <- commandArgs(trailingOnly = TRUE)
n_iter <- if (length(args)) as.integer(args[1]) else 1000L
d <- readRDS("data/derived/analysis_data.rds")
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)

set.seed(20260923)
nulls <- lapply(setNames(nm = names(d$scenarios)), function(s) {
  message("Scenario ", s, " ...")
  run_nulls(d$visits, d$watches, s, n_iter = n_iter)
})
set.seed(20260924)
nulls_span <- lapply(setNames(nm = c("mono", "poly_fm")), function(s)
  run_nulls(d$visits, d$watches, s, n_iter = n_iter, duration = "arrival_span_h"))

saveRDS(list(main = nulls, span = nulls_span), "output/results/nulls.rds")
