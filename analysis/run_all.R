# Run the whole pipeline in order, then render the manuscript
# Run from the repo root: Rscript analysis/run_all.R
steps <- c("02_prepare_data.R", "03_null_models.R", "04_models.R", "04a_power.R",
           "04b_pseudo_pairs.R", "05_figures.R",
           "06_build_results.R")
for (s in steps) {
  message("== ", s)
  if (system2("Rscript", file.path("analysis", s)) != 0) stop(s, " failed")
}
system2("quarto", c("render", "manuscript/paper.qmd"))
