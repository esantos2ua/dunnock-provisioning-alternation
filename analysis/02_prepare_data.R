# Build the analysis-ready visit and watch tables
#
# Run from the repo root: Rscript analysis/02_prepare_data.R

source("R/alternation.R")
source("R/data.R")

visits <- read_visits("data/raw")
watches <- summarise_watches(visits)
scenarios <- scenario_watches(watches)

saveRDS(list(visits = visits, watches = watches, scenarios = scenarios),
        "data/derived/analysis_data.rds")

message("Visits: ", nrow(visits), " | watches: ", nrow(watches),
        " | scenario sizes: ", paste(names(scenarios), lengths(scenarios), collapse = ", "))
