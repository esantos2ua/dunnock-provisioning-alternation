# Power of the observed-vs-null comparison to detect turn-taking
#
# For each scenario and each injected turn-taking probability p, simulate data
# sets with no coordination (null D applied to the observed watches), inject the
# rule of inject_turn_taking(), and run nulls D, E and F with the primary
# inference rule (nest-level cluster bootstrap). Detection = the bootstrap
# interval of the observed-vs-null odds ratio lies above one.
# Run from the repo root: Rscript analysis/04a_power.R [n_rep]  (about 1 h on 8 cores; resumes from output/results/power_parts)

suppressMessages(library(parallel))
for (f in c("alternation", "data", "nulls", "models", "power")) source(file.path("R", paste0(f, ".R")))

args <- commandArgs(trailingOnly = TRUE)
n_rep <- if (length(args)) as.integer(args[1]) else 200L
probs <- c(0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5)
d <- readRDS("data/derived/analysis_data.rds")
nests <- setNames(d$watches$nest, d$watches$watch)

RNGkind("L'Ecuyer-CMRG")
set.seed(20261002)
longs <- lapply(setNames(nm = names(d$scenarios)), function(s) make_long(d$visits, d$scenarios[[s]], s))
dir.create("output/results/power_parts", recursive = TRUE, showWarnings = FALSE)
# One scenario at a time, saved as it finishes, so an interrupted run can resume
sims <- do.call(rbind, lapply(names(longs), function(s) {
  part <- file.path("output/results/power_parts", paste0(s, "_", n_rep, ".rds"))
  if (file.exists(part)) return(readRDS(part))
  grid <- expand.grid(rep = seq_len(n_rep), p = probs)
  res <- mclapply(seq_len(nrow(grid)), function(i)
    cbind(scenario = s, rep = grid$rep[i],
          power_replicate(longs[[s]], nests, grid$p[i], n_iter = 100, B = 500)),
    mc.cores = max(1, detectCores() - 2), mc.set.seed = TRUE)
  out <- do.call(rbind, res)
  saveRDS(out, part)
  message("Scenario ", s, " done")
  out
}))

power <- aggregate(cbind(obs_rate, base_rate, null_rate, log_or,
                         detect_boot = boot_lo > 0, detect_wald = wald_lo > 0,
                         below_boot = boot_hi < 0) ~ scenario + null + p,
                   data = sims, FUN = mean)
power$pp <- 100 * (power$obs_rate - power$base_rate)
power <- power[order(power$scenario, power$null, power$p), ]

dir.create("output/results", recursive = TRUE, showWarnings = FALSE)
saveRDS(list(sims = sims, power = power, n_rep = n_rep, probs = probs), "output/results/power.rds")
print(power, digits = 3, row.names = FALSE)
