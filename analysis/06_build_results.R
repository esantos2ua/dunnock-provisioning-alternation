# Single build step: every number the manuscript reports, as display strings
#
# Writes build/results.json. manuscript/paper.qmd reads it (R$...) so that no
# statistic in the text is typed by hand.
# Run from the repo root: Rscript analysis/06_build_results.R

suppressMessages(library(jsonlite))

d <- readRDS("data/derived/analysis_data.rds")
m <- readRDS("output/results/models.rds")
w <- d$watches

MINUS <- "−"
DASH <- "–"
num <- function(x, digits = 2) {
  s <- formatC(x, format = "f", digits = digits, big.mark = ",")
  sub("^-", MINUS, s)
}
int <- function(x) formatC(x, format = "d", big.mark = ",")
pct <- function(x, digits = 1) paste0(num(100 * x, digits), "%")
rng <- function(lo, hi, f = num, ...) paste0(f(lo, ...), DASH, f(hi, ...))
# Ratio intervals: a bound that would print as 1.00 gets a third decimal, so the
# text never shows "1.00" for an interval that does (or does not) exclude one
num_ratio <- function(x) if (num(x) == "1.00" && abs(x - 1) > 1e-9) num(x, 3) else num(x)
rng_ratio <- function(lo, hi) paste0(num_ratio(lo), DASH, num_ratio(hi))
pval <- function(p) if (is.na(p)) "n/a" else if (p < 0.001) "< 0.001" else paste("=", num(p, 3))
or_ci <- function(est, lo, hi) list(or = num(exp(est)), ci = rng_ratio(exp(lo), exp(hi)),
                                    excludes_one = (lo > 0 | hi < 0),
                                    above_one = lo > 0, below_one = hi < 0)

# ---- Sample ----------------------------------------------------------------------
cat_n <- function(expr) sum(with(w, eval(expr)))
sample <- list(
  visits = int(nrow(d$visits)),
  watches = int(nrow(w)),
  nests = int(length(unique(w$nest))),
  nests_mono = int(length(unique(w$nest[w$system == "monogamous"]))),
  nests_poly = int(length(unique(w$nest[w$system == "polyandrous"]))),
  seasons = int(length(unique(w$season))),
  watches_mono = int(sum(w$system == "monogamous")),
  watches_poly = int(sum(w$system == "polyandrous")),
  females = int(length(unique(w$female))),
  males = int(length(unique(c(w$male_a, w$male_b[!is.na(w$male_b)])))),
  single_visitor = int(cat_n(quote(n_birds_seen < 2))),
  poly_one_male = int(cat_n(quote(system == "polyandrous" & n_f > 0 & (n_a > 0) + (n_b > 0) == 1))),
  poly_males_only = int(cat_n(quote(system == "polyandrous" & n_f == 0 & n_birds_seen > 1))),
  duration_median_min = num(60 * median(w$length_h), 1),
  duration_range_min = rng(60 * min(w$length_h), 60 * max(w$length_h), num, 1),
  overflow_watches = int(sum(w$span_h > w$length_h + 1e-9)),
  flagged = int(sum(w$flagged)),
  flagged_roster = int(sum(w$flagged & w$watch %in% c("focal9", "focal585"))),
  flagged_order = int(sum(w$flagged & w$watch == "focal150")),
  n_iter = int(m$n_iter),
  among_pooled_pct = pct(min(m$share_among_pooled), 0),
  chick_age = paste(sort(unique(w$chick_age)), collapse = " and "),
  watches_age8 = int(sum(w$chick_age == 8)),
  watches_age9 = int(sum(w$chick_age == 9)),
  brood_range = rng(min(w$brood_size), max(w$brood_size), int),
  brood_mode = int(as.integer(names(which.max(table(w$brood_size))))),
  chicks_range = rng(min(w$n_chicks), max(w$n_chicks), int),
  same_day_watches = int(sum(as.numeric(ave(w$nest, w$date, FUN = function(x) length(unique(x)))) > 1))
)
sample$scenario_n <- lapply(d$scenarios, function(x) int(length(x)))

# ---- Observed vs nulls, per scenario ------------------------------------------
scen <- lapply(setNames(nm = names(d$scenarios)), function(s) {
  v <- m$vs_null[m$vs_null$scenario == s, ]
  o <- m$observed_pooled[m$observed_pooled$scenario == s, ]
  nulls <- lapply(setNames(nm = v$null), function(k) {
    r <- v[v$null == k, ]
    c(or_ci(r$log_or, r$lo, r$hi),
      list(boot_ci = rng_ratio(exp(r$boot_lo), exp(r$boot_hi)), boot_excludes_one = r$boot_lo > 0 | r$boot_hi < 0,
           rate = pct(r$null_mean), interval = rng(100 * r$null_lo, 100 * r$null_hi, num, 1),
           p = pval(r$p_upper), p_two = pval(r$p_two), fit_ok = r$ok,
           pp = num(100 * (r$observed - r$null_mean), 1)))
  })
  list(n = int(nrow(m$watch_tables[[s]])), rate = pct(o$rate), ci = rng(100 * o$lo, 100 * o$hi, num, 1),
       nulls = nulls)
})

# ---- Mating systems -----------------------------------------------------------------
sm <- m$system_models
system <- list(
  n_mono = int(m$n_system[["monogamous"]]), n_poly = int(m$n_system[["polyandrous"]]),
  rate_mono = pct(m$system_rate[["monogamous"]]), rate_poly = pct(m$system_rate[["polyandrous"]]),
  observed = or_ci(sm["observed", "est"], sm["observed", "lo"], sm["observed", "hi"]),
  excess = or_ci(sm["excess_within", "est"], sm["excess_within", "lo"], sm["excess_within", "hi"]),
  seen_observed = or_ci(sm["trio_seen_observed", "est"], sm["trio_seen_observed", "lo"], sm["trio_seen_observed", "hi"]),
  seen_excess = or_ci(sm["trio_seen_excess", "est"], sm["trio_seen_excess", "lo"], sm["trio_seen_excess", "hi"])
)
tc <- m$carer_class
carers <- list(
  n_mono = int(tc$n[["monogamous"]]), n_poly = int(tc$n[["polyandrous"]]),
  observed = or_ci(tc$observed[["est"]], tc$observed[["lo"]], tc$observed[["hi"]]),
  excess = or_ci(tc$excess_within[["est"]], tc$excess_within[["lo"]], tc$excess_within[["hi"]]),
  vs_within = lapply(split(tc$vs_within, tc$vs_within$system), function(r) or_ci(r$est, r$lo, r$hi))
)

# ---- Sensitivity analyses -------------------------------------------------------------
no_flag <- lapply(split(m$no_flag, m$no_flag$scenario), function(r)
  c(or_ci(r$est, r$lo, r$hi), n = int(r$n_watches)))
span_b <- lapply(split(m$span_b, m$span_b$scenario), function(r)
  c(or_ci(r$est, r$lo, r$hi), rate = pct(r$null_mean)))

# ---- Synchrony ------------------------------------------------------------------------
sync_fmt <- function(r) list(
  observed = int(r$observed), switches = int(r$switches),
  share = pct(r$observed / r$switches), null_mean = num(r$null_mean, 1),
  interval = rng(r$null_lo, r$null_hi, num, 1), ratio = num(r$observed / r$null_mean),
  p_lower = pval(r$p_lower), p_upper = pval(r$p_upper),
  above = r$p_upper < 0.05, below = r$p_lower < 0.05)
sync <- lapply(split(m$synchrony[m$synchrony$measure == "S", ], m$synchrony$scenario[m$synchrony$measure == "S"]), sync_fmt)
sync2 <- lapply(split(m$synchrony[m$synchrony$measure == "S2", ], m$synchrony$scenario[m$synchrony$measure == "S2"]), sync_fmt)

# ---- Models of Ihle et al. (2019) ---------------------------------------------
ihle <- lapply(split(m$ihle_type, m$ihle_type$scenario), function(x)
  lapply(split(x, x$null), function(r) list(rr = num(r$rr), ci = rng_ratio(r$lo, r$hi),
                                           excludes_one = r$lo > 1 | r$hi < 1, above_one = r$lo > 1)))
ihle_rr_range <- function(nulls) {
  x <- m$ihle_type[m$ihle_type$null %in% nulls, ]
  list(lo = num(min(x$rr)), hi = num(max(x$rr)), all_include_one = all(x$lo < 1 & x$hi > 1),
       all_above_one = all(x$lo > 1))
}
ihle_iv <- lapply(setNames(nm = rownames(m$ihle_iv)), function(k) {
  r <- m$ihle_iv[k, ]
  list(rr = num(r[["rr"]]), ci = rng_ratio(r[["lo"]], r[["hi"]]), excludes_one = r[["lo"]] > 1 | r[["hi"]] < 1)
})

# ---- Provisioning rates -----------------------------------------------------------------
rt <- m$rate_tab
key <- c("Female, pair" = "f_pair", "Male, pair" = "m_pair", "Female, trio" = "f_trio",
         "Alpha male, trio" = "a_trio", "Beta male, trio" = "b_trio")
rates <- setNames(lapply(seq_len(nrow(rt)), function(i)
  list(rate = num(rt$rate[i]), ci = rng(rt$lo[i], rt$hi[i]))), key[rt$group])
rc <- m$rate_contrasts
contrasts <- lapply(setNames(nm = rownames(rc)), function(k)
  list(ratio = num(rc[k, 1]), ci = rng_ratio(rc[k, 2], rc[k, 3]),
       excludes_one = rc[k, 2] > 1 | rc[k, 3] < 1))

# ---- Directed dyads in trios and regularity (CV) ----------------------------------
dg <- m$dyads_glmm
dyads <- lapply(split(dg, dg$null), function(x) lapply(split(x, x$dyad), function(r)
  list(rr = num(r$rr), ci = rng_ratio(r$lo, r$hi), excludes_one = r$lo > 1 | r$hi < 1)))
dr <- m$dyads_rand
dyads_obs <- lapply(split(dr[dr$null == "within", ], dr$dyad[dr$null == "within"]), function(r)
  list(observed = int(r$observed), expected = num(r$expected, 1), p = pval(r$p_two),
       rand_below = r$p_two < 0.05 && r$observed < r$expected,
       rand_above = r$p_two < 0.05 && r$observed > r$expected))
dyad_rr_range <- function(nulls) {
  x <- dg[dg$null %in% nulls, ]
  list(lo = num(min(x$rr)), hi = num(max(x$rr)), all_include_one = all(x$lo < 1 & x$hi > 1))
}
cvr <- m$cv_by_role
cv_role <- setNames(lapply(seq_len(nrow(cvr)), function(i)
  list(median = num(cvr$cv[i, "median"]), iqr = rng(cvr$cv[i, "q1"], cvr$cv[i, "q3"]),
       n = int(cvr$cv[i, "n"]))), paste(cvr$system, cvr$role, sep = "_"))
cv_models <- lapply(split(m$cv_models, m$cv_models$scenario), function(r)
  c(or_ci(r$est, r$lo, r$hi), n = int(r$n)))

# ---- Brood size ------------------------------------------------------------------
b <- m$brood
brood <- lapply(split(b, b$scenario), function(r) list(
  obs = or_ci(r$obs_est, r$obs_lo, r$obs_hi), excess = or_ci(r$exc_est, r$exc_lo, r$exc_hi),
  iv = list(rr = num(r$iv_rr), ci = rng_ratio(r$iv_lo, r$iv_hi),
            excludes_one = r$iv_lo > 1 | r$iv_hi < 1)))
rc_ch <- m$rate_chicks
brood_sum <- list(
  rate_1 = num(m$rate_by_chicks[["1"]], 1), rate_max = num(m$rate_by_chicks[[length(m$rate_by_chicks)]], 1),
  max_chicks = names(m$rate_by_chicks)[length(m$rate_by_chicks)],
  rate_rr = num(exp(rc_ch[["est"]])), rate_ci = rng_ratio(exp(rc_ch[["lo"]]), exp(rc_ch[["hi"]])),
  rate_above = rc_ch[["lo"]] > 0,
  cor_brood = num(m$cor_chicks_brood),
  none_excludes = !any(b$obs_lo > 0 | b$obs_hi < 0 | b$exc_lo > 0 | b$exc_hi < 0 | b$iv_lo > 1 | b$iv_hi < 1),
  or_lo = num(exp(min(c(b$obs_lo, b$exc_lo)))), or_hi = num(exp(max(c(b$obs_hi, b$exc_hi)))))

# ---- Review additions (2026-10-02): equivalence, sensitivity, synchrony, power, pseudo-pairs
# Pairs vs polyandrous groups: the OR interval as percentage points at the pairs rate
pp_at <- function(rate, log_or) 100 * (plogis(qlogis(rate) + log_or) - rate)
system$pp_ci <- rng(pp_at(m$system_rate[["monogamous"]], sm["observed", "lo"]),
                    pp_at(m$system_rate[["monogamous"]], sm["observed", "hi"]), num, 1)
system$pp_ci_lo <- num(pp_at(m$system_rate[["monogamous"]], sm["observed", "lo"]), 1)
system$pp_ci_hi <- num(pp_at(m$system_rate[["monogamous"]], sm["observed", "hi"]), 1)
system$pp_excess_ci <- rng(pp_at(m$system_rate[["monogamous"]], sm["excess_within", "lo"]),
                           pp_at(m$system_rate[["monogamous"]], sm["excess_within", "hi"]), num, 1)
no_overflow <- lapply(split(m$no_overflow, m$no_overflow$scenario), function(r)
  c(or_ci(r$est, r$lo, r$hi), n = int(r$n_watches),
    list(boot_ci = rng_ratio(exp(r$boot_lo), exp(r$boot_hi)), boot_above_one = r$boot_lo > 0,
         boot_excludes_one = r$boot_lo > 0 | r$boot_hi < 0)))
sm_s <- m$sync_models
sync_boot <- lapply(split(sm_s, sm_s$measure), function(x) lapply(split(x, x$model), function(r)
  list(ratio = num(r$est), ci = rng_ratio(r$lo, r$hi), below = r$hi < 1, above = r$lo > 1,
       excludes_one = r$lo > 1 | r$hi < 1)))
sy <- m$synchrony
sy$p_two <- pmin(1, 2 * pmin(sy$p_lower, sy$p_upper))
sy$p_holm <- p.adjust(sy$p_two, "holm")
sync_holm <- lapply(split(sy, sy$measure), function(x) lapply(split(x, x$scenario), function(r)
  list(p = pval(r$p_holm), sig = r$p_holm < 0.05)))

pw <- readRDS("output/results/power.rds")
pt <- pw$power
# Injected effect (percentage points over the uncoordinated base) at which the
# primary rule reaches 80% power, by linear interpolation
pp_for_power <- function(x, target = 0.8) {
  x <- x[order(x$pp), ]
  k <- which(x$detect_boot >= target)[1]
  if (is.na(k)) return(NA_real_)
  if (k == 1) return(x$pp[1])
  x$pp[k - 1] + (target - x$detect_boot[k - 1]) * (x$pp[k] - x$pp[k - 1]) /
    (x$detect_boot[k] - x$detect_boot[k - 1])
}
power <- lapply(split(pt, pt$scenario), function(x) lapply(split(x, x$null), function(y) list(
  pp80 = if (is.na(pp_for_power(y))) paste0("> ", num(max(y$pp), 1)) else num(pp_for_power(y), 1),
  pp80_raw = if (is.na(pp_for_power(y))) 999 else pp_for_power(y), pp_max = max(y$pp),
  type1 = pct(y$detect_boot[y$p == 0], 1), type1_raw = y$detect_boot[y$p == 0],
  base_rate = pct(y$base_rate[y$p == 0]),
  rows = lapply(seq_len(nrow(y)), function(i) list(p = y$p[i], pp = num(y$pp[i], 1),
                                                   power = pct(y$detect_boot[i], 0))))))
power_meta <- list(n_rep = int(pw$n_rep), p_max = num(max(pw$probs), 2))

ps <- readRDS("output/results/pseudo_pairs.rds")
ps_fmt <- function(tab) lapply(split(tab, tab$null), function(r)
  c(or_ci(r$log_or, r$wald_lo, r$wald_hi),
    list(boot_ci = rng_ratio(exp(r$boot_lo), exp(r$boot_hi)), boot_excludes_one = r$boot_lo > 0 | r$boot_hi < 0,
         observed = pct(r$observed), null_rate = pct(r$null_rate),
         pp = num(100 * (r$observed - r$null_rate), 1))))
pseudo <- list(n = int(ps$n_pseudo), n_real = int(ps$n_real), n_dates = int(ps$n_dates),
               n_watches = int(ps$n_watches_used), pseudo = ps_fmt(ps$pseudo), real = ps_fmt(ps$real),
               difference = c(or_ci(ps$difference[["est"]], ps$difference[["lo"]], ps$difference[["hi"]])))

R <- list(sample = sample, scen = scen, system = system, carers = carers,
          brood = brood, brood_sum = brood_sum,
          dyads = dyads, dyads_obs = dyads_obs,
          dyad_def = dyad_rr_range(c("within", "switched")),
          cv_role = cv_role, cv_models = cv_models,
          cv_all_include_zero = all(m$cv_models$lo < 0 & m$cv_models$hi > 0),
          no_flag = no_flag, no_overflow = no_overflow, sync_boot = sync_boot,
          sync_holm = sync_holm, power = power, power_meta = power_meta, pseudo = pseudo,
          span_b = span_b, sync = sync, sync2 = sync2, rates = rates,
          ihle = ihle, ihle_iv = ihle_iv,
          ihle_def = ihle_rr_range(c("within", "switched", "sorted")),
          ihle_b = ihle_rr_range("simulated"),
          rate_contrasts = contrasts,
          all_fits_ok = all(m$vs_null$ok) && all(m$fits_ok),
          n_reduced_random = int(sum(!m$vs_null$full_random)),
          dispersion_range = rng(min(m$diagnostics$pearson), max(m$diagnostics$pearson)),
          dharma_range = rng(min(m$diagnostics$dharma_ratio), max(m$diagnostics$dharma_ratio)),
          dharma_p_max = pval(max(m$diagnostics$dharma_p)),
          rate_family = m$rate_family,
          rate_delta_aic = num(abs(diff(m$rate_aic)), 1),
          re_sd_max = num(max(m$vs_null$re_sd_max), 2),
          built = format(Sys.time(), "%Y-%m-%d"))
dir.create("build", showWarnings = FALSE)
write_json(R, "build/results.json", auto_unbox = TRUE, pretty = TRUE)
message("Wrote build/results.json")
