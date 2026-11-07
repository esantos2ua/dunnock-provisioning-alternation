# Fit the models and the sensitivity analyses
#
# Run from the repo root: Rscript analysis/04_models.R

source("R/alternation.R")
source("R/data.R")
source("R/nulls.R")
source("R/models.R")

d <- readRDS("data/derived/analysis_data.rds")
nl <- readRDS("output/results/nulls.rds")
set.seed(20260925)

wt <- lapply(nl$main, watch_table, watches = d$watches)
nulls <- names(null_labels)

# 1. Observed vs each null, per scenario --------------------------------------
vs_null <- do.call(rbind, lapply(names(wt), function(s) {
  do.call(rbind, lapply(nulls, function(k) {
    m <- compare_to_null(wt[[s]], k)
    rt <- randomization_test(nl$main[[s]], k)
    bt <- cluster_boot_or(wt[[s]], k)
    data.frame(scenario = s, null = k, log_or = m[["est"]], lo = m[["lo"]], hi = m[["hi"]],
               boot_lo = bt[["lo"]], boot_hi = bt[["hi"]],
               re_sd_max = m[["re_sd_max"]], ok = as.logical(m[["ok"]]),
               full_random = as.logical(m[["full_random"]]), observed = rt$observed,
               null_mean = rt$null_mean, null_lo = rt$null_lo, null_hi = rt$null_hi,
               p_upper = rt$p_upper, p_two = rt$p_two, n_watches = nrow(wt[[s]]))
  }))
}))
observed_pooled <- do.call(rbind, lapply(names(nl$main), function(s)
  data.frame(scenario = s, t(pooled_observed(nl$main[[s]])))))

# Dispersion and random-effect SDs of the observed-vs-within models
diagnostics <- do.call(rbind, lapply(names(wt), function(s) {
  x <- wt[[s]]; x$off <- x$off_within
  f <- fit_binomial("1 + offset(off)", x)
  data.frame(scenario = s, t(dispersion(f)), max_sd = re_sd_max(f))
}))

# 2. Monogamous pairs vs polyandrous groups (female vs male alternation) ------
pm <- rbind(wt$mono, wt$poly_fm)
fit_sys_obs <- fit_binomial("system", pm)
pm$off <- pm$off_within
fit_sys_excess <- fit_binomial("system + offset(off)", pm)
poly <- wt$poly_fm
poly$trio_seen <- as.integer(poly$n_birds_seen == 3)
poly$off <- poly$off_within
fit_seen_excess <- fit_binomial("trio_seen + offset(off)", poly)
fit_seen_obs <- fit_binomial("trio_seen", poly)
system_models <- rbind(
  observed = coef_ci(fit_sys_obs, "systempolyandrous"),
  excess_within = coef_ci(fit_sys_excess, "systempolyandrous"),
  trio_seen_observed = coef_ci(fit_seen_obs, "trio_seen"),
  trio_seen_excess = coef_ci(fit_seen_excess, "trio_seen")
)
# Model-based observed rate for each system (inverse logit of the fixed effects)
sys_rate <- plogis(c(monogamous = fixef(fit_sys_obs)$cond[[1]],
                     polyandrous = sum(fixef(fit_sys_obs)$cond)))

# 3. Sensitivity: classify watches by the number of carers that visited -----
pm_carers <- pm
pm_carers$system <- factor(pm_carers$system_carers, levels = c("monogamous", "polyandrous"))
fit_carers_obs <- fit_binomial("system", pm_carers)
fit_carers_excess <- fit_binomial("system + offset(off)", pm_carers)
carer_class <- list(
  n = table(pm_carers$system),
  observed = coef_ci(fit_carers_obs, "systempolyandrous"),
  excess_within = coef_ci(fit_carers_excess, "systempolyandrous"),
  vs_within = do.call(rbind, lapply(c("monogamous", "polyandrous"), function(g) {
    x <- pm_carers[pm_carers$system == g, ]
    data.frame(system = g, t(compare_to_null(x, "within")))
  }))
)

# 4. Sensitivity: drop flagged watches ----------------------------------------
no_flag <- do.call(rbind, lapply(names(wt), function(s) {
  x <- wt[[s]][!wt[[s]]$flagged, ]
  data.frame(scenario = s, n_watches = nrow(x), t(compare_to_null(x, "within")))
}))

# 4b. Sensitivity: drop the watches whose visits run past the recorded length
#     with the primary bootstrap and the Wald interval ----------------
overflow_ids <- d$watches$watch[d$watches$span_h > d$watches$length_h + 1e-9]
no_overflow <- do.call(rbind, lapply(names(wt), function(s) {
  x <- wt[[s]][!wt[[s]]$watch %in% overflow_ids, ]
  bt <- cluster_boot_or(x, "within")
  data.frame(scenario = s, n_watches = nrow(x), t(compare_to_null(x, "within")),
             boot_lo = bt[["lo"]], boot_hi = bt[["hi"]])
}))

# 5. Sensitivity: dataset B with the arrival span as watch duration ------
span_b <- do.call(rbind, lapply(names(nl$span), function(s) {
  x <- watch_table(nl$span[[s]], d$watches)
  rt <- randomization_test(nl$span[[s]], "simulated")
  data.frame(scenario = s, t(compare_to_null(x, "simulated")), null_mean = rt$null_mean)
}))

# 6. Synchrony, against the within-watch null: presence overlap (S) and the
#    2-min arrival window of Ihle et al. 2019 (S2) --------------------------------
synchrony <- do.call(rbind, lapply(names(nl$main), function(s) rbind(
  data.frame(scenario = s, synchrony_test(nl$main[[s]], "within", "S")),
  data.frame(scenario = s, synchrony_test(nl$main[[s]], "within", "S2")))))

# 6a. Synchrony with the primary inference rule: pooled ratio of
#     observed to expected (under D) synchrony with a nest-level cluster
#     bootstrap, per scenario and as the ratio polyandrous / pairs (female vs
#     male). Holm-adjusted p-values across the eight pooled randomization tests.
sync_watch_table <- function(res, watches, measure) {
  x <- watches[match(res$watch, watches$watch), c("watch", "nest", "system")]
  x$y <- res$observed[[measure]]
  x$e <- rowMeans(res$nulls$within[[measure]])
  x
}
boot_by_nest <- function(x, stat, B = 2000) {
  by_nest <- split(seq_len(nrow(x)), x$nest)
  est <- replicate(B, stat(x[unlist(by_nest[sample.int(length(by_nest), replace = TRUE)],
                                    use.names = FALSE), ]))
  c(est = stat(x), lo = unname(quantile(est, 0.025)), hi = unname(quantile(est, 0.975)))
}
sync_ratio <- function(x) sum(x$y) / sum(x$e)
sync_system <- function(x) sync_ratio(x[x$system == "polyandrous", ]) /
  sync_ratio(x[x$system == "monogamous", ])
sync_models <- do.call(rbind, lapply(c("S", "S2"), function(ms) {
  tabs <- lapply(nl$main, sync_watch_table, watches = d$watches, measure = ms)
  rbind(do.call(rbind, lapply(names(tabs), function(s)
          data.frame(measure = ms, model = s, t(boot_by_nest(tabs[[s]], sync_ratio))))),
        data.frame(measure = ms, model = "system",
                   t(boot_by_nest(rbind(tabs$mono, tabs$poly_fm), sync_system))))
}))
synchrony$p_two <- pmin(1, 2 * pmin(synchrony$p_lower, synchrony$p_upper))
synchrony$p_holm <- p.adjust(synchrony$p_two, "holm")

# 6b. The models of Ihle et al. (2019), for comparison with ours ---------------
stacked <- lapply(nl$main, stack_datasets, watches = d$watches)
ihle_type <- do.call(rbind, lapply(names(stacked), function(s)
  data.frame(scenario = s, fit_ihle_type(stacked[[s]]))))
pm_stack <- rbind(stack_datasets(nl$main$mono, d$watches, "within"),
                  stack_datasets(nl$main$poly_fm, d$watches, "within"))
poly_stack <- stack_datasets(nl$main$poly_fm, d$watches, "within")
poly_stack$trio_seen <- as.integer(poly_stack$n_birds_seen == 3)
ihle_iv <- rbind(system = fit_ihle_iv(pm_stack, "system"),
                 trio_seen = fit_ihle_iv(poly_stack, "trio_seen"))

# 7. Provisioning rates by role and group type --------------------------------
birds <- bird_visit_table(d$watches)
fit_rate <- fit_rates(birds)
rate_tab <- data.frame(group = levels(birds$group),
                       est = fixef(fit_rate)$cond, se = sqrt(diag(vcov(fit_rate)$cond)))
rate_tab$rate <- exp(rate_tab$est)
rate_tab$lo <- exp(rate_tab$est - 1.96 * rate_tab$se)
rate_tab$hi <- exp(rate_tab$est + 1.96 * rate_tab$se)
rate_contrast <- function(g1, g2) {
  b <- fixef(fit_rate)$cond; V <- vcov(fit_rate)$cond
  i <- match(paste0("group", c(g1, g2)), names(b))
  est <- b[i[1]] - b[i[2]]
  se <- sqrt(V[i[1], i[1]] + V[i[2], i[2]] - 2 * V[i[1], i[2]])
  exp(c(ratio = est, lo = est - 1.96 * se, hi = est + 1.96 * se))
}
rate_contrasts <- rbind(
  female_trio_vs_pair = rate_contrast("Female, trio", "Female, pair"),
  alpha_vs_pair_male = rate_contrast("Alpha male, trio", "Male, pair"),
  beta_vs_alpha = rate_contrast("Beta male, trio", "Alpha male, trio")
)

# 8. Directed dyads in trios (who follows whom), against the nulls ----------------
trio <- nl$main$poly_fab
dyads_rand <- do.call(rbind, lapply(names(null_labels), function(k) dyad_test(trio, k)))
dyads_glmm <- do.call(rbind, lapply(c("within", "switched", "sorted"), function(k)
  do.call(rbind, lapply(dyad_levels, function(dy)
    data.frame(null = k, dyad = dy, t(dyad_glmm(trio, d$watches, k, dy)))))))

# 9. Regularity of visits (CV of intervals) and excess alternation over D ----------
cv <- interval_cv(d$visits)
cv$system <- d$watches$system[match(cv$watch, d$watches$watch)]
cv_by_role <- aggregate(cv ~ system + role, data = cv[!is.na(cv$cv), ],
                        FUN = function(x) c(median = median(x), q1 = unname(quantile(x, 0.25)),
                                            q3 = unname(quantile(x, 0.75)), n = length(x)))
watch_cv <- tapply(cv$cv, cv$watch, mean, na.rm = TRUE)
cv_models <- do.call(rbind, lapply(names(wt), function(s) {
  x <- wt[[s]]
  x$cv <- watch_cv[x$watch]
  x <- x[is.finite(x$cv), ]
  x$cv_z <- as.vector(scale(x$cv))
  x$off <- x$off_within
  f <- fit_binomial("cv_z + offset(off)", x)
  data.frame(scenario = s, n = nrow(x), t(coef_ci(f, "cv_z")), ok = fit_ok(f))
}))

# 10. Brood size (chicks present at the watch): workload vs alternation ----------
#     Ihle et al. (2019) showed that brood size produces spurious effects on raw
#     alternation through provisioning rate; we test observed and excess
#     alternation and their model iv (type x brood size).
w_all <- d$watches
rate_by_chicks <- tapply(w_all$n_visits / w_all$duration_h, w_all$n_chicks, mean)
fit_rate_chicks <- glmmTMB(n_visits ~ scale(n_chicks) + offset(log(duration_h)) + (1 | nest) + (1 | season),
                           data = w_all, family = nbinom2)
brood <- do.call(rbind, lapply(names(wt), function(s) {
  x <- wt[[s]]
  mu <- mean(x$n_chicks); sdv <- sd(x$n_chicks)
  x$chicks_z <- (x$n_chicks - mu) / sdv
  x$off <- x$off_within
  f_obs <- fit_binomial("chicks_z", x)
  f_exc <- fit_binomial("chicks_z + offset(off)", x)
  st <- stack_datasets(nl$main[[s]], d$watches, "within")
  st$chicks_z <- (w_all$n_chicks[match(st$watch, w_all$watch)] - mu) / sdv
  iv <- fit_ihle_iv(st, "chicks_z")
  o <- coef_ci(f_obs, "chicks_z"); e <- coef_ci(f_exc, "chicks_z")
  data.frame(scenario = s, n = nrow(x), sd_chicks = sdv,
             obs_est = o[["est"]], obs_lo = o[["lo"]], obs_hi = o[["hi"]],
             exc_est = e[["est"]], exc_lo = e[["lo"]], exc_hi = e[["hi"]],
             iv_rr = iv[["rr"]], iv_lo = iv[["lo"]], iv_hi = iv[["hi"]])
}))

saveRDS(list(vs_null = vs_null, observed_pooled = observed_pooled,
             system_models = system_models, system_rate = sys_rate,
             n_system = table(pm$system), carer_class = carer_class,
             no_flag = no_flag, no_overflow = no_overflow,
             n_overflow = length(overflow_ids), sync_models = sync_models, span_b = span_b, synchrony = synchrony,
             rate_tab = rate_tab, rate_contrasts = rate_contrasts,
             rate_params = lapply(nl$main, `[[`, "rate_params"),
             share_among_pooled = sapply(nl$main, `[[`, "share_among_pooled"),
             n_iter = nl$main$mono$n_iter, birds_n = nrow(birds),
             fits_ok = c(sys_obs = fit_ok(fit_sys_obs), sys_excess = fit_ok(fit_sys_excess),
                         rate = fit_ok(fit_rate)),
             rate_family = family(fit_rate)$family, rate_aic = attr(fit_rate, "aic"),
             diagnostics = diagnostics, ihle_type = ihle_type, ihle_iv = ihle_iv,
             dyads_rand = dyads_rand, dyads_glmm = dyads_glmm,
             cv_by_role = cv_by_role, cv_models = cv_models,
             brood = brood, rate_by_chicks = rate_by_chicks,
             rate_chicks = coef_ci(fit_rate_chicks, "scale(n_chicks)"),
             chicks_table = table(w_all$n_chicks),
             cor_chicks_brood = cor(w_all$n_chicks, w_all$brood_size),
             watch_tables = wt),
        "output/results/models.rds")
print(vs_null[, c("scenario", "null", "log_or", "lo", "hi", "observed", "null_mean", "p_upper", "p_two")], digits = 3)
print(system_models, digits = 3)
print(rate_tab[, c("group", "rate", "lo", "hi")], digits = 3)
print(synchrony, digits = 3)
print(sync_models, digits = 3)
print(no_overflow, digits = 3)
print(ihle_type, digits = 3)
print(ihle_iv, digits = 3)
print(dyads_rand[dyads_rand$null %in% c("within", "switched"), ], digits = 3)
print(dyads_glmm, digits = 3)
print(cv_by_role)
print(cv_models, digits = 3)
print(brood, digits = 3)
