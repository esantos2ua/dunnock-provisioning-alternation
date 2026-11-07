# Statistical comparison of observed and null alternation

suppressMessages({ library(glmmTMB); library(DHARMa) })

null_labels <- c(simulated = "Simulated (B)", among = "Among watches (C)",
                 within = "Within watch (D)", switched = "Switched (E)",
                 sorted = "Sorted (F)")
scenario_labels <- c(mono = "Monogamous pairs, female vs male",
                     poly_fm = "Polyandrous groups, female vs males",
                     poly_fab = "Polyandrous trios, among individuals",
                     poly_ab = "Polyandrous groups, alpha vs beta male")

random_part <- "(1 | nest) + (1 | female) + (1 | male_a) + (1 | season)"

#' Expected null alternation rate for each watch, with a 0.5/1 continuity
#' correction so the logit offset is finite
null_rate <- function(n) {
  amax <- if (is.matrix(n$Amax)) rowSums(n$Amax) else n$Amax * ncol(n$A)
  (rowSums(n$A) + 0.5) / (amax + 1)
}

#' Watch-level table for one scenario: observed counts plus each null's offset
watch_table <- function(res, watches) {
  w <- watches[match(res$watch, watches$watch), ]
  w$A <- res$observed$A
  w$Amax <- res$observed$Amax
  w$S <- res$observed$S
  for (k in names(res$nulls)) w[[paste0("off_", k)]] <- qlogis(null_rate(res$nulls[[k]]))
  w
}

#' Binomial GLMM of observed alternation. If the full random-effect structure
#' does not give a positive-definite Hessian (variances at the boundary), the
#' random part is reduced to nest only, then dropped; the structure used is
#' stored in attr(fit, "random").
fit_binomial <- function(formula_rhs, data) {
  for (re in c(random_part, "(1 | nest)", "")) {
    f <- as.formula(paste("cbind(A, Amax - A) ~", formula_rhs, if (nzchar(re)) paste("+", re)))
    fit <- suppressWarnings(glmmTMB(f, data = data, family = binomial))
    if (fit_ok(fit)) break
  }
  attr(fit, "random") <- if (nzchar(re)) re else "none"
  fit
}

# Wald estimate and 95% CI for one fixed effect (conditional model, link scale)
coef_ci <- function(fit, term) {
  est <- fixef(fit)$cond[[term]]
  se <- sqrt(diag(vcov(fit)$cond))[[term]]
  c(est = est, lo = est - 1.96 * se, hi = est + 1.96 * se)
}

# Largest random-effect SD, and whether the fit converged with a positive-definite Hessian
re_sd_max <- function(fit) {
  sds <- unlist(lapply(VarCorr(fit)$cond, attr, "stddev"))
  if (length(sds)) max(sds) else 0
}
fit_ok <- function(fit) {
  fit$fit$convergence == 0 && isTRUE(fit$sdr$pdHess) &&
    all(is.finite(suppressWarnings(sqrt(diag(vcov(fit)$cond)))))
}

#' Observed vs one null: intercept of a model with offset logit(null rate) is
#' the log odds ratio of observed to expected alternation
compare_to_null <- function(wt, null) {
  wt$off <- wt[[paste0("off_", null)]]
  fit <- fit_binomial("1 + offset(off)", wt)
  out <- coef_ci(fit, "(Intercept)")
  c(out, re_sd_max = re_sd_max(fit), ok = fit_ok(fit),
    full_random = attr(fit, "random") == random_part)
}

#' Population-level randomization test on the pooled alternation rate
#' (sum A / sum A_max over watches), one value per null iteration
randomization_test <- function(res, null) {
  n <- res$nulls[[null]]
  obs <- sum(res$observed$A) / sum(res$observed$Amax)
  amax <- if (is.matrix(n$Amax)) colSums(n$Amax) else rep(sum(n$Amax), ncol(n$A))
  sims <- colSums(n$A) / amax
  R <- length(sims)
  p_upper <- (1 + sum(sims >= obs)) / (R + 1)
  p_lower <- (1 + sum(sims <= obs)) / (R + 1)
  list(observed = obs, null_mean = mean(sims),
       null_lo = unname(quantile(sims, 0.025)), null_hi = unname(quantile(sims, 0.975)),
       p_upper = if (R > 1) p_upper else NA, p_two = if (R > 1) min(1, 2 * min(p_upper, p_lower)) else NA,
       iterations = R)
}

#' Pooled observed rate with a watch-level bootstrap interval
pooled_observed <- function(res, B = 2000) {
  A <- res$observed$A
  M <- res$observed$Amax
  boots <- replicate(B, { i <- sample.int(length(A), replace = TRUE); sum(A[i]) / sum(M[i]) })
  c(rate = sum(A) / sum(M), lo = unname(quantile(boots, 0.025)), hi = unname(quantile(boots, 0.975)))
}

#' Synchrony against one null; `measure` is "S" (presence overlap) or "S2"
#' (alternated arrival within 2 min)
synchrony_test <- function(res, null = "within", measure = "S") {
  s <- colSums(res$nulls[[null]][[measure]])
  obs <- sum(res$observed[[measure]])
  R <- length(s)
  list(measure = measure, observed = obs, switches = sum(res$observed$A), null_mean = mean(s),
       null_lo = unname(quantile(s, 0.025)), null_hi = unname(quantile(s, 0.975)),
       p_lower = (1 + sum(s <= obs)) / (R + 1), p_upper = (1 + sum(s >= obs)) / (R + 1))
}

# ---- Models of Ihle et al. (2019) ------------------------------------------------

#' Per-watch median of the null iterations as an integer; halves are rounded up
#' or down at random (Ihle et al.'s median_integer)
median_integer <- function(x) {
  m <- median(x)
  if (m %% 1 == 0.5) m + sample(c(-0.5, 0.5), 1) else m
}

#' Observed and null alternation stacked in long format, one row per watch and
#' data set, with each null summarized by its per-watch median (Ihle et al. 2019)
stack_datasets <- function(res, watches, nulls = names(res$nulls)) {
  w <- watches[match(res$watch, watches$watch), c("watch", "nest", "season", "system",
                                                  "female", "male_a", "n_birds_seen")]
  obs <- cbind(w, type = "observed", A = res$observed$A, Amax = res$observed$Amax)
  sims <- lapply(nulls, function(k) {
    n <- res$nulls[[k]]
    A <- apply(n$A, 1, median_integer)
    Amax <- if (is.matrix(n$Amax)) apply(n$Amax, 1, median_integer) else n$Amax
    cbind(w, type = k, A = A, Amax = Amax)
  })
  d <- do.call(rbind, c(list(obs), sims))
  d$type <- factor(d$type, levels = c("observed", nulls))
  d$obs_id <- factor(seq_len(nrow(d)))
  d
}

#' Ihle et al.'s comparison of observed and null alternation: Poisson GLMM of A
#' on data-set type with a random intercept for watch. Returns, for each null,
#' the rate ratio observed/null (exp of minus the type coefficient).
fit_ihle_type <- function(stacked) {
  fit <- glmmTMB(A ~ type + (1 | watch), data = stacked, family = poisson)
  b <- fixef(fit)$cond[-1]
  se <- sqrt(diag(vcov(fit)$cond))[-1]
  data.frame(null = sub("^type", "", names(b)), rr = exp(-b),
             lo = exp(-b - 1.96 * se), hi = exp(-b + 1.96 * se), ok = fit_ok(fit))
}

#' Ihle et al.'s recommended model (their model iv) for a predictor of
#' coordination: observed and within-watch null stacked, Poisson GLMM with
#' type x predictor. The interaction is the difference between groups in
#' observed alternation relative to the null (rate ratio).
fit_ihle_iv <- function(stacked, predictor) {
  f <- as.formula(paste("A ~ type *", predictor,
                        "+ (1 | watch) + (1 | nest) + (1 | female) + (1 | male_a) + (1 | season) + (1 | obs_id)"))
  fit <- suppressWarnings(glmmTMB(f, data = stacked, family = poisson))
  if (!fit_ok(fit)) {
    f <- update(f, . ~ . - (1 | female) - (1 | male_a) - (1 | season))
    fit <- suppressWarnings(glmmTMB(f, data = stacked, family = poisson))
  }
  b <- fixef(fit)$cond
  V <- vcov(fit)$cond
  term <- grep(":", names(b), value = TRUE)
  # sign: 'observed' is the reference type, so observed-minus-null = -interaction
  est <- -b[[term]]
  se <- sqrt(V[term, term])
  c(rr = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se), ok = fit_ok(fit))
}

#' Visits per bird per watch, with zeros for roster birds that never came
bird_visit_table <- function(watches) {
  mk <- function(role, n, id) data.frame(watch = watches$watch, nest = watches$nest,
                                          season = watches$season, system = watches$system,
                                          duration_h = watches$duration_h, role = role,
                                          visits = n, bird = id, stringsAsFactors = FALSE)
  b <- rbind(mk("f", watches$n_f, watches$female), mk("a", watches$n_a, watches$male_a),
             mk("b", watches$n_b, watches$male_b))
  b <- b[!(b$system == "monogamous" & b$role == "b"), ]
  b$group <- factor(paste(b$system, b$role),
                    levels = c("monogamous f", "monogamous a", "polyandrous f",
                               "polyandrous a", "polyandrous b"),
                    labels = c("Female, pair", "Male, pair", "Female, trio",
                               "Alpha male, trio", "Beta male, trio"))
  b
}

#' Provisioning-rate GLMM; the count family is chosen by AIC (Poisson vs
#' negative binomial, nbinom2)
fit_rates <- function(birds) {
  f <- visits ~ 0 + group + offset(log(duration_h)) + (1 | bird) + (1 | watch) +
    (1 | nest) + (1 | season)
  fits <- list(poisson = glmmTMB(f, data = birds, family = poisson),
               nbinom2 = glmmTMB(f, data = birds, family = nbinom2))
  aic <- sapply(fits, AIC)
  best <- fits[[which.min(aic)]]
  attr(best, "aic") <- aic
  best
}

#' Nest-level cluster bootstrap of the observed-vs-null log odds ratio
#'
#' Robustness check for compare_to_null(): alternation counts are
#' underdispersed relative to the binomial and every random-effect variance was
#' estimated at zero, so the GLMM Wald intervals are conservative. Resampling
#' nests keeps the non-independence of watches from the same nest.
cluster_boot_or <- function(wt, null, B = 2000) {
  wt$off <- wt[[paste0("off_", null)]]
  by_nest <- split(seq_len(nrow(wt)), wt$nest)
  est <- replicate(B, {
    i <- unlist(by_nest[sample.int(length(by_nest), replace = TRUE)], use.names = FALSE)
    x <- wt[i, ]
    coef(suppressWarnings(glm(cbind(A, Amax - A) ~ 1 + offset(off), family = binomial, data = x)))[[1]]
  })
  c(lo = unname(quantile(est, 0.025)), hi = unname(quantile(est, 0.975)))
}

#' Dispersion of a fitted GLMM: Pearson ratio and the DHARMa simulation test
dispersion <- function(fit) {
  pearson <- sum(residuals(fit, type = "pearson")^2) / df.residual(fit)
  sim <- simulateResiduals(fit, n = 1000, plot = FALSE)
  t <- testDispersion(sim, plot = FALSE)
  c(pearson = pearson, dharma_ratio = unname(t$statistic), dharma_p = t$p.value)
}


# ---- Directed dyads in trios -----------------------------------------------------

#' Observed vs null totals of each directed transition (pooled over watches),
#' with the 2.5-97.5% range across iterations and a two-sided randomization p
dyad_test <- function(res, null) {
  obs <- colSums(res$observed$D)
  sims <- apply(res$nulls[[null]]$D, c(2, 3), sum)          # dyad x iteration
  R <- ncol(sims)
  data.frame(dyad = names(obs), null = null, observed = obs, expected = rowMeans(sims),
             lo = apply(sims, 1, quantile, 0.025), hi = apply(sims, 1, quantile, 0.975),
             p_two = if (R > 1) pmin(1, 2 * pmin((1 + rowSums(sims >= obs)) / (R + 1),
                                                 (1 + rowSums(sims <= obs)) / (R + 1))) else NA,
             row.names = NULL)
}

#' Watch-level Poisson GLMM for one directed transition: observed count with
#' offset log(expected under the null); exp(intercept) is the rate ratio
dyad_glmm <- function(res, watches, null, dyad) {
  w <- watches[match(res$watch, watches$watch), c("watch", "nest")]
  w$y <- res$observed$D[, dyad]
  w$e <- apply(res$nulls[[null]]$D[, dyad, , drop = FALSE], 1, mean)
  w <- w[w$e > 0, ]
  fit <- suppressWarnings(glmmTMB(y ~ 1 + offset(log(e)) + (1 | nest), data = w, family = poisson))
  if (!fit_ok(fit)) fit <- glmmTMB(y ~ 1 + offset(log(e)), data = w, family = poisson)
  est <- fixef(fit)$cond[[1]]
  se <- sqrt(vcov(fit)$cond[1, 1])
  c(rr = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))
}
