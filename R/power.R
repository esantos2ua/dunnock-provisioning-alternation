# Power of the observed-vs-null comparison to detect turn-taking
#
# Each simulated data set starts from the observed watches with every bird's
# intervals reshuffled within the watch (the generating process of null D), so
# visit numbers, interval distributions and refractory periods are realistic and
# there is no response between carers. A known turn-taking rule is then
# injected, and the same pipeline as for the real data (nulls D, E and F, the
# logit-offset model and the nest-level cluster bootstrap) is run on the result.

#' Inject turn-taking into one scenario's visits
#'
#' Visits are replayed in time order, each bird keeping its own sequence of
#' intervals. When the next visit would repeat the class of the previous visit
#' and a carer of another class still has visits to make, the bird waits with
#' probability `p`: its arrival is moved to just after the next arrival of a
#' carer of another class (`lag`, in hours, drawn uniformly), and its later
#' visits follow from there with their original intervals. With p = 0 the
#' visits are unchanged. Visit numbers are kept, so A_max does not change.
#'
#' @param x Data frame with watch, role, bird, class, arrive, dur.
inject_turn_taking <- function(x, p, lag = c(0.5, 2) / 60) {
  if (p == 0) return(x[, c("watch", "role", "bird", "class", "arrive", "dur")])
  do.call(rbind, lapply(split(x, x$watch), function(v) {
    v <- v[order(v$bird, v$arrive), ]
    birds <- split(seq_len(nrow(v)), v$bird)
    gaps <- lapply(birds, function(i) diff(c(0, v$arrive[i])))
    cls <- vapply(birds, function(i) v$class[i[1]], "")
    n <- lengths(birds)
    k <- setNames(rep(1L, length(birds)), names(birds))
    nxt <- vapply(gaps, `[`, 0, 1)
    t_out <- numeric(nrow(v))
    row_out <- integer(nrow(v))
    j <- 0L
    last <- NA_character_
    while (j < nrow(v)) {
      alive <- k <= n
      b <- names(birds)[which.min(ifelse(alive, nxt, Inf))]
      others <- alive & cls != cls[[b]]
      if (!is.na(last) && cls[[b]] == last && any(others) && runif(1) < p) {
        nxt[[b]] <- min(nxt[others]) + runif(1, lag[1], lag[2])
        next
      }
      j <- j + 1L
      t_out[j] <- nxt[[b]]
      row_out[j] <- birds[[b]][k[[b]]]
      last <- cls[[b]]
      k[[b]] <- k[[b]] + 1L
      if (k[[b]] <= n[[b]]) nxt[[b]] <- nxt[[b]] + gaps[[b]][k[[b]]]
    }
    out <- v[row_out, c("watch", "role", "bird", "class", "dur")]
    out$arrive <- t_out
    out
  }))
}

#' Long table (as make_long) from a data frame of visits with watch, role,
#' arrive and dur
relong <- function(x) {
  x <- x[order(x$watch, x$role, x$arrive), ]
  bird <- paste(x$watch, x$role, sep = ":")
  data.frame(watch = x$watch, role = x$role, bird = bird, class = x$class,
             slot = ave(seq_along(bird), bird, FUN = seq_along), arrive = x$arrive,
             dur = x$dur, gap = ave(x$arrive, bird, FUN = function(t) c(t[1], diff(t))),
             stringsAsFactors = FALSE)
}

#' Observed scores and the D, E and F null expectations for a long table
#'
#' Returns a watch table (watch, A, Amax, off_within, off_switched, off_sorted)
#' with the logit offsets used by compare_to_null(), and the pooled null rates.
def_offsets <- function(L, n_iter = 100) {
  lev <- sort(unique(L$watch))
  obs <- observed_scores(L)
  rescore <- function(gap) {
    a <- arrivals_from_gaps(gap, L$bird)
    score_long(L$watch, a, L$dur, L$class, lev)$A
  }
  null_A <- list(
    within = replicate(n_iter, rescore(null_within(L))),
    switched = replicate(n_iter, rescore(null_switched(L))),
    sorted = matrix(rescore(null_sorted(L)), ncol = 1)
  )
  wt <- data.frame(watch = lev, A = obs$A, Amax = obs$Amax, stringsAsFactors = FALSE)
  for (k in names(null_A))
    wt[[paste0("off_", k)]] <- qlogis((rowSums(null_A[[k]]) + 0.5) /
                                        (obs$Amax * ncol(null_A[[k]]) + 1))
  list(wt = wt, null_rate = vapply(null_A, function(a) sum(a) / (sum(obs$Amax) * ncol(a)), 0))
}

# Pooled log odds ratio of observed to null alternation, with its Wald interval
glm_or <- function(wt, null) {
  wt$off <- wt[[paste0("off_", null)]]
  fit <- suppressWarnings(glm(cbind(A, Amax - A) ~ 1 + offset(off), family = binomial, data = wt))
  est <- coef(fit)[[1]]
  se <- sqrt(vcov(fit)[1, 1])
  c(log_or = est, wald_lo = est - 1.96 * se, wald_hi = est + 1.96 * se)
}

#' One power replicate: simulate, inject, run nulls D, E and F, and estimate
#' the observed-vs-null log odds ratio with its cluster-bootstrap interval
power_replicate <- function(long, nests, p, n_iter = 100, B = 500) {
  base <- long
  base$arrive <- arrivals_from_gaps(null_within(long), long$bird)
  L <- relong(inject_turn_taking(base, p))
  base_obs <- score_long(base$watch, base$arrive, base$dur, base$class, sort(unique(L$watch)))
  o <- def_offsets(L, n_iter)
  wt <- o$wt
  wt$nest <- nests[wt$watch]
  do.call(rbind, lapply(c("within", "switched", "sorted"), function(k) {
    bt <- cluster_boot_or(wt, k, B = B)
    data.frame(null = k, p = p, base_rate = sum(base_obs$A) / sum(base_obs$Amax),
               obs_rate = sum(wt$A) / sum(wt$Amax), null_rate = o$null_rate[[k]],
               t(glm_or(wt, k)), boot_lo = bt[["lo"]], boot_hi = bt[["hi"]])
  }))
}
