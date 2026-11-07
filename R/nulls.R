# Null models for provisioning alternation (after Ihle et al. 2019)
#
# The among-watch null pools intervals across watches, and in the switched
# null the bird whose intervals are swapped is picked at random. Everything
# works on a long table with one row per visit.
#
# Columns of the long table built by make_long():
#   watch, bird (watch:role), class (what alternation is counted between),
#   slot (visit number within bird), arrive, dur (visit duration, h),
#   gap (interval since the bird's previous arrival; first arrival for slot 1)

make_long <- function(visits, watches_in, scenario) {
  v <- visits[visits$watch %in% watches_in, ]
  v <- do.call(rbind, lapply(split(v, v$watch), scenario_visits, scenario = scenario))
  v <- v[order(v$watch, v$role, v$arrive), ]
  bird <- paste(v$watch, v$role, sep = ":")
  slot <- ave(seq_along(bird), bird, FUN = seq_along)
  gap <- ave(v$arrive, bird, FUN = function(t) c(t[1], diff(t)))
  data.frame(watch = v$watch, role = v$role, bird = bird,
             class = scenario_class(v$role, scenario), slot = slot,
             arrive = v$arrive, dur = v$depart - v$arrive, gap = gap,
             stringsAsFactors = FALSE)
}

# Arrival times from gaps: cumulative sum within bird (rows grouped by bird,
# in slot order)
arrivals_from_gaps <- function(gap, bird) ave(gap, bird, FUN = cumsum)

#' Alternation, maximum alternation and synchrony for every watch at once
#'
#' Two synchrony measures, the same for observed and null data:
#' - S (presence overlap): an arrival while a bird of another class is still at
#'   the nest, judged against the preceding visit;
#' - S2 (Ihle et al. 2019; Bebbington & Hatchwell 2016): an alternated visit
#'   whose arrival is within `sync_window` hours (2 min) of the previous arrival.
#' With `dyads = TRUE`, also returns D: a watch x transition matrix counting
#' each directed transition between classes ("f>a" = a visit by a that follows
#' a visit by f), including repeats ("a>a").
dyad_levels <- as.vector(outer(c("f", "a", "b"), c("f", "a", "b"), paste, sep = ">"))
score_long <- function(watch, arrive, dur, class, watch_levels, sync_window = 2 / 60,
                       dyads = FALSE) {
  o <- order(watch, arrive)
  w <- factor(watch[o], levels = watch_levels)
  cl <- class[o]
  a <- arrive[o]
  d <- a + dur[o]
  n <- length(w)
  same_watch <- w[-1] == w[-n]
  switch_ <- same_watch & cl[-1] != cl[-n]
  overlap <- switch_ & a[-1] < d[-n]
  window <- switch_ & (a[-1] - a[-n]) <= sync_window
  per_watch <- function(x) {
    v <- as.vector(rowsum(as.integer(x), w[-1], reorder = TRUE)[watch_levels, 1])
    replace(v, is.na(v), 0L)
  }
  alt <- per_watch(switch_)
  syn <- per_watch(overlap)
  syn2 <- per_watch(window)
  counts <- table(w, cl)
  amax <- apply(counts, 1, max_alternations)
  out <- list(A = alt, Amax = as.integer(amax), S = syn, S2 = syn2)
  if (dyads) {
    key <- factor(paste(cl[-n], cl[-1], sep = ">")[same_watch], levels = dyad_levels)
    out$D <- unclass(table(w[-1][same_watch], key))[watch_levels, , drop = FALSE]
  }
  out
}

observed_scores <- function(long, dyads = FALSE) {
  lev <- sort(unique(long$watch))
  score_long(long$watch, long$arrive, long$dur, long$class, lev, dyads = dyads)
}

# ---- Randomization nulls (C, D, E, F): permute each bird's gaps -------------

# Shuffle x within groups g
shuffle_within <- function(x, g) {
  idx <- split(seq_along(x), g)
  out <- x
  for (i in idx) if (length(i) > 1) out[i] <- x[i[sample.int(length(i))]]
  out
}

#' Dataset D: gaps shuffled within each bird within its watch (Johnstone et al. 2014)
null_within <- function(long) {
  g <- long$gap
  inner <- long$slot > 1
  g[inner] <- shuffle_within(long$gap[inner], long$bird[inner])
  g
}

#' Dataset C: gaps shuffled among birds of the same role and the same whole
#' number of visits per hour, across watches (Bebbington & Hatchwell 2016).
#' `rate_bin` is one value per row (constant within bird).
null_among <- function(long, rate_bin) {
  g <- long$gap
  inner <- long$slot > 1
  key <- paste(long$role, rate_bin, sep = ":")
  g[inner] <- shuffle_within(long$gap[inner], key[inner])
  g
}

#' Dataset E: one randomly chosen bird per watch has consecutive gaps swapped
#' pairwise (2<->3, 4<->5, ...); the other birds are unchanged
null_switched <- function(long) {
  g <- long$gap
  birds_by_watch <- split(unique(long$bird[long$slot > 2]),
                          sub(":[^:]*$", "", unique(long$bird[long$slot > 2])))
  pick <- vapply(birds_by_watch, function(b) b[sample.int(length(b), 1)], "")
  for (b in pick) {
    i <- which(long$bird == b & long$slot > 1)
    k <- length(i)
    if (k < 2) next
    pairs <- seq(1, k - 1, by = 2)
    swapped <- i
    swapped[pairs] <- i[pairs + 1]
    swapped[pairs + 1] <- i[pairs]
    g[i] <- long$gap[swapped]
  }
  g
}

#' Dataset F: every bird's gaps sorted in increasing order (deterministic)
null_sorted <- function(long) {
  g <- long$gap
  inner <- long$slot > 1
  g[inner] <- ave(long$gap[inner], long$bird[inner], FUN = sort)
  g
}

# ---- Simulated null (B) -------------------------------------------------------

#' Lognormal parameters for per-bird visit rates, on the expected scale
#'
#' The observed rate r = n / D carries Poisson noise with variance lambda / D,
#' so var(lambda) is estimated as var(r) - mean(r / D) (Ihle et al. 2019).
#' Moment matching then gives meanlog and sdlog (meanlog = log(m) - sdlog^2 / 2).
rate_lognormal <- function(rate, duration) {
  m <- mean(rate)
  v <- max(var(rate) - mean(rate / duration), 1e-6 * m^2)
  s2 <- log(1 + v / m^2)
  c(meanlog = log(m) - s2 / 2, sdlog = sqrt(s2))
}

#' Dataset B: for each observed watch, simulate the same birds (roles) with
#' independent Poisson arrivals over the watch's own duration. Each bird's
#' expected rate is drawn from the lognormal for its role; counts of zero are
#' redrawn because every bird in the observed data visited at least once.
#' Visit durations are drawn from the observed durations.
simulate_b <- function(long, duration_h, params) {
  birds <- unique(long[, c("watch", "role", "class", "bird")])
  out <- vector("list", nrow(birds))
  for (i in seq_len(nrow(birds))) {
    p <- params[[birds$role[i]]]
    D <- duration_h[[birds$watch[i]]]
    n <- 0
    while (n == 0) n <- rpois(1, rlnorm(1, p["meanlog"], p["sdlog"]) * D)
    out[[i]] <- data.frame(watch = birds$watch[i], class = birds$class[i],
                           arrive = runif(n, 0, D),
                           dur = sample(long$dur, n, replace = TRUE))
  }
  do.call(rbind, out)
}

# ---- Driver -----------------------------------------------------------------

#' Run all nulls for one scenario
#'
#' Returns observed scores and, for each null, matrices (watch x iteration) of
#' A and S, plus A_max (a matrix for B, a vector otherwise).
run_nulls <- function(visits, watches, scenario, n_iter = 1000,
                      duration = c("duration_h", "arrival_span_h")) {
  duration <- match.arg(duration)
  ids <- scenario_watches(watches)[[scenario]]
  long <- make_long(visits, ids, scenario)
  lev <- sort(unique(long$watch))
  dyads <- scenario == "poly_fab"
  obs <- observed_scores(long, dyads)

  wt <- watches[match(lev, watches$watch), ]
  dur_h <- setNames(pmax(wt[[duration]], 1 / 60), wt$watch)   # guard: >= 1 min
  bird_n <- table(long$bird)
  bird_watch <- sub(":[^:]*$", "", names(bird_n))
  bird_role <- sub("^.*:", "", names(bird_n))
  bird_rate <- as.vector(bird_n) / dur_h[bird_watch]
  params <- lapply(split(seq_along(bird_rate), bird_role), function(i)
    rate_lognormal(bird_rate[i], dur_h[bird_watch[i]]))
  rate_bin <- round(bird_rate)[match(long$bird, names(bird_n))]

  rescore <- function(gap) {
    a <- arrivals_from_gaps(gap, long$bird)
    s <- score_long(long$watch, a, long$dur, long$class, lev, dyads = dyads)
    s[intersect(c("A", "S", "S2", "D"), names(s))]
  }
  collect <- function(res, amax) {
    out <- list(A = sapply(res, `[[`, "A"), S = sapply(res, `[[`, "S"),
                S2 = sapply(res, `[[`, "S2"), Amax = amax)
    if (dyads) out$D <- simplify2array(lapply(res, `[[`, "D"))   # watch x dyad x iter
    out
  }
  iterate <- function(fun) collect(replicate(n_iter, rescore(fun()), simplify = FALSE), obs$Amax)

  sim <- replicate(n_iter, {
    b <- simulate_b(long, dur_h, params)
    score_long(b$watch, b$arrive, b$dur, b$class, lev, dyads = dyads)
  }, simplify = FALSE)
  sorted <- rescore(null_sorted(long))

  list(
    scenario = scenario, watch = lev, n_iter = n_iter, duration = duration,
    rate_params = params,
    share_among_pooled = mean(ave(rep(1, length(bird_rate)),
                                  paste(bird_role, round(bird_rate)), FUN = length) > 1),
    observed = obs,
    nulls = list(
      simulated = collect(sim, sapply(sim, `[[`, "Amax")),
      among = iterate(function() null_among(long, rate_bin)),
      within = iterate(function() null_within(long)),
      switched = iterate(function() null_switched(long)),
      sorted = collect(list(sorted), obs$Amax)
    )
  )
}
