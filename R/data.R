# Data preparation: from the visit-level sheets to analysis-ready tables

# Watches flagged during data auditing (a visit out of chronological order, or
# a colour-ring combination that does not match the nest roster); kept in the
# main analysis, dropped in a sensitivity analysis
flagged_watches <- c("focal9", "focal150", "focal585")

#' Read both sheets and keep the columns the analysis needs
#'
#' Visits are ordered by arrival time within each watch. Birds are
#' identified within a watch by their role (`f`, `a`, `b`), not by ring code,
#' so ring-reading errors do not split one bird into two.
read_visits <- function(dir = "data/raw") {
  files <- file.path(dir, c("mono_sheet_coordination.csv", "poli_sheet_coordination.csv"))
  raw <- do.call(rbind, lapply(files, function(f) {
    d <- read.csv(f, stringsAsFactors = FALSE)
    d[, c("est.rep", "nest.id", "sms", "female", "male.a", "male.b", "focal.id",
          "focal.date", "focal.length.h", "n.chicks", "brood.size", "chick.age",
          "visit.id", "parental.status", "parental.sex", "start.time", "end.time")]
  }))
  visits <- data.frame(
    season = factor(raw$est.rep),
    nest = paste0("s", raw$est.rep, "_n", raw$nest.id),
    system = factor(ifelse(raw$sms == "mono", "monogamous", "polyandrous"),
                    levels = c("monogamous", "polyandrous")),
    female = raw$female,
    male_a = raw$male.a,
    male_b = raw$male.b,
    watch = raw$focal.id,
    date = as.Date(raw$focal.date),
    length_h = raw$focal.length.h,
    n_chicks = raw$n.chicks,
    brood_size = raw$brood.size,
    chick_age = raw$chick.age,
    visit = raw$visit.id,
    role = raw$parental.status,
    sex = raw$parental.sex,
    arrive = raw$start.time,
    depart = raw$end.time,
    stringsAsFactors = FALSE
  )
  visits <- visits[order(visits$watch, visits$arrive), ]
  rownames(visits) <- NULL
  visits
}

#' One row per nest watch, with its observation window
summarise_watches <- function(visits) {
  by_watch <- split(visits, visits$watch)
  out <- do.call(rbind, lapply(by_watch, function(v) {
    span <- max(v$depart) - min(v$arrive)
    data.frame(
      watch = v$watch[1], season = v$season[1], nest = v$nest[1],
      system = v$system[1], female = v$female[1], male_a = v$male_a[1],
      male_b = v$male_b[1], date = v$date[1], n_chicks = v$n_chicks[1],
      brood_size = v$brood_size[1], chick_age = v$chick_age[1],
      length_h = v$length_h[1],
      span_h = span,
      duration_h = max(v$length_h[1], span),
      arrival_span_h = max(v$arrive) - min(v$arrive),
      n_visits = nrow(v),
      n_f = sum(v$role == "f"), n_a = sum(v$role == "a"), n_b = sum(v$role == "b"),
      flagged = v$watch[1] %in% flagged_watches,
      stringsAsFactors = FALSE
    )
  }))
  out$n_birds_seen <- (out$n_f > 0) + (out$n_a > 0) + (out$n_b > 0)
  # Classification by the number of carers that visited (sensitivity analysis)
  out$system_carers <- ifelse(out$n_f > 0 & (out$n_a > 0) + (out$n_b > 0) == 1,
                              "monogamous", "polyandrous")
  rownames(out) <- NULL
  out
}

#' Which watches enter each alternation scenario
#'
#' - `mono`: socially monogamous, both parents visited; female vs male
#' - `poly_fm`: socially polyandrous, female and at least one male visited;
#'   female vs males
#' - `poly_fab`: socially polyandrous, all three visited; among individuals
#' - `poly_ab`: socially polyandrous, both males visited; alpha vs beta in the
#'   male-only sequence (computed the same way for observed and null data)
scenario_watches <- function(watches) {
  w <- watches
  list(
    mono = w$watch[w$system == "monogamous" & w$n_f > 0 & w$n_a > 0],
    poly_fm = w$watch[w$system == "polyandrous" & w$n_f > 0 & (w$n_a + w$n_b) > 0],
    poly_fab = w$watch[w$system == "polyandrous" & w$n_f > 0 & w$n_a > 0 & w$n_b > 0],
    poly_ab = w$watch[w$system == "polyandrous" & w$n_a > 0 & w$n_b > 0]
  )
}

#' The class used to count alternation in each scenario
scenario_class <- function(role, scenario) {
  switch(scenario,
         mono = , poly_fm = ifelse(role == "f", "f", "m"),
         poly_fab = , poly_ab = role)
}

#' The visits of one watch that enter a scenario (poly_ab drops female visits)
scenario_visits <- function(v, scenario) {
  if (scenario == "poly_ab") v[v$role != "f", , drop = FALSE] else v
}

#' Regularity of each bird's visits: coefficient of variation (SD / mean) of its
#' intervals between arrivals, for birds with at least two intervals in a watch
#' (Halliwell et al. 2024 use the mean CV across carers as a covariate)
interval_cv <- function(visits) {
  v <- visits[order(visits$watch, visits$role, visits$arrive), ]
  birds <- split(v$arrive, paste(v$watch, v$role, sep = ":"))
  out <- do.call(rbind, lapply(names(birds), function(k) {
    g <- diff(birds[[k]])
    data.frame(watch = sub(":[^:]*$", "", k), role = sub("^.*:", "", k),
               n_intervals = length(g),
               cv = if (length(g) >= 2 && mean(g) > 0) sd(g) / mean(g) else NA_real_)
  }))
  out
}
