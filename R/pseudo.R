# Pseudo-pairs: a female and a male from different nests filmed on the same day
# (after Baldan, Hinde & Lessells 2019)
#
# The two birds cannot respond to each other, but they share the weather and
# the time of day. Times are aligned on the time since the start of filming
# (clock times were not kept), and each pseudo-pair is cut to the shorter of
# the two watches. If shared environmental variation produces the small excess
# over the within-watch null seen in real pairs, pseudo-pairs should show it too.

#' Visits of every pseudo-pair: female of watch i with the alpha (or only) male
#' of watch j, for all ordered combinations of different nests filmed on the
#' same date. Returns a visits table with watch (pseudo-pair id), date, role,
#' class, arrive, dur.
pseudo_pair_visits <- function(visits, watches, ids) {
  w <- watches[watches$watch %in% ids, ]
  combos <- do.call(rbind, lapply(split(w, w$date), function(x) {
    if (length(unique(x$nest)) < 2) return(NULL)
    g <- expand.grid(i = seq_len(nrow(x)), j = seq_len(nrow(x)))
    g <- g[x$nest[g$i] != x$nest[g$j], ]
    data.frame(f_watch = x$watch[g$i], m_watch = x$watch[g$j], date = x$date[g$i],
               end = pmin(x$duration_h[g$i], x$duration_h[g$j]), stringsAsFactors = FALSE)
  }))
  out <- do.call(rbind, lapply(seq_len(nrow(combos)), function(k) {
    cb <- combos[k, ]
    v <- visits[(visits$watch == cb$f_watch & visits$role == "f") |
                  (visits$watch == cb$m_watch & visits$role == "a"), ]
    v <- v[v$arrive < cb$end, ]
    if (length(unique(v$role)) < 2) return(NULL)
    data.frame(watch = paste(cb$f_watch, cb$m_watch, sep = "~"), f_watch = cb$f_watch,
               m_watch = cb$m_watch, date = cb$date, role = v$role,
               class = ifelse(v$role == "f", "f", "m"), arrive = v$arrive,
               dur = v$depart - v$arrive, stringsAsFactors = FALSE)
  }))
  out
}

#' Real pairs on the same watches: female vs alpha male only
real_pair_visits <- function(visits, ids, watches) {
  v <- visits[visits$watch %in% ids & visits$role %in% c("f", "a"), ]
  data.frame(watch = v$watch, f_watch = v$watch, m_watch = v$watch,
             date = watches$date[match(v$watch, watches$watch)], role = v$role,
             class = ifelse(v$role == "f", "f", "m"), arrive = v$arrive,
             dur = v$depart - v$arrive, stringsAsFactors = FALSE)
}

#' Pooled log odds ratio of observed to null alternation, with a bootstrap
#' that resamples filming dates (pseudo-pairs share birds within a date)
date_boot_or <- function(wt, null, B = 2000) {
  wt$off <- wt[[paste0("off_", null)]]
  by_date <- split(seq_len(nrow(wt)), wt$date)
  est <- replicate(B, {
    i <- unlist(by_date[sample.int(length(by_date), replace = TRUE)], use.names = FALSE)
    coef(suppressWarnings(glm(cbind(A, Amax - A) ~ 1 + offset(off), family = binomial,
                              data = wt[i, ])))[[1]]
  })
  c(lo = unname(quantile(est, 0.025)), hi = unname(quantile(est, 0.975)))
}
