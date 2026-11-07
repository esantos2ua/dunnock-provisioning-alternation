# Pseudo-pairs of a female and a male from different nests filmed on the same
# day, against nulls D, E and F
# Run from the repo root: Rscript analysis/04b_pseudo_pairs.R

for (f in c("alternation", "data", "nulls", "models", "power", "pseudo"))
  source(file.path("R", paste0(f, ".R")))
d <- readRDS("data/derived/analysis_data.rds")
w <- d$watches
set.seed(20261003)

# Watches in which the female and the alpha (or only) male both visited
ids <- w$watch[w$n_f > 0 & w$n_a > 0]
pp <- pseudo_pair_visits(d$visits, w, ids)
real_ids <- intersect(ids, c(pp$f_watch, pp$m_watch))
rp <- real_pair_visits(d$visits, real_ids, w)

score <- function(v) {
  o <- def_offsets(relong(v), n_iter = 500)
  wt <- o$wt
  key <- v[!duplicated(v$watch), c("watch", "date", "f_watch", "m_watch")]
  wt <- merge(wt, key, by = "watch")
  tab <- do.call(rbind, lapply(c("within", "switched", "sorted"), function(k) {
    bt <- date_boot_or(wt, k)
    data.frame(null = k, observed = sum(wt$A) / sum(wt$Amax), null_rate = o$null_rate[[k]],
               t(glm_or(wt, k)), boot_lo = bt[["lo"]], boot_hi = bt[["hi"]])
  }))
  list(wt = wt, tab = tab)
}
pseudo <- score(pp)
real <- score(rp)

# Difference between real and pseudo-pairs in the excess over D (date bootstrap)
both <- rbind(cbind(pseudo$wt, real = 0L), cbind(real$wt, real = 1L))
fit_diff <- function(x) coef(suppressWarnings(glm(cbind(A, Amax - A) ~ real + offset(off_within),
                                                  family = binomial, data = x)))[["real"]]
by_date <- split(seq_len(nrow(both)), both$date)
diff_boot <- replicate(2000, fit_diff(both[unlist(by_date[sample.int(length(by_date), replace = TRUE)]), ]))
difference <- c(est = fit_diff(both), lo = unname(quantile(diff_boot, 0.025)),
                hi = unname(quantile(diff_boot, 0.975)))

out <- list(pseudo = pseudo$tab, real = real$tab, difference = difference,
            n_pseudo = nrow(pseudo$wt), n_real = nrow(real$wt),
            n_dates = length(unique(pseudo$wt$date)),
            n_watches_used = length(unique(c(pseudo$wt$f_watch, pseudo$wt$m_watch))))
saveRDS(out, "output/results/pseudo_pairs.rds")
print(out, digits = 3)
