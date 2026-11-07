source(file.path("..", "..", "R", "alternation.R"))
source(file.path("..", "..", "R", "data.R"))
source(file.path("..", "..", "R", "nulls.R"))

# Two toy watches: a pair and a trio
toy_visits <- function() {
  mk <- function(watch, role, t) data.frame(watch = watch, role = role, arrive = t,
                                             depart = t + 0.01, stringsAsFactors = FALSE)
  rbind(
    mk("w1", "f", c(0.1, 0.5, 0.9, 1.3)), mk("w1", "a", c(0.3, 0.7, 1.1)),
    mk("w2", "f", c(0.2, 0.6)), mk("w2", "a", c(0.1, 0.4, 0.8, 1.6, 1.7)),
    mk("w2", "b", c(1.0, 1.2, 1.5))
  )
}

test_that("observed scores match direct counts", {
  long <- make_long(toy_visits(), c("w1", "w2"), "poly_fab")
  s <- observed_scores(long)
  # w1 in time order: f a f a f a f; w2: a f a f a b b b a a
  expect_equal(s$A, c(6, count_alternations(c("a","f","a","f","a","b","b","b","a","a"))))
  expect_equal(s$Amax, c(max_alternations(c(4, 3)), max_alternations(c(2, 5, 3))))
  expect_equal(s$S, c(0, 0))
})

test_that("randomization nulls keep each bird's visits and gap multiset", {
  set.seed(1)
  long <- make_long(toy_visits(), c("w1", "w2"), "poly_fab")
  inner <- long$slot > 1
  for (g in list(null_within(long), null_switched(long), null_sorted(long))) {
    expect_equal(g[!inner], long$gap[!inner])          # first arrivals unchanged
    expect_equal(tapply(g[inner], long$bird[inner], sort),
                 tapply(long$gap[inner], long$bird[inner], sort))
  }
  expect_true(all(tapply(null_sorted(long)[inner], long$bird[inner],
                         function(x) !is.unsorted(x))))
})

test_that("the among null only mixes gaps within role and rate bin", {
  set.seed(2)
  long <- make_long(toy_visits(), c("w1", "w2"), "poly_fab")
  bin <- ifelse(long$role == "a", 1, long$slot)   # make every 'a' share one bin
  g <- null_among(long, bin)
  inner <- long$slot > 1
  keep <- inner & long$role != "a"
  expect_equal(g[keep], long$gap[keep])
  a <- inner & long$role == "a"
  expect_equal(sort(g[a]), sort(long$gap[a]))
})

test_that("male-only scenario drops female visits", {
  long <- make_long(toy_visits(), "w2", "poly_ab")
  expect_false(any(long$role == "f"))
  expect_equal(observed_scores(long)$A, count_alternations(c("a","a","a","b","b","b","a","a")))
})

test_that("lognormal parameters reproduce the target mean", {
  p <- rate_lognormal(c(2, 4, 6, 8), c(2, 2, 2, 2))
  expect_equal(exp(p[["meanlog"]] + p[["sdlog"]]^2 / 2), 5)
})

test_that("simulated watches have every bird visiting at least once", {
  set.seed(3)
  long <- make_long(toy_visits(), c("w1", "w2"), "poly_fab")
  params <- list(f = c(meanlog = 0, sdlog = 0.5), a = c(meanlog = 0, sdlog = 0.5),
                 b = c(meanlog = 0, sdlog = 0.5))
  b <- simulate_b(long, c(w1 = 1.5, w2 = 1.8), params)
  expect_equal(nrow(unique(b[, c("watch", "class")])), 5)
  expect_true(all(b$arrive >= 0 & b$arrive <= 1.8))
})

test_that("directed dyad counts add up to the transitions in each watch", {
  long <- make_long(toy_visits(), c("w1", "w2"), "poly_fab")
  s <- observed_scores(long, dyads = TRUE)
  expect_equal(unname(rowSums(s$D)), c(6, 9))                 # n visits - 1
  switches <- s$D[, !colnames(s$D) %in% c("f>f", "a>a", "b>b"), drop = FALSE]
  expect_equal(unname(rowSums(switches)), s$A)
  # w2 in time order: a f a f a b b b a a
  expect_equal(unname(s$D["w2", c("a>f", "f>a", "a>b", "b>b", "b>a", "a>a")]),
               c(2, 2, 1, 2, 1, 1))
})
