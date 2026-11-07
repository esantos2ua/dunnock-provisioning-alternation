source(file.path("..", "..", "R", "alternation.R"))
source(file.path("..", "..", "R", "data.R"))
source(file.path("..", "..", "R", "nulls.R"))
source(file.path("..", "..", "R", "power.R"))

# Tests for the power simulation and pseudo-pairs

toy <- function() {
  data.frame(watch = "w1", role = c("f", "f", "f", "a", "a", "a"),
             bird = paste0("w1:", c("f", "f", "f", "a", "a", "a")),
             class = c("f", "f", "f", "m", "m", "m"),
             arrive = c(0.1, 0.2, 0.3, 0.5, 0.6, 0.7), dur = 0.01,
             stringsAsFactors = FALSE)
}

test_that("p = 0 leaves the visits unchanged", {
  x <- toy()
  y <- inject_turn_taking(x, 0)
  expect_equal(y$arrive, x$arrive)
})

test_that("injection keeps every bird's number of visits", {
  set.seed(1)
  y <- inject_turn_taking(toy(), 0.5)
  expect_equal(as.vector(table(y$role)[c("f", "a")]), c(3, 3))
})

test_that("p = 1 produces complete alternation when counts allow it", {
  y <- inject_turn_taking(toy(), 1)
  y <- y[order(y$arrive), ]
  expect_equal(count_alternations(y$class), max_alternations(c(3, 3)))
})

test_that("relong rebuilds gaps that sum to the arrival times", {
  L <- relong(toy())
  expect_equal(arrivals_from_gaps(L$gap, L$bird), L$arrive)
})
