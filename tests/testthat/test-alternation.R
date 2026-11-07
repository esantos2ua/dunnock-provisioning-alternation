source(file.path("..", "..", "R", "alternation.R"))

# Exhaustive search: most alternations achievable by any ordering of the visits
brute_max <- function(counts) {
  memo <- new.env()
  best <- function(c, last) {
    if (sum(c) == 0) return(0)
    key <- paste(c(c, last), collapse = ",")
    if (!is.null(memo[[key]])) return(memo[[key]])
    out <- -Inf
    for (i in which(c > 0)) {
      cc <- c
      cc[i] <- cc[i] - 1
      out <- max(out, best(cc, i) + (last != 0 && last != i))
    }
    memo[[key]] <- out
    out
  }
  best(counts, 0)
}

test_that("count_alternations counts switches of visitor or class", {
  expect_equal(count_alternations(c("F", "M", "F", "M")), 3)
  expect_equal(count_alternations(c("M", "M", "F", "F")), 1)
  expect_equal(count_alternations("F"), 0)
  # trio seen through sex: F-a-b-F has 2 sex switches but 3 individual switches
  v <- c("F", "a", "b", "F")
  expect_equal(count_alternations(v), 3)
  expect_equal(count_alternations(v, class = c("F", "M", "M", "F")), 2)
})

test_that("max_alternations matches Ihle et al. (2019) for pairs", {
  ihle <- function(pm, pf) {
    pt <- pm + pf
    if (pm == pf) pt - abs(pm - pf) - 1 else pt - abs(pm - pf)
  }
  for (pm in 1:12) for (pf in 1:12) {
    expect_equal(max_alternations(c(pm, pf)), ihle(pm, pf))
  }
})

test_that("max_alternations equals exhaustive search for pairs and trios", {
  grid <- expand.grid(a = 0:8, b = 0:8, c = 0:8)
  for (i in seq_len(nrow(grid))) {
    counts <- unlist(grid[i, ])
    if (sum(counts) < 2) next
    expect_equal(max_alternations(counts), brute_max(counts), info = toString(counts))
  }
})

test_that("A_max is right when the top bird made just over half the visits", {
  # top bird made 12 of 22 visits (54.5%): the maximum is 20, not 21
  expect_equal(brute_max(c(12, 10, 0)), 20)
  expect_equal(max_alternations(c(12, 10, 0)), 20)
})
