# Alternation metrics for provisioning sequences
#
# A nest watch is an ordered vector of visitor identities (one element per
# provisioning visit). Alternation is a visit by a different caregiver (or
# caregiver class, e.g. sex) from the one that made the previous visit.

#' Count observed alternations (A_obs) in a visit sequence
#'
#' @param visitor Vector of visitor identities, in time order.
#' @param class Optional vector of the same length mapping each visit to a
#'   caregiver class (e.g. sex). When given, alternation is counted between
#'   classes rather than between individuals.
count_alternations <- function(visitor, class = NULL) {
  x <- if (is.null(class)) visitor else class
  stopifnot(length(x) == length(visitor))
  if (length(x) < 2) return(0L)
  sum(x[-1] != x[-length(x)])
}

#' Maximum possible alternations (A_max) given visit counts per caregiver
#'
#' Works for any number of caregivers (or classes). If the most active
#' caregiver made no more than one visit beyond all the others combined,
#' every consecutive pair can alternate (A_max = P_T - 1); otherwise the
#' surplus visits must run consecutively and A_max = 2 * (P_T - P_max).
#' For two caregivers this reduces to Ihle et al. (2019):
#' P_T - |P_m - P_f| - 1 if P_m == P_f, else P_T - |P_m - P_f|.
#'
#' @param counts Integer vector of visits per caregiver.
max_alternations <- function(counts) {
  counts <- counts[counts > 0]
  total <- sum(counts)
  if (total < 2) return(0L)
  top <- max(counts)
  rest <- total - top
  as.integer(if (top <= rest + 1) total - 1 else 2 * rest)
}

#' Alternation rate (A_obs / A_max) for one nest watch
alternation_rate <- function(visitor, class = NULL) {
  units <- if (is.null(class)) visitor else class
  a_max <- max_alternations(as.vector(table(units)))
  if (a_max == 0) return(NA_real_)
  count_alternations(visitor, class) / a_max
}
