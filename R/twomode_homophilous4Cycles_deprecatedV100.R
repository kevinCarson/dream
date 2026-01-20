## The Creation of Homophliy Four Cyles
## see Fujimoto, Kayo, Tom Snijders, and Thomas Valente. 2018.
##   "Multivariate dynamics of one-mode and two-mode networks: Explaining similarity
##   in sports participation among friends." Network Science, Vol 6(3), pp. 370-395
## Code written by Kevin Carson (kacarson@arizona.edu) and Deigo Leal (https://www.diegoleal.info/)
## Last Updated: 09-07-24
#' @title Compute Fujimoto, Snijders, and Valente's (2018) Homophilous Four-Cycles for Two-Mode Networks
#' @name computeHomFourCycles
#' @param net The two-mode adjacency matrix.
#' @param mem The vector of membership values that the homophilous four-cycles will be based on.
#' @return The vector of counts of homophilous four-cycles for the two-mode network.
#' @export
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `computeHomFourCycles()` has been deprecated starting on version 1.0.0 of the `dream` package. Please use the `netstats_tm_homfourcycles()` function and see the `NEWS.md` file for more details.
#'
#'
#' This function computes the number of homophilous four-cycles in
#' a two-mode network as proposed by Fujimoto, Snijders, and Valente (2018: 380).
#' See Fujimoto, Snijders, and Valente (2018) for more details about this
#' measure.
#' @details
#' Following Fujimoto, Snijders, and Valente (2018: 380), the number
#' of homophilous four-cycles for actor *i* is:
#' \deqn{ \sum_{j} \sum_{a\neq b} y_{ia}y_{ib}y_{ja}y_{jb}I{{v_{i} = v_{j}}}}
#' where *y* is the two-mode adjacency matrix, *v* is the vector of
#' membership scores (e.g., sports/club membership), *a* and *b* represent
#' the level two groups, and \eqn{I{v_i = v_j}} is the indicator function that
#' is 1 if the values are the same and 0 if not.
#'
#'
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#' @references
#' Fujimoto, Kayo, Tom A.B. Snijders, and Thomas W. Valente. 2018. "Multivariate
#' dynamics of one-mode and two-mode networks: Explaining similarity in sports
#' participation among friends." *Network Science* 6(3): 370-395.
#'
#' @examples
#'
#' # For this example, we use the Davis Southern Women's Dataset.
#' data("southern.women")
#' #creating a random binary membership vector
#' set.seed(9999)
#' membership <- sample(0:1, nrow(southern.women), replace = TRUE)
#' #the homophilous four-cycle values
#' computeHomFourCycles(southern.women, mem = membership)


#### Formula for Calculating Two-Mode Four-Cycle Homophily Effect
##
##     Σj  Σa != b   Yia Yib Yja Yjb  I{Vi = Vj}
##
##     In which I is the indicator function defined as 1 if i and j share the same
##     categorical membership, and 0 elsewhere. Y represents the two-mode affliation
##     matrix, and V represents a vector of categorical membership values.
##
####
computeHomFourCycles <- function(net, # affilation matrix
                            mem # membership homophily
                            ) {
  lifecycle::deprecate_warn("1.0.0", " computeHomFourCycles()", "netstats_tm_homfourcycles()")
  pbypmat <- (net %*% t(net)) # getting the nxn matrix of shared affilation ties
  diag(pbypmat) <- 0 # changing the diagonal to 0
  pbypmat[pbypmat < 2] <- 0 # if the value is less than 2, then 0, i.e., no four cycles possible
  same_matrix <- matrix(0, nrow(net), nrow(net)) # creating a nxn matrix of same membership ties
  for (i in 1:nrow(same_matrix)) { # for all actors
    same_matrix[i, ] <- ifelse(mem == mem[i], 1, 0) # make the network ties to be 1 if same membership and 0 if not
  }
  pbypmat2 <- pbypmat - 1 # (n - 1) possible network ties
  pbypmat2[pbypmat2 < 0] <- 0 # if the value is negative, turn back to 0
  pbypmat3 <- pbypmat * pbypmat2 # n*(n-1) essentially all possible four cycle combinations
  pbypmat3 <- pbypmat3 / 2
  homophilycycles <- pbypmat3 * same_matrix # multiplying the cycles by the same-shared homophily value
  numcycles <- rowSums(homophilycycles) # the number of cycles for each network actor

  ifelse(is.null(rownames(net)),
         names( numcycles)  <- NULL,
         names( numcycles) <- rownames(net))

  return(numcycles) # return the number of homophily four cycles
}
