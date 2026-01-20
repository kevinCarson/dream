#######
# Computing the number of paths of length k between two actors in a graph



#' @title Compute the Number of Walks of Length K in a One-Mode Network
#' @name computeNPaths
#' @param net An unweighted one-mode network adjacency matrix.
#' @param k A numerical value that corresponds to the length of the paths to be computed.
#' @return An *n* x *n* matrix of counts of walks.
#' @export
#'
#'
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `computeNPaths()` has been deprecated starting on version 1.0.0 of the `dream` package. Please use the `netstats_om_nwalks()` function and see the `NEWS.md` file for more details.
#'
#'
#'
#' This function calculates the number of walks of length *k*
#' between any two vertices in an unweighted one-mode network.
#'
#' @details
#' A nice result from graph theory is that the number of walks of length *k* between
#' vertices *i* and *j* can be found by:
#' \deqn{ A_{ij}^k }
#'
#' This function is similar to the functions provided in *igraph* that provide the walks between
#' two vertices. The main difference is that this function provides the counts
#' of paths between all vertices in the network. In addition, this function assumes that
#' there are no self-loops (i.e., the diagonal of the matrix is 0).
#'
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#'
#'
#' @examples
#'
#' # For this example, we generate a random one-mode graph with the sna package.
#' #creating the random network with 10 actors
#' set.seed(9999)
#' rnet <- matrix(sample(c(0,1), 10*10, replace = TRUE, prob = c(0.8,0.2)),
#'                nrow = 10, ncol = 10, byrow = TRUE)
#' diag(rnet) <- 0 #setting self ties to 0
#' #counting the paths of length 2
#' computeNPaths(rnet, k = 2)
#' #counting the paths of length 5
#' computeNPaths(rnet, k = 5)

computeNPaths <- function(net, # the network adjacency matrix
                   k){ #the requested path length
  lifecycle::deprecate_warn("1.0.0", "computeNPaths()", "netstats_om_nwalks()")
  diag(net) <- 0 #ensuring there are no self-loops encase the user misspecificed
  net1 <- net  #making a copy of the original matrix
  for(i in 1:(k-1)){ #for i in the requested path length (minus 1)
                     #if k = 2, for example, matrix mulitiplication should
                     #only be done once
    net <- net %*% net1 #multiplying the matrix by the original
  }
  return(net) #returning the n x n
}











