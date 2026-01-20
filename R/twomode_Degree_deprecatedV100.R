## Density of Two-Mode Networks
## see Wasserman, Stanley and Katherine Faust. 1995.
## Code written by Kevin Carson (kacarson@arizona.edu) and Deigo Leal (https://www.diegoleal.info/)
## Last Updated: 09-07-24

#' @title Compute Degree Centrality Values for Two-Mode Networks
#' @name computeTMDegree
#' @param net A two-mode adjacency matrix
#' @param level1 TRUE/FALSE. TRUE indicates that the degree centrality will be computed for level 1 nodes. FALSE indicates that the degree centrality will be computed for level 2 nodes. Set to TRUE by default.
#' @return The vector of two-mode level-specific degree centrality values.
#' @export
#'
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `computeTMDegree()` has been deprecated starting on version 1.0.0 of the `dream` package. Please use the `netstats_tm_degreecent()` function and see the `NEWS.md` file for more details.
#'
#'
#'
#' This function computes the degree centrality values for two-mode
#' networks following Knoke and Yang (2020). The computed degree centrality is
#' based on the specified level. That is, in an affiliation matrix, the density
#' can be computed on the symmetric *g x g* co-membership matrix of
#' level 1 actors (e.g., medical doctors) or on the symmetric *h x h* shared actors matrix for level 2
#' groups (e.g., hospitals) based on their shared members.
#'
#' @details
#' Following Knoke and Yang (2020), the computation of degree for two-mode affiliation networks
#' is level specific. A two-mode affiliation matrix *X* with dimensions *g x h*, where *g* is
#' the number of level 1 nodes (e.g., medical doctors) and *h* is the number of level 2 nodes
#' (i.e., hospitals). If the function is defined on the level 1 nodes,
#' the degree centrality of an actor *i* is computed as:
#' \deqn{ X^{G} = XX^{T} }
#' \deqn{ C_{D}^{G}(g_{i}) = \sum_{i = 1}^{g}  x_{ij}^{g} \quad (i \neq j)}
#' In contrast, if it is defined on the level 2 nodes, the degree centrality of
#' an actor *i* is computed as:
#' \deqn{ X^{H} = X^{T}X }
#' \deqn{ C_{D}^{H}(h_{i}) = \sum_{i = 1}^{h}  x_{ij}^{h} \quad (i \neq j)}
#'
#'
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#' @references
#' Knoke, David and Song Yang. 2020. *Social Network Analysis*. Sage: Quantitative
#' Applications in the Social Sciences (154)

#' @examples
#' #Replicating the biparitate graph presented in Knoke and Yang (2020: 109)
#' knoke_yang_PC <- matrix(c(1,1,0,0, 1,1,0,0,
#'                           1,1,1,0, 0,0,1,1,
#'                           0,0,1,1), byrow = TRUE,
#'                           nrow = 5, ncol = 4)
#'colnames(knoke_yang_PC) <- c("Rubio-R","McConnell-R", "Reid-D", "Sanders-D")
#'rownames(knoke_yang_PC) <- c("UPS", "MS", "HD", "SEU", "ANA")
#'computeTMDegree(knoke_yang_PC, level1 = TRUE) #this value matches the book
#'computeTMDegree(knoke_yang_PC, level1 = FALSE) #this value matches the book

computeTMDegree <- function(net,#a two-mode network adjancency matrix
                      level1 = TRUE){#Boolean: TRUE indicating if the density should be
  #computed on level 1 nodes, FALSE computes graph density for the level 2 nodes
  lifecycle::deprecate_warn("1.0.0", " computeTMDegree()", "netstats_tm_degreecent()")
  if(level1 == TRUE){ #should level 1 density be computed?
    xA <- net %*% t(net) #if so, get the level 1 x level 1 matrix
  }else{ #if not
    xA <- t(net) %*% net #get the level 2 x level 2 matrix
  }
  diag(xA) <- 0 #making the diagonal zero
  xA[xA > 1] <- 1
  degree <- rowSums(xA) #the degree of each actor in the specific network

  if(level1 == TRUE){
    ifelse(is.null(rownames(net)),
           names(degree) <- NULL,
           names(degree) <- rownames(net))
  }else{
    ifelse(is.null(colnames(net)),
           names(degree) <- NULL,
           names(degree) <-colnames(net))
  }
  return(degree) #returning the network
}
