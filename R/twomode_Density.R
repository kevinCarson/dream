## Density of Two-Mode Networks
## see Wasserman, Stanley and Katherine Faust. 1995.
## Code written by Kevin Carson (kacarson@arizona.edu) and Deigo Leal (https://www.diegoleal.info/)
## Last Updated: 09-07-24

#' @title Compute Level-Specific Graph Density for Two-Mode Networks
#' @name netstats_tm_density
#' @param net A two-mode adjacency matrix.
#' @param binary TRUE/FALSE. TRUE indicates that the transposed matrices will be binarized (see Wasserman and Faust 1995: 316). FALSE indicates that the transposed matrices will not be binarized. Set to FALSE by default.
#' @param level1 TRUE/FALSE. TRUE indicates that the graph density will be computed for level 1 nodes. FALSE indicates that the graph density will be computed for level 2 nodes. Set to FALSE by default.
#' @return The level-specific network density for the two-mode graph.
#' @export
#'
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' This function computes the density of a two-mode network following
#' Wasserman and Faust (1994) and Knoke and Yang (2020). The density is computed
#' based on the specified level. That is, in an affiliation matrix, density
#' can be computed on the symmetric *g x g* matrix of co-membership for the
#' level 1 actors or on the symmetric *h x h* matrix of shared actors for level 2
#' groups.
#'
#' @details
#' Following Wasserman and Faust (1994) and Knoke and Yang (2020), the computation
#' of density for two-mode networks is level specific. A two-mode matrix *X* with
#' dimensions *g x h*, where *g* is the number of level 1 nodes (e.g., medical doctors)
#' and *h* is the number of level 2 nodes (i.e., hospitals). If
#' the function is defined on the level 1 nodes, the density is computed as:
#'
#' \deqn{ X^{g} = XX^{T} }
#' \deqn{ D^{g} = \frac{\sum_{i = 1}^{g}\sum_{j = 1}^{g} x_{ij}^{g} }{g(g-1)}}
#' In contrast, if it is defined on the level 2 nodes, the density is:
#' \deqn{ X^{h} = X^{T}X }
#' \deqn{ D^{h} = \frac{\sum_{i = 1}^{h}\sum_{j = 1}^{h} x_{ij}^{h} }{h(h-1)}}
#'
#' Moreover, as discussed in Wasserman and Faust (1994: 316), the density can be
#' based on the dichotomous relations instead of the shared membership values.
#' This can be specified by *binary = TRUE*.
#'
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#' @references
#' Wasserman, Stanley and Katherine Faust. 1994. *Social Network Analysis: Methods
#' and Applications*. Cambridge University Press.
#'
#' Knoke, David and Song Yang. 2020. *Social Network Analysis*. Sage: Quantitative
#' Applications in the Social Sciences (154).
#'
#' @examples
#' #Replicating the biparitate graph presented in Knoke and Yang (2020: 109)
#' knoke_yang_PC <- matrix(c(1,1,0,0, 1,1,0,0,
#'                           1,1,1,0, 0,0,1,1,
#'                           0,0,1,1), byrow = TRUE,
#'                           nrow = 5, ncol = 4)
#'colnames(knoke_yang_PC) <- c("Rubio-R","McConnell-R", "Reid-D", "Sanders-D")
#'rownames(knoke_yang_PC) <- c("UPS", "MS", "HD", "SEU", "ANA")
#'#compute two-mode density for level 1
#'#note: this value does not match that of Knoke and Yang (which we believe
#'#is a typo in that book), but does match that of Wasserman and
#'#Faust (1995: 317) for the ceo dataset.
#'netstats_tm_density(knoke_yang_PC, level1 = TRUE)
#'#compute two-mode density for level 2.
#'#note: this value matches that of the book
#'netstats_tm_density(knoke_yang_PC, level1 = FALSE)
#'

netstats_tm_density <- function(net,#a two-mode network adjancency matrix
                    binary = FALSE, #Boolean: TRUE indicating if the tranposed matrices
                                    #should be binarized (see Wasserman and Faust 199*: 316)
                    level1 = TRUE){ #Boolean: TRUE indicating if the density should be
                                    #computed on level 1 nodes, FALSE computes
                                    #graph density for the level 2 nodes

  if(level1 == TRUE){ #should level 1 density be computed?
    xA <- net %*% t(net) #if so, get the level 1 x level 1 matrix
  }else{ #if not
    xA <- t(net) %*% net #get the level 2 x level 2 matrix
  }
  diag(xA) <- 0 #making the diagonal 0
  denom <- nrow(xA)*(nrow(xA) - 1) #the number of possible ties
  if(binary == TRUE){ #if the density should be binarized
    xA[xA > 1] <- 1 #if the values are greater than 1, make that value 1
  }
  num <- sum(xA) #the number of co-membership ties
  gdens <- num / denom #the density formula
  return(gdens) #return the density values
}

