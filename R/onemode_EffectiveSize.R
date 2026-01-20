##### Code Written to Compute Burt's Effective Size and Constraint
##### Code written by Kevin Carson and Diego Leal
##### Last Updated: 04-18-24

#' @title Compute Burt's (1992) Effective Size for Ego Networks from a Sociomatrix
#' @name netstats_om_effective
#' @param net The one-mode sociomatrix with network ties.
#' @param isolates The numerical value that represents what value will isolates be given. Set to NA by default.
#' @param pendants The numerical value that represents what value will pendant vertices be given. Set to 1 by default. Pendant vertices are those nodes who have one outgoing tie.
#' @param inParallel TRUE/FALSE. TRUE indicates that parallel processing will be used to compute the statistic with the *foreach* package. FALSE indicates that parallel processing will not be used. Set to FALSE by default.
#' @param nCores If inParallel = TRUE, the number of computing cores for parallel processing. If this value is not specified, then the function internally provides it by dividing the number of available cores in half.
#' @return The vector of ego network effective size values.
#' @import foreach
#' @import doParallel
#' @import parallel
#' @import Rcpp
#' @export
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' This function computes Burt's (1992) one-mode ego effective size based upon a sociomatrix (see details).
#' @details The formula for Burt's (1992; see also Borgatti 1997) one-mode ego effective size is:
#' \deqn{ E_{i} = \sum_{j} 1 - \sum_{q}p_{iq}m_{jq} ; q \neq i \neq j}
#' where \eqn{E_{i}} is the ego effective size for an ego *i*.
#' \eqn{p_{iq}} is formulated as:
#' \deqn{\frac{(z_{iq} + z_{qi}) }{\sum_{j}(z_{ij} + z_{ji})} ;  i \neq j}
#' and \eqn{m_{jq}} is:
#' \deqn{m_{jq} = \frac{(z_{jq} + z_{qj})}{max(z_{jk} + z_{kj})}}
#'
#' While this function internally locates isolates (i.e., nodes
#' who have no ties) and pendants (i.e., nodes who only have
#' one tie), the user should specify what values for constraint are returned for them via the *isolates* and
#' *pendants* options. In particular, pendant vertices are those nodes who have one outgoing tie.
#'
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#' @references
#' Burt, Ronald. 1992. *Structural Holes: The Social Structure of Competition*.
#' Harvard University Press.
#'
#' Borgatti, Stephen. 1997. "Structural Holes: Unpacking Burt's Redundancy Measures." *Connections* 20(1): 35-38.
#'
#' @examples
#' # For this example, we recreate the ego network provided in Borgatti (1997):
#' BorgattiEgoNet <- matrix(
#'  c(0,1,0,0,0,0,0,0,1,
#'    1,0,0,0,0,0,0,0,1,
#'    0,0,0,1,0,0,0,0,1,
#'    0,0,1,0,0,0,0,0,1,
#'    0,0,0,0,0,1,0,0,1,
#'   0,0,0,0,1,0,0,0,1,
#'   0,0,0,0,0,0,0,1,1,
#'    0,0,0,0,0,0,1,0,1,
#'    1,1,1,1,1,1,1,1,0),
#'  nrow = 9, ncol = 9, byrow = TRUE)
#'colnames(BorgattiEgoNet) <- rownames(BorgattiEgoNet) <- c("A", "B", "C",
#'                                                          "D", "E", "F",
#'                                                         "G", "H", "ego")
#'#the effective size value for the ego replicates that provided in Borgatti (1997)
#'netstats_om_effective(BorgattiEgoNet)
#'
#' # For this example, we recreate the ego network provided in Burt (1992: 56):
#'BurtEgoNet <- matrix(c(
#'   0,1,0,0,1,1,1,
#'  1,0,0,1,0,0,1,
#'  0,0,0,0,0,0,1,
#'  0,1,0,0,0,0,1,
#'  1,0,0,0,0,0,1,
#'  1,0,0,0,0,0,1,
#'  1,1,1,1,1,1,0),
#'  nrow = 7, ncol = 7)
#' colnames(BurtEgoNet) <- rownames(BurtEgoNet) <- c("A", "B", "C", "D", "E",
#'                                                  "F", "ego")
#'#the effective size value for the ego replicates that provided in Burt (1992: 56)
#'netstats_om_effective(BurtEgoNet)


#################################################################################
###     Burt's Effective Size Measure for Ego Networks
#################################################################################

netstats_om_effective <- function(net, # the full sociomatrix
                  inParallel = FALSE, # should this be computed in parallel?
                  nCores = NULL, # the number of cores for computing in parallel?
                  isolates = NA, # what value should isolates get?
                  pendants = 1) { # what value should be given to pendant vertices?
  #-----------------#
  # Checking if colnames exist
  #-----------------#
  if(is.null(colnames(net))){ #if colnames do not exist name the matrix
    colnames(net)<-rownames(net)<-paste("ego",1:nrow(net),sep = "") #naming the egos9
  }

  if (inParallel == FALSE) {
    effectivesize <- rep(0, nrow(net)) # creating an empty vector to store constraint values for any ego i

    for(i in 1:nrow(net)){
      ###---------------------------------------------------------------###
      #       Extracting the ego networks
      ###---------------------------------------------------------------###
      egos_out <- which(net[i, ] == 1) # All Alters were Netij = 1  ( outdegree for ego i )
      egos_in <- which(net[, i] == 1) # All Alters were Netji = 1 ( indegree for ego i )
      alters <- c(i, unique(c(egos_out, egos_in))) # combining the unique cases
      egonet <- net[alters, alters] # creating the ego network

      ###---------------------------------------------------------------###
      #       Computing the effective size for the network
      ###---------------------------------------------------------------###
      if(length(egonet) == 1){ # if they are a isolate (i.e., no ties)
        effectivesize[i] <- isolates  # return the value for isolates
        next #skip to the next person
      }
      if(nrow(egonet) < 3){ # if they are a pendant (i.e., only one tie)
        effectivesize[i] <- pendants  # return the value for pendants
        next #skip to the next person
      }
      effectivesize[i] <- burteffective(net = egonet,
                                           nactors = nrow(egonet)) #compute the value for effective size
    }
    ifelse(is.null(rownames(net)),
           names(effectivesize)  <- NULL,
           names(effectivesize) <- rownames(net))

   return(effectivesize) # Return the effective size values for the network

  }else{
      if (is.null(nCores)) {
        nCores <- round(parallel::detectCores()/2)
      }
      myCluster <- parallel::makeForkCluster(nnodes = nCores)
      doParallel::registerDoParallel(myCluster)
      effectiveSize <- foreach::foreach(i = 1:nrow(net), .combine = c) %dopar%
        ({
          egos_out <- which(net[i, ] == 1)
          egos_in <- which(net[, i] == 1)
          ego_centric <- c(i, unique(c(egos_out, egos_in)))
          Z <- net[unique(ego_centric), unique(ego_centric)]
          if (length(Z) == 1) {
            return(isolates)
          }
          if (nrow(Z) < 3) {
            return(pendants)
          }
          burteffective(net = Z, nactors = nrow(Z))
          })
      parallel::stopCluster(myCluster)
      ifelse(is.null(rownames(net)), names(effectivesize) <- NULL,
             names(effectivesize) <- rownames(net))
      return(effectiveSize)
    }
}
