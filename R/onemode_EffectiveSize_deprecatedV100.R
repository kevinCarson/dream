##### Code Written to Compute Burt's Effective Size and Constraint
##### Code written by Kevin Carson and Diego Leal
##### Last Updated: 04-18-24

#' @title Compute Burt's (1992) Effective Size for Ego Networks from a Sociomatrix
#' @name computeBurtsES
#' @param net The one-mode sociomatrix with network ties.
#' @param isolates The numerical value that represents what value will isolates be given. Set to NA by default.
#' @param pendants The numerical value that represents what value will pendant vertices be given. Set to 1 by default. Pendant vertices are those nodes who have one outgoing tie.
#' @param inParallel TRUE/FALSE. TRUE indicates that parallel processing will be used to compute the statistic with the *foreach* package. FALSE indicates that parallel processing will not be used. Set to FALSE by default.
#' @param nCores If inParallel = TRUE, the number of computing cores for parallel processing. If this value is not specified, then the function internally provides it by dividing the number of available cores in half.
#' @return The vector of ego network effective size values.
#' @import foreach
#' @import doParallel
#' @import parallel
#' @export
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `computeBurtsES()` has been deprecated starting on version 1.0.0 of the `dream` package. Please use the `netstats_om_effective()` function and see the `NEWS.md` file for more details.
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
#' *pendants* options. Pendant vertices are those nodes who have one outgoing tie.
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
#'computeBurtsES(BorgattiEgoNet)
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
#'computeBurtsES(BurtEgoNet)


#################################################################################
###     Burt's Effective Size Measure for Ego Networks
#################################################################################

computeBurtsES <- function(net, # the full sociomatrix
                  inParallel = FALSE, # should this be computed in parallel?
                  nCores = NULL, # the number of cores for computing in parallel?
                  isolates = NA, # what value should isolates get?
                  pendants = 1) { # what value should be given to pendant vertices?

  lifecycle::deprecate_warn("1.0.0", "computeBurtsES()", "netstats_om_effective()")
  #-----------------#
  # Checking if colnames exist
  #-----------------#
  if(is.null(colnames(net))){ #if colnames do not exist name the matrix
    colnames(net)<-rownames(net)<-paste("ego",1:nrow(net),sep = "") #naming the egos9
  }
  n_actors <- nrow(net) # the number of actors in the network
  if (inParallel == FALSE) {
    effectivesize <- rep(0, nrow(net)) # creating an empty vector to store constraint values for any ego i

    egoNets <- list() # creating an empty list to store all ego Networks
    for (i in 1:nrow(net)) { # for all actors in the network (i.e., all egos)
      egos_out <- which(net[i, ] == 1) # All Alters were Netij = 1  ( outdegree for ego i )
      egos_in <- which(net[, i] == 1) # All Alters were Netji = 1 ( indegree for ego i )
      ego_centric <- c(i, unique(c(egos_out, egos_in))) # combining the unique cases
      egoNets[[i]] <- net[unique(ego_centric), unique(ego_centric)] # creating the ego network
    }

    for (i in 1:length(effectivesize)) { # for all actors in the network
      Z <- egoNets[[i]] # set the Z matrix to the current ego i network
      if (length(Z) == 1) {
        effectivesize[i] <- isolates
        next
      } # if the ego i is an isolate: effective size = 0,
      # then go to the next ego
      if (nrow(Z) < 3) {
        effectivesize[i] <- pendants
        next
      } # if the ego i has only one alter: effective size = 1,
      # then go to the next ego [alters have no redundancy]
      pij1 <- Z + t(Z) # connection zij is made symmetric prior to computing pij [Zij + Zji]
      pij <- pij1 / rowSums(pij1) # the redundancy for all actors: Sum(Zij + Zji) / Total Outdegree
      mjq1 <- Z + t(Z) # connection zij is made symmetric prior to computing mjq [Zij + Zji]
      mjq2 <- apply(mjq1, 1, max) # getting the maximum energy that any alter spends with any other alter max(Zij + Zji)
      mjq <- mjq1 / mjq2 # Zij + Zji / max(Zij + Zji)  ; Burt's marginal strength measure
      sumq <- pij %*% t(mjq) #  pij * mjq   (total energy spent by ego i) * marginal energy expenditure for alter pairs
      sumqego <- sumq[1, 2:ncol(sumq)] # The sum of the redundancy in ego i's network
      effectivesize[i] <- (sum(1 - sumqego)) # Following Burt, Effective size for any ego i, is
      # Sum[ 1 - rij] for all alters j in ego i's network
    }
    ifelse(is.null(rownames(net)),
           names(effectivesize)  <- NULL,
           names(effectivesize) <- rownames(net))

   return(effectivesize) # Return the effective size values for the network
  }

  if (inParallel == TRUE) {
    if (is.null(nCores)) { # do this user provide the number of cores to use?
      nCores <- round(parallel::detectCores() / 2) # if not, detect the number of cores, and use half of them
    }
    myCluster <- parallel::makeForkCluster(nnodes = nCores) # creating the cluster
    doParallel::registerDoParallel(myCluster) # registering the cluster

    # using the foreach package to do the computations in parallel
    effectiveSize <- foreach::foreach(i = 1:n_actors, .combine = c) %dopar%

      ({
        egos_out <- which(net[i, ] == 1) # All Alters were Netij = 1  ( outdegree for ego i )
        egos_in <- which(net[, i] == 1) # All Alters were Netji = 1 ( indegree for ego i )
        ego_centric <- c(i, unique(c(egos_out, egos_in))) # combining the unique cases
        Z <- net[unique(ego_centric), unique(ego_centric)] # creating the ego network
        if (length(Z) == 1) {
          return(isolates)
        } # if the ego i is an isolate: effective size = 0,
        # then go to the next ego
        if (nrow(Z) < 3) {
          return(pendants)
        } # if the ego i has only one alter: effective size = 1,
        # then go to the next ego [alters have no redundancy]
        pij1 <- Z + t(Z) # connection zij is made symmetric prior to computing pij [Zij + Zji]
        pij <- pij1 / rowSums(pij1) # the redundancy for all actors: Sum(Zij + Zji) / Total Outdegree
        mjq1 <- Z + t(Z) # connection zij is made symmetric prior to computing mjq [Zij + Zji]
        mjq2 <- apply(mjq1, 1, max) # getting the maximum energy that any alter spends with any other alter max(Zij + Zji)
        mjq <- mjq1 / mjq2 # Zij + Zji / max(Zij + Zji)  ; Burt's marginal strength measure
        sumq <- pij %*% t(mjq) #  pij * mjq   (total energy spent by ego i) * marginal energy expenditure for alter pairs
        sumqego <- sumq[1, 2:ncol(sumq)] # The sum of the redundancy in ego i's network
        (sum(1 - sumqego)) # Following Burt, Effective size for any ego i, is
        # Sum[ 1 - rij] for all alters j in ego i's network
      })

    parallel::stopCluster(myCluster) # closing the cluster
    ifelse(is.null(rownames(net)),
           names(effectivesize)  <- NULL,
           names(effectivesize) <- rownames(net))

    return(effectiveSize)
  }
}
