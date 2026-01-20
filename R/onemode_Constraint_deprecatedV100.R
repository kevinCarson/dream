##### Code Written to Compute Burt's Effective Size and Constraint
##### Code written by Kevin Carson and Diego Leal
##### Last Updated: 09-09-24

#' @title Compute Burt's (1992) Constraint for Ego Networks from a Sociomatrix
#' @name computeBurtsConstraint
#' @param net A one-mode sociomatrix with network ties.
#' @param isolates What value should isolates be given? Set to NA by default.
#' @param pendants What value should be given to pendant vertices? Pendant vertices are those nodes who have one outgoing tie. Set to 1 by default.
#' @param inParallel TRUE/FALSE. TRUE indicates that parallel processing will be used to compute the statistic with the *foreach* package. FALSE indicates that parallel processing will not be used. Set to FALSE by default.
#' @param nCores If inParallel = TRUE, the number of computing cores for parallel processing. If this value is not specified, then the function internally provides it by dividing the number of available cores in half.
#' @return The vector of ego network constraint values.
#' @import foreach
#' @import doParallel
#' @import parallel
#' @export

#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `computeBurtsConstraint()` has been deprecated starting on version 1.0.0 of the `dream` package. Please use the `netstats_om_constraint()` function and see the `NEWS.md` file for more details.
#'
#'
#' This function computes Burt's (1992) one-mode ego constraint based upon a sociomatrix.
#' @details The formula for Burt's (1992) one-mode ego constraint is:
#' \deqn{c_{ij} = \left(p_{ij} + \sum_{q} p_{iq} p_{qj}\right)^2 \quad ; \; q \neq i \neq j}
#' where:
#' \itemize{
#'   \item \eqn{p_{iq}} is formulated as: \eqn{p_{iq} = \frac{z_{iq} + z_{qi}}{\sum_{j}(z_{ij} + z_{ji})} \quad ; \; i \neq j}
#' }
#' Finally, the aggregate constraint of an ego *i* is:
#' \deqn{C_{i} = \sum_{j} c_{ij}}
#' While this function internally locates isolates (i.e., nodes
#' who have no ties) and pendants, that is, nodes who only have
#' one outgoing tie, the user should specify what values for constraint are returned for them via the *isolates* and
#' *pendants* options.
#'
#' Lastly, this function allows users to compute the values in parallel via the
#' *foreach*, *doParallel*, and *parallel* R packages.
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#' @references
#' Burt, Ronald. 1992. *Structural Holes: The Social Structure of Competition*.
#' Harvard University Press.

#' @examples
#'
#' # For this example, we recreate the ego network provided in Burt (1992: 56):
#' BurtEgoNet <- matrix(c(
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
#' #the constraint value for the ego replicates that provided in Burt (1992: 56)
#' computeBurtsConstraint(BurtEgoNet)
#'
#'

#################################################################################
###     Burt's Constraint Measure for Ego Networks
#################################################################################

computeBurtsConstraint <- function(net, # the full sociomatrix
                  inParallel = FALSE, # should this be computed in parallel?
                  nCores = NULL, # the number of cores for computing in parallel?
                  isolates = NA, # what value should isolates get?
                  pendants = 1) { # what value should be given to pendant vertices?

  lifecycle::deprecate_warn("1.0.0", "computeBurtsConstraint()", "netstats_om_constraint()")
  #-----------------#
  # Checking if colnames exist
  #-----------------#
  if(is.null(colnames(net))){ #if colnames do not exist name the matrix
    colnames(net)<-rownames(net)<-paste("ego",1:nrow(net),sep = "") #naming the egos9
  }

  n_actors <- nrow(net) # the number of actors in the network
  if (inParallel == FALSE) {
    constraint <- rep(0, nrow(net)) # creating an empty vector to store constraint values for any ego i

    egoNets <- list() # creating an empty list to store all ego Networks

    for (i in 1:nrow(net)) {
      egos_out <- which(net[i, ] == 1) # All Alters were Netij = 1  ( outdegree for ego i )
      egos_in <- which(net[, i] == 1) # All Alters were Netji = 1 ( indegree for ego i )
      ego_centric <- c(i, unique(c(egos_out, egos_in))) # combining the unique cases
      egoNets[[i]] <- net[unique(ego_centric), unique(ego_centric)] # creating the ego network
    }

    for (i in 1:length(constraint)) { # for all actors in the network
      neti <- egoNets[[i]] # set the network to the current ego i network
      if (length(neti) == 1) {
        constraint[i] <- isolates
        next
      } # if the ego i is an isolate: constraint = 0,
      # then go to the next ego
      if (nrow(neti) < 3) {
        constraint[i] <- pendants
        next
      } # if the ego i has only one alter: constraint = 0, ,
      # then go to the next ego, ego with one alter cannot be constrained
      Z <- neti + t(neti) # connection zij is made symmetric prior to computing pij
      sumzq <- rowSums(Z) # the outdegree of ego i
      pij <- Z / sumzq # pij: proportion of time ego i spends with any alter j
      pij[is.na(pij)] <- 0 # values divided by 0 go to zero [this acts as a check for isolates in the network with 0 degree]
      psumq <- pij %*% pij # Sum Q: that is, the proportion i and j spend with any alter q
      diag(psumq) <- 0 # Making the diagonal zero, that is an ego or alter cannot be constrained in relation to themselves
      Ci <- rowSums((pij + psumq)^2) # Following Burt: Constraint C = Sum(pij + Sum pij * pij)^2
      constraint[i] <- (round(as.numeric(Ci), 4))[1] # rounding the constraint for the ego and adding it to the constraint value vector
    }

    ifelse(is.null(rownames(net)),
           names(constraint)  <- NULL,
           names(constraint) <- rownames(net))

    return(constraint) # returning the constraint vector for the full network
  }


  if (inParallel == TRUE) {
    if (is.null(nCores)) { # do this user provide the number of cores to use?
      nCores <- round(parallel::detectCores() / 2) # if not, detect the number of cores, and use half of them
    }
    myCluster <- parallel::makeForkCluster(nnodes = nCores) # creating the cluster
    doParallel::registerDoParallel(myCluster) # registering the cluster

    # using the foreach package to do the computations in parallel
    constraint <- foreach::foreach(i = 1:n_actors, .combine = c) %dopar%

      ({
        egos_out <- which(net[i, ] == 1) # All Alters were Netij = 1  ( outdegree for ego i )
        egos_in <- which(net[, i] == 1) # All Alters were Netji = 1 ( indegree for ego i )
        ego_centric <- c(i, unique(c(egos_out, egos_in))) # combining the unique cases
        neti <- net[unique(ego_centric), unique(ego_centric)] # creating the ego network
        if (length(neti) == 1) {
          return(isolates)
        } # if the ego i is an isolate: constraint = 0,
        # then go to the next ego
        if (nrow(neti) < 3) {
          return(pendants)
        } # if the ego i has only one alter: constraint = 0, ,
        # then go to the next ego, ego with one alter cannot be constrained
        Z <- neti + t(neti) # connection zij is made symmetric prior to computing pij
        sumzq <- rowSums(Z) # the outdegree of ego i
        pij <- Z / sumzq # pij: proportion of time ego i spends with any alter j
        pij[is.na(pij)] <- 0 # values divided by 0 go to zero [this acts as a check for isolates in the network with 0 degree]
        psumq <- pij %*% pij # Sum Q: that is, the proportion i and j spend with any alter q
        diag(psumq) <- 0 # Making the diagonal zero, that is an ego or alter cannot be constrained in relation to themselves
        Ci <- rowSums((pij + psumq)^2) # Following Burt: Constraint C = Sum(pij + Sum pij * pij)^2
        (round(as.numeric(Ci), 4))[1] # rounding the constraint for the ego and adding it to the constraint value vector
      })

    parallel::stopCluster(myCluster) # closing the cluster

    ifelse(is.null(rownames(net)),
           names(constraint)  <- NULL,
           names(constraint) <- rownames(net))

    return(constraint)
  }
}
