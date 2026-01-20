## The Computation of Two-Mode Redundancy in Two-Mode Static Network
## See Bruchard and Cornwell 2018 for the two-mode network
## Code written by Kevin Carson (kacarson@arizona.edu) and Diego Leal (https://www.diegoleal.info/)
## Last Updated: 10-29-24

#' @title Compute Burchard and Cornwell's (2018) Two-Mode Redundancy
#' @name netstats_tm_redundancy
#' @param net A two-mode adjacency matrix or affiliation matrix.
#' @param isolates What value should isolates be given? Preset to be NA.
#' @param weighted TRUE/FALSE. TRUE indicates the resulting statistic will be based on the weighted formula (see the details section). FALSE indicates the statistic will be based on the original non-weighted formula. Set to FALSE by default.
#' @param inParallel TRUE/FALSE. TRUE indicates that parallel processing will be used to compute the statistic with the *foreach* package. FALSE indicates that parallel processing will not be used. Set to FALSE by default.
#' @param nCores If inParallel = TRUE, the number of computing cores for parallel processing. If this value is not specified, then the function internally provides it by dividing the number of available cores in half.
#' @return An *n x n* matrix with level 1 redundancy scores for actors in a two-mode network.
#' @import foreach
#' @import doParallel
#' @import parallel
#' @export


#' @description
#' `r lifecycle::badge("stable")`
#'
#' This function calculates the values for two mode redundancy
#' for weighted and unweighted two-mode networks based on Burchard and Cornwell (2018).
#' @details The formula for two-mode redundancy is:
#' \deqn{r_{ij} = \frac{|\sigma(j) \cap \sigma(i)|}{|\sigma(i)|}}
#' where:
#' \itemize{
#'   \item \eqn{r_{ij}} is the redundancy of ego *i* with respect to actor *j*.
#'   \item \eqn{|\sigma(j) \cap \sigma(i)|} is the number of same-class contacts (e.g., medical doctors in a hospital) that *i* and *j* both share.
#'   \item \eqn{|\sigma(i)|} is the number of same-class contacts of ego *i*.
#' }
#' The two-mode redundancy is ego-bound, that is, the redundancy is only based on the
#' two-mode ego network of *i*. Put differently, \eqn{r_{ij}} only considers the perspective of the ego.
#' This function allows the user to compute the scores in parallel through the *foreach* and *doParallel* R packages.
#' If the matrix is weighted, the user should specify *weighted = TRUE*. Following Burchard and Cornwell (2018),
#' the formula for two-mode weighted redundancy is:
#' \deqn{r_{ij} = \frac{|\sigma(j) \cap \sigma(i)|}{|\sigma(i)| \times w_t}}
#' where \eqn{w_t} is the average of the tie weights that *i* and *j* send
#' to their shared opposite class contacts.
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#' @references
#' Burchard, Jake and Benjamin Cornwell. 2018. "Structural Holes and bridging
#' in two-mode networks." *Social Networks* 55:11-20.
#'
#' @examples
#'
#' # For this example, we recreate Figure 2 in Burchard and Cornwell (2018: 13)
#'BCNet <- matrix(
#'  c(1,1,0,0,
#'    1,0,1,0,
#'    1,0,0,1,
#'    0,1,1,1),
#'  nrow = 4, ncol = 4, byrow = TRUE)
#'colnames(BCNet) <- c("1", "2", "3", "4")
#'rownames(BCNet) <- c("i", "j", "k", "m")
#'#this values replicate those reported by Burchard and Cornwell (2018: 14)
#'netstats_tm_redundancy(BCNet)
#'
#'
#'#For this example, we recreate Figure 9 in Burchard and Cornwell (2018:18)
#'#for weighted two mode networks.
#'BCweighted <- matrix(c(1,2,1, 1,0,0,
#'                       0,2,1,0,0,1),
#'                       nrow = 4, ncol = 3,
#'                       byrow = TRUE)
#'rownames(BCweighted) <- c("i", "j", "k", "l")
#'netstats_tm_redundancy(BCweighted, weighted = TRUE)
#'
#'

### Note formula for two-mode redundancy
#
#           | group membership of i n group membership of j |
#   Rij = ------------------------------------------------------*
#                    | group membership of i |
#
###

netstats_tm_redundancy <- function(net, # the two mode network
                       inParallel = FALSE, # should this be computed in parallel?
                       nCores = NULL, # the number of cores for computing in parallel?
                       isolates = NA,
                       weighted = FALSE) { # what value should isolates get?

  n_actors <- nrow(net) # the number of actors in the network

  # rij = shared level 1 nodes between actors i and j    /    all level 1 nodes that actor i is tied to
  ################### 3
  ###
  ###    We define redundancy as follows. The extent to which j is a redundant contact of an ego, i’s,
  ###    other same-class contacts is measured by ..... above forumla ... where sigma(i) is the set of
  ###    all same-class contacts of the node i. In otherwords, j is a redundant contact of i’s other
  ###    same-class contacts to the extent that a large proportion of its total same-class contacts are shared
  ###    with j’ s.  (Burchard and Cornwell 2018: pp 14)
  ###
  ###################

  if(weighted == FALSE){ #Are the two mode networks weighted? No

  if (inParallel == FALSE) { # if the user does not want to use parallel

    redund.net <- matrix(0, nrow = nrow(net), ncol = nrow(net)) # an empty nxn matrix to store the rij scores
    for (i in 1:n_actors) { # for all actors

      egoneti <- net[, net[i, ] == 1] # extracting the ego 2 mode network for an actor i
      pnp <- (egoneti %*% t(egoneti)) # getting the shared group membership for any two level 1 actors
      diag(pnp) <- 0 # setting the diagonal to 0
      pnp[pnp > 1] <- 1 # setting all true co-membership ties back to 1
      sigmai <- rowSums(pnp) # getting the shared group membership of i
      same_nodes <- pnp %*% t(pnp) # getting the shared alters
      diag(same_nodes) <- 0 # setting the diagonal to zero
      testing_file <- same_nodes / sigmai # dividing the number of shared alters by the total number of alters
      redund.net[i, ] <- testing_file[i, ] # adding these values to the empty matrix
    }

    redund.net[is.nan(redund.net)] <- isolates # the value specified to isolates

    # checking if the user provided names to the level-1 actors
    ifelse(is.null(rownames(net)), # if it is null
      {
        rownames(redund.net) <- rownames(redund.net) # keep it null
        colnames(redund.net) <- rownames(redund.net)
      }, # keep it null
      # if not null
      {
        rownames(redund.net) <- rownames(net) # change names to the one's provided by user
        colnames(redund.net) <- rownames(net)
      }
    ) # change names to the one's provided by user

    return(redund.net) # return the redundany matrix
  }

  if (inParallel == TRUE) { # if the user does want to use parallel
    if (is.null(nCores)) { # do this user provide the number of cores to use?
      nCores <- round(parallel::detectCores() / 2) # if not, detect the number of cores, and use half of them
    }
    myCluster <- parallel::makeForkCluster(nnodes = nCores) # creating the cluster
    doParallel::registerDoParallel(myCluster) # registering the cluster

    # using the foreach package to do the computations in parallel
    redund.net <- foreach::foreach(i = 1:n_actors, .combine = rbind) %dopar%

      ({
        egoneti <- net[, net[i, ] == 1] # extracting the ego 2 mode network for an actor i
        pnp <- (egoneti %*% t(egoneti)) # getting the shared group membership for any two level 1 actors
        diag(pnp) <- 0  # setting the diagonal to 0
        pnp[pnp > 1] <- 1 # setting all true co-membership ties back to 1
        sigmai <- rowSums(pnp) # getting the shared group membership of i
        same_nodes <- pnp %*% t(pnp) # getting the shared alters
        diag(same_nodes) <- 0 # setting the diagonal to zero
        testing_file <- same_nodes / sigmai # dividing the number of shared alters by the total number of alters
        testing_file[i, ] # adding these values to the empty matrix
      })
    parallel::stopCluster(myCluster) # closing the cluster

    redund.net[is.nan(redund.net)] <- isolates # the value specified to isolates

    # checking if the user provided names to the level-1 actors
    ifelse(is.null(rownames(net)), # if it is null
      {
        rownames(redund.net) <- rownames(redund.net) # keep it null
        colnames(redund.net) <- rownames(redund.net)
      }, # keep it null
      # if not null
      {
        rownames(redund.net) <- rownames(net) # change names to the one's provided by user
        colnames(redund.net) <- rownames(net)
      }
    ) # change names to the one's provided by user

    return(redund.net) # return the redundany matrix
  }

  }

  if(weighted == TRUE){ #Are the two mode networks weighted? Yes

    if (inParallel == FALSE) { # if the user does not want to use parallel

      redund.net <- matrix(0, nrow = nrow(net), ncol = nrow(net)) # an empty nxn matrix to store the rij scores
      for (i in 1:n_actors) { # for all actors

        egoneti <- net[, net[i, ] > 0] # extracting the ego 2 mode network for an actor i
        pnp <- (egoneti %*% t(egoneti)) # getting the shared group membership for any two level 1 actors
        diag(pnp) <- 0 # setting the diagonal to 0
        pnp[pnp > 1] <- 1 # setting all true co-membership ties back to 1
        sigmai <- rowSums(pnp) # getting the shared group membership of i

        #weight correction
        wt <- matrix(0, nrow = nrow(net ), ncol = nrow(net ))
        for(k in 1:nrow(net )){
          for(j in 1:nrow( net )){
            if(k != j & k < j){ #the matrix is symmetric, so only go thorugh the upper
              imem <- net[k,]
              jmem <- net[j,]
              samemem <- which(imem != 0 & jmem != 0)
              wt[k,j] <- mean(c(imem[samemem], jmem[samemem]))
              wt[k,j] <- wt[k,j] * sigmai[k]
            }}
        }
        wt1 <- wt + t(wt)
        wt1[is.nan(wt1)] <- 0
        same_nodes <- pnp %*% t(pnp) # getting the shared alters
        diag(same_nodes) <- 0 # setting the diagonal to zero
        testing_file <- same_nodes / wt1 # dividing the number of shared alters by the total number of alters
        redund.net[i, ] <- testing_file[i, ] # adding these values to the empty matrix
      }

      redund.net[is.nan(redund.net)] <- isolates # the value specified to isolates
      diag(redund.net) <- 0

      # checking if the user provided names to the level-1 actors
      if(is.null(rownames(net))){ # if it is null

               rownames(redund.net) <- rownames(redund.net) # keep it null
               colnames(redund.net) <- rownames(redund.net)# keep it null
      }else{ # if not null
               rownames(redund.net) <- rownames(net) # change names to the one's provided by user
               colnames(redund.net) <- rownames(net) # change names to the one's provided by user
       }


      return(redund.net) # return the redundany matrix
    }

    if (inParallel == TRUE) { # if the user does want to use parallel
      if (is.null(nCores)) { # do this user provide the number of cores to use?
        nCores <- round(parallel::detectCores() / 2) # if not, detect the number of cores, and use half of them
      }
      myCluster <- parallel::makeForkCluster(nnodes = nCores) # creating the cluster
      doParallel::registerDoParallel(myCluster) # registering the cluster

      # using the foreach package to do the computations in parallel
      redund.net <- foreach::foreach(i = 1:n_actors, .combine = rbind) %dopar%

        ({
          egoneti <- net[, net[i, ] == 1] # extracting the ego 2 mode network for an actor i
          pnp <- (egoneti %*% t(egoneti)) # getting the shared group membership for any two level 1 actors
          diag(pnp) <- 0  # setting the diagonal to 0
          pnp[pnp > 1] <- 1 # setting all true co-membership ties back to 1
          sigmai <- rowSums(pnp) # getting the shared group membership of i
          same_nodes <- pnp %*% t(pnp) # getting the shared alters
          diag(same_nodes) <- 0 # setting the diagonal to zero
          #weight correction
          wt <- matrix(0, nrow = nrow(net ), ncol = nrow(net ))
          for(k in 1:nrow(net )){
            for(j in 1:nrow( net )){
              if(k != j & k < j){ #the matrix is symmetric, so only go thorugh the upper
                imem <- net[k,]
                jmem <- net[j,]
                samemem <- which(imem != 0 & jmem != 0)
                wt[k,j] <- mean(c(imem[samemem], jmem[samemem]))
                wt[k,j] <- wt[k,j] * sigmai[k]
              }}
          }
          wt1 <- wt + t(wt)
          wt1[is.nan(wt1)] <- 0
          testing_file <- same_nodes / wt1 # dividing the number of shared alters by the total number of alters
          testing_file[i, ] # adding these values to the empty matrix
        })
      parallel::stopCluster(myCluster) # closing the cluster

      redund.net[is.nan(redund.net)] <- isolates # the value specified to isolates
      diag(redund.net) <- 0


      # checking if the user provided names to the level-1 actors
      ifelse(is.null(rownames(net)), # if it is null
             {
               rownames(redund.net) <- rownames(redund.net) # keep it null
               colnames(redund.net) <- rownames(redund.net)
             }, # keep it null
             # if not null
             {
               rownames(redund.net) <- rownames(net) # change names to the one's provided by user
               colnames(redund.net) <- rownames(net)
             }
      ) # change names to the one's provided by user

      return(redund.net) # return the redundany matrix
    }
  }
}
