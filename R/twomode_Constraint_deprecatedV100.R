## The Computation of Two-Mode Constraint in Two-Mode Static Network
## See Bruchard and Cornwell 2018 for the two-mode network
## Code written by Kevin Carson (kacarson@arizona.edu) and Diego Leal (https://www.diegoleal.info/)
## Last Updated: 10-28-24

#' @title Compute Burchard and Cornwell's (2018) Two-Mode Constraint
#' @name computeBCConstraint
#' @param net A two-mode adjacency matrix or affiliation matrix.
#' @param isolates What value should isolates be given? Preset to be NA.
#' @param returnCIJmat TRUE/FALSE. TRUE indicates that the full constraint matrix, that is, the network constraint from an alter j on node i, will be returned to the user. FALSE indicates that the total constraint will be returned. Set to FALSE by default.
#' @param weighted TRUE/FALSE. TRUE indicates the resulting statistic will be based on the weighted formula (see the details section). FALSE indicates the statistic will be based on the original non-weighted formula. Set to FALSE by default.
#' @return The vector of two-mode constraint scores for level 1 actors in a two-mode network.
#' @export
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `computeBCConstraint()` has been deprecated starting on version 1.0.0 of the `dream` package. Please use the `netstats_tm_constraint()` function and see the `NEWS.md` file for more details.
#'
#'
#' This function calculates the values for two-mode network constraint
#' for weighted and unweighted two-mode networks based on Burchard and Cornwell (2018).
#' @details Following Burchard and Cornwell (2018), the formula for two-mode constraint is:
#' \deqn{c_{ij} = \left(\frac{|\zeta(j) \cap \zeta(i)|}{|\zeta^{(i*)}|}\right)^2}
#' where:
#' \itemize{
#'   \item \eqn{c_{ij}} is the constraint of ego *i* with respect to actor *j*.
#'   \item \eqn{|\zeta(j) \cap \zeta(i)|} is the number of opposite-class contacts that *i* and *j* both share.
#'   \item The denominator, \eqn{|\zeta^{(i*)}|}, represents the total number of opposite-class contacts of ego *i* excluding pendants, that is, level 2 groups that only have one member.
#' }
#' The total constraint for ego *i* is given by:
#' \deqn{C_{i} = \sum_{j \in \sigma(i)} c_{ij}}
#' The function returns the aggregate constraint for each actor; however, the user can specify the function to return the constraint matrix by setting *returnCIJmat* to TRUE.
#'
#' The function can also compute constraint for weighted two-mode networks by setting *weighted* to TRUE. The formula for two-mode weighted constraint is:
#' \deqn{c_{ij} = \left(\frac{|\zeta(j) \cap \zeta(i)|}{|\zeta^{(i*)}|}\right)^2 \times w_{t}}
#' where \eqn{w_{t}} is the average of the tie weights that *i* and *j* send to their shared opposite-class contacts.
#'
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#' @references
#' Burchard, Jake and Benjamin Cornwell. 2018. "Structural Holes and Bridging
#' in Two-Mode Networks." *Social Networks* 55:11-20.
#'
#' @examples
#'
#'# For this example, we recreate Figure 2 in Burchard and Cornwell (2018: 13)
#'BCNet <- matrix(
#'  c(1,1,0,0,
#'    1,0,1,0,
#'    1,0,0,1,
#'    0,1,1,1),
#'  nrow = 4, ncol = 4, byrow = TRUE)
#'colnames(BCNet) <- c("1", "2", "3", "4")
#'rownames(BCNet) <- c("i", "j", "k", "m")
#'#library(sna) #To plot the two mode network, we use the sna R package
#'#gplot(BCNet, usearrows = FALSE,
#'#      gmode = "twomode", displaylabels = TRUE)
#'computeBCConstraint(BCNet)
#'
#'#For this example, we recreate Figure 9 in Burchard and Cornwell (2018:18) for
#'#weighted two mode networks.
#'BCweighted <- matrix(c(1,2,1, 1,0,0,
#'                       0,2,1,0,0,1),
#'                     nrow = 4, ncol = 3,
#'                     byrow = TRUE)
#'rownames(BCweighted) <- c("i", "j", "k", "l")
#'computeBCConstraint(BCweighted, weighted = TRUE)
#'
#'
#'
#'


# cij = shared level 2 nodes between actors i and j    /    all level 2 nodes that actor i is tied to minus pendants
# pendants are defined as level 2 nodes that have only 1 contact
################### 3
###
###            (shared level 2 contacts between actor i and j)
###   Cij = --------------------------------------------------
###            (level 2 contacts of actor i minus pendants)
###
###   Cij <-   (Cij) ^ 2
###
################### (Bruchard and Cornwell 2018: pp. 15)

computeBCConstraint <- function(net, # the two mode network
                      isolates = NA,
                      returnCIJmat = FALSE,
                      weighted = FALSE) { # what value should isolates get?
  lifecycle::deprecate_warn("1.0.0", " computeBCConstraint()", "netstats_tm_constraint()")
  if(weighted == FALSE){
  shared_groups <- net %*% t(net)
  diag(shared_groups) <- 0
  #### Need to compute the number of pendant groups
  lev2_degrees <- base::colSums(net) # the membership list of each group
  lev2_pendant <- base::which(lev2_degrees < 2) # finding groups where they have less than 2 level 1 nodes, thus they are pendants
  # Restructing the matrix to account for this:
  numberofpendants <- base::length(lev2_pendant)
  if (numberofpendants != 0) { # if there is pendant groups
    for (i in 1:numberofpendants) { # for all the level 2 pendants
      net[, lev2_pendant[i]] <- 0 # changing to no membership
    }
  }
  level2contacts <- base::rowSums(net) # the number of level 2  contacts of actor i minus pendants
  CIJ <- (shared_groups / level2contacts)^2 # See the above formula (Bruchard and Cornwell 2018: pp. 15)
  if (returnCIJmat == TRUE) { # if the user wants the matrix
    CIJ[base::is.nan(CIJ)] <- isolates # changing the value given to isolates from the user
    return(CIJ) # return the matrix
  }
  total_constraint <- base::rowSums(CIJ) # the total sum of constraint for each ego
  ifelse(base::is.null(base::rownames(net)),
         names(total_constraint)  <- NULL,
         names(total_constraint) <- rownames(net))
  return(total_constraint) # return the value
  }

  if(weighted == TRUE){
    net1 <- net
    net[net > 0] <- 1
    shared_groups <- net %*% t(net)
    diag(shared_groups) <- 0
    #### Need to compute the number of pendant groups
    lev2_degrees <- colSums(net) # the membership list of each group
    lev2_pendant <- which(lev2_degrees < 2) # finding groups where they have less than 2 level 1 nodes, thus they are pendants
    # Restructing the matrix to account for this:
    numberofpendants <- length(lev2_pendant)
    if (numberofpendants != 0) { # if there is pendant groups
      for (i in 1:numberofpendants) { # for all the level 2 pendants
        net[, lev2_pendant[i]] <- 0 # changing to no membership
      }
    }
    level2contacts <- rowSums(net) # the number of level 2  contacts of actor i minus pendants
    CIJ <- (shared_groups / level2contacts)^2 # See the above formula (Bruchard and Cornwell 2018: pp. 15)
    #now getting the weighted network
    #weight correction
    wt <- matrix(0, nrow = nrow(net ), ncol = nrow(net ))
    for(k in 1:nrow(net1 )){
      for(j in 1:nrow( net1 )){
        if(k != j & k < j){ #the matrix is symmetric, so only go thorugh the upper
          imem <- net1[k,]
          jmem <- net1[j,]
          samemem <- which(imem != 0 & jmem != 0)
          wt[k,j] <- base::mean(c(imem[samemem], jmem[samemem]))
          wt[k,j] <- wt[k,j]
        }}
    }
    wt1 <- wt + t(wt)
    wt1[is.nan(wt1)] <- 0
    CIJ <- CIJ * wt1
    if (returnCIJmat == TRUE) { # if the user wants the matrix
      CIJ[is.nan(CIJ)] <- isolates # changing the value given to isolates from the user
      return(CIJ) # return the matrix
    }
    total_constraint <- rowSums(CIJ) # the total sum of constraint for each ego
    ifelse(is.null(rownames(net)),
           names(total_constraint)  <- NULL,
           names(total_constraint) <- rownames(net))
    return(total_constraint) # return the value
  }

}
