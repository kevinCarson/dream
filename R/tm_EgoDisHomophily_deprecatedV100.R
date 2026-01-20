## The Creation of Two Mode Homophily Distance
## see Fujimoto, Kayo, Tom Snijders, and Thomas Valente. 2018.
##   "Multivariate dynamics of one-mode and two-mode networks: Explaining similarity
##   in sports participation among friends." Network Science, Vol 6(3), pp. 370-395
## Code written by Kevin Carson (kacarson@arizona.edu) and Deigo Leal (https://www.diegoleal.info/)
## Last Updated: 09-07-24
#' @title Compute Fujimoto, Snijders, and Valente's (2018) Ego Homophily Distance for Two-Mode Networks
#' @name computeTMEgoDis
#' @param net The two-mode adjacency matrix.
#' @param mem The vector of membership values that the homophilous four cycles will be based on.
#' @param standardize TRUE/FALSE. TRUE indicates that the sores will be standardized by the number of level 2 nodes the level 1 node is connected to. FALSE indicates that the scores will not be standardized. Set to FALSE by default.
#' @return The vector of two-mode ego homophily distance.
#' @export
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `computeTMEgoDis()` has been deprecated starting on version 1.0.0 of the `dream` package. Please use the `netstats_tm_egodistance()` function and see the `NEWS.md` file for more details.
#'
#'
#' This function computes the ego homophily distance in two-mode
#' networks as  proposed by Fujimoto, Snijders, and Valente (2018: 380).
#' See Fujimoto, Snijders, and Valente (2018) for more details about this
#' measure.
#'
#' @details
#' The formula for ego homophily distance in two-mode networks is:
#'\deqn{Ego2Dist_{i} = \sum_{a}y_{ia}{1 - |v_i - p_ia |}     }
#' where:
#' \itemize{
#'   \item \eqn{\sum_a} sums across all level 2 nodes in the network
#'   \item \eqn{y_{ia}} is the 1 if node i is tied to node a and 0 else.
#'   \item \eqn{v_i} is the value of the respondent. Within the function this is
#'   predefined to be 1 if there are multiple categories.
#'   \item \eqn{p_ia} is the proportion of same-category actors that are tied to
#'   node a not including the ego itself.
#'   \item \eqn{|v_i - p_ia|} is equal to 1 if all the level 1 nodes that are tied
#'   to the level 2 node share the same categorical membership and 0 if all
#'   level 1 nodes are a different category.
#'
#' }
#'
#' If the ego is a level 2 isolate or a level 2 pendant, that is, only one level 1
#' node (e.g., patient) is connected to that specific level 2 node (e.g., medical doctor),
#' then they are given a value of 0. In particular, the contribution to
#' the ego distance for a pendant is 0. The ego distance value can be standardized
#' by the number of groups which would provide the average ego distance as a
#' proportion between 0 and 1.
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
#'#the ego 2 mode distance non-standardized
#'computeTMEgoDis(southern.women, mem = membership)
#'#the ego 2 mode distance standardized
#'computeTMEgoDis(southern.women, mem = membership, standardize = TRUE)
#'
computeTMEgoDis <- function(net, #the two-mode adjacency matrix
                      mem,#the vector of membership scores
                      standardize = FALSE){ #to standardize the scores for all non_pendant groups
  lifecycle::deprecate_warn("1.0.0", " computeTMEgoDis()", "netstats_tm_egodistance()")
  dist <- rep(0, nrow(net)) #creating an empty vector to store the homophily distance variables
  for(i in 1:nrow(net)){ #for all level 1 actors in the network
    #recreate the membership vector to 1 for all i == Vi and 0 if not
    new <- rep(0, length(mem)) #the membership values for each actor
    new[which(mem == mem[i])] <- 1 #the actors who share the same values
    vi <- new[i] #the current actors value
    mode2tie <- as.numeric(which(net[i,] == 1)) #the actors who the node is tied to
    if(length(mode2tie) != 0){ #if the actor is not an isolate
      for(j in 1:length(mode2tie)){ #for all groups in network Y
        a <- mode2tie[j] #the current group
        if(net[i,a] == 1){
          if(sum(net[,a]) == 1){  #if there is only one actor in the group
            dist[i] <- dist[i] + 0 #since, pia is 0/0 and illdefined
          }else{
            pia <- (sum(new[net[,a] == 1])-1) / (sum((net[,a] == 1)) - 1)
            dist[i] <- dist[i] + (net[i,a] * ( 1 - abs((1 - pia)))    )
            }
      }

    }

    }
  }
  if(standardize == TRUE){
    pendants <- which(colSums(net) == 1) #which groups only have 1 actor
    net[,pendants] <- 0 #making there ties now to 0
    ngroups <- rowSums(net) #the non-pendenant group membership scores
    dist <- dist / ngroups #standardizing the scores to be within 0 to 1.
  }

  return(dist)
}










