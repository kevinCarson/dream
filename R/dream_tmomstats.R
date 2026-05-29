# the master file for dream R package functions related to the computation
# of one-mode and two-mode temporal/static network statistics

#' @title Compute Potential for Intercultural Brokerage (PIB) Based on Leal (2025)
#' @name netstats_om_pib
#' @param net The one-mode adjacency matrix.
#' @param g.mem The vector of membership values that the brokerage scores will be based on.
#' @param symmetric TRUE/FALSE. TRUE indicates that network matrix will be treated as symmetric. FALSE indicates that the network matrix will be treated as asymmetric. Set to TRUE by default.
#' @param triad.type The string value (or vector) that indicates what specific triadic (star) structures the potential for cultural brokerage will be computed for. Possible values are "ANY", "OTS", "ITS", "MTS" (see the details section). The function defaults to “ANY”.
#' @param count TRUE/FALSE. TRUE indicates that the number of culturally brokered open triangles will be returned. FALSE indicates that the proportion of culturally brokered open triangles to all open triangles will be returned (see the details section). Set to TRUE by default.
#' @param isolate If count = FALSE, the numerical value that will be given to isolates. This value is set to NA by default, as 0/0 is undefined. The user can specify this value!
#' @import Rcpp
#' @return The vector of interpersonal cultural brokerage values for the one-mode network.
#' @export
#'
#'
#'
#'
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' Following Leal (2025), this function calculates node’s Potential for
#' Intercultural Brokerage (PIB) in a one-mode network, that is, brokerage
#' based on nodes’ distinct group memberships. For example, users can examine
#' PIB based on actors’ gender. The option count determines what is returned
#' by the function. If count is TRUE, then the count of ‘culturally’ dissimilar
#' pairs brokered by ego is included (i.e., ego’s total count of brokered open
#' triangles where the alters at the two endpoints of said open triangles
#' are ‘culturally’ dissimilar from one another). If count is FALSE, the
#' proportion of ego’s brokered open triangles where the endpoints
#' are ‘culturally’ dissimilar out of all of ego’s brokered open
#' triangles (regardless of the cultural identity of the alters) is returned. The
#' formula for computing interpersonal brokerage is presented in the details section.
#'
#' @details
#' Following Leal (2025), the formula for interpersonal brokerage is:
#'
#' \deqn{ \text{PIB}_i = \sum_{j < k} \frac{S_{jik}}{S_{jk}} m_{jk}, \quad S_{jik} \neq 0 \text{ and } i \neq j \neq k }
#'
#' where:
#' \itemize{
#'   \item \eqn{S_{jik} = 1} if there is an (un)directed two-path connecting actors *j* and *k* through actor *i*; 0 otherwise.
#'   \item \eqn{m_{jk} = 1} if actors *j* and *k* are on different sides of a symbolic boundary; 0 otherwise.
#'   \item Following Gould (1989), \eqn{S_{jik}} represents the total number of two-paths between actors *j* and *k*.
#' }
#' If the network is non-symmetric (i.e., the user specified symmetric = FALSE), then
#' the function can compute the cultural brokerage scores for different star structures. The possible
#' values are: "ANY", which computes the scores for all structures, where a tie exists
#' between *i* and *j*, *j* and *k*, and one does not exist between *i* and *k*. "OTS" computes the
#' values for outgoing two-stars (i<-j->k or the 021D triad according to the M.A.N. notation; see Wasserman and Faust 1994), where j is the broker. "ITS" computes the
#' values for incoming two-stars (i->j<-k or the 021U triad according to the M.A.N. notation; see Wasserman and Faust 1994), where j is the broker. "MTS" computes PIB for
#' mixed triadic structures (i<-j<-k or i->j->k or the 021C triad according to the M.A.N. notation; see Wasserman and Faust 1994). If not specified, the function defaults to the
#' "ANY" category. This function can also compute all of the formations at once.
#'
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#' @references
#' Gould, Roger. 1989. "Power and Social Structure in Community Elites." *Social Forces* 68(2): 531-552.
#'
#' Leal, Diego F. 2025. "Locating Cultural Holes Brokers in Diffusion Dynamics Across Bright Symbolic Boundaries." *Sociological Methods & Research* \doi{10.1177/00491241251322517}
#'
#' Wasserman, Stanley and Katherine Faust. 1994. *Social Network Analysis: Methods and Applications*. Cambridge: Cambridge University Press.
#'
#' @examples
#'
#' # For this example, we recreate Figure 3 in Leal (2025)
#' LealNet <- matrix( c(
#'  0,1,0,0,0,0,0,
#'  1,0,1,1,0,0,0,
#'  0,1,0,0,1,1,0,
#'  0,1,0,0,1,0,0,
#'  0,0,1,1,0,0,0,
#'  0,0,1,0,0,0,1,
#'  0,0,0,0,0,1,0),
#'  nrow = 7, ncol = 7, byrow = TRUE)
#'
#' colnames(LealNet) <- rownames(LealNet) <- c("A", "B", "C","D",
#'                                            "E", "F", "G")
#' categorical_variable <- c(0,0,1,0,0,0,0)
#' #These values are exactly the same as reported by Leal (2025)
#' netstats_om_pib(LealNet,
#'    symmetric = TRUE,
#'    g.mem = categorical_variable)
#'
#'
#'
netstats_om_pib <- function(net,
                            g.mem,
                            symmetric = TRUE,
                            triad.type = NULL, #can be of c("ANY", "OTS", "ITS", "MTS")
                            count = TRUE,
                            isolate = NA){

  if(symmetric == FALSE){
    message("Note: the network is specified as symmetric by default. However, you
           specified the network as asymmetric. Therefore, if traid.type was
           not specified, the function defaults to 'ANY'.")}
  if(symmetric == FALSE & is.null(triad.type)){ #if the triad type was not specifed,
    #simply compute all of the possible starst
    triad.type <- "ANY"
  }
  if(sum(!(triad.type %in% c("ANY", "OTS", "ITS", "MTS"))) != 0){
    base::stop("an unknown type was specifed in the triad.type argument. Please
               update this and retry!")
  }
  if(symmetric){
    triad.type <- "ANY" # this does not matter, as there are no structural differences in symmetric networks.
  }

  ntypes <- length(triad.type)
  results <- matrix(0,
                    nrow = nrow(net),
                    ncol = length(ntypes))
  colnames(results) <- triad.type
  for(i in 1:ntypes){
    results[,i] <-  pibcpp(net = net,
                           gmem = g.mem,
                           symmetric = symmetric,
                           traidtype = triad.type[i],
                           count=count,
                           isolate=isolate)
  }
  if(is.null(rownames(net))){
    rownames(results) <- paste0("Node_", 1:nrow(net))
  }else{
    rownames(results) <- rownames(net)

  }
  return(results)

}





#' @title Compute Burchard and Cornwell's (2018) Two-Mode Constraint
#' @name netstats_tm_constraint
#' @param net A two-mode adjacency matrix or affiliation matrix.
#' @param isolates What value should isolates be given? Preset to be NA.
#' @param returnCIJmat TRUE/FALSE. TRUE indicates that the full constraint matrix, that is, the network constraint from an alter j on node i, will be returned to the user. FALSE indicates that the total constraint will be returned. Set to FALSE by default.
#' @param weighted TRUE/FALSE. TRUE indicates the resulting statistic will be based on the weighted formula (see the details section). FALSE indicates the statistic will be based on the original non-weighted formula. Set to FALSE by default.
#' @return The vector of two-mode constraint scores for level 1 actors in a two-mode network.
#' @export
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' This function calculates the values for two-mode network constraint
#' for weighted and unweighted two-mode networks based on Burchard and Cornwell (2018).
#' @details Following Burchard and Cornwell (2018), the formula for two-mode constraint is:
#' \deqn{c_{ij} = \left(\frac{|\zeta(j) \cap \zeta(i)|}{|\zeta^{(i*)}|}\right)^2}
#' where:
#' \itemize{
#'   \item \eqn{c_{ij}} is the constraint of ego *i* with respect to actor *j*.
#'   \item \eqn{|\zeta(j) \cap \zeta(i)|} is the number of opposite-class contacts that *i* and *j* both share.
#'   \item The denominator, \eqn{|\zeta^{(i*)}|}, represents the total number of opposite-class contacts of ego *i* excluding pendants. Pendants are level 2 groups that only have one member (i.e., incoming tie).
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
#'netstats_tm_constraint(BCNet)
#'
#'#For this example, we recreate Figure 9 in Burchard and Cornwell (2018:18) for
#'#weighted two mode networks.
#'BCweighted <- matrix(c(1,2,1, 1,0,0,
#'                       0,2,1,0,0,1),
#'                     nrow = 4, ncol = 3,
#'                     byrow = TRUE)
#'rownames(BCweighted) <- c("i", "j", "k", "l")
#'netstats_tm_constraint(BCweighted, weighted = TRUE)
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

netstats_tm_constraint <- function(net, # the two mode network
                                   isolates = NA,
                                   returnCIJmat = FALSE,
                                   weighted = FALSE) { # what value should isolates get?
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




#' @title Compute Burchard and Cornwell's (2018) Two-Mode Effective Size
#' @name netstats_tm_effective
#' @param net A two-mode adjacency matrix or affiliation matrix
#' @param isolates What value should isolates be given? Preset to be NA.
#' @param weighted TRUE/FALSE. TRUE indicates the resulting statistic will be based on the weighted formula (see the details section). FALSE indicates the statistic will be based on the original non-weighted formula. Set to FALSE by default.
#' @param inParallel TRUE/FALSE. TRUE indicates that parallel processing will be used to compute the statistic with the *foreach* package. FALSE indicates that parallel processing will not be used. Set to FALSE by default.
#' @param nCores If inParallel = TRUE, the number of computing cores for parallel processing. If this value is not specified, then the function internally provides it by dividing the number of available cores in half.
#' @return The vector of two-mode effective size values for level 1 actors in a two-mode network.
#' @import foreach
#' @import doParallel
#' @import parallel
#' @export

#' @description
#' `r lifecycle::badge("stable")`
#'
#' This function calculates the values for two-mode effective size for
#' weighted and unweighted two-mode networks based on Burchard and Cornwell (2018).
#' @details The formula for two-mode effective size is:
#' \deqn{ES_{i} = |\sigma(i)| - \sum_{j \in \sigma(i)} r_{ij}}
#' where:
#' \itemize{
#'   \item \eqn{ES_{i}} is the effective size of ego *i*.
#'   \item \eqn{|\sigma(i)|} is the number of same-class contacts of ego *i*.
#'   \item \eqn{\sum_{j \in \sigma(i)} r_{ij}} is the summation of the redundancy
#'   for each alter *j* in the two-mode ego network of *i*.
#' }
#' This function allows the user to compute the scores in parallel through the *foreach* and *doParallel* R packages.
#' If the matrix is weighted, the user should specify *weighted = TRUE*. If the matrix is weighted, following
#' Burchard and Cornwell (2018), the formula for two-mode weighted redundancy is:
#' \deqn{r_{ij} = \frac{|\sigma(j) \cap \sigma(i)|}{|\sigma(i)| \times w_t}}
#' where \eqn{w_t} is the average of the tie weights that *i* and *j* send
#' to their shared opposite class contacts.

#'
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#' @references
#' Burchard, Jake and Benjamin Cornwell. 2018. "Structural Holes and Bridging
#' in Two-Mode Networks." *Social Networks* 55:11-20.

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
#'#library(sna) #To plot the two mode network, we use the sna R package
#'#gplot(BCNet, usearrows = FALSE,
#'#      gmode = "twomode", displaylabels = TRUE)
#'netstats_tm_effective(BCNet)
#'
#'#In this example, we recreate Figure 9 in Burchard and Cornwell (2018:18)
#'#for weighted two mode networks.
#'BCweighted <- matrix(c(1,2,1, 1,0,0,
#'                       0,2,1,0,0,1),
#'                       nrow = 4, ncol = 3,
#'                       byrow = TRUE)
#'rownames(BCweighted) <- c("i", "j", "k", "l")
#'netstats_tm_effective(BCweighted, weighted = TRUE)
#'



netstats_tm_effective <- function(net, # the two mode network
                                  inParallel = FALSE, # should this be computed in parallel?
                                  nCores = NULL, # the number of cores for computing in parallel?
                                  isolates = NA,# what value should isolates get?
                                  weighted = FALSE) {  #is the matrix weighted?

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

  if(weighted == FALSE){ #Is the matrix weighted? No

    if (inParallel == FALSE) { # if the user does not want to use parallel

      redund.net <- matrix(0, nrow = nrow(net), ncol = nrow(net)) # an empty nxn matrix to store the rij scores
      sigma <- rep(0, n_actors) # the sigma vector (i.e., the number of level 1 alters)
      for (i in 1:n_actors) { # for all actors

        egoneti <- net[, net[i, ] == 1] # extracting the ego 2 mode network for an actor i
        pnp <- (egoneti %*% t(egoneti)) # getting the shared group membership for any two level 1 actors
        diag(pnp) <- 0 # setting the diagonal to 0
        pnp[pnp > 1] <- 1 # setting all true co-membership ties back to 1
        sigmai <- rowSums(pnp) # getting the shared group membership of i
        sigma[i] <- sigmai[i]
        same_nodes <- pnp %*% t(pnp) # getting the shared alters
        diag(same_nodes) <- 0 # setting the diagonal to zero
        testing_file <- same_nodes / sigmai # dividing the number of shared alters by the total number of alters
        redund.net[i, ] <- testing_file[i, ] # adding these values to the empty matrix
      }

      redund.net[is.nan(redund.net)] <- isolates # the value specified to isolates
      ES <- rep(0, n_actors) # pre-creating a vector to hold the effective size values
      Sij <- rowSums(redund.net, na.rm = TRUE) # the total sum of redundancy for all egos
      ##
      ##
      ##        ESi = sigma(i) - sumRij
      ##
      ##       #sigma(i) is the number of same level contacts ego i has
      ##       #sumRij is the sum of redundnacy between ego i and all same level alters j
      ##       #A special note for level 1 isolates and level 1 isolates with a tie to level 2
      ##       #isolates (i.e., level 2 groups with only 1 node as a member) (pp. 15)
      ##
      ES <- as.numeric(sigma - Sij) # Number of Alter - Sum of all Sij  (Burchard and Cornwell 2018: pp 14)
      ES[is.nan(ES)] <- isolates # updating the user-specified value for isolates
      ifelse(is.null(rownames(net)),
             names(ES)  <- NULL,
             names(ES) <- rownames(net))
      return(ES)
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

      sigma <- foreach::foreach(i = 1:n_actors, .combine = rbind) %dopar%

        ({
          egoneti <- net[, net[i, ] == 1] # extracting the ego 2 mode network for an actor i
          pnp <- (egoneti %*% t(egoneti)) # getting the shared group membership for any two level 1 actors
          diag(pnp) <- 0 # setting the diagonal to 0
          pnp[pnp > 1] <- 1 # setting all true co-membership ties back to 1
          sigmai <- rowSums(pnp) # getting the shared group membership of i
          sigmai[i]
        })

      parallel::stopCluster(myCluster) # closing the cluster

      redund.net[is.nan(redund.net)] <- isolates # the value specified to isolates
      ES <- rep(0, n_actors) # pre-creating a vector to hold the effective size values
      Sij <- rowSums(redund.net, na.rm = TRUE) # the total sum of redundancy for all egos
      ##
      ##
      ##        ESi = sigma(i) - sumRij
      ##
      ##       #sigma(i) is the number of same level contacts ego i has
      ##       #sumRij is the sum of redundnacy between ego i and all same level alters j
      ##       #A special note for level 1 isolates and level 1 isolates with a tie to level 2
      ##       #isolates (i.e., level 2 groups with only 1 node as a member) (pp. 15)
      ##
      ES <- as.numeric(sigma - Sij) # Number of Alter - Sum of all Sij  (Burchard and Cornwell 2018: pp 14)
      ES[is.nan(ES)] <- isolates # updating the user-specified value for isolates
      ifelse(is.null(rownames(net)),
             names(ES)  <- NULL,
             names(ES) <- rownames(net))
      return(ES)
    }
  }
  if(weighted == TRUE){ #Is the matrix weighted? Yes

    if (inParallel == FALSE) { # if the user does not want to use parallel

      redund.net <- matrix(0, nrow = nrow(net), ncol = nrow(net)) # an empty nxn matrix to store the rij scores
      sigma <- rep(0, n_actors) # the sigma vector (i.e., the number of level 1 alters)
      for (i in 1:n_actors) { # for all actors

        egoneti <- net[, net[i, ] > 0] # extracting the ego 2 mode network for an actor i
        pnp <- (egoneti %*% t(egoneti)) # getting the shared group membership for any two level 1 actors
        diag(pnp) <- 0 # setting the diagonal to 0
        pnp[pnp > 1] <- 1 # setting all true co-membership ties back to 1
        sigmai <- rowSums(pnp) # getting the shared group membership of i
        sigma[i] <- sigmai[i]
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
      ES <- rep(0, n_actors) # pre-creating a vector to hold the effective size values
      Sij <- rowSums(redund.net, na.rm = TRUE) # the total sum of redundancy for all egos
      ##
      ##
      ##        ESi = sigma(i) - sumRij
      ##
      ##       #sigma(i) is the number of same level contacts ego i has
      ##       #sumRij is the sum of redundnacy between ego i and all same level alters j
      ##       #A special note for level 1 isolates and level 1 isolates with a tie to level 2
      ##       #isolates (i.e., level 2 groups with only 1 node as a member) (pp. 15)
      ##
      ES <- as.numeric(sigma - Sij) # Number of Alter - Sum of all Sij  (Burchard and Cornwell 2018: pp 14)
      ES[is.nan(ES)] <- isolates # updating the user-specified value for isolates
      ifelse(is.null(rownames(net)),
             names(ES)  <- NULL,
             names(ES) <- rownames(net))
      return(ES)
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
          egoneti <- net[, net[i, ] > 0] # extracting the ego 2 mode network for an actor i
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

      sigma <- foreach::foreach(i = 1:n_actors, .combine = rbind) %dopar%

        ({
          egoneti <- net[, net[i, ] > 0] # extracting the ego 2 mode network for an actor i
          pnp <- (egoneti %*% t(egoneti)) # getting the shared group membership for any two level 1 actors
          diag(pnp) <- 0 # setting the diagonal to 0
          pnp[pnp > 1] <- 1 # setting all true co-membership ties back to 1
          sigmai <- rowSums(pnp) # getting the shared group membership of i
          sigmai[i]
        })


      parallel::stopCluster(myCluster) # closing the cluster

      redund.net[is.nan(redund.net)] <- isolates # the value specified to isolates
      ES <- rep(0, n_actors) # pre-creating a vector to hold the effective size values
      Sij <- rowSums(redund.net, na.rm = TRUE) # the total sum of redundancy for all egos
      ##
      ##
      ##        ESi = sigma(i) - sumRij
      ##
      ##       #sigma(i) is the number of same level contacts ego i has
      ##       #sumRij is the sum of redundnacy between ego i and all same level alters j
      ##       #A special note for level 1 isolates and level 1 isolates with a tie to level 2
      ##       #isolates (i.e., level 2 groups with only 1 node as a member) (pp. 15)
      ##
      ES <- as.numeric(sigma - Sij) # Number of Alter - Sum of all Sij  (Burchard and Cornwell 2018: pp 14)
      ES[is.nan(ES)] <- isolates # updating the user-specified value for isolates
      ifelse(is.null(rownames(net)),
             names(ES)  <- NULL,
             names(ES) <- rownames(net))
      return(ES)
    }

  }

}



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




#' @title Compute Burt's (1992) Constraint for Ego Networks from a Sociomatrix
#' @name netstats_om_constraint
#' @param net A one-mode sociomatrix with network ties.
#' @param isolates What value should isolates be given? Set to NA by default.
#' @param pendants What value should be given to pendant vertices? Set to 1 by default. Pendant vertices are those nodes who have one outgoing tie.
#' @param inParallel TRUE/FALSE. TRUE indicates that parallel processing will be used to compute the statistic with the *foreach* package. FALSE indicates that parallel processing will not be used. Set to FALSE by default.
#' @param nCores If inParallel = TRUE, the number of computing cores for parallel processing. If this value is not specified, then the function internally provides it by dividing the number of available cores in half.
#' @return The vector of ego network constraint values.
#' @import foreach
#' @import doParallel
#' @import parallel
#' @import Rcpp
#' @export

#' @description
#' `r lifecycle::badge("stable")`
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
#' who have no ties) and pendants (i.e., nodes who only have
#' one tie), the user should specify what values for constraint are returned for them via the *isolates* and
#' *pendants* options. In particular, pendant vertices are those nodes who have one outgoing tie.
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
#' netstats_om_constraint(BurtEgoNet)
#'
#'

#################################################################################
###     Burt's Constraint Measure for Ego Networks
#################################################################################

netstats_om_constraint <- function(net, # the full sociomatrix
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
      constraint[i] <- burtconstraint(net = neti, nactors = nrow(neti)) #compute the value for constraint for the actor in c++
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
        burtconstraint(net = neti, nactors = nrow(neti)) #compute the value for constraint for the actor in c++
      })

    parallel::stopCluster(myCluster) # closing the cluster

    ifelse(is.null(rownames(net)),
           names(constraint)  <- NULL,
           names(constraint) <- rownames(net))

    return(constraint)
  }
}





#' @title Compute Fujimoto, Snijders, and Valente's (2018) Homophilous Four-Cycles for Two-Mode Networks
#' @name netstats_tm_homfourcycles
#' @param net The two-mode adjacency matrix.
#' @param mem The vector of membership values that the homophilous four-cycles will be based on.
#' @return The vector of counts of homophilous four-cycles for the two-mode network.
#' @export
#'
#' @description
#' `r lifecycle::badge("stable")`
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
#' netstats_tm_homfourcycles(southern.women, mem = membership)


#### Formula for Calculating Two-Mode Four-Cycle Homophily Effect
##
##     Σj  Σa != b   Yia Yib Yja Yjb  I{Vi = Vj}
##
##     In which I is the indicator function defined as 1 if i and j share the same
##     categorical membership, and 0 elsewhere. Y represents the two-mode affliation
##     matrix, and V represents a vector of categorical membership values.
##
####
netstats_tm_homfourcycles <- function(net, # affilation matrix
                                      mem # membership homophily
) {

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






#' @title Compute Fujimoto, Snijders, and Valente's (2018) Ego Homophily Distance for Two-Mode Networks
#' @name netstats_tm_egodistance
#' @param net The two-mode adjacency matrix.
#' @param mem The vector of membership values that the homophilous four cycles will be based on.
#' @param standardize TRUE/FALSE. TRUE indicates that the sores will be standardized by the number of level 2 nodes the level 1 node is connected to. FALSE indicates that the scores will not be standardized. Set to FALSE by default.
#' @import Rcpp
#' @return The vector of two-mode ego homophily distance.
#' @export
#'
#' @description
#' `r lifecycle::badge("stable")`
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
#'netstats_tm_egodistance(southern.women, mem = membership)
#'#the ego 2 mode distance standardized
#'netstats_tm_egodistance(southern.women, mem = membership, standardize = TRUE)
#'
netstats_tm_egodistance <- function(net, #the two-mode adjacency matrix
                                    mem,#the vector of membership scores
                                    standardize = FALSE){ #to standardize the scores for all non_pendant groups

  dist <- tmegodist(net,mem) #computing the distance measures in c++
  if(standardize){ #if we need to standardize the scores
    ngroups <- rowSums(net) #computing the rowsums of the network
    ngroups[ngroups == 1] <- 0 #doing a simple check for pendants in the network
    dist <- dist/ngroups #standardizing by the number of groups in the network
  }
  return(dist) #returning the distance values
}









#' @title Compute Degree Centrality Values for Two-Mode Networks
#' @name netstats_tm_degreecent
#' @param net A two-mode adjacency matrix
#' @param level1 TRUE/FALSE. TRUE indicates that the degree centrality will be computed for level 1 nodes. FALSE indicates that the degree centrality will be computed for level 2 nodes. Set to TRUE by default.
#' @import Rcpp
#' @return The vector of two-mode level-specific degree centrality values.
#' @export
#'
#'
#' @description
#' `r lifecycle::badge("stable")`
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
#'netstats_tm_degreecent(knoke_yang_PC, level1 = TRUE) #this value matches the book
#'netstats_tm_degreecent(knoke_yang_PC, level1 = FALSE) #this value matches the book

netstats_tm_degreecent <- function(net,#a two-mode network adjancency matrix
                                   level1 = TRUE){#Boolean: TRUE indicating if the density should be
  #computed on level 1 nodes, FALSE computes graph density for the level 2 nodes

  if(level1 == TRUE){ #should level 1 density be computed?
    degree <- tmdegcentraility(net) #computing level 1 degree in c++
    ifelse(is.null(rownames(net)),
           names(degree) <- NULL,
           names(degree) <- rownames(net))
  }else{ #if not
    degree <- tmdegcentraility(t(net))#computing level 2 degree in c++
    ifelse(is.null(colnames(net)),
           names(degree) <- NULL,
           names(degree) <-colnames(net))
  }
  return(degree) #returning the network
}








#' @title Compute the Number of Walks of Length K in a One-Mode Network
#' @name netstats_om_nwalks
#' @param net An unweighted one-mode network adjacency matrix.
#' @param k A numerical value that corresponds to the length of the paths to be computed.
#' @return An *n* x *n* matrix of counts of paths.
#' @export
#'
#'
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' This function calculates the number of walks of length *k*
#' between any two vertices in an unweighted one-mode network.
#'
#' @details
#' A nice result from graph theory is that the number of walks of length *k* between
#' vertices *i* and *j* can be found by:
#' \deqn{ A_{ij}^k }
#'
#' This function assumes that there are no self-loops (i.e., the diagonal of the matrix is 0).
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
#' #counting the walks of length 2
#' netstats_om_nwalks(rnet, k = 2)
#' #counting the walks of length 5
#' netstats_om_nwalks(rnet, k = 5)

netstats_om_nwalks <- function(net, # the network adjacency matrix
                               k){ #the requested path length
  diag(net) <- 0 #ensuring there are no self-loops encase the user misspecificed
  net1 <- net  #making a copy of the original matrix
  for(i in 1:(k-1)){ #for i in the requested path length (minus 1)
    #if k = 2, for example, matrix mulitiplication should
    #only be done once
    net <- net %*% net1 #multiplying the matrix by the original
  }
  return(net) #returning the n x n
}











