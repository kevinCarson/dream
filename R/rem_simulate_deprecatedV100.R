## Simulating Random Relational Event Sequences
## Code written by Kevin Carson (kacarson@arizona.edu) and Deigo Leal (https://www.diegoleal.info/)
## Last Updated: 11-22-24
#' @title Simulate a Random One-Mode Relational Event Sequence
#' @name simulateRESeq
#' @param n_actors The number of potential actors in the event sequence.
#' @param n_events The number of simulated events for the relational event sequence.
#' @param inertia TRUE/FALSE. True indicates the effect will be included (see the details section). FALSE indicates the effect will not be included.
#' @param inertia_p If *inertia* = TRUE, the numerical value that corresponds to the parameter weight for the inertia statistic.
#' @param recip TRUE/FALSE. True indicates the effect will be included (see the details section). FALSE indicates the effect will not be included.
#' @param recip_p If *recip* = TRUE, the numerical value that corresponds to the parameter weight for the reciprocity statistic.
#' @param sender_outdegree TRUE/FALSE. True indicates the effect will be included (see the details section). FALSE indicates the effect will not be included.
#' @param sender_outdegree_p If *sender_outdegree* = TRUE, the numerical value that corresponds to the parameter weight for the outdegree statistic.
#' @param sender_indegree TRUE/FALSE. True indicates the effect will be included (see the details section). FALSE indicates the effect will not be included.
#' @param sender_indegree_p If *sender_indegree* = TRUE, the numerical value that corresponds to the parameter weight for the indegree statistic.
#' @param target_outdegree TRUE/FALSE. True indicates the effect will be included (see the details section). FALSE indicates the effect will not be included.
#' @param target_outdegree_p If *target_outdegree* = TRUE, the numerical value that corresponds to the parameter weight for the outdegree statistic.
#' @param target_indegree TRUE/FALSE. True indicates the effect will be included (see the details section). FALSE indicates the effect will not be included.
#' @param target_indegree_p If *target_indegree* = TRUE, the numerical value that corresponds to the parameter weight for the indegree statistic.
#' @param assort Boolean. TRUE/FALSE. True indicates the effect will be included (see the details section). FALSE indicates the effect will not be included.
#' @param assort_p If *assort* = TRUE, the numerical value that corresponds to the parameter weight for the assortativity statistic.
#' @param trans_trips TRUE/FALSE. True indicates the effect will be included (see the details section). FALSE indicates the effect will not be included.
#' @param trans_trips_p If *trans_trips* = TRUE, the numerical value that corresponds to the parameter weight for the transitive triplets statistic.
#' @param three_cycles TRUE/FALSE. True indicates the effect will be included (see the details section). FALSE indicates the effect will not be included.
#' @param three_cycles_p If *three_cycles* = TRUE, the numerical value that corresponds to the parameter weight for the three cycles statistic.
#' @param starting_events A *n* x 2 dataframe with *n* starting events and 2 columns. The first column should be the sender and the second should be the target.
#' @param returnStats TRUE/FALSE. TRUE indicates that the requested network statistics will be returned alongside the simulated relational event sequence. FALSE indicates that only the simulated relational event sequence will be returned. Set to FALSE by default.
#' @return A data frame that contains the simulated relational event sequence with the sufficient statistics (if requested).
#' @importFrom collapse whichv
#' @import data.table
#' @export


#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `simulateRESeq()` has been deprecated starting on version 1.0.0 of the `dream` package. Instead, please use the `simulate_rem_seq()` function and see the `NEWS.md` file for more details.
#'
#'
#'
#' The function allows users to simulate a random one-mode relational event sequence
#' between *n* actors for *k* events. Importantly, this function follows the methods
#' discussed in Butts (2008), Amati, Lomi, and Snijders (2024), and Scheter and
#' Quintane (2021). See the details for more information on this algorithm. Critically,
#' this function can be used to simulate a random event sequence, to assess the goodness of
#' fit for ordinal timing relational event models (see Amati, Lomi, and Snijders 2024), and simulate
#' random outcomes for relational outcome models.
#'
#' @details
#' Following the authors listed in the descriptions section, the probability of
#' selecting a new event for *t+1* based on the past relational history, *\eqn{H_{t}}*,  from *\eqn{0<t<t+1}*
#' is given by:
#'
#' \deqn{ p(e_{t}) = \frac{\lambda{ij}(t;\theta)}{\sum_{(u,v)\in R_{t}} \lambda_{uv}(t;\theta)} }
#'
#' where *(i,j,t)* is the triplet that corresponds to the dyadic pair with sender *i*
#' and target *j* at time *t* contained in the full risk set, *\eqn{R_{t}}*, based on the
#' past relational history. \eqn{\lambda_{ij}(t;\theta)} is formulated as:
#'
#' \deqn{ \lambda_{ij}(t;\theta) = e^{\sum_{p}\theta_{p} X_{ijp}(H_{t})} }
#'
#' where \eqn{\theta_{p}} corresponds to the specific parameter weight given by the
#' user, and  \eqn{X_{ijp}} represents the value of the specific statistic based on the
#' current past relational history *\eqn{H_{t}}*.
#'
#' Following Scheter and Quintane (2021) and Amati, Lomi, and Snijders (2024), the
#' algorithm for simulating the random relational sequence for *k* events is:
#' \itemize{
#'   \item 1. Initialize the full risk set, *\eqn{R_{t}}*, which is the full Cartesian plot of actors.
#'   \item 2. Randomly sample the first event \eqn{ e_{1} } and add that event into the relational history, *\eqn{H_{t}}*.
#'   \item 3. Until *i* = *k*, compute the sufficient statistics for each event in the risk set,
#'   sample a new event \eqn{ e_{i} } based on the probability function specified above, and add that element into the relational history.
#'   \item 4. End when *i* > *k*.
#' }
#'
#' Currently, the function supports 6 statistics for one-mode networks. These are:
#' \itemize{
#' \item Inertia: \eqn{{n}_{ijt}}
#' \item Reciprocity: \eqn{{n}_{jit}}
#' \item Target Indegree: \eqn{\sum_{k} {n}_{kjt}}
#' \item Target Outdegree: \eqn{\sum_{k} {n}_{jkt}}
#' \item Sender Outdegree: \eqn{\sum_{k} {n}_{ikt}}
#' \item Sender Indegree: \eqn{\sum_{k} {n}_{kit}}
#' \item Assortativity: \eqn{\sum_{k} {n}_{kit} \cdot \sum_{k} {n}_{ikt}}
#' \item Transitive Triplets: \eqn{\sum_{k} {n}_{ikt} \cdot {n}_{kjt}}
#' \item Three Cycles: \eqn{\sum_{k} {n}_{jkt} \cdot {n}_{kit}}
#' }
#' Where *n* represents the counts of past events, *i* is the event sender, and *j* is the event target. See Scheter and Quintane (2021)
#' and Butts (2008) for a further discussion of these statistics.
#'
#' Users are allowed to insert a starting event sequence to base the simulation on. A few things are worth nothing. The starting
#' event sequence should be a matrix with *n* rows indicating the number of starting events and 2 columns, with the
#' first representing the event senders and the second column representing the event targets. Internally, the number
#' of actors is ignored, as the number of possible actors in the risk set is based only on the actors present in the
#' starting event sequence. Finally, the sender and target actor IDs should be numerical values.
#'
#'
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#' @references
#' Amati, Viviana, Alessandro Lomi, and Tom A.B. Snijders. 2024. "A goodness of fit
#' framework for relational event models." *Journal of the Royal Statistical Society Series A: Statistics in Society*
#' 187(4): 967-988.
#'
#' Butts, Carter T. "A Relational Framework for Social Action." *Sociological Methodology*
#' 38: 155-200.
#'
#' Schecter, Aaron and Eric Quintane. 2021 "The Power, Accuracy, and Precision of the Relational
#' Event Model." *Organizational Research Methods* 24(4): 802-829.
#'
#'
#' @examples
#'#Creating a random relational sequence with 5 actors and 25 events
#'rem1<- simulateRESeq(n_actors = 25,
#'                      n_events = 1000,
#'                      inertia = TRUE,
#'                      inertia_p = 0.12,
#'                      recip = TRUE,
#'                      recip_p = 0.08,
#'                      sender_outdegree = TRUE,
#'                      sender_outdegree_p = 0.09,
#'                      target_indegree = TRUE,
#'                      target_indegree_p = 0.05,
#'                      assort = TRUE,
#'                      assort_p = -0.01,
#'                      trans_trips = TRUE,
#'                      trans_trips_p = 0.09,
#'                      three_cycles = TRUE,
#'                      three_cycles_p = 0.04,
#'                      starting_events = NULL,
#'                      returnStats = TRUE)
#'rem1
#'
#' #Creating a random relational sequence with 100 actors and 1000 events with
#' #only inertia and reciprocity
#'rem2 <- simulateRESeq(n_actors = 100,
#'                      n_events = 1000,
#'                      inertia = TRUE,
#'                      inertia_p = 0.12,
#'                      recip = TRUE,
#'                      recip_p = 0.08,
#'                      returnStats = TRUE)
#'rem2
#'
#'#Creating a random relational sequence based on the starting sequence with
#'#only inertia and reciprocity
#'rem3 <- simulateRESeq(n_actors = 100, #does not matter can be any value, this is
#'                                     #overridden by the starting event sequence
#'                     n_events = 100,
#'                     inertia = TRUE,
#'                     inertia_p = 0.12,
#'                     recip = TRUE,
#'                     recip_p = 0.08,
#'                     #a random starting event sequence
#'                     starting_events = matrix(c(1:10, 10:1),
#'                     nrow = 10, ncol = 2,  byrow = FALSE),
#'                     returnStats = TRUE)
#'rem3


simulateRESeq <- function(n_actors,
                         n_events,
                         inertia = FALSE,
                         inertia_p = 0,
                         recip = FALSE,
                         recip_p = 0,
                         sender_outdegree = FALSE,
                         sender_outdegree_p = 0,
                         sender_indegree = FALSE,
                         sender_indegree_p = 0,
                         target_outdegree = FALSE,
                         target_outdegree_p = 0,
                         target_indegree = FALSE,
                         target_indegree_p = 0,
                         assort = FALSE,
                         assort_p = 0,
                         trans_trips = FALSE,
                         trans_trips_p = 0,
                         three_cycles = FALSE,
                         three_cycles_p = 0,
                         starting_events = NULL,
                         returnStats = FALSE){
  lifecycle::deprecate_warn("1.0.0", " simulateRESeq()", "simulate_rem_seq()")
  ##################################################################################
  #
  #         Check User Inputs
  #
  ##################################################################################
  if(inertia == FALSE &
     sender_outdegree == FALSE &
     target_outdegree == FALSE &
     recip == FALSE &
     sender_indegree == FALSE &
     target_indegree == FALSE &
     assort == FALSE &
     trans_trips == FALSE &
     three_cycles == FALSE){
    base::stop("No effects were requested. Therefore, stoping the simulation.")
  }
  if(trans_trips== TRUE & trans_trips_p == 0){
    base::stop("You requested the transitive triplet effect be included, however, specificed
               a parameter of 0. This indicates that transitive triplet has no effect on the
               event probabilties. Please either update the trans_trips_p parameter or
               make trans_trips = FALSE.")
  }
  if(three_cycles == TRUE & three_cycles_p == 0){
    base::stop("You requested the three cycles effect be included, however, specificed
               a parameter of 0. This indicates that transitive triplet has no effect on the
               event probabilties. Please either update the three_cycles_p parameter or
               make three_cycles = FALSE.")
  }
  if(assort== TRUE & assort_p == 0){
    base::stop("You requested the assortativity effect be included, however, specificed
               a parameter of 0. This indicates that assortativity  has no effect on the
               event probabilties. Please either update the assort_p parameter or
               make assort = FALSE.")
  }

  if(sender_indegree == TRUE & sender_indegree_p == 0){
    base::stop("You requested the sender indegree effect be included, however, specificed
               a parameter of 0. This indicates that indegree has no effect on the
               event probabilties. Please either update the sender_indegree_p parameter or
               make sender_indegree = FALSE.")
  }
  if(sender_outdegree == TRUE & sender_outdegree_p == 0){
    base::stop("You requested the sender outdegree effect be included, however, specificed
               a parameter of 0. This indicates that outdegree has no effect on the
               event probabilties. Please either update the sender_outdegree_p parameter or
               make sender_outdegree = FALSE.")
  }

  if(target_indegree == TRUE & target_indegree_p == 0){
    base::stop("You requested the target indegree effect be included, however, specificed
               a parameter of 0. This indicates that indegree has no effect on the
               event probabilties. Please either update the target_indegree_p parameter or
               make target_indegree = FALSE.")
  }
  if(target_outdegree == TRUE & target_outdegree_p == 0){
    base::stop("You requested the target outdegree effect be included, however, specificed
               a parameter of 0. This indicates that outdegree has no effect on the
               event probabilties. Please either update the target_outdegree_p parameter or
               make target_outdegree = FALSE.")
  }

  if(recip == TRUE & recip_p == 0){
    base::stop("You requested the reciprocity effect be included, however, specificed
               a parameter of 0. This indicates that reciprocity has no effect on the
               event probabilties. Please either update the recip_p parameter or
               make recip = FALSE.")
  }
  if(inertia == TRUE & inertia_p == 0){
    base::stop("You requested the inertia effect be included, however, specificed
               a parameter of 0. This indicates that inertia has no effect on the
               event probabilties. Please either update the inertia_p parameter or
               make inertia = FALSE.")
  }
  if(!is.null(starting_events)){
    if(nrow(starting_events) >= n_events){
      base::stop("You included a starting event list, however, this event list is
               greater than or equal to the number of simulated events that you
               requested. Please update this to either increase the number of
               simulated events or reduce the size of the starting events.")
    }
    if(ncol(starting_events) != 2){
      #there should only be two columns: the event senders (column 1) and the event targets (column 2)
      base::stop("You included a starting event list, however, the number of columns in this
                 event list is not equal to 2. Please update this such that the first
                 column is the starting event senders and the second column is the starting
                 event targets.")
    }
    if(!is.numeric(starting_events[,1]) | !is.numeric(starting_events[,2])){
      #to assist with the simulation procedure, both of the event targets and senders id
      #should be numeric (given we have to create a full possible risk set with new actors)
      base::stop("The starting events senders and targets should be numeric vectors, that is,
                 the respective actor ids should be numeric values. This is requested since
                 we have to create new actor ids for the full risk set. One solution would be to
                 make the non-numeric values as factor variable then force it to be numeric.")
    }

  }

  if(is.null(starting_events)){

    #intializing the risk set, namely, the relational events that could occur
    risk_set <- data.table::CJ(sender = 1:n_actors,
                               target = 1:n_actors)
    #removing self loops from the events, by assumption, the target and sender cannot be the same
    risk_set <- risk_set[risk_set$sender != risk_set$target]
    possible_events <- nrow(risk_set)

    #intializing the first event
    ei <- sample(1:possible_events, 1) #randomly sample the first event
    #intializing the first event
    events <- data.table::data.table(
      eventID = 1,
      sender = risk_set$sender[ei],
      target = risk_set$target[ei])
    #vector of sufficient statistics
    inertiaSS <- numeric(possible_events)
    reciprocitySS <- numeric(possible_events)
    sender_indegreeSS <- numeric(possible_events)
    sender_outdegreeSS <- numeric(possible_events)
    target_indegreeSS <- numeric(possible_events)
    target_outdegreeSS <- numeric(possible_events)
    assortSS <- numeric(possible_events)
    trans_tripsSS <- numeric(possible_events)
    three_cyclesSS <- numeric(possible_events)
    #vector of sufficient statistics for selected events
    inertiaEF <- numeric(n_events)
    reciprocityEF <-numeric(n_events)
    sender_indegreeEF <-numeric(n_events)
    sender_outdegreeEF <-numeric(n_events)
    target_indegreeEF <-numeric(n_events)
    target_outdegreeEF <-numeric(n_events)
    assortEF <- numeric(n_events)
    trans_tripsEF <- numeric(n_events)
    three_cyclesEF <- numeric(n_events)
    if(trans_trips == TRUE){ #precreating the count matrix
      trans_mat <- matrix(0, #a matrix of 0s
                          ncol = n_actors, #with n rows
                          nrow = n_actors) #and n columns
    }
    if(three_cycles == TRUE){ #precreating the count matrix
      cycles_mat <- matrix(0, #a matrix of 0s
                           ncol = n_actors, #with n rows
                           nrow = n_actors) #and n columns
    }


    #initializing the past relational history
    History <- events
    #starting the simulation
    for(i in 2:n_events){ #starting at the 2nd event to the n event
      A_past <- History #the past relational history
      #
      #   No need to go through each possible event, just simply look at
      #   what values need to be updated!
      #
      #

      #most recent event
      mrevent <- History[(i-1)] #the most recent updated event

      ################################################################################
      #     Compute past inertia values
      ################################################################################
      if(inertia == TRUE){
        senderUpdate <- collapse::whichv(risk_set$sender,   mrevent$sender   )
        targetUpdate <- collapse::whichv(risk_set$target,   mrevent$target   )
        eventsUpdate <- base::intersect(senderUpdate, targetUpdate)
        if(length(eventsUpdate) != 0){
          inertiaSS[eventsUpdate] <- (inertiaSS[eventsUpdate]  + 1)
        }

      }
      ################################################################################
      #     Compute past reciprocity values
      ################################################################################
      if(recip == TRUE){
        targetUpdate <- collapse::whichv(risk_set$sender,   mrevent$target   )
        senderUpdate <- collapse::whichv(risk_set$target,   mrevent$sender   )
        eventsUpdate <- base::intersect(senderUpdate, targetUpdate)
        if(length(eventsUpdate) != 0){
          reciprocitySS[eventsUpdate] <- ( reciprocitySS[eventsUpdate]  + 1)
        }

      }
      ################################################################################
      #     Compute past target indegree values
      ################################################################################
      if(target_indegree == TRUE){
        targetUpdate <- collapse::whichv(risk_set$target,   mrevent$target   )
        if(length(targetUpdate) != 0){
          target_indegreeSS[ targetUpdate ] <- ( target_indegreeSS[ targetUpdate ]  + 1)
        }
      }
      ################################################################################
      #     Compute past sender indegree values
      ################################################################################
      if(sender_indegree == TRUE){
        targetUpdate <- collapse::whichv(risk_set$sender,   mrevent$target   )
        if(length(targetUpdate) != 0){
          sender_indegreeSS[ targetUpdate ] <- ( sender_indegreeSS[ targetUpdate ]  + 1)
        }
      }
      ################################################################################
      #     Compute past sender outdegree values
      ################################################################################
      if(sender_outdegree == TRUE){
        senderUpdate <- collapse::whichv(risk_set$sender,   mrevent$sender   )
        if(length(senderUpdate ) != 0){
          sender_outdegreeSS[ senderUpdate  ] <- ( sender_outdegreeSS[ senderUpdate  ]  + 1)
        }

      }
      ################################################################################
      #     Compute past target outdegree values
      ################################################################################
      if(target_outdegree == TRUE){
        senderUpdate <- collapse::whichv(risk_set$target,   mrevent$sender   )
        if(length(senderUpdate ) != 0){
          target_outdegreeSS[ senderUpdate  ] <- ( target_outdegreeSS[ senderUpdate  ]  + 1)
        }

      }
      ################################################################################
      #     Compute past assortativity values
      ################################################################################
      if(assort == TRUE){

        if(target_indegree == TRUE & sender_outdegree == TRUE){
          assortSS <- target_indegreeSS * sender_outdegreeSS
        }

        if(target_indegree == FALSE & sender_outdegree == TRUE){

          targetUpdate <- collapse::whichv(risk_set$target,   mrevent$target   )
          if(length(targetUpdate) != 0){
            target_indegreeSS[ targetUpdate ] <- ( target_indegreeSS[ targetUpdate ]  + 1)
          }
          assortSS <- target_indegreeSS * sender_outdegreeSS
        }

        if(target_indegree == TRUE & sender_outdegree == FALSE){

          senderUpdate <- collapse::whichv(risk_set$sender,   mrevent$sender   )
          if(length(senderUpdate ) != 0){
            sender_outdegreeSS[ senderUpdate  ] <- ( sender_outdegreeSS[ senderUpdate  ]  + 1)
          }
          assortSS <- target_indegreeSS * sender_outdegreeSS
        }

        if(target_indegree == FALSE & sender_outdegree == FALSE){
          targetUpdate <- collapse::whichv(risk_set$target,   mrevent$target   )
          if(length(targetUpdate) != 0){
            target_indegreeSS[ targetUpdate ] <- ( target_indegreeSS[ targetUpdate ]  + 1)
          }
          senderUpdate <- collapse::whichv(risk_set$sender,   mrevent$sender   )
          if(length(senderUpdate ) != 0){
            sender_outdegreeSS[ senderUpdate  ] <- ( sender_outdegreeSS[ senderUpdate  ]  + 1)
          }
          assortSS <- target_indegreeSS * sender_outdegreeSS
        }

      }
      ################################################################################
      #     Compute past transitive triplets value
      ################################################################################
      if(trans_trips == TRUE){
        trans_mat[mrevent$sender, mrevent$target] <- 1 +  trans_mat[mrevent$sender, mrevent$target]
        updateT <- trans_mat %*% trans_mat #getting the full product of the matrices
        #(note these are directed ties and should be
        # treated as such)
        diag(updateT) <- NA #making the diagonal NA, using this for later subsetting
        values <- base::as.numeric(base::t(updateT)) #transposing the matrix and then extracting it's values
        trans_tripsSS <- values[!is.na(values)] #removing the missing values from the matrix
        #making this the new network matrix
      }

      if(three_cycles == TRUE){
        cycles_mat[mrevent$sender, mrevent$target] <- 1 +  cycles_mat[mrevent$sender, mrevent$target]
        updateT <- trans_mat %*% trans_mat #getting the full product of the matrices
        updateT <- t(updateT) #transposing the matrix, since currently, the values are
        #j ->i and need to be i -> j
        #(note these are directed ties and should be
        # treated as such)
        diag(updateT) <- NA #making the diagonal NA, using this for later subsetting
        values <- base::as.numeric(base::t(updateT)) #transposing the matrix and then extracting it's values
        three_cyclesSS <- values[!is.na(values)] #removing the missing values from the matrix
        #making this the new network matrix
      }



      ################################################################################
      #     Randomly sample the next event
      ################################################################################
      SS_ij <- matrix(c(inertiaSS,
                        reciprocitySS,
                        sender_indegreeSS,
                        sender_outdegreeSS,
                        target_indegreeSS,
                        target_outdegreeSS,
                        assortSS,
                        trans_tripsSS,
                        three_cyclesSS),
                      nrow = possible_events,
                      ncol = 9)
      rem_effects <- c(inertia_p,
                       recip_p,
                       sender_indegree_p,
                       sender_outdegree_p,
                       target_indegree_p,
                       target_outdegree_p,
                       assort_p,
                       trans_trips_p,
                       three_cycles_p)
      hazard <- exp(SS_ij%*%rem_effects)   #choosing the dyadic event probabilisticly
      pij <- hazard / sum(hazard) #probability of selection
      ei <- sample(1:nrow(risk_set),  #sample from the selection of events
                   1,  #one event
                   prob = pij) #with probability pij
      events <-  data.table::rbindlist(list(events,
                                            data.table::data.table(
                                              eventID = i,
                                              sender = risk_set$sender[ei],
                                              target = risk_set$target[ei])))
      inertiaEF[i] <- inertiaSS[ei]
      reciprocityEF[i] <- reciprocitySS[ei]
      sender_indegreeEF[i] <- sender_indegreeSS[ei]
      sender_outdegreeEF[i] <- sender_outdegreeSS[ei]
      target_indegreeEF[i] <-  target_indegreeSS[ei]
      target_outdegreeEF[i] <-  target_outdegreeSS[ei]
      assortEF[i] <- assortSS[ei]
      trans_tripsEF[i] <- trans_tripsSS[ei]
      three_cyclesEF[i] <- three_cyclesSS[ei]


      History <- events
    }
    if(returnStats == TRUE){
      if(inertia == TRUE){
        events$inertia <- inertiaEF
      }
      if(recip == TRUE){
        events$reciprocity <- reciprocityEF
      }
      if(sender_indegree == TRUE){
        events$sender_indegree <- sender_indegreeEF
      }
      if(sender_outdegree == TRUE){
        events$sender_outdegree <- sender_outdegreeEF
      }
      if(target_indegree == TRUE){
        events$target_indegree <- target_indegreeEF
      }
      if(target_outdegree == TRUE){
        events$target_outdegree <- target_outdegreeEF
      }
      if(assort == TRUE){
        events$assort <- assortEF
      }
      if(trans_trips == TRUE){
        events$trans_trip <- trans_tripsEF
      }
      if(three_cycles == TRUE){
        events$three_cycles <-  three_cyclesEF
      }
    }
  }

  if(!is.null(starting_events)){

    #checking the number of actors within the dataset, if the starting events are used
    #we only use those actors that are involved in the starting event dataset, that is,
    #we do not generate new actors.

    n_actors <- length(unique(c(starting_events[,1],starting_events[,2])))
    #intializing the risk set, namely, the relational events that could occur
    risk_set <- data.table::CJ(sender = 1:n_actors,
                               target = 1:n_actors)
    #removing self loops from the events, by assumption, the target and sender cannot be the same
    risk_set <- risk_set[risk_set$sender != risk_set$target]
    possible_events <- nrow(risk_set)

    #intializing the first event
    events <- data.table::data.table(
      eventID = 1:nrow(starting_events),
      sender = starting_events[,1],
      target = starting_events[,2])
    #vector of sufficient statistics
    inertiaSS <- numeric(possible_events)
    reciprocitySS <- numeric(possible_events)
    sender_indegreeSS <- numeric(possible_events)
    sender_outdegreeSS <- numeric(possible_events)
    target_indegreeSS <- numeric(possible_events)
    target_outdegreeSS <- numeric(possible_events)
    assortSS <- numeric(possible_events)
    trans_tripsSS <- numeric(possible_events)
    three_cyclesSS <- numeric(possible_events)
    #vector of sufficient statistics for selected events
    inertiaEF <- numeric(n_events)
    reciprocityEF <-numeric(n_events)
    sender_indegreeEF <-numeric(n_events)
    sender_outdegreeEF <-numeric(n_events)
    target_indegreeEF <-numeric(n_events)
    target_outdegreeEF <-numeric(n_events)
    assortEF <- numeric(n_events)
    trans_tripsEF <- numeric(n_events)
    three_cyclesEF <- numeric(n_events)
    if(trans_trips == TRUE){ #precreating the count matrix
      trans_mat <- matrix(0, #a matrix of 0s
                          ncol = n_actors, #with n rows
                          nrow = n_actors) #and n columns
    }
    if(three_cycles == TRUE){ #precreating the count matrix
      cycles_mat <- matrix(0, #a matrix of 0s
                           ncol = n_actors, #with n rows
                           nrow = n_actors) #and n columns
    }


    #initializing the past relational history
    History <- starting_events
    #creating a new dataset to compute the statistics from
    starting_eventsC <- data.table::data.table(eventID = 1:nrow(starting_events),
                                               sender = starting_events[,1],
                                               target = starting_events[,2])
    # Need to compute the relevant values for the statistics
    #computing the values for the past statistics
    for(i in 2:nrow(starting_eventsC)){ #starting at the 2nd event to the n event
      #
      #   No need to go through each possible event, just simply look at
      #   what values need to be updated!
      #
      #

      #most recent event
      mrevent <- starting_eventsC[(i-1)] #the most recent updated event

      ################################################################################
      #     Compute past inertia values
      ################################################################################
      if(inertia == TRUE){
        senderUpdate <- collapse::whichv(risk_set$sender,   mrevent$sender   )
        targetUpdate <- collapse::whichv(risk_set$target,   mrevent$target   )
        eventsUpdate <- base::intersect(senderUpdate, targetUpdate)
        if(length(eventsUpdate) != 0){
          inertiaSS[eventsUpdate] <- (inertiaSS[eventsUpdate]  + 1)
        }

      }
      ################################################################################
      #     Compute past reciprocity values
      ################################################################################
      if(recip == TRUE){
        targetUpdate <- collapse::whichv(risk_set$sender,   mrevent$target   )
        senderUpdate <- collapse::whichv(risk_set$target,   mrevent$sender   )
        eventsUpdate <- base::intersect(senderUpdate, targetUpdate)
        if(length(eventsUpdate) != 0){
          reciprocitySS[eventsUpdate] <- ( reciprocitySS[eventsUpdate]  + 1)
        }

      }
      ################################################################################
      #     Compute past target indegree values
      ################################################################################
      if(target_indegree == TRUE){
        targetUpdate <- collapse::whichv(risk_set$target,   mrevent$target   )
        if(length(targetUpdate) != 0){
          target_indegreeSS[ targetUpdate ] <- ( target_indegreeSS[ targetUpdate ]  + 1)
        }
      }
      ################################################################################
      #     Compute past sender indegree values
      ################################################################################
      if(sender_indegree == TRUE){
        targetUpdate <- collapse::whichv(risk_set$sender,   mrevent$target   )
        if(length(targetUpdate) != 0){
          sender_indegreeSS[ targetUpdate ] <- ( sender_indegreeSS[ targetUpdate ]  + 1)
        }
      }
      ################################################################################
      #     Compute past sender outdegree values
      ################################################################################
      if(sender_outdegree == TRUE){
        senderUpdate <- collapse::whichv(risk_set$sender,   mrevent$sender   )
        if(length(senderUpdate ) != 0){
          sender_outdegreeSS[ senderUpdate  ] <- ( sender_outdegreeSS[ senderUpdate  ]  + 1)
        }

      }
      ################################################################################
      #     Compute past target outdegree values
      ################################################################################
      if(target_outdegree == TRUE){
        senderUpdate <- collapse::whichv(risk_set$target,   mrevent$sender   )
        if(length(senderUpdate ) != 0){
          target_outdegreeSS[ senderUpdate  ] <- ( target_outdegreeSS[ senderUpdate  ]  + 1)
        }

      }
      ################################################################################
      #     Compute past assortativity values
      ################################################################################
      if(assort == TRUE){

        if(target_indegree == TRUE & sender_outdegree == TRUE){
          assortSS <- target_indegreeSS * sender_outdegreeSS
        }

        if(target_indegree == FALSE & sender_outdegree == TRUE){

          targetUpdate <- collapse::whichv(risk_set$target,   mrevent$target   )
          if(length(targetUpdate) != 0){
            target_indegreeSS[ targetUpdate ] <- ( target_indegreeSS[ targetUpdate ]  + 1)
          }
          assortSS <- target_indegreeSS * sender_outdegreeSS
        }

        if(target_indegree == TRUE & sender_outdegree == FALSE){

          senderUpdate <- collapse::whichv(risk_set$sender,   mrevent$sender   )
          if(length(senderUpdate ) != 0){
            sender_outdegreeSS[ senderUpdate  ] <- ( sender_outdegreeSS[ senderUpdate  ]  + 1)
          }
          assortSS <- target_indegreeSS * sender_outdegreeSS
        }

        if(target_indegree == FALSE & sender_outdegree == FALSE){
          targetUpdate <- collapse::whichv(risk_set$target,   mrevent$target   )
          if(length(targetUpdate) != 0){
            target_indegreeSS[ targetUpdate ] <- ( target_indegreeSS[ targetUpdate ]  + 1)
          }
          senderUpdate <- collapse::whichv(risk_set$sender,   mrevent$sender   )
          if(length(senderUpdate ) != 0){
            sender_outdegreeSS[ senderUpdate  ] <- ( sender_outdegreeSS[ senderUpdate  ]  + 1)
          }
          assortSS <- target_indegreeSS * sender_outdegreeSS
        }

      }
      ################################################################################
      #     Compute past transitive triplets value
      ################################################################################
      if(trans_trips == TRUE){
        trans_mat[mrevent$sender, mrevent$target] <- 1 +  trans_mat[mrevent$sender, mrevent$target]
        updateT <- trans_mat %*% trans_mat #getting the full product of the matrices
        #(note these are directed ties and should be
        # treated as such)
        diag(updateT) <- NA #making the diagonal NA, using this for later subsetting
        values <- base::as.numeric(base::t(updateT)) #transposing the matrix and then extracting it's values
        trans_tripsSS <- values[!is.na(values)] #removing the missing values from the matrix
        #making this the new network matrix
      }

      if(three_cycles == TRUE){
        cycles_mat[mrevent$sender, mrevent$target] <- 1 +  cycles_mat[mrevent$sender, mrevent$target]
        updateT <- cycles_mat %*% cycles_mat #getting the full product of the matrices
        updateT <- t(updateT) #transposing the matrix, since currently, the values are
        #j ->i and need to be i -> j
        #(note these are directed ties and should be
        # treated as such)
        diag(updateT) <- NA #making the diagonal NA, using this for later subsetting
        values <- base::as.numeric(base::t(updateT)) #transposing the matrix and then extracting it's values
        three_cyclesSS <- values[!is.na(values)] #removing the missing values from the matrix
        #making this the new network matrix
      }

      ei <- which(risk_set$sender == starting_eventsC[i]$sender &
                    risk_set$target == starting_eventsC[i]$target)
      inertiaEF[i] <- inertiaSS[ei]
      reciprocityEF[i] <- reciprocitySS[ei]
      sender_indegreeEF[i] <- sender_indegreeSS[ei]
      sender_outdegreeEF[i] <- sender_outdegreeSS[ei]
      target_indegreeEF[i] <-  target_indegreeSS[ei]
      target_outdegreeEF[i] <-  target_outdegreeSS[ei]
      assortEF[i] <- assortSS[ei]
      trans_tripsEF[i] <- trans_tripsSS[ei]
      three_cyclesEF[i] <- three_cyclesSS[ei]

    }

    ### Now we need to start sampling new events based on probabilites
    starthere <- nrow(starting_eventsC) + 1
    #initializing the past relational history
    History <- events
    for(i in starthere:n_events){
      #starting at the 2nd event to the n event
      A_past <- History #the past relational history
      #
      #   No need to go through each possible event, just simply look at
      #   what values need to be updated!
      #
      #

      #most recent event
      mrevent <- History[(i-1)] #the most recent updated event

      ################################################################################
      #     Compute past inertia values
      ################################################################################
      if(inertia == TRUE){
        senderUpdate <- collapse::whichv(risk_set$sender,   mrevent$sender   )
        targetUpdate <- collapse::whichv(risk_set$target,   mrevent$target   )
        eventsUpdate <- base::intersect(senderUpdate, targetUpdate)
        if(length(eventsUpdate) != 0){
          inertiaSS[eventsUpdate] <- (inertiaSS[eventsUpdate]  + 1)
        }

      }
      ################################################################################
      #     Compute past reciprocity values
      ################################################################################
      if(recip == TRUE){
        targetUpdate <- collapse::whichv(risk_set$sender,   mrevent$target   )
        senderUpdate <- collapse::whichv(risk_set$target,   mrevent$sender   )
        eventsUpdate <- base::intersect(senderUpdate, targetUpdate)
        if(length(eventsUpdate) != 0){
          reciprocitySS[eventsUpdate] <- ( reciprocitySS[eventsUpdate]  + 1)
        }

      }
      ################################################################################
      #     Compute past target indegree values
      ################################################################################
      if(target_indegree == TRUE){
        targetUpdate <- collapse::whichv(risk_set$target,   mrevent$target   )
        if(length(targetUpdate) != 0){
          target_indegreeSS[ targetUpdate ] <- ( target_indegreeSS[ targetUpdate ]  + 1)
        }
      }
      ################################################################################
      #     Compute past sender indegree values
      ################################################################################
      if(sender_indegree == TRUE){
        targetUpdate <- collapse::whichv(risk_set$sender,   mrevent$target   )
        if(length(targetUpdate) != 0){
          sender_indegreeSS[ targetUpdate ] <- ( sender_indegreeSS[ targetUpdate ]  + 1)
        }
      }
      ################################################################################
      #     Compute past sender outdegree values
      ################################################################################
      if(sender_outdegree == TRUE){
        senderUpdate <- collapse::whichv(risk_set$sender,   mrevent$sender   )
        if(length(senderUpdate ) != 0){
          sender_outdegreeSS[ senderUpdate  ] <- ( sender_outdegreeSS[ senderUpdate  ]  + 1)
        }

      }
      ################################################################################
      #     Compute past target outdegree values
      ################################################################################
      if(target_outdegree == TRUE){
        senderUpdate <- collapse::whichv(risk_set$target,   mrevent$sender   )
        if(length(senderUpdate ) != 0){
          target_outdegreeSS[ senderUpdate  ] <- ( target_outdegreeSS[ senderUpdate  ]  + 1)
        }

      }
      ################################################################################
      #     Compute past assortativity values
      ################################################################################
      if(assort == TRUE){

        if(target_indegree == TRUE & sender_outdegree == TRUE){
          assortSS <- target_indegreeSS * sender_outdegreeSS
        }

        if(target_indegree == FALSE & sender_outdegree == TRUE){

          targetUpdate <- collapse::whichv(risk_set$target,   mrevent$target   )
          if(length(targetUpdate) != 0){
            target_indegreeSS[ targetUpdate ] <- ( target_indegreeSS[ targetUpdate ]  + 1)
          }
          assortSS <- target_indegreeSS * sender_outdegreeSS
        }

        if(target_indegree == TRUE & sender_outdegree == FALSE){

          senderUpdate <- collapse::whichv(risk_set$sender,   mrevent$sender   )
          if(length(senderUpdate ) != 0){
            sender_outdegreeSS[ senderUpdate  ] <- ( sender_outdegreeSS[ senderUpdate  ]  + 1)
          }
          assortSS <- target_indegreeSS * sender_outdegreeSS
        }

        if(target_indegree == FALSE & sender_outdegree == FALSE){
          targetUpdate <- collapse::whichv(risk_set$target,   mrevent$target   )
          if(length(targetUpdate) != 0){
            target_indegreeSS[ targetUpdate ] <- ( target_indegreeSS[ targetUpdate ]  + 1)
          }
          senderUpdate <- collapse::whichv(risk_set$sender,   mrevent$sender   )
          if(length(senderUpdate ) != 0){
            sender_outdegreeSS[ senderUpdate  ] <- ( sender_outdegreeSS[ senderUpdate  ]  + 1)
          }
          assortSS <- target_indegreeSS * sender_outdegreeSS
        }

      }
      ################################################################################
      #     Compute past transitive triplets value
      ################################################################################
      if(trans_trips == TRUE){
        trans_mat[mrevent$sender, mrevent$target] <- 1 +  trans_mat[mrevent$sender, mrevent$target]
        updateT <- trans_mat %*% trans_mat #getting the full product of the matrices
        #(note these are directed ties and should be
        # treated as such)
        diag(updateT) <- NA #making the diagonal NA, using this for later subsetting
        values <- base::as.numeric(base::t(updateT)) #transposing the matrix and then extracting it's values
        trans_tripsSS <- values[!is.na(values)] #removing the missing values from the matrix
        #making this the new network matrix
      }

      if(three_cycles == TRUE){
        cycles_mat[mrevent$sender, mrevent$target] <- 1 +  cycles_mat[mrevent$sender, mrevent$target]
        updateT <- trans_mat %*% trans_mat #getting the full product of the matrices
        updateT <- t(updateT) #transposing the matrix, since currently, the values are
        #j ->i and need to be i -> j
        #(note these are directed ties and should be
        # treated as such)
        diag(updateT) <- NA #making the diagonal NA, using this for later subsetting
        values <- base::as.numeric(base::t(updateT)) #transposing the matrix and then extracting it's values
        three_cyclesSS <- values[!is.na(values)] #removing the missing values from the matrix
        #making this the new network matrix
      }



      ################################################################################
      #     Randomly sample the next event
      ################################################################################
      SS_ij <- matrix(c(inertiaSS,
                        reciprocitySS,
                        sender_indegreeSS,
                        sender_outdegreeSS,
                        target_indegreeSS,
                        target_outdegreeSS,
                        assortSS,
                        trans_tripsSS,
                        three_cyclesSS),
                      nrow = possible_events,
                      ncol = 9)
      rem_effects <- c(inertia_p,
                       recip_p,
                       sender_indegree_p,
                       sender_outdegree_p,
                       target_indegree_p,
                       target_outdegree_p,
                       assort_p,
                       trans_trips_p,
                       three_cycles_p)
      hazard <- exp(SS_ij%*%rem_effects)   #choosing the dyadic event probabilisticly
      pij <- hazard / sum(hazard) #probability of selection
      ei <- sample(1:nrow(risk_set),  #sample from the selection of events
                   1,  #one event
                   prob = pij) #with probability pij
      events <-  data.table::rbindlist(list(events,
                                            data.table::data.table(
                                              eventID = i,
                                              sender = risk_set$sender[ei],
                                              target = risk_set$target[ei])))
      inertiaEF[i] <- inertiaSS[ei]
      reciprocityEF[i] <- reciprocitySS[ei]
      sender_indegreeEF[i] <- sender_indegreeSS[ei]
      sender_outdegreeEF[i] <- sender_outdegreeSS[ei]
      target_indegreeEF[i] <-  target_indegreeSS[ei]
      target_outdegreeEF[i] <-  target_outdegreeSS[ei]
      assortEF[i] <- assortSS[ei]
      trans_tripsEF[i] <- trans_tripsSS[ei]
      three_cyclesEF[i] <- three_cyclesSS[ei]


      History <- events
    }
    if(returnStats == TRUE){
      if(inertia == TRUE){
        events$inertia <- inertiaEF
      }
      if(recip == TRUE){
        events$reciprocity <- reciprocityEF
      }
      if(sender_indegree == TRUE){
        events$sender_indegree <- sender_indegreeEF
      }
      if(sender_outdegree == TRUE){
        events$sender_outdegree <- sender_outdegreeEF
      }
      if(target_indegree == TRUE){
        events$target_indegree <- target_indegreeEF
      }
      if(target_outdegree == TRUE){
        events$target_outdegree <- target_outdegreeEF
      }
      if(assort == TRUE){
        events$assort <- assortEF
      }
      if(trans_trips == TRUE){
        events$trans_trip <- trans_tripsEF
      }
      if(three_cycles == TRUE){
        events$three_cycles <-  three_cyclesEF
      }
    }
  }
  return(events)
}





