## Code written by Kevin Carson (kacarson@arizona.edu) and Deigo Leal (https://www.diegoleal.info/)
## Last Updated: 07-15-24

#' @title Compute Butts' (2008) Recency Network Statistic for Event Dyads in a Relational Event Sequence
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `computeRecency()` has been deprecated starting on version 1.0.0 of the `dream` package. Please use the `remstats_recency()` function and see the `NEWS.md` file for more details.
#'
#'
#' This function computes the recency network sufficient statistic for
#' a relational event sequence (see Butts 2008; Vu et al. 2015; Meijerink-Bosman et al. 2022). The recency statistic
#' captures the tendency in which more recent events (i.e., an exchange between two medical doctors) are more
#' likely to reoccur in comparison to events that happened in the distant past (see Butts 2008 for a discussion). This measure allows for recency scores to be only
#' computed for the sampled events, while creating the statistics based on the full event
#' sequence. Moreover, the function allows users to specify relational relevancy for the statistic and
#' employ a sliding windows framework for large relational sequences.
#' @name computeRecency
#' @param observed_time The vector of event times from the pre-processing event sequence.
#' @param observed_sender The vector of event senders from the pre-processing event sequence.
#' @param observed_receiver The vector of event receivers from the pre-processing event sequence
#' @param processed_time The vector of event times from the post-processing event sequence (i.e., the event sequence that contains the observed and null events).
#' @param processed_sender The vector of event senders from the post-processing event sequence (i.e., the event sequence that contains the observed and null events).
#' @param processed_receiver The vector of event receivers from the post-processing event sequence (i.e., the event sequence that contains the observed and null events).
#' @param type A string value that specifies which recency formula will be used to compute the statistics. The options are "raw.diff", "inv.diff.plus1", "rank.ordered.count" (see the details section).
#' @param i_neighborhood TRUE/FALSE. TRUE indicates that the recency statistic will be computed in reference to the sender’s past relational history (see details section). FALSE indicates that the persistence statistic will be computed in reference to the target’s past relational history (see details section). Set to TRUE by default.
#' @param nopastEvents The numerical value that specifies what value should be given to events in which the sender has sent not past ties (i's neighborhood when i_neighborhood = TRUE) or has not received any past ties (j's neighborhood when i_neighborhood = FALSE). Set to NA by default.
#' @param dependency TRUE/FALSE. TRUE indicates that temporal relevancy will be modeled (see the details section). FALSE indicates that temporal relevancy will not be modeled, that is, all past events are relevant (see the details section). Set to FALSE by default.
#' @param relationalTimeSpan If dependency = TRUE, a numerical value that corresponds to the temporal span for relational relevancy, which must be the same measurement unit as the observed_time and processed_time objects. When dependency = TRUE, the relevant events are events that have occurred between current event time, *t*, and *t - relationalTimeSpan*. For example, if the time measurement is the number of days since the first event and the value for relationalTimeSpan is set to 10, then only those events which occurred in the past 10 days are included in the computation of the statistic.
#' @param sliding_windows TRUE/FALSE. TRUE indicates that the sliding windows computational approach will
#' be used to compute the resulting network statistic, while FALSE indicates the approach will not be used. Set
#' to FALSE by default. It’s important to note that the sliding windows framework should only be used
#' when the pre-processed event sequence is ‘big’, such as the 360 million pre-processed event sequence
#' used in Lerner and Lomi (2020), as it aims to reduce the computational burden of sorting ‘big’ datasets. In general,
#' most pre-processed event sequences will not need to use the sliding windows
#' approach. There is not a strict cutoff for ‘big’ dataset. This definition depends on both the
#' size of the observed event sequence and the post-processing sampling dataset. For instance,
#' according to our internal tests, when the event sequence is relatively large (i.e., 100,000
#' observed events) with probability of sampling from the observed event sequence set to 0.05
#' and using 10 controls per sampled event, the sliding windows framework for computing repetition
#' is about 11% faster than the non-sliding windows framework. Yet, in a smaller dataset
#' (i.e., 10,000 observed events) the sliding windows framework is about 25% slower than the
#' non-sliding framework with the same conditions as before.
#' @param window_size If sliding_windows is set to TRUE, the sizes of the windows that are used for the sliding windows computational framework. If NA, the function internally divides the dataset into ten slices (may not be optimal).
#' @param processed_seqIDs If sliding_windows is set to TRUE, the vector of event sequence IDs from the post-processing event sequence. The event sequence IDs represents the index for when the event occurred in the observed event sequence (e.g., the 5th event in the sequence will have a value of 5 in this vector).
#' @import data.table
#' @importFrom collapse whichv
#' @return The vector of recency network statistics for the relational event sequence.
#' @export
#' @details This function calculates the recency network sufficient statistic for a relational event based on
#' Butts (2008), Vu et al. (2015), or Meijerink-Bosman et al. (2022).
#' Depending on the type and neighborhood requested, different formulas will be used.
#'
#' In the below equations, when *i_neighborhood* is TRUE:
#'  \deqn{t^{*} = max(t \in \left\{(s',r',t') \in E : s'= s \land r'= r  \land t'<t \right\}) }
#'
#' When *i_neighborhood* is FALSE, the following formula is used:
#'  \deqn{t^{*} = max(t \in \left\{(s',r',t') \in E : s'= r \land r'= s  \land t'<t \right\}) }
#'
#' The formula for recency for event \eqn{e_i} with type set to "raw.diff" and *i_neighborhood* is TRUE (Vu et al. 2015):
#' \deqn{recency_{e_i} = t_{e_i} - t^{*} }
#' where \eqn{t^{*}}, is the most recent time in
#' which the past event has the same receiver and sender as the current event. If there are no past events within the current dyad, then
#' the value defaults to the *nopastEvents* argument.
#'
#' The formula for recency for event \eqn{e_i} with type set to "raw.diff" and *i_neighborhood* is FALSE (Vu et al. 2015):
#' \deqn{recency_{e_i} = t_{e_i} - t^{*}    }
#' where \eqn{t^{*}}, is the most recent time in
#' which the past event's sender is the current event receiver and the past event receiver is the current event sender.  If there are no past events within the current dyad, then
#' the value defaults to the *nopastEvents* argument.
#'
#' The formula for recency for event \eqn{e_i} with type set to "inv.diff.plus1" and *i_neighborhood* is TRUE (Meijerink-Bosman et al. 2022):
#' \deqn{recency_{e_i} =\frac{1}{t_{e_i} - t^{*} + 1} }
#' where \eqn{t^{*}}, is the most recent time in
#' which the past event has the same receiver and sender as the current event. If there are no past events within the current dyad, then
#' the value defaults to the *nopastEvents* argument.
#'
#' The formula for recency for event \eqn{e_i} with type set to "inv.diff.plus1" and *i_neighborhood* is FALSE (Meijerink-Bosman et al. 2022):
#' \deqn{recency_{e_i} = \frac{1}{t_{e_i} - t^{*} + 1}         }
#' where \eqn{t^{*}}, is the most recent time in
#' which the past event's sender is the current event receiver and the past event receiver is the current event sender.  If there are no past events within the current dyad, then
#' the value defaults to the *nopastEvents* argument.
#'
#' The formula for recency for event \eqn{e_i} with type set to "rank.ordered.count" and *i_neighborhood* is TRUE (Butts 2008):
#' \deqn{recency_{e_i} = \rho(s(e_i), r(e_i), A_t)^{-1}}
#' where \eqn{\rho(s(e_i), r(e_i), A_t) }, is the current event receiver's rank amongst the current sender's recent relational events. That is, as Butts (2008: 174) argues,
#' "\eqn{\rho(s(e_i), r(e_i), A_t) } is j’s recency rank among i’s in-neighborhood. Thus, if j is the last person to have called i, then \eqn{\rho(s(e_i), r(e_i), A_t)^{-1}} = 1. This falls to 1/2 if j is the second
#' most recent person to call i, 1/3 if j is the third most recent person, etc." Moreover, if j is not in i's neighborhood, the value defaults to infinity. If there are no past events with the current sender, then
#' the value defaults to the *nopastEvents* argument.
#'
#' The formula for recency for event \eqn{e_i} with type set to "rank.ordered.count" and *i_neighborhood* is FALSE (Butts 2008):
#' \deqn{recency_{e_i} =  \rho(r(e_i), s(e_i), A_t)^{-1}}
#' where \eqn{\rho(r(e_i), s(e_i), A_t) }, is the current event sender's rank amongst the current receiver's recent relational events. That is, this measure is the same as above
#' where the dyadic pair is flipped for the past relational events. Moreover, if j is not in i's neighborhood, the value defaults to infinity. If there are no past events with the current sender, then
#' the value defaults to the *nopastEvents* argument.
#'
#' Finally, researchers interested in modeling temporal relevancy (see Quintane, Mood, Dunn, and Falzone 2022) can specify the relational time span, that is, length of time for which events are considered
#' relationally relevant. This should be specified via the option *relationalTimeSpan* with *dependency* set to TRUE.
#'
#'
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#' @references
#' Butts, Carter T. 2008. "A relational event framework for social action." *Sociological Methodology* 38(1): 155-200.
#'
#' Meijerink-Bosman, Marlyne, Roger Leenders, and Joris Mulder. 2022. "Dynamic relational event modeling: Testing, exploring,
#' and applying." *PLOS One* 17(8): e0272309.
#'
#' Quintane, Eric, Martin Wood, John Dunn, and Lucia Falzon. 2022. “Temporal
#' Brokering: A Measure of Brokerage as a Behavioral Process.” *Organizational Research Methods*
#' 25(3): 459-489.
#'
#' Vu, Duy, Philippa Pattison, and Garry Robbins. 2015. "Relational event models for social learning in MOOCs." *Social Networks* 43: 121-135.
#'
#'@examples
#'
#'
#'# A Dummy One-Mode Event Dataset
#'events <- data.frame(time = 1:18,
#'                                 eventID = 1:18,
#'                                 sender = c("A", "B", "C",
#'                                            "A", "D", "E",
#'                                            "F", "B", "A",
#'                                            "F", "D", "B",
#'                                            "G", "B", "D",
#'                                            "H", "A", "D"),
#'                                 target = c("B", "C", "D",
#'                                            "E", "A", "F",
#'                                            "D", "A", "C",
#'                                            "G", "B", "C",
#'                                            "H", "J", "A",
#'                                            "F", "C", "B"))
#'
#'# Creating the Post-Processing Event Dataset with Null Events
#'eventSet <- processOMEventSeq(data = events,
#'                           time = events$time,
#'                           eventID = events$eventID,
#'                           sender = events$sender,
#'                           receiver = events$target,
#'                           p_samplingobserved = 1.00,
#'                          n_controls = 6,
#'                          seed = 9999)
#'
#'# Compute Recency Statistic without Sliding Windows Framework and
#'# No Temporal Dependency
#'eventSet$recency_rawdiff <- computeRecency(
#'  observed_time = events$time,
#'  observed_receiver = events$target,
#'  observed_sender = events$sender,
#'  processed_time = eventSet$time,
#'  processed_receiver = eventSet$receiver,
#'  processed_sender = eventSet$sender,
#'  type = "raw.diff",
#'  dependency = FALSE,
#'  i_neighborhood = TRUE,
#'  nopastEvents = 0)
#'
#'# Compute Recency Statistic without Sliding Windows Framework and
#'# No Temporal Dependency
#'eventSet$recency_inv <- computeRecency(
#'  observed_time = events$time,
#'  observed_receiver = events$target,
#'  observed_sender = events$sender,
#'  processed_time = eventSet$time,
#'  processed_receiver = eventSet$receiver,
#'  processed_sender = eventSet$sender,
#'  type = "inv.diff.plus1",
#'  dependency = FALSE,
#'  i_neighborhood = TRUE,
#'  nopastEvents = 0)
#'
#'
#'# Compute Recency Statistic without Sliding Windows Framework and
#'# No Temporal Dependency
#'eventSet$recency_rank <- computeRecency(
#'  observed_time = events$time,
#'  observed_receiver = events$target,
#'  observed_sender = events$sender,
#'  processed_time = eventSet$time,
#'  processed_receiver = eventSet$receiver,
#'  processed_sender = eventSet$sender,
#'  type = "rank.ordered.count",
#'  dependency = FALSE,
#'  i_neighborhood = TRUE,
#'  nopastEvents = 0)
#'
#'# Compute Recency Statistic with Sliding Windows Framework and No Temporal Dependency
#'eventSet$recency_rawdiffSW <- computeRecency(
#'  observed_time = events$time,
#'  observed_receiver = events$target,
#'  observed_sender = events$sender,
#'  processed_time = eventSet$time,
#'  processed_receiver = eventSet$receiver,
#'  processed_sender = eventSet$sender,
#'  type = "raw.diff",
#'  dependency = FALSE,
#'  i_neighborhood = TRUE,
#'  sliding_windows = TRUE,
#'  processed_seqIDs = eventSet$sequenceID,
#'  nopastEvents = 0)
#'
#'
#'# Compute Recency Statistic with Sliding Windows Framework and No Temporal Dependency
#'eventSet$recency_invSW <- computeRecency(
#'  observed_time = events$time,
#'  observed_receiver = events$target,
#'  observed_sender = events$sender,
#'  processed_time = eventSet$time,
#'  processed_receiver = eventSet$receiver,
#'  processed_sender = eventSet$sender,
#'  type = "inv.diff.plus1",
#'  dependency = FALSE,
#'  i_neighborhood = TRUE,
#'  sliding_windows = TRUE,
#'  processed_seqIDs = eventSet$sequenceID,
#'  nopastEvents = 0)
#'
#'
#'# Compute Recency Statistic with Sliding Windows Framework and No Temporal Dependency
#'eventSet$recency_rankSW <- computeRecency(
#'  observed_time = events$time,
#'  observed_receiver = events$target,
#'  observed_sender = events$sender,
#'  processed_time = eventSet$time,
#'  processed_receiver = eventSet$receiver,
#'  processed_sender = eventSet$sender,
#'  type = "rank.ordered.count",
#'  dependency = FALSE,
#'  i_neighborhood = TRUE,
#'  sliding_windows = TRUE,
#'  processed_seqIDs = eventSet$sequenceID,
#'  nopastEvents = 0)
#'


computeRecency <-   function(    observed_time, # variable (column) name that contains the time variable
                                   observed_sender, # variable (column) name that contains the sender variable
                                   observed_receiver, # variable (column) name that contains the receiver variable
                                   processed_time,
                                   processed_sender, # variable (column) name that contains the sender variable
                                   processed_receiver, # variable (column) name that contains the receiver variable
                                   type = c("raw.diff", "inv.diff.plus1", "rank.ordered.count"),
                                   i_neighborhood = TRUE, #should the recency be computed on the i's neighborhood or j's neighborhood
                                   dependency = FALSE, #Boolean for temporal dependency
                                   relationalTimeSpan = NULL, #if dependency == TRUE, this should specific the associated temporal time span
                                   nopastEvents = NA, #the value given to events where there are no past events (defaults to NA)
                                   sliding_windows = FALSE, # TRUE = we want to use the sliding windows framework
                                   processed_seqIDs = NULL, #If true, the user should insert the placement of sampled events in the original file
                                   window_size = NA # the sizes of the windows that we will use, if NA, we will compute it internally
) {

  #base::cat("Checking Data Structure and User Inputs.......") # outputting status to user
  lifecycle::deprecate_warn("1.0.0", " computeRecency()", "remstats_recency()")
  ########################################################
  #
  #   Checking for Errors in User Inputs
  #
  ########################################################


  n_forrealevents <- base::length(observed_time) #the number of real events provided by user
  n_forsampledevents <- base::length(processed_time)#the number of sampled events provided by user

  if (length(observed_time) != n_forrealevents) { # if the length of the time vector is not the same as the full dataset
    base::stop("Error: The length of the provided events time is not the same length as the events dataset") # stop computation and tell the user
  }
  if (length(observed_sender) != n_forrealevents) { # if the length of the time vector is not the same as the full dataset
    base::stop("Error: The length of the provided events sender is not the same length as the events dataset") # stop computation and tell the user
  }
  if (length(observed_receiver) != n_forrealevents) { # if the length of the time vector is not the same as the full dataset
    base::stop("Error: The length of the provided events receiver is not the same length as the events dataset") # stop computation and tell the user
  }
  if (length(processed_time) != n_forsampledevents) { # if the length of the time vector is not the same as the full dataset
    base::stop("Error: The length of the provided eventSet time is not the same length as the sampld events dataset") # stop computation and tell the user
  }
  if (length(processed_sender) != n_forsampledevents) { # if the length of the time vector is not the same as the full dataset
    base::stop("Error: The length of the provided eventSet sender is not the same length as the sampld events dataset") # stop computation and tell the user
  }
  if (length(processed_receiver) != n_forsampledevents) { # if the length of the time vector is not the same as the full dataset
    base::stop("Error: The length of the provided eventSet receiver is not the same length as the sampld events dataset") # stop computation and tell the user
  }
  if(sliding_windows == TRUE & is.null(processed_seqIDs)){
    base::stop("Error: Sliding windows was specified to be true, however, the processed_seqIDs argument is missing. Please add this and restart the function!") # stop computation and tell the user
  }
  if(sliding_windows == TRUE & length(processed_seqIDs) != n_forsampledevents){
    base::stop("Error: The processed_seqIDs argument is not the same length as the processed_sender vector.") # stop computation and tell the user
  }
  if(((type %in% c("raw.diff", "inv.diff.plus1", "rank.ordered.count")) == FALSE) |
     length(type) != 1) {
    # stop computation and tell the user
    base::stop("Error: The argument for type was misspecified. The argument must be of type: raw.diff, inv.diff.plus1, or rank.ordered.count. Also, only one type can be specified at a time.")
  }
  if(dependency == TRUE & is.null(relationalTimeSpan)){
    base::stop("Error: Temporal dependency was requested, however, the relationalTimeSpan is missing. Please add this and restart the function.") # stop computation and tell the user

  }
  ########################################################
  # Renaming the columns to match the user inputs
  ########################################################
  events <- data.table::data.table(sender = observed_sender,# renaming the sender column
                                   receiver = observed_receiver, # renaming the receiver column
                                   time = observed_time)# renaming the time column

  eventSet <- data.table::data.table(sender = processed_sender,# renaming the sender column
                                     receiver = processed_receiver, # renaming the receiver column
                                     time = processed_time)# renaming the time column

  ########################################################
  #### Clearing User Inputs for Memory
  ########################################################
  rm(list = c("observed_time", "observed_sender", "observed_receiver",
              "processed_time", "processed_sender", "processed_receiver"))

  ########################################################
  #
  #    If the User Did not want to use the sliding windows framework
  #
  ########################################################

  if (sliding_windows == FALSE) {
    #base::cat("Starting computation of repetition scores without the sliding windows framework.......") # outputting status to user
    n_observed_events <- nrow(eventSet) # number of observed events in dataset
    recency <- numeric(n_observed_events) # empty vector to store computed statistics
    for (i in 1:n_observed_events) { # for all observed events in the dataset
      senderi <- eventSet$sender[i] ####  current user for event i
      receiveri <- eventSet$receiver[i] ####  current article for event i
      timei <- eventSet$time[i] ####  current time for event i

      ########################################################
      #
      #    Checking the neighborhood type
      #
      ########################################################
      ########################################################
      # all ties sent from the current sender (checking the senders outdegree)
      ########################################################
      if(i_neighborhood == TRUE){


        ####### For repetition!
        fullhistoryi <- events[ events$sender == senderi]
        if(base::nrow(fullhistoryi) == 0){recency[i] <- nopastEvents; next}
        fullhistoryi  <- fullhistoryi[fullhistoryi$time < timei]
        if(base::nrow(fullhistoryi) == 0){recency[i] <- nopastEvents; next}

        if(base::nrow(fullhistoryi[fullhistoryi$receiver == receiveri]) == 0){
          recency[i] <- Inf; next} #if j is not in i's neighborhood, value becomes 0 automatically

        if(dependency == TRUE){ #if the user wants temporal dependency to be modeled
          timeSpani <- timei - relationalTimeSpan #creating the temporal time span (the minimum relevant time)
          fullhistoryi <- fullhistoryi[fullhistoryi$time >= timeSpani]#getting only relevant events
          if(base::nrow(fullhistoryi) == 0){recency[i] <- nopastEvents; next} #if there are no remaining events
          if(base::nrow(fullhistoryi[fullhistoryi$receiver == receiveri]) == 0){
            recency[i] <- Inf; next} #if j is not in i's neighborhood, value becomes 0 automatically
        }
        ##### raw time difference, (inverse of time difference), inverse + 1, or rank difference (amongst recency of events)
        if(type == "raw.diff"){ #https://www-sciencedirect-com.ezproxy1.library.arizona.edu/science/article/pii/S0378873315000477
          recency[i]<- timei - base::max(fullhistoryi[fullhistoryi$receiver == receiveri]$time)
        }
        if(type == "inv.diff.plus1"){#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0272309
          recency[i]<- 1/((timei - base::max(fullhistoryi[fullhistoryi$receiver == receiveri]$time)) + 1)
        }
        if(type == "rank.ordered.count"){#Butts 2008
          reodered <- fullhistoryi[.N:1] #reflipping the dataset
          reodered <- unique(reodered, by = c("receiver")) #uniquing the data.table so that we can find the rank of the response
          #finding the specific rank of the individual
          recency[i] <- 1/collapse::whichv(reodered$receiver,receiveri)#the inverse of the rank!
        }


      }
      ########################################################
      # all ties recieved by the current sender (checking the senders indegree)
      ########################################################
      if(i_neighborhood == FALSE){
        ####### For repetition!
        fullhistoryi <- events[ events$receiver == senderi]
        if(base::nrow(fullhistoryi) == 0){recency[i] <- nopastEvents; next}
        fullhistoryi  <- fullhistoryi[fullhistoryi$time < timei]
        if(base::nrow(fullhistoryi) == 0){recency[i] <- nopastEvents; next}

        if(base::nrow(fullhistoryi[fullhistoryi$sender == receiveri]) == 0){
          recency[i] <- Inf; next} #if j is not in i's neighborhood, value becomes 0 automatically

        if(dependency == TRUE){ #if the user wants temporal dependency to be modeled
          timeSpani <- timei - relationalTimeSpan #creating the temporal time span (the minimum relevant time)
          fullhistoryi <- fullhistoryi[fullhistoryi$time >= timeSpani]#getting only relevant events
          if(base::nrow(fullhistoryi) == 0){recency[i] <- nopastEvents; next} #if there are no remaining events
          if(base::nrow(fullhistoryi[fullhistoryi$sender == receiveri]) == 0){
            recency[i] <- Inf; next} #if j is not in i's neighborhood, value becomes 0 automatically
        }
        ##### raw time difference, (inverse of time difference), inverse + 1, or rank difference (amongst recency of events)
        if(type == "raw.diff"){ #https://www-sciencedirect-com.ezproxy1.library.arizona.edu/science/article/pii/S0378873315000477
          recency[i]<- timei - base::max(fullhistoryi[fullhistoryi$sender == receiveri]$time)
        }
        if(type == "inv.diff.plus1"){#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0272309
          recency[i]<- 1/((timei - base::max(fullhistoryi[fullhistoryi$sender == receiveri]$time)) + 1)
        }
        if(type == "rank.ordered.count"){#Butts 2008
          reodered <- fullhistoryi[.N:1] #reflipping the dataset
          reodered <- unique(reodered, by = c("sender")) #uniquing the data.table so that we can find the rank of the response
          #finding the specific rank of the individual
          recency[i] <- 1/collapse::whichv(reodered$sender,receiveri)#the inverse of the rank!
        }

      }
      rm(list = c("fullhistoryi", "senderi","receiveri","timei" ))
    }
  }

  ########################################################
  #
  #    If the User Did Want To Use the Sliding Windows Framework
  #
  ########################################################
  if (sliding_windows == TRUE) {
    n_observed_events <- nrow(eventSet) # number of observed events in dataset
    recency <- rep(0, n_observed_events) # empty vector to store computed statistics
    #base::message("Setting up data structure for the sliding windows framework.......") # outputting status to user
    if (is.na(window_size)) { # if the user did not specify the size of the windows
      window_size <- round(nrow(eventSet) / 10) # an estimate of the window size of the dataset is simply dividing it into 10 datasets
    }
    ##### Given how R handles odd and even sequences, we have to internally check if the number of events is even odd
    starting_blocks <- seq(from = 1, to = nrow(eventSet), by = window_size) # creating a sequence from 1 to n by frames of window size
    ending_blocks <- seq(from = window_size, to = nrow(eventSet), by = window_size) # creating a sequence from windowsize to n by frames of window size
    if (length(starting_blocks) != length(ending_blocks)) { # if the vectors are not of the same length
      ending_blocks <- c(ending_blocks, nrow(eventSet)) # append the vector so they are the same length
    }
    # creating a sliding window framework
    sliding_window <- data.table::data.table(
      start_block = starting_blocks,
      stop_block = ending_blocks
    )
    eventSet$EVENT <- processed_seqIDs # the sampled eventnet IDs
    sliding_window$eventENDS <- eventSet[ending_blocks]$EVENT # getting the ending event from the block [the index in which the block ends]
    sliding_window$eventSTARTS <- rep(1, length(starting_blocks)) # creating a dummy variable for the starting blocks for each events (in other words, the starting of relationally relevancy)
    sliding_window$time <- eventSet[starting_blocks]$time # the eventtimes of the starting blocks (of the eventtime of the first event in each block)
    if(dependency == TRUE){
    sliding_window$minimum_times <- sliding_window$time - relationalTimeSpan #the relevant relational time span
    for (c in 2:nrow(sliding_window)) { # from 2 to all other windows
      sliding_window$eventSTARTS[c] <- base::min(which(events$time >= sliding_window$minimum_times[c])) # getting the minimum event that provides relevancy for each block
    }
    }else{
      sliding_window$minimum_times <- 0
      sliding_window$eventSTARTS <- 1
    }

    ######################################################################
    ####### Creating the Vector to Store Event Weights
    ######################################################################
    n_blocks <- nrow(sliding_window) # number of observed events in dataset

    for (j in 1:n_blocks) { # for all blocks from the sliding windows framework
      eventsj <- events[sliding_window$eventSTARTS[j]:sliding_window$eventENDS[j]] # get the current events information for the jth block
      new_start <- sliding_window$start_block[j] # the starting sampled event for block j
      new_stop <- sliding_window$stop_block[j] # the ending sampled event for block j

      for (i in new_start:new_stop) { # for all sampled events in block j

        senderi <- eventSet$sender[i] ####  current user for event i
        receiveri <- eventSet$receiver[i] ####  current article for event i
        timei <- eventSet$time[i] ####  current time for event i

        ########################################################
        #
        #    Checking the neighborhood type
        #
        ########################################################
        ########################################################
        # all ties sent from the current sender (checking the senders outdegree)
        ########################################################
        if(i_neighborhood == TRUE){

          ####### For repetition!
          fullhistoryi <- eventsj[ eventsj$sender == senderi]
          if(base::nrow(fullhistoryi) == 0){recency[i] <- nopastEvents; next}
          fullhistoryi  <- fullhistoryi[fullhistoryi$time < timei]
          if(base::nrow(fullhistoryi) == 0){recency[i] <- nopastEvents; next}

          if(base::nrow(fullhistoryi[fullhistoryi$receiver == receiveri]) == 0){
            recency[i] <- Inf; next} #if j is not in i's neighborhood, value becomes 0 automatically

          if(dependency == TRUE){ #if the user wants temporal dependency to be modeled
            timeSpani <- timei - relationalTimeSpan #creating the temporal time span (the minimum relevant time)
            fullhistoryi <- fullhistoryi[fullhistoryi$time >= timeSpani]#getting only relevant events
            if(base::nrow(fullhistoryi) == 0){recency[i] <- nopastEvents; next} #if there are no remaining events
            if(base::nrow(fullhistoryi[fullhistoryi$receiver == receiveri]) == 0){
              recency[i] <- Inf; next} #if j is not in i's neighborhood, value becomes 0 automatically
          }
          ##### raw time difference, (inverse of time difference), inverse + 1, or rank difference (amongst recency of events)
          if(type == "raw.diff"){ #https://www-sciencedirect-com.ezproxy1.library.arizona.edu/science/article/pii/S0378873315000477
            recency[i]<- timei - base::max(fullhistoryi[fullhistoryi$receiver == receiveri]$time)
          }
          if(type == "inv.diff.plus1"){#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0272309
            recency[i]<- 1/((timei - base::max(fullhistoryi[fullhistoryi$receiver == receiveri]$time)) + 1)
          }
          if(type == "rank.ordered.count"){#Butts 2008
            reodered <- fullhistoryi[.N:1] #reflipping the dataset
            reodered <- unique(reodered, by = c("receiver")) #uniquing the data.table so that we can find the rank of the response
            #finding the specific rank of the individual
            recency[i] <- 1/collapse::whichv(reodered$receiver,receiveri)#the inverse of the rank!
          }


        }
        ########################################################
        # all ties recieved by the current sender (checking the senders indegree)
        ########################################################
        if(i_neighborhood == FALSE){
          ####### For repetition!
          fullhistoryi <- eventsj[ eventsj$receiver == senderi]
          if(base::nrow(fullhistoryi) == 0){recency[i] <- nopastEvents; next}
          fullhistoryi  <- fullhistoryi[fullhistoryi$time < timei]
          if(base::nrow(fullhistoryi) == 0){recency[i] <- nopastEvents; next}

          if(base::nrow(fullhistoryi[fullhistoryi$sender == receiveri]) == 0){
            recency[i] <- Inf; next} #if j is not in i's neighborhood, value becomes 0 automatically

          if(dependency == TRUE){ #if the user wants temporal dependency to be modeled
            timeSpani <- timei - relationalTimeSpan #creating the temporal time span (the minimum relevant time)
            fullhistoryi <- fullhistoryi[fullhistoryi$time >= timeSpani]#getting only relevant events
            if(base::nrow(fullhistoryi) == 0){recency[i] <- nopastEvents; next} #if there are no remaining events
            if(base::nrow(fullhistoryi[fullhistoryi$sender == receiveri]) == 0){
              recency[i] <- Inf; next} #if j is not in i's neighborhood, value becomes 0 automatically
          }
          ##### raw time difference, (inverse of time difference), inverse + 1, or rank difference (amongst recency of events)
          if(type == "raw.diff"){ #https://www-sciencedirect-com.ezproxy1.library.arizona.edu/science/article/pii/S0378873315000477
            recency[i]<- timei - base::max(fullhistoryi[fullhistoryi$sender == receiveri]$time)
          }
          if(type == "inv.diff.plus1"){#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0272309
            recency[i]<- 1/((timei - base::max(fullhistoryi[fullhistoryi$sender == receiveri]$time)) + 1)
          }
          if(type == "rank.ordered.count"){#Butts 2008
            reodered <- fullhistoryi[.N:1] #reflipping the dataset
            reodered <- unique(reodered, by = c("sender")) #uniquing the data.table so that we can find the rank of the response
            #finding the specific rank of the individual
            recency[i] <- 1/collapse::whichv(reodered$sender,receiveri)#the inverse of the rank!
          }

        }
        rm(list = c("fullhistoryi", "senderi","receiveri","timei" ))    }
    }
  }
  ########################################################
  #
  #   Returning the values back to the user
  #
  ########################################################
  return(recency)# return the vector of values
}

