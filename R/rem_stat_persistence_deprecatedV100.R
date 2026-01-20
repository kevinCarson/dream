## The Creation of Repetition Scores for Relational Event Models (Can be Used in Both one-mode and two-mode networks)
## See Butts 2008; Lerner and Lomi 2020; REM R Package (Brandenberger 2018)
## Code written by Kevin Carson (kacarson@arizona.edu) and Deigo Leal (https://www.diegoleal.info/)
## Last Updated: 07-15-24

#' @title Compute Butts' (2008) Persistence Network Statistic for Event Dyads in a Relational Event Sequence
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `computePersistence()` has been deprecated starting on version 1.0.0 of the `dream` package. Please use the `remstats_persistence()` function and see the `NEWS.md` file for more details.
#'
#'
#' This function computes the persistence network sufficient statistic for
#' a relational event sequence (see Butts 2008). Persistence measures the proportion of past ties sent from the event sender that went to the current event receiver.
#' Furthermore, this measure allows for persistence scores to be only
#' computed for the sampled events, while creating the weights based on the full event
#' sequence. Moreover, the function allows users to specify relational relevancy for the statistic and
#' employ a sliding windows framework for large relational sequences.

#' @name computePersistence
#' @param observed_time The vector of event times from the pre-processing event sequence.
#' @param observed_sender The vector of event senders from the pre-processing event sequence.
#' @param observed_receiver The vector of event receivers from the pre-processing event sequence
#' @param processed_time The vector of event times from the post-processing event sequence (i.e., the event sequence that contains the observed and null events).
#' @param processed_sender The vector of event senders from the post-processing event sequence (i.e., the event sequence that contains the observed and null events).
#' @param processed_receiver The vector of event receivers from the post-processing event sequence (i.e., the event sequence that contains the observed and null events).
#' @param sender TRUE/FALSE. TRUE indicates that the persistence statistic will be computed in reference to the sender’s past relational history (see details section). FALSE indicates that the persistence statistic will be computed in reference to the target’s past relational history (see details section). Set to TRUE by default.
#' @param nopastEvents The numerical value that specifies what value should be given to events in which the sender has sent not past ties (i's neighborhood when sender = TRUE) or has not received any past ties (j's neighborhood when sender = FALSE). Set to NA by default.
#' @param dependency TRUE/FALSE. TRUE indicates that temporal relevancy will be modeled (see the details section). FALSE indicates that temporal relevancy will not be modeled, that is, all past events are relevant (see the details section). Set to FALSE by default.
#' @param relationalTimeSpan If dependency = TRUE, a numerical value that corresponds to the temporal span for relational relevancy, which must be the same measurement unit as the observed_time and processed_time objects. When dependency = TRUE, the relevant events are events that have occurred between current event time, *t*, and *t-relationalTimeSpan*. For example, if the time measurement is the number of days since the first event and the value for relationalTimeSpan is set to 10, then only those events which occurred in the past 10 days are included in the computation of the statistic.
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
#' @return The vector of persistence network statistics for the relational event sequence.
#' @export
#'
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#' @references
#' Butts, Carter T. 2008. "A relational event framework for social action." *Sociological Methodology* 38(1): 155-200.
#'
#' Quintane, Eric, Martin Wood, John Dunn, and Lucia Falzon. 2022. “Temporal
#' Brokering: A Measure of Brokerage as a Behavioral Process.” *Organizational Research Methods*
#' 25(3): 459-489.
#'
#'@details The function calculates the persistence network sufficient statistic for a relational event sequence based on Butts (2008).
#'
#'The formula for persistence for event \eqn{e_i} with reference to the sender's past relational history is:
#'\deqn{Persistence_{e_{i}} = \frac{d(s(e_{i}),r(e_{i}), A_t)}{d(s(e_{i}), A_t)} }
#'
#'where  \eqn{d(s(e_{i}),r(e_{i}), A_t)} is the number of past events where the current event sender sent a tie to the current event receiver, and \eqn{d(s(e_{i}), A_t)} is the number of past events where the current sender sent a tie.
#'
#'The formula for persistence for event \eqn{e_i} with reference to the target's past relational history is:
#'\deqn{Persistence_{e_{i}} = \frac{d(s(e_{i}),r(e_{i}), A_t)}{d(r(e_{i}), A_t)} }
#'
#'where  \eqn{d(s(e_{i}),r(e_{i}), A_t)} is the number of past events where the current event sender sent a tie to the current event receiver, and \eqn{d(r(e_{i}), A_t)} is the number of past events where the current receiver recieved a tie.
#'
#'Moreover, researchers interested in modeling temporal relevancy (see Quintane, Mood, Dunn, and Falzone 2022) can specify the relational time span, that is, length of time for which events are considered
#'relationally relevant. This should be specified via the option *relationalTimeSpan* with *dependency* set to TRUE.
#'
#'
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
#'                           n_controls = 6,
#'                           seed = 9999)
#'
#'#Compute Persistence with respect to the sender's past relational history without
#'#the sliding windows framework and no temporal dependency
#'eventSet$persist <- computePersistence(observed_time = events$time,
#'                                         observed_receiver = events$target,
#'                                         observed_sender = events$sender,
#'                                         processed_time = eventSet$time,
#'                                         processed_receiver = eventSet$receiver,
#'                                         processed_sender = eventSet$sender,
#'                                         sender = TRUE,
#'                                         nopastEvents = 0)
#'
#'#Compute Persistence with respect to the sender's past relational history with
#'#the sliding windows framework and no temporal dependency
#'eventSet$persistSW <- computePersistence(observed_time = events$time,
#'                                         observed_receiver = events$target,
#'                                         observed_sender = events$sender,
#'                                         processed_time = eventSet$time,
#'                                         processed_receiver = eventSet$receiver,
#'                                         processed_sender = eventSet$sender,
#'                                         sender = TRUE,
#'                                         sliding_windows = TRUE,
#'                                         processed_seqIDs = eventSet$sequenceID,
#'                                         nopastEvents = 0)
#'
#'#The results with and without the sliding windows are the same (see correlation
#'#below). Using the sliding windows method is recommended when the data are
#'#big' so that memory allotment is more efficient.
#'cor(eventSet$persist,eventSet$persistSW)
#'
#'
#'#Compute Persistence with respect to the sender's past relational history without
#'#the sliding windows framework and temporal dependency
#'eventSet$persistDep <- computePersistence(observed_time = events$time,
#'                                         observed_receiver = events$target,
#'                                         observed_sender = events$sender,
#'                                         processed_time = eventSet$time,
#'                                         processed_receiver = eventSet$receiver,
#'                                         processed_sender = eventSet$sender,
#'                                         sender = TRUE,
#'                                         dependency = TRUE,
#'                                         relationalTimeSpan = 5, #the past 5 events
#'                                         nopastEvents = 0)
#'
#'#Compute Persistence with respect to the receiver's past relational history without
#'#the sliding windows framework and no temporal dependency
#'eventSet$persistT <- computePersistence(observed_time = events$time,
#'                                         observed_receiver = events$target,
#'                                         observed_sender = events$sender,
#'                                         processed_time = eventSet$time,
#'                                         processed_receiver = eventSet$receiver,
#'                                         processed_sender = eventSet$sender,
#'                                         sender = FALSE,
#'                                         nopastEvents = 0)
#'
#'#Compute Persistence with respect to the receiver's past relational history with
#'#the sliding windows framework and no temporal dependency
#'eventSet$persistSWT <- computePersistence(observed_time = events$time,
#'                                         observed_receiver = events$target,
#'                                         observed_sender = events$sender,
#'                                         processed_time = eventSet$time,
#'                                         processed_receiver = eventSet$receiver,
#'                                         processed_sender = eventSet$sender,
#'                                         sender = FALSE,
#'                                         sliding_windows = TRUE,
#'                                         processed_seqIDs = eventSet$sequenceID,
#'                                         nopastEvents = 0)
#'
#'#The results with and without the sliding windows are the same (see correlation
#'#below). Using the sliding windows method is recommended when the data are
#'#big' so that memory allotment is more efficient.
#'cor(eventSet$persistT,eventSet$persistSWT)
#'
#'
#'#Compute Persistence with respect to the receiver's past relational history without
#'#the sliding windows framework and temporal dependency
#'eventSet$persistDepT <- computePersistence(observed_time = events$time,
#'                                         observed_receiver = events$target,
#'                                         observed_sender = events$sender,
#'                                         processed_time = eventSet$time,
#'                                         processed_receiver = eventSet$receiver,
#'                                         processed_sender = eventSet$sender,
#'                                         sender = FALSE,
#'                                         dependency = TRUE,
#'                                         relationalTimeSpan = 5, #the past 5 events
#'                                         nopastEvents = 0)
#'


computePersistence <-   function(observed_time, # variable (column) name that contains the time variable
                                   observed_sender, # variable (column) name that contains the sender variable
                                   observed_receiver, # variable (column) name that contains the receiver variable
                                   processed_time,
                                   processed_sender, # variable (column) name that contains the sender variable
                                   processed_receiver, # variable (column) name that contains the receiver variable
                                   sender = TRUE, #if persistence should be computed with reference for the sender or the receiver
                                   dependency = FALSE, #Boolean for temporal dependency
                                   relationalTimeSpan = NULL, #if dependency == TRUE, this should specific the associated temporal time span
                                   nopastEvents = NA, #the value given to events where there are no past events (defaults to NA)
                                   sliding_windows = FALSE, # TRUE = we want to use the sliding windows framework
                                   processed_seqIDs = NULL, #If true, the user should insert the placement of sampled events in the original file
                                   window_size = NA # the sizes of the windows that we will use, if NA, we will compute it internally
) {

  #base::cat("Checking Data Structure and User Inputs.......") # outputting status to user
  lifecycle::deprecate_warn("1.0.0", " computePersistence()", "remstats_persistence()")
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
    persistence <- numeric(n_observed_events) # empty vector to store computed statistics
    for (i in 1:n_observed_events) { # for all observed events in the dataset
      senderi <- eventSet$sender[i] ####  current user for event i
      receiveri <- eventSet$receiver[i] ####  current article for event i
      timei <- eventSet$time[i] ####  current time for event i
      ####### For persistence!


      ############################################################################
      #
      #   Persistence with respect to the sender: (shared) / outdegree
      #
      ############################################################################
      if(sender == TRUE){
      fullhistoryi <- events[ events$sender == senderi]
      if(base::nrow(fullhistoryi) == 0){persistence[i] <- nopastEvents; next}
      fullhistoryi  <- fullhistoryi [fullhistoryi$time < timei]
      if(base::nrow(fullhistoryi) == 0){persistence[i] <- nopastEvents; next}
      if(dependency == TRUE){ #if the user wants temporal dependency to be modeled
      timeSpani <- timei - relationalTimeSpan #creating the temporal time span (the minimum relevant time)
      fullhistoryi <- fullhistoryi[fullhistoryi$time >= timeSpani]#getting only relevant events
      if(base::nrow(fullhistoryi) == 0){persistence[i] <- nopastEvents; next}
      }
      fullN <- base::nrow(fullhistoryi) #the total number of past sent ties
      sharedN <- base::nrow(fullhistoryi[fullhistoryi$receiver == receiveri]) #n ties sent to current reciever
      persistence[i] <- sharedN/fullN #add the shared proportion
      }
      ############################################################################
      #
      #   Persistence with respect to the target: (shared) / indegree
      #
      ############################################################################
      if(sender == FALSE){
        fullhistoryi <- events[ events$receiver == receiveri]
        if(base::nrow(fullhistoryi) == 0){persistence[i] <- nopastEvents; next}
        fullhistoryi  <- fullhistoryi[fullhistoryi$time < timei]
        if(base::nrow(fullhistoryi) == 0){persistence[i] <- nopastEvents; next}
        if(dependency == TRUE){ #if the user wants temporal dependency to be modeled
          timeSpani <- timei - relationalTimeSpan #creating the temporal time span (the minimum relevant time)
          fullhistoryi <- fullhistoryi[fullhistoryi$time >= timeSpani]#getting only relevant events
          if(base::nrow(fullhistoryi) == 0){persistence[i] <- nopastEvents; next}
        }
        fullN <- base::nrow(fullhistoryi) #the total number of past received ties
        sharedN <- base::nrow(fullhistoryi[fullhistoryi$sender == senderi]) #n ties sent to current reciever
        persistence[i] <- sharedN/fullN #add the shared proportion
      }

      rm(list = c("fullhistoryi", "fullN", "sharedN", "senderi",
                 "receiveri", "timei" ))
    }
  }

  ########################################################
  #
  #    If the User Did Want To Use the Sliding Windows Framework
  #
  ########################################################
  if (sliding_windows == TRUE) {
    n_observed_events <- nrow(eventSet) # number of observed events in dataset
    persistence <- rep(0, n_observed_events) # empty vector to store computed statistics
    base::message("Setting up data structure for the sliding windows framework.......") # outputting status to user
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
        ############################################################################
        #
        #   Persistence Attachment with respect to the sender: (shared) / outdegree
        #
        ############################################################################
        if(sender == TRUE){
          fullhistoryi <- eventsj[eventsj$sender == senderi]
          if(base::nrow(fullhistoryi) == 0){persistence[i] <- nopastEvents; next}
          fullhistoryi  <- fullhistoryi[fullhistoryi$time < timei]
          if(base::nrow(fullhistoryi) == 0){persistence[i] <- nopastEvents; next}
          if(dependency == TRUE){ #if the user wants temporal dependency to be modeled
            timeSpani <- timei - relationalTimeSpan #creating the temporal time span (the minimum relevant time)
            fullhistoryi <- fullhistoryi[fullhistoryi$time >= timeSpani]#getting only relevant events
            if(base::nrow(fullhistoryi) == 0){persistence[i] <- nopastEvents; next}
          }
          fullN <- base::nrow(fullhistoryi) #the total number of past sent ties
          sharedN <- base::nrow(fullhistoryi[fullhistoryi$receiver == receiveri]) #n ties sent to current reciever
          persistence[i] <- sharedN/fullN #add the shared proportion
        }
        ############################################################################
        #
        #   Persistence with respect to the target: (shared) / indegree
        #
        ############################################################################
        if(sender == FALSE){
          fullhistoryi <- eventsj[eventsj$receiver == receiveri]
          if(base::nrow(fullhistoryi) == 0){persistence[i] <- nopastEvents; next}
          fullhistoryi  <- fullhistoryi[fullhistoryi$time < timei]
          if(base::nrow(fullhistoryi) == 0){persistence[i] <- nopastEvents; next}
          if(dependency == TRUE){ #if the user wants temporal dependency to be modeled
            timeSpani <- timei - relationalTimeSpan #creating the temporal time span (the minimum relevant time)
            fullhistoryi <- fullhistoryi[fullhistoryi$time >= timeSpani]#getting only relevant events
            if(base::nrow(fullhistoryi) == 0){persistence[i] <- nopastEvents; next}
          }
          fullN <- base::nrow(fullhistoryi) #the total number of past received ties
          sharedN <- base::nrow(fullhistoryi[fullhistoryi$sender == senderi]) #n ties sent to current reciever
          persistence[i] <- sharedN/fullN #add the shared proportion
        }

        rm(list = c("fullhistoryi", "fullN", "sharedN", "senderi",
                    "receiveri", "timei" ))
        }
    }
  }
  ########################################################
  #
  #   Returning the values back to the user
  #
  ########################################################
  return(persistence)# return the vector of values
}
