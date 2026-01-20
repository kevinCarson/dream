## The Creation of Indegree Statistic Scores for Relational Event Models
## See Butts 2008; Lerner and Lomi 2020; REM R Package (Brandenberger 2018)
## Code written by Kevin Carson (kacarson@arizona.edu) and Deigo Leal (https://www.diegoleal.info/)
## Last Updated: 07-15-24

#' @title Compute the Indegree Network Statistic for Event Receivers in a Relational Event Sequence
#' @name computeReceiverIndegree
#' @param observed_time The vector of event times from the pre-processing event sequence.
#' @param observed_receiver The vector of event receivers from the pre-processing event sequence
#' @param processed_time The vector of event times from the post-processing event sequence (i.e., the event sequence that contains the observed and null events).
#' @param processed_receiver The vector of event receivers from the post-processing event sequence (i.e., the event sequence that contains the observed and null events).
#' @param counts TRUE/FALSE. TRUE indicates that the counts of past events should be computed (see the details section). FALSE indicates that the temporal exponential weighting function should be used to downweigh past events (see the details section). Set to FALSE by default.
#' @param halflife A numerical value that is the halflife value to be used in the exponential weighting function (see details section). Preset to 2 (should be updated by the user based on substantive context).
#' @param dyadic_weight A numerical value for the dyadic cutoff weight that represents the numerical cutoff value for temporal relevancy based on the exponential weighting function. For example, a numerical value of 0.01, indicates that an exponential weight less than 0.01 will become 0 and that events with such value (or smaller values) will not be included in the sum of the past event weights (see the details section). Set to 0 by default.
#' @param Lerneretal_2013 TRUE/FALSE. TRUE indicates that the Lerner et al. (2013) exponential weighting function will be used (see the details section). FALSE indicates that the Lerner and Lomi (2020) exponential weighting function will be used (see the details section). Set to FALSE by default
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
#' @return The vector of receiver indegree statistics for the relational event sequence.
#' @export

#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `computeReceiverIndegree()` has been deprecated starting on version 1.0.0 of the `dream` package. Please use the `remstats_degree()` function and see the `NEWS.md` file for more details.
#'
#'
#' The function computes the indegree network sufficient statistic for event receivers
#' in a relational event sequence (see Lerner and Lomi 2020; Butts 2008).
#' This measure allows for the indegree scores to be computed only for the sampled
#' events, while creating the weights based on the full event sequence (see
#' Lerner and Lomi 2020; Vu et al. 2015). The function allows users to use two different weighting functions,
#' reduce computational runtime, employ a sliding windows framework for large relational sequences, and
#' specify a dyadic cutoff for relational relevancy.
#'
#'
#'@details The function calculates receiver indegree scores for relational event
#'sequences based on the exponential weighting function used in either Lerner
#'and Lomi (2020) or Lerner et al. (2013).
#'
#'Following Lerner and Lomi (2020), the exponential weighting function in
#'relational event models is:
#'\deqn{w(s, r, t) = e^{-(t-t') \cdot \frac{ln(2)}{T_{1/2}} }}
#'
#'Following Lerner et al. (2013), the exponential weighting function in
#'relational event models is:
#'\deqn{w(s, r, t) = e^{-(t-t') \cdot \frac{ln(2)}{T_{1/2}} } \cdot \frac{ln(2)}{T_{1/2}}}
#'
#'In both of the above equations, *s* is the current event sender, *r* is the
#'current event receiver (target), *t* is the current event time, *t'* is the
#'past event times that meet the weight subset, and \eqn{T_{1/2}} is the halflife parameter.
#'
#'The formula for receiver indegree for event \eqn{e_i} is:
#'\deqn{reciever indegree_{e_{i}} = w(s', r, t) }
#'
#'That is, all past events in which the event receiver is the current receiver.
#'
#'Moreover, researchers interested in modeling temporal relevancy (see Quintane,
#'Mood, Dunn, and Falzone 2022; Lerner and Lomi 2020) can specify the dyadic
#'weight cutoff, that is, the minimum value for which the weight is considered
#'relationally relevant. Users who do not know the specific dyadic cutoff value to use, can use the
#'\code{\link{computeRemDyadCut}} function.
#'
#'Following Butts (2008), if the counts of the past events are requested, the formula for receiver indegree for
#'event \eqn{e_i} is:
#'\deqn{repetition_{e_{i}} = d(r' = r, t') }
#'where, \eqn{d()} is the number of past events where the past event receiver, *r'*, is the
#'current event receiver (target). Moreover, the counting equation
#'can be used in tandem with relational relevancy, by specifying the halflife parameter, exponential
#'weighting function, and the dyadic cut off weight values. If the user is not interested in modeling
#'relational relevancy, then those value should be left at their defaults.
#'

#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#' @references
#'
#' Butts, Carter T. 2008. "A Relational Event Framework for Social Action." *Sociological Methodology* 38(1): 155-200.
#'
#'Quintane, Eric, Martin Wood, John Dunn, and Lucia Falzon. 2022. “Temporal
#'Brokering: A Measure of Brokerage as a Behavioral Process.” *Organizational Research Methods*
#'25(3): 459-489.
#'
#'Lerner, Jürgen and Alessandro Lomi. 2020. “Reliability of relational event
#'model estimates under sampling: How to fit a relational event model to 360
#'million dyadic events.” *Network Science* 8(1): 97-135.
#'
#'Lerner, Jürgen, Margit Bussman, Tom A.B. Snijders, and Ulrik Brandes. 2013. "
#'Modeling Frequency and Type of Interaction in Event Networks."
#'*The Corvinus Journal of Sociology and Social Policy* 4(1): 3-32.
#'
#' Vu, Duy, Philippa Pattison, and Garry Robins. 2015. "Relational event models for social learning in MOOCs." *Social Networks* 43: 121-135.



#'@examples
#'events <- data.frame(time = 1:18,
#'                                 eventID = 1:18,
#'                                 sender = c("A", "B", "C",
#'                                            "A", "D", "E",
#'                                            "F", "B", "A",
#'                                            "F", "D", "B",
#'                                            "G", "B", "D",
#'                                           "H", "A", "D"),
#'                                target = c("B", "C", "D",
#'                                           "E", "A", "F",
#'                                           "D", "A", "C",
#'                                           "G", "B", "C",
#'                                           "H", "J", "A",
#'                                           "F", "C", "B"))
#'
#'eventSet <- processOMEventSeq(data = events,
#'                       time = events$time,
#'                       eventID = events$eventID,
#'                       sender = events$sender,
#'                       receiver = events$target,
#'                       p_samplingobserved = 1.00,
#'                       n_controls = 1,
#'                       seed = 9999)
#'
#'# Computing Target Indegree Statistics without the sliding windows framework
#'eventSet$target_indegree <- computeReceiverIndegree(
#'    observed_time = events$time,
#'    observed_receiver = events$target,
#'    processed_time = eventSet$time,
#'    processed_receiver = eventSet$receiver,
#'    halflife = 2, #halflife parameter
#'    dyadic_weight = 0,
#'    Lerneretal_2013 = FALSE)
#'
#'# Computing Target Indegree Statistics with the sliding windows framework
#'eventSet$target_indegreeSW <- computeReceiverIndegree(
#'    observed_time = events$time,
#'    observed_receiver = events$target,
#'    processed_time = eventSet$time,
#'    processed_receiver = eventSet$receiver,
#'    halflife = 2, #halflife parameter
#'    processed_seqIDs = eventSet$sequenceID,
#'    dyadic_weight = 0,
#'    sliding_window = TRUE,
#'    Lerneretal_2013 = FALSE)
#'
#'#The results with and without the sliding windows are the same (see correlation
#'#below). Using the sliding windows method is recommended when the data are
#'#big' so that memory allotment is more efficient.
#'cor(eventSet$target_indegree , eventSet$target_indegreeSW )
#'
#'# Computing Target Indegree Statistics with the counts of events being returned
#'eventSet$target_indegreeC <- computeReceiverIndegree(
#'    observed_time = events$time,
#'    observed_receiver = events$target,
#'    processed_time = eventSet$time,
#'    processed_receiver = eventSet$receiver,
#'    halflife = 2, #halflife parameter
#'    processed_seqIDs = eventSet$sequenceID,
#'    dyadic_weight = 0,
#'    sliding_window = TRUE,
#'    Lerneretal_2013 = FALSE,
#'    counts = TRUE)
#'
#'cbind(eventSet$target_indegree,
#'      eventSet$target_indegreeSW,
#'      eventSet$target_indegreeC)


########################################################################################################
#  Events = the full event sequence
#  eventSet = the sampled event sequence
#  Time = the name for the time variable
#  receiver = the name for the receiver variable
#  eventID = the name of the event sequence variable
#  sliding_windows = logical value, should the sliding windows framework be used (TRUE = yes; FALSE = no)
#  halflife = halflife parameter for exponential weighting function
#  dyadic_weight = numerical value for relational relevance cutoff for events (set to 0 for all events to count)
#  window_size = size of windows for the sliding windows framework (if NA, value will be internally computed)
#  Lerneretal_2013 = which version of the exponential weighting version should be used (see weighting function)
#  returnOnlyValues = FALSE indicates return the eventSet dataframe with the repetition values added, if TRUE only values are returned in a vector
#  countsofevents = A logical value that indicates if we want the raw number of events or the exponential weighting function used: (see Butts 2008: 195 d(i,j,Ak))
#  note: that the mef can still be used if we only want the counts
########################################################################################################

computeReceiverIndegree <- function(observed_time, # variable (column) name that contains the time variable
                               observed_receiver, # variable (column) name that contains the receiver variable
                               processed_time,
                               processed_receiver, # variable (column) name that contains the receiver variable
                               sliding_windows = FALSE, # TRUE = we want to use the sliding windows framework
                               processed_seqIDs = NULL, #If true, the user should insert the placement of sampled events in the original file
                               counts = FALSE, #Logical indicating if the raw counts of events should be returned or the exponential weighting function should be used (TRUE = counts; FALSE = exponential weighting)
                               halflife=2, # the half life value for the weighting function
                               dyadic_weight=0.00, # dyadic cutoff weight for events that no longer matter
                               window_size = NA, # the sizes of the windows that we will use, if NA, we will compute it internally
                               Lerneretal_2013 = FALSE
) {
  #base::cat("Checking Data Structure and User Inputs.......") # outputting status to user

  ########################################################
  #
  #   Checking for Errors in User Inputs
  #
  ########################################################
  lifecycle::deprecate_warn("1.0.0", " computeReceiverIndegree()", "remstats_degree()")
  n_forrealevents <- base::length(observed_time) #the number of real events provided by user
  n_forsampledevents <- base::length(processed_time)#the number of sampled events provided by user

  if (length(observed_time) != n_forrealevents) { # if the length of the time vector is not the same as the full dataset
    base::stop("Error: The length of the provided events time is not the same length as the events dataset") # stop computation and tell the user
  }
  if (length(observed_receiver) != n_forrealevents) { # if the length of the time vector is not the same as the full dataset
    base::stop("Error: The length of the provided events receiver is not the same length as the events dataset") # stop computation and tell the user
  }
  if (length(processed_time) != n_forsampledevents) { # if the length of the time vector is not the same as the full dataset
    base::stop("Error: The length of the provided eventSet time is not the same length as the sampld events dataset") # stop computation and tell the user
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
  ########################################################
  # Renaming the columns to match the user inputs
  ########################################################
  events <- data.table::data.table(receiver = observed_receiver, # renaming the receiver column
                                   time = observed_time)# renaming the time column

  eventSet <- data.table::data.table(receiver = processed_receiver, # renaming the receiver column
                                     time = processed_time)# renaming the time column

  ########################################################
  #### Clearing User Inputs for Memory
  ########################################################
  rm(list = c("observed_time", "observed_receiver", "processed_time", "processed_receiver"))

  ########################################################
  #
  #    If the User Did not want to use the sliding windows framework
  #
  ########################################################
  if (sliding_windows == FALSE) {
    #base::cat("Starting computation of indegree scores without the sliding windows framework.......") # outputting status to user
    n_observed_events <- nrow(eventSet) # number of observed events in dataset
    indegree <- rep(0, n_observed_events) # empty vector to store computed statistics
    for (i in 1:n_observed_events) { # for all observed events in the dataset
      receiveri <- eventSet$receiver[i] ####  current article for event i
      timei <- eventSet$time[i] ####  current time for event i
      ####### For repetition!
      indegreei <- events[events$receiver == receiveri] # all events that have the current article and user
      if(nrow(indegreei) == 0){next}
      ### for all events that occurred before the current, compute the weight with a dyadic cutoff at dyadic_weighti
      indegreei <- indegreei[indegreei$time < timei]
      if(nrow(indegreei) == 0){next}
        indegreei$weight <- remExpWeights(
          current = timei,
          past = indegreei$time, halflife = halflife,
          dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013
        )
        if(counts == TRUE){ #if the user wants the counts of events to be returned
          indegree[i] <- length(indegreei$weight)  # add the sum of the past event weights to the ith spot in the inertia vector
        }else{
          indegree[i] <- sum(indegreei$weight)  # add the sum of the past event weights to the ith spot in the inertia vector
        }
        rm(list = c("indegreei"))
    }
  }

  ########################################################
  #
  #    If the User Did Want To Use the Sliding Windows Framework
  #
  ########################################################
  if (sliding_windows == TRUE) {
    #base::cat("Setting up data structure for the sliding windows framework.......") # outputting status to user
    if (is.na(window_size)) {
      window_size <- round(nrow(events) / 10) # an estimate of the window size of the dataset is simply dividing it into 10 datasets
    }

    ##### Given how R handles odd and even sequences, we have to internally check if the number of events is even odd
    starting_blocks <- seq(from = 1, to = nrow(events), by = window_size)
    ending_blocks <- seq(from = window_size, to = nrow(events), by = window_size)
    if (length(starting_blocks) != length(ending_blocks)) {
      ending_blocks <- c(ending_blocks, nrow(events))
    }
    sliding_window <- data.table::data.table(
      start_block = starting_blocks,
      stop_block = ending_blocks
    )

    #### Computing the minimum effective event time for each sampled event (event time that exceeds relational relevancy)
    #### see minimum_effective_time()
    # eventSet[, minimum_event := redr:::minimum_effective_time(
    #   eventtime = time, dyadicweight = dyadic_weight,
    #   halflife = halflife, Lerneretal_2013 = Lerneretal_2013
    # )]
    eventSet$minimum_event <-  remExpWeights(
      current = eventSet$time, dyadic_weight = dyadic_weight,
      halflife = halflife, Lerneretal_2013 = Lerneretal_2013,
      exp.weights = FALSE)

    ###### Doing a quick test: it should be noted that too large of a halflife parameter in the Lerner et al. 2013
    ###### weighting function results in minimum effective time greater than the eventTime. Therefore, it cannot be
    ###### properly approximated in this case.
    test_vec <- eventSet$minimum_event - eventSet$time # getting the difference between the values, these should all be negative
    areanygreater <- sum((test_vec) > 0) # checking if we have any positive values,
    if (areanygreater != 0) { # if we do, stop the function, and tell the user this
      base::stop("Error: Unfortunately the combination of the provided halflife parameter and the weighting function (i.e., using the \n
                 Lerner et al. 2013 specification) resulted in a minimum effective time that is greater than the eventTimes, therefore, \n
                 the sliding windows framework cannot be used. Please restart with sliding_windows == FALSE. Please see the documentation \n
                 for this function and the minimum effective time documentation. As always, we hope you know what you're doing.....") # stop computation and tell the user
    }


    ######################################################################
    ####### Creating the Vector to Store Event Weights
    ######################################################################
    n_observed_events <- nrow(eventSet) # number of observed events in dataset
    indegree <- rep(0, n_observed_events) # empty vector to store computed statistics
    eventIDS <- processed_seqIDs # the sampled eventnet IDS

    ######################################################################
    ####### Computing Statistics on a Sliding Frame Method
    ######################################################################
    #base::cat("Starting computation of indegree scores with sliding windows framework.......") # outputting status to user
    n_blocks <- nrow(sliding_window) # the number of sliding windows
    for (j in 1:n_blocks) { # go through each of the blocks (i.e., subsetting the datamatrix to only handle
      # an arbitrary number of rows at a time)

      ######################################################################
      ####### Subsetting full event sequence based on sliding window sizes
      ######################################################################
      eventsj <- events[sliding_window$start_block[j]:sliding_window$stop_block[j]]


      ########################################################
      #
      #   Comment: Eventually allow the user to add the working directory and file name encase the WikiEvent dataset is being used
      #
      ########################################################


      ### for jumping ahead in a moving frame
      ### the sampled event has occured before the current blocks, then skip
      ### this is because event statistics are not computed for future events
      ### only the network of past events
      ifelse(j == 1, new_start <- 1,
        new_start <- base::min(which(eventIDS >= sliding_window$start_block[j]))
      )
      if(new_start == Inf){ #if there is no longer a need to continue moving through the dataset
        break #break from the sliding windows inner for loop
      }

      ifelse(j == n_blocks, n_sampled <- nrow(eventSet),
        n_sampled <- base::suppressWarnings({base::min(which(eventSet$minimum_event >= max(eventsj$time))) - 1}) # max possible event is the minimum minus 1, since minimum will have a weight of zero
      )

      ifelse(n_sampled == Inf, n_sampled <- nrow(eventSet), n_sampled)

      for (i in new_start:n_sampled) { # for i in all sampled events
        receiveri <- eventSet$receiver[i] ####  current article for event i
        timei <- eventSet$time[i] ####  current time for event i
        ####### For repetition!
        indegreei <- eventsj[eventsj$receiver == receiveri] # all events that have the current article and user
        if(nrow(indegreei) == 0){next}
        ### for all events that occurred before the current, compute the weight with a dyadic cutoff at dyadic_weighti
        indegreei <- indegreei[indegreei$time < timei]
        if(nrow(indegreei) == 0){next}
          indegreei$weight <- remExpWeights(
            current = timei,
            past = indegreei$time, halflife = halflife,
            dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013
          )
          if(counts == TRUE){ #if the user wants the counts of events to be returned
            indegree[i] <- length(indegreei$weight) + indegree[i]  # add the sum of the past event weights to the ith spot in the inertia vector
          }else{
            indegree[i] <- sum(indegreei$weight) + indegree[i] # add the sum of the past event weights to the ith spot in the inertia vector
          }
        ###### clearing memory for consumption purposes
        rm(list = c("receiveri", "timei", "indegreei"))
      }
      rm(list = c("eventsj"))
    }
  }
  ########################################################
  #
  #   Returning the values back to the user
  #
  ########################################################
  return(indegree) # return the vector of values
}
