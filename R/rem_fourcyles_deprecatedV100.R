## The Creation of Four-Cycle Statistics for Large Relational Events into R
## Purpose: Compute Four Cycle Statistics
## See Butts 2008; Lerner and Lomi 2020; REM R Package (Brandenberger 2018)
## Code written by Kevin Carson (kacarson@arizona.edu) and Deigo Leal (https://www.diegoleal.info/)
## Last Updated: 07-31-24
#' @title Compute the Four-Cycles Network Statistic for Event Dyads in a Relational Event Sequence
#' @name  computeFourCycles
#' @param observed_time The vector of event times from the pre-processing event sequence.
#' @param observed_sender The vector of event senders from the pre-processing event sequence.
#' @param observed_receiver The vector of event receivers from the pre-processing event sequence
#' @param processed_time The vector of event times from the post-processing event sequence (i.e., the event sequence that contains the observed and null events).
#' @param processed_sender The vector of event senders from the post-processing event sequence (i.e., the event sequence that contains the observed and null events).
#' @param processed_receiver The vector of event receivers from the post-processing event sequence (i.e., the event sequence that contains the observed and null events).
#' @param counts TRUE/FALSE. TRUE indicates that the counts of past events should be computed (see the details section). FALSE indicates that the temporal exponential weighting function should be used to downweigh past events (see the details section). Set to FALSE by default.
#' @param halflife A numerical value that is the halflife value to be used in the exponential weighting function (see details section). Preset to 2 (should be updated by the user based on substantive context).
#' @param dyadic_weight A numerical value for the dyadic cutoff weight that represents the numerical cutoff value for temporal relevancy based on the exponential weighting function. For example, a numerical value of 0.01, indicates that an exponential weight less than 0.01 will become 0 and that events with such value (or smaller values) will not be included in the sum of the past event weights (see the details section). Set to 0 by default.
#' @param Lerneretal_2013 TRUE/FALSE. TRUE indicates that the Lerner et al. (2013) exponential weighting function will be used (see the details section). FALSE indicates that the Lerner and Lomi (2020) exponential weighting function will be used (see the details section). Set to FALSE by default.
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
#' @param priorStats TRUE/FALSE. Set to FALSE by default. TRUE indicates that the user has previously computed the sender outdegree and target indegree network statistics. Set to FALSE by default. The four-cycles network statistics is computationally burdensome. If priorStats =TRUE, the function speeds things up by setting the statistic for an event dyad to 0 if either a) the current event sender was not a sender in a previous event or b) the current event receiver was not a receiver in a past event, then the four-cycles statistics for that event dyad will be 0.
#' @param sender_OutDeg If priorStats = TRUE, the vector of previously computed sender outdegree scores.
#' @param receiver_InDeg If priorStats = TRUE, the vector of previously computed receiver indegree scores.
#' @import data.table
#' @importFrom collapse fnrow
#' @importFrom collapse funique
#' @importFrom collapse whichv
#' @return The vector of four-cycle statistics for the two-mode relational event sequence.
#' @export
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `computeFourCycles()` has been deprecated starting on version 1.0.0 of the `dream` package. Please use the `remstats_fourcycles()` function and see the `NEWS.md` file for more details.
#'
#'
#' The function computes the four-cycles network sufficient statistic for a two-mode relational
#' sequence with the exponential weighting function (Lerner and Lomi 2020). In essence, the
#' four-cycles measure captures the tendency for clustering to occur in the network of past
#' events, whereby an event is more likely to occur between a sender node *a* and receiver
#' node *b* given that *a* has interacted with other receivers in past events who have
#' received events from other senders that interacted with *b* (e.g., Duxbury and Haynie 2021, Lerner and Lomi 2020). The function
#' allows users to use two different weighting functions, reduce computational runtime, employ a
#' sliding windows framework for large relational sequences, and specify a dyadic cutoff for relational relevancy.
#'
#'@details The function calculates the four-cycles network statistic for two-mode relational event models
#'based on the exponential weighting function used in either Lerner and Lomi
#'(2020) or Lerner et al. (2013).
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
#'past event times that meet the weight subset (in this case, all events that
#'have the same sender and receiver), and \eqn{T_{1/2}} is the halflife parameter.
#'
#'The formula for four-cycles for event \eqn{e_i} is:
#'\deqn{four cycles_{e_{i}} = \sqrt[3]{\sum_{s' and r'} w(s', r, t) \cdot w(s, r', t) \cdot w(s', r', t)}}
#'
#'That is, the four-cycle measure captures all the past event structures in which the
#'current event pair, sender *s* and target *r* close a four-cycle. In particular, it
#'finds all events in which: a past sender *s'* had a relational event with
#'target *r*, a past target *r'* had a relational event with current sender *s*, and finally,
#'a relational event occurred between sender *s'* and target *r'*.
#'
#'Four-cycles are computationally expensive, especially for large relational event
#'sequences (see Lerner and Lomi 2020 for a discussion on this), therefore this
#'function allows the user to input previously computed target indegree and sender
#'outdegree scores to reduce the runtime. Relational events where
#'either the event target or event sender were not involved in any prior relational
#'events (i.e., a target indegree or sender outdegree score of 0) will close no-four
#'cycles. This function exploits this feature.
#'
#'
#'Moreover, researchers interested in modeling temporal relevancy (see Quintane,
#'Mood, Dunn, and Falzone 2022; Lerner and Lomi 2020) can specify the dyadic
#'weight cutoff, that is, the minimum value for which the weight is considered
#'relationally relevant. Users who do not know the specific dyadic cutoff value to use, can use the
#'\code{\link{computeRemDyadCut}} function.
#'
#'
#'Following Lerner and Lomi (2020), if the counts of the past events are requested, the formula for four-cycles formation for
#'event \eqn{e_i} is:
#'\deqn{four cycles_{e_{i}} = \sum_{i=1}^{|S'|} \sum_{j=1}^{|R'|} \min\left[d(s'_{i}, r, t),\ d(s, r'_{j}, t),\ d(s'_{i}, r'_{j}, t)\right]}
#'where, \eqn{d()} is the number of past events that meet the specific set operations, \eqn{d(s'_{i},r,t)} is the number
#'of past events where the current event receiver received a tie from another sender \eqn{s'_{i}}, \eqn{d(s,r'_{j},t)} is the number
#'of past events where the current event sender sent a tie to another receiver \eqn{r'_{j}}, and \eqn{d(s'_{i},r'_{j},t)} is the
#'number of past events where the sender \eqn{s'_{i}} sent a tie to the receiver \eqn{r'_{j}}. Moreover, the counting
#'equation can leverage relational relevancy, by specifying the halflife parameter, exponential
#'weighting function, and the dyadic cut off weight values (see the above sections for help with this). If the user is not interested in modeling
#'relational relevancy, then those value should be left at their default values.


#'@examples
#'data("WikiEvent2018.first100k")
#'WikiEvent2018 <- WikiEvent2018.first100k[1:1000,] #the first one thousand events
#'WikiEvent2018$time <- as.numeric(WikiEvent2018$time) #making the variable numeric
#'### Creating the EventSet By Employing Case-Control Sampling With M = 5 and
#'### Sampling from the Observed Event Sequence with P = 0.01
#'EventSet <- processTMEventSeq(
#'  data = WikiEvent2018, # The Event Dataset
#'  time = WikiEvent2018$time, # The Time Variable
#'  eventID = WikiEvent2018$eventID, # The Event Sequence Variable
#'  sender = WikiEvent2018$user, # The Sender Variable
#'  receiver = WikiEvent2018$article, # The Receiver Variable
#'  p_samplingobserved = 0.01, # The Probability of Selection
#'  n_controls = 8, # The Number of Controls to Sample from the Full Risk Set
#'  seed = 9999) # The Seed for Replication
#'
#'#### Estimating the Four-Cycle Statistic Without the Sliding Windows Framework
#'EventSet$fourcycle <- computeFourCycles(
#'    observed_time = WikiEvent2018$time,
#'    observed_sender = WikiEvent2018$user,
#'    observed_receiver = WikiEvent2018$article,
#'    processed_time = EventSet$time,
#'    processed_sender = EventSet$sender,
#'    processed_receiver = EventSet$receiver,
#'    halflife = 2.592e+09, #halflife parameter
#'    dyadic_weight = 0,
#'    Lerneretal_2013 = FALSE)
#'
#'#### Estimating the Four-Cycle Statistic With the Sliding Windows Framework
#'EventSet$cycle4SW <- computeFourCycles(
#'    observed_time = WikiEvent2018$time,
#'    observed_sender = WikiEvent2018$user,
#'    observed_receiver = WikiEvent2018$article,
#'    processed_time = EventSet$time,
#'    processed_sender = EventSet$sender,
#'    processed_receiver = EventSet$receiver,
#'    processed_seqIDs = EventSet$sequenceID,
#'    halflife = 2.592e+09, #halflife parameter
#'    dyadic_weight = 0,
#'    sliding_window = TRUE,
#'    Lerneretal_2013 = FALSE)
#'
#'#The results with and without the sliding windows are the same (see correlation
#'#below). Using the sliding windows method is recommended when the data are
#'#big' so that memory allotment is more efficient.
#'cor(EventSet$fourcycle, EventSet$cycle4SW)
#'
#'#### Estimating the Four-Cycle Statistic  with the Counts of Events Returned
#'EventSet$cycle4C <- computeFourCycles(
#'    observed_time = WikiEvent2018$time,
#'    observed_sender = WikiEvent2018$user,
#'    observed_receiver = WikiEvent2018$article,
#'    processed_time = EventSet$time,
#'    processed_sender = EventSet$sender,
#'    processed_receiver = EventSet$receiver,
#'    processed_seqIDs = EventSet$sequenceID,
#'    halflife = 2.592e+09, #halflife parameter
#'    dyadic_weight = 0,
#'    sliding_window = FALSE,
#'    counts = TRUE,
#'    Lerneretal_2013 = FALSE)
#'
#'cbind(EventSet$fourcycle,
#'      EventSet$cycle4SW,
#'      EventSet$cycle4C)
#'
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#' @references
#'Duxbury, Scott and Dana Haynie. 2021. "Shining a Light on the Shadows: Endogenous Trade
#'Structure and the Growth of an Online Illegal Market." *American Journal of Sociology* 127(3): 787-827.
#'
#'
#'Quintane, Eric, Martin Wood, John Dunn, and Lucia Falzon. 2022. “Temporal
#'Brokering: A Measure of Brokerage as a Behavioral Process.” *Organizational Research Methods*
#'25(3): 459-489.
#'
#'Lerner, Jürgen and Alessandro Lomi. 2020. “Reliability of relational event
#'model estimates under sampling: How to fit a relational event model to 360
#'million dyadic events.” *Network Science* 8(1): 97-135.
#'
#'Lerner, Jürgen, Margit Bussman, Tom A.B. Snijders, and Ulrik Brandes. 2013. "Modeling
#'Frequency and Type of Interaction in Event Networks." *The Corvinus Journal of Sociology and Social Policy* 4(1): 3-32.
#'
#'
#
#'
#'
#'

########################################################################################################
#  Events = the full event sequence
#  eventSet = the sampled event sequence
#  Time = the name for the time variable
#  sender = the name for the sender variable
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
#  Prior Computations: did they already compute indegree or outdegree?
########################################################################################################

computeFourCycles <- function(observed_time, # variable (column) name that contains the time variable
                            observed_sender, # variable (column) name that contains the sender variable
                            observed_receiver, # variable (column) name that contains the receiver variable
                            processed_time,
                            processed_sender, # variable (column) name that contains the sender variable
                            processed_receiver, # variable (column) name that contains the receiver variable
                            sliding_windows = FALSE, # TRUE = we want to use the sliding windows framework
                            processed_seqIDs = NULL, #If true, the user should insert the placement of sampled events in the original file
                            counts = FALSE, #Logical indicating if the raw counts of events should be returned or the exponential weighting function should be used (TRUE = counts; FALSE = exponential weighting)
                            halflife=2, # the half life value for the weighting function
                            dyadic_weight=0.00, # dyadic cutoff weight for events that no longer matter
                            window_size = NA, # the sizes of the windows that we will use, if NA, we will compute it internally
                            Lerneretal_2013 = FALSE,
                            priorStats = FALSE, # Did the user already compute the indegree or outdegree?
                            sender_OutDeg = NULL,
                            receiver_InDeg = NULL
) {
  ################################################################################
  #       Note: Four Cycle in Relational Event Model is a very expensive statistic
  #       to compute. Therefore, this function provides the user the possibility to
  #       include the outdegree and indegree statistic. Since, four-cycles is
  #       mathematically a lower order function of the event users, indegree and outdegree,
  #       then, if the values are 0 for either one of the statisitcs, then it's not possible
  #       for a four-cycle to be closed.
  ################################################################################

  #base::cat("Checking Data Structure and User Inputs.......") # outputting status to user

  ########################################################
  #
  #   Checking for Errors in User Inputs
  #
  ########################################################
  lifecycle::deprecate_warn("1.0.0", "computeFourCycles()", "remstats_fourcycles()")



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
  ########################################################
  # Renaming the columns to match the user inputs
  ########################################################
  events <- data.table::data.table(sender = observed_sender,# renaming the sender column
                                   receiver = observed_receiver, # renaming the receiver column
                                   time = observed_time)# renaming the time column

  eventSet <- data.table::data.table(sender = processed_sender,# renaming the sender column
                                     receiver = processed_receiver, # renaming the receiver column
                                     time = processed_time)# renaming the time column



  if (priorStats == TRUE) { # if the user entered prior computations
    eventSet$outdegree <- sender_OutDeg # creating a new variable for the outdegree object based on user entry
    eventSet$indegree <-  receiver_InDeg # creating a new variable for the indegree object based on user entry
  }

  ########################################################
  #### Clearing User Inputs for Memory
  ########################################################
  rm(list = c("observed_time", "observed_sender", "observed_receiver",
              "processed_time", "processed_sender", "processed_receiver",
              "sender_OutDeg", "receiver_InDeg"))


  fourcycle <- rep(0, nrow(eventSet)) # Creating a Empty Vector To Store the Computed Statistics
  ########################################################
  #
  #    If the User Did not want to use the sliding windows framework
  #
  ########################################################
  if (sliding_windows == FALSE) {
    ########################################################
    # If there were prior computations
    ########################################################
    if (priorStats == TRUE) {
      for (i in 1:nrow(eventSet)) { # for all sampled events that we need to compute values for

        #  If outdegree or indegree is equal to 0, then there will be no possible four-cycle values,
        #  therefore we skip this event
        if (eventSet[i]$outdegree == 0 | eventSet[i]$indegree == 0) {
          next
        } else {
          senderi <- eventSet$sender[i] # current sender
          timei <- eventSet$time[i] # current time
          receiveri <- eventSet$receiver[i] # current receiver

          #### getting all level 2 ties!
          targetallevents <- events[events$receiver == receiveri] # past events that the receiver (receiver) is involved in
          if(collapse::fnrow(targetallevents) == 0){next}
          senderallevents <- events[events$sender == senderi] # past events that the sender (sender) is involved in
          if(collapse::fnrow(senderallevents) == 0){next}
          #targetallevents <- targetallevents[time < timei, weightbj := remExpWeights(current = timei, past = time, halflife = halflife, dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013)] # compute the weight of prior events
          senderallevents <- senderallevents[senderallevents$time < timei]
          if(collapse::fnrow( senderallevents) == 0){next}
          senderallevents$weightia <- remExpWeights(current = timei,
                                                    past = senderallevents$time,
                                                    halflife = halflife,
                                                    dyadic_weight = dyadic_weight,
                                                    Lerneretal_2013 = Lerneretal_2013)
          senderallevents <- senderallevents[senderallevents$weightia != 0] # keep if only weight is not equal to 0 (this is important as we are removing events before computation)
          if(collapse::fnrow(senderallevents) == 0){next}
          targetallevents <- targetallevents[targetallevents$time < timei]
          if(collapse::fnrow(targetallevents) == 0){next}
          # targetallevents <- targetallevents[time < timei, weightbj := remExpWeights(current = timei, past = time, halflife = halflife, dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013)] # compute the weight of prior events
          targetallevents$weightbj <- remExpWeights(current = timei,
                                                    past = targetallevents$time,
                                                    halflife = halflife,
                                                    dyadic_weight = dyadic_weight,
                                                    Lerneretal_2013 = Lerneretal_2013)

          targetallevents <- targetallevents[targetallevents$weightbj != 0] # keep if only weight is not equal to 0 (this is important as we are removing events before computation)
          if(collapse::fnrow(targetallevents) == 0){next}
          #### getting only the past events, and then computing the event weight
          ### if event weight is not equal to 0, then keep them! (this lessens the overall computation time by removing potential event orders that would sum to 0)
          #### getting only the past events, and then computing the event weight
          ### if event weight is not equal to 0, then keep them! (this lessens the overall computation time by removing potential event orders that would sum to 0)
          if (collapse::fnrow(senderallevents) == 0 | collapse::fnrow(targetallevents) == 0) {
            next
          } # if there are no events with anyone involved (i.e, the same events)
          ## getting only unique cases
          pastsenderloop <- collapse::funique(targetallevents$sender)
          ## getting only unique cases
          pasttargetloop <- collapse::funique(senderallevents$receiver)
          ##### removing the current sender from the unique list
          inloop <- collapse::whichv(pastsenderloop, senderi)
          if (length(inloop) != 0) {
            pastsenderloop <- pastsenderloop[-(inloop)] # remove sender
          }
          ##### removing the current target from the unique list
          inloop <- collapse::whichv(pasttargetloop, receiveri)
          if (length(inloop) != 0) {
            pasttargetloop <- pasttargetloop[-(inloop)] # remove target
          }
          if (length(pastsenderloop) == 0 | length(pasttargetloop) == 0) {
            next
          } # if there are no events with anyone involved (i.e, the same events)
         # data.table::setkey(events, receiver) # using the data.table set key function [no need to set key anymore]
          targets_with_past_senders <- events[events$receiver %in% pasttargetloop]#DT original code: events[J(pasttargetloop)]
          if(collapse::fnrow( targets_with_past_senders) == 0){next}
          targets_with_past_senders <- targets_with_past_senders[targets_with_past_senders$time < timei] # give me only the prior time events
          if(collapse::fnrow( targets_with_past_senders) == 0){next}
          #targets_with_past_senders[, weight := remExpWeights(current = timei, past = time, halflife = halflife, dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013)] # compute the weight of all events
          targets_with_past_senders$weight <- remExpWeights(current = timei,
                                                            past = targets_with_past_senders$time,
                                                            halflife = halflife,
                                                            dyadic_weight = dyadic_weight,
                                                            Lerneretal_2013 = Lerneretal_2013)
          #targets_with_past_senders[, weight := remExpWeights(current = timei, past = time, halflife = halflife, dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013)] # compute the weight of all events
          targets_with_past_senders <- targets_with_past_senders[targets_with_past_senders$weight != 0] # keep if only weight is not equal to 0 (this is important as we are removing events before computation)
          if(collapse::fnrow( targets_with_past_senders) == 0){next}
          # all events that have the unique targets within for the creation of four cycles
          ### if event weight is not equal to 0, then keep them! (this lessens the overall computation time by removing potential event orders that would sum to 0)
          senders_with_past_targets <- events[events$sender %in% pastsenderloop] # all events that have the unique targets within for the creation of four cycles
          if(collapse::fnrow( senders_with_past_targets ) == 0){next}
          senders_with_past_targets <- senders_with_past_targets[senders_with_past_targets$time < timei] # give me only the prior time events
          if(collapse::fnrow( senders_with_past_targets ) == 0){next}
          #senders_with_past_targets <- #senders_with_past_targets[, weight := remExpWeights(current = timei, past = time, halflife = halflife, dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013)] # compute the weight of all events
          senders_with_past_targets$weight <-  remExpWeights(current = timei,
                                                             past = senders_with_past_targets$time,
                                                             halflife = halflife,
                                                             dyadic_weight = dyadic_weight,
                                                             Lerneretal_2013 = Lerneretal_2013)
          senders_with_past_targets <- senders_with_past_targets[senders_with_past_targets$weight != 0] # keep if only weight is not equal to 0 (this is important as we are removing events before computation)
          if(collapse::fnrow( senders_with_past_targets ) == 0){next}
          ### if event weight is not equal to 0, then keep them! (this lessens the overall computation time by removing potential event orders that would sum to 0)
          # if there are no events with anyone involved (i.e, the same events)
          #### combining the two data sets to give me all prior events with senders that share the same receiver as the current and receivers that share the same sender as past receivers
          open_events <- senders_with_past_targets[targets_with_past_senders,
                                                   on = colnames(senders_with_past_targets)[-which(colnames(senders_with_past_targets) == "weight")],
                                                   nomatch = 0
          ]
          if(collapse::fnrow( open_events ) == 0){next}
          #### getting only the past events, and then computing the event weight
          #open_events <- open_events[time < timei, weightij := remExpWeights(current = timei, past = time, halflife = halflife, dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013)]
          open_events <- open_events[open_events$time < timei]
          if(collapse::fnrow( open_events ) == 0){next}
          open_events$weightij <-  remExpWeights(current = timei,
                                                 past = open_events$time,
                                                 halflife = halflife,
                                                 dyadic_weight = dyadic_weight,
                                                 Lerneretal_2013 = Lerneretal_2013)
          open_events <- open_events[open_events$weightij != 0] # keep if only weight is not equal to 0 (this is important as we are removing events before computation)
          ### if event weight is not equal to 0, then keep them! (this lessens the overall computation time by removing potential event orders that would sum to 0)
          if (collapse::fnrow(open_events) == 0) {
            next
          } # if there are no events in which the current event is being closed (i.e, the same events)
          ### merging the event cases from w(i,j) and w(i,A)
          if(counts == TRUE){
            counti <- 0
            pairs <- unique(open_events,by = c("sender","receiver")) #getting the unique pairs for the summation count
            for(z in 1:nrow(pairs)){
              counti <- counti + base::min(sum(senderallevents$receiver == pairs$receiver[z]), #d(s,r')
                                           sum(targetallevents$sender == pairs$sender[z]),#d(r,s')
                                           sum(open_events$sender == pairs$sender[z] &
                                                 open_events$receiver  == pairs$receiver[z]))#d(s',r')
            }
            fourcycle[i]  <- counti
            ###### clearing memory for consumption purposes
            rm(list = c(
              "open_events", "senders_with_past_targets", "targets_with_past_senders",
              "targetallevents", "pastsenderloop", "pasttargetloop", "senderallevents"
            ))
          }else{

            ### merging the event cases from w(i,j) and w(i,A)
            d31 <- data.table::merge.data.table(open_events, senderallevents, by = "receiver", all.x = T, all.y = T, allow.cartesian = T)
            if(collapse::fnrow( d31) == 0){next}
            d31t <- d31[stats::complete.cases(d31)] # complete cases, removing the cases that dont fully merge (i.e., are missing one of the weights)
            if(collapse::fnrow( d31t) == 0){next}
            ### merging the event cases from w(i,j) and w(i,A) to now include w(U, j)
            d31t3 <- data.table::merge.data.table(d31t, targetallevents, by.x = "sender.x", by.y = "sender", all.x = T, all.y = T, allow.cartesian = T)
            if(collapse::fnrow( d31t3) == 0){next}
            d31t4 <- d31t3[stats::complete.cases(d31t3)] # complete cases, removing the cases that dont fully merge (i.e., are missing one of the weights)
            if(collapse::fnrow(d31t4) == 0){next}
            d31t5 <- unique(d31t4) # getting only non-redudant unique cases
            if(collapse::fnrow( d31t5) == 0){next}
            # d31t5[, product := with(d31t5, weightij * weightia * weightbj)] # now computing the product of all three weights, but removing the instances whereby events are recycled from the past 1 million event dataset
            #fourcycle[i] <- sum(d31t5$product)^(1 / 3) # sum all values, take the cubic root, and then finally, add that value to the computed vector ith spot
            fourcycle[i] <- sum(d31t5$weightij * d31t5$weightia * d31t5$weightbj )^(1 / 3) # sum all values, take the cubic root, and then finally, add that value to the computed vector ith spot
            ###### clearing memory for consumption purposes
            rm(list = c(
              "open_events", "senders_with_past_targets", "targets_with_past_senders",
              "targetallevents", "pastsenderloop", "pasttargetloop", "senderallevents",
              "d31", "d31t", "d31t3", "d31t4", "d31t5"
            ))
            }                ###### clearing memory for consumption purposes

        }
      }
    }

    ########################################################
    # If there were not prior computations
    ########################################################
    if (priorStats == FALSE) {
      for (i in 1:nrow(eventSet)) { # for all sampled events that we need to compute values for
        senderi <- eventSet$sender[i] # current sender
        timei <- eventSet$time[i] # current time
        receiveri <- eventSet$receiver[i] # current receiver

        #### getting all level 2 ties!
        targetallevents <- events[events$receiver == receiveri] # past events that the receiver (receiver) is involved in
        if(collapse::fnrow(targetallevents) == 0){next}
        senderallevents <- events[events$sender == senderi] # past events that the sender (sender) is involved in
        if(collapse::fnrow(senderallevents) == 0){next}
        #targetallevents <- targetallevents[time < timei, weightbj := remExpWeights(current = timei, past = time, halflife = halflife, dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013)] # compute the weight of prior events
        senderallevents <- senderallevents[senderallevents$time < timei]
        if(collapse::fnrow( senderallevents) == 0){next}
        senderallevents$weightia <- remExpWeights(current = timei,
                                                  past = senderallevents$time,
                                                  halflife = halflife,
                                                  dyadic_weight = dyadic_weight,
                                                  Lerneretal_2013 = Lerneretal_2013)
        senderallevents <- senderallevents[senderallevents$weightia != 0] # keep if only weight is not equal to 0 (this is important as we are removing events before computation)
        if(collapse::fnrow(senderallevents) == 0){next}
        targetallevents <- targetallevents[targetallevents$time < timei]
        if(collapse::fnrow(targetallevents) == 0){next}
        # targetallevents <- targetallevents[time < timei, weightbj := remExpWeights(current = timei, past = time, halflife = halflife, dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013)] # compute the weight of prior events
        targetallevents$weightbj <- remExpWeights(current = timei,
                                                  past = targetallevents$time,
                                                  halflife = halflife,
                                                  dyadic_weight = dyadic_weight,
                                                  Lerneretal_2013 = Lerneretal_2013)

        targetallevents <- targetallevents[targetallevents$weightbj != 0] # keep if only weight is not equal to 0 (this is important as we are removing events before computation)
        if(collapse::fnrow(targetallevents) == 0){next}
        #### getting only the past events, and then computing the event weight
        ### if event weight is not equal to 0, then keep them! (this lessens the overall computation time by removing potential event orders that would sum to 0)
        #### getting only the past events, and then computing the event weight
        ### if event weight is not equal to 0, then keep them! (this lessens the overall computation time by removing potential event orders that would sum to 0)
        if (collapse::fnrow(senderallevents) == 0 | collapse::fnrow(targetallevents) == 0) {
          next
        } # if there are no events with anyone involved (i.e, the same events)
        ## getting only unique cases
        pastsenderloop <- collapse::funique(targetallevents$sender)
        ## getting only unique cases
        pasttargetloop <- collapse::funique(senderallevents$receiver)
        ##### removing the current sender from the unique list
        inloop <- collapse::whichv(pastsenderloop, senderi)
        if (length(inloop) != 0) {
          pastsenderloop <- pastsenderloop[-(inloop)] # remove sender
        }
        ##### removing the current target from the unique list
        inloop <- collapse::whichv(pasttargetloop, receiveri)
        if (length(inloop) != 0) {
          pasttargetloop <- pasttargetloop[-(inloop)] # remove target
        }
        if (length(pastsenderloop) == 0 | length(pasttargetloop) == 0) {
          next
        } # if there are no events with anyone involved (i.e, the same events)
        # data.table::setkey(events, receiver) # using the data.table set key function [no need to set key anymore]
        targets_with_past_senders <- events[events$receiver %in% pasttargetloop]#DT original code: events[J(pasttargetloop)]
        if(collapse::fnrow( targets_with_past_senders) == 0){next}
        targets_with_past_senders <- targets_with_past_senders[targets_with_past_senders$time < timei] # give me only the prior time events
        if(collapse::fnrow( targets_with_past_senders) == 0){next}
        #targets_with_past_senders[, weight := remExpWeights(current = timei, past = time, halflife = halflife, dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013)] # compute the weight of all events
        targets_with_past_senders$weight <- remExpWeights(current = timei,
                                                          past = targets_with_past_senders$time,
                                                          halflife = halflife,
                                                          dyadic_weight = dyadic_weight,
                                                          Lerneretal_2013 = Lerneretal_2013)
        #targets_with_past_senders[, weight := remExpWeights(current = timei, past = time, halflife = halflife, dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013)] # compute the weight of all events
        targets_with_past_senders <- targets_with_past_senders[targets_with_past_senders$weight != 0] # keep if only weight is not equal to 0 (this is important as we are removing events before computation)
        if(collapse::fnrow( targets_with_past_senders) == 0){next}
        # all events that have the unique targets within for the creation of four cycles
        ### if event weight is not equal to 0, then keep them! (this lessens the overall computation time by removing potential event orders that would sum to 0)
        senders_with_past_targets <- events[events$sender %in% pastsenderloop] # all events that have the unique targets within for the creation of four cycles
        if(collapse::fnrow( senders_with_past_targets ) == 0){next}
        senders_with_past_targets <- senders_with_past_targets[senders_with_past_targets$time < timei] # give me only the prior time events
        if(collapse::fnrow( senders_with_past_targets ) == 0){next}
        #senders_with_past_targets <- #senders_with_past_targets[, weight := remExpWeights(current = timei, past = time, halflife = halflife, dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013)] # compute the weight of all events
        senders_with_past_targets$weight <-  remExpWeights(current = timei,
                                                           past = senders_with_past_targets$time,
                                                           halflife = halflife,
                                                           dyadic_weight = dyadic_weight,
                                                           Lerneretal_2013 = Lerneretal_2013)
        senders_with_past_targets <- senders_with_past_targets[senders_with_past_targets$weight != 0] # keep if only weight is not equal to 0 (this is important as we are removing events before computation)
        if(collapse::fnrow( senders_with_past_targets ) == 0){next}
        ### if event weight is not equal to 0, then keep them! (this lessens the overall computation time by removing potential event orders that would sum to 0)
        # if there are no events with anyone involved (i.e, the same events)
        #### combining the two data sets to give me all prior events with senders that share the same receiver as the current and receivers that share the same sender as past receivers
        open_events <- senders_with_past_targets[targets_with_past_senders,
                                                 on = colnames(senders_with_past_targets)[-which(colnames(senders_with_past_targets) == "weight")],
                                                 nomatch = 0
        ]
        if(collapse::fnrow( open_events ) == 0){next}
        #### getting only the past events, and then computing the event weight
        #open_events <- open_events[time < timei, weightij := remExpWeights(current = timei, past = time, halflife = halflife, dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013)]
        open_events <- open_events[open_events$time < timei]
        if(collapse::fnrow( open_events ) == 0){next}
        open_events$weightij <-  remExpWeights(current = timei,
                                               past = open_events$time,
                                               halflife = halflife,
                                               dyadic_weight = dyadic_weight,
                                               Lerneretal_2013 = Lerneretal_2013)
        open_events <- open_events[open_events$weightij != 0] # keep if only weight is not equal to 0 (this is important as we are removing events before computation)
        ### if event weight is not equal to 0, then keep them! (this lessens the overall computation time by removing potential event orders that would sum to 0)
        if (collapse::fnrow(open_events) == 0) {
          next
        } # if there are no events in which the current event is being closed (i.e, the same events)
        ### merging the event cases from w(i,j) and w(i,A)
        if(counts == TRUE){
          counti <- 0
          pairs <- unique(open_events,by = c("sender","receiver")) #getting the unique pairs for the summation count
          for(z in 1:nrow(pairs)){
            counti <- counti + base::min(sum(senderallevents$receiver == pairs$receiver[z]), #d(s,r')
                                         sum(targetallevents$sender == pairs$sender[z]),#d(r,s')
                                         sum(open_events$sender == pairs$sender[z] &
                                               open_events$receiver  == pairs$receiver[z]))#d(s',r')
          }
          fourcycle[i]  <- counti
          ###### clearing memory for consumption purposes
          rm(list = c(
            "open_events", "senders_with_past_targets", "targets_with_past_senders",
            "targetallevents", "pastsenderloop", "pasttargetloop", "senderallevents"
          ))
        }else{

          ### merging the event cases from w(i,j) and w(i,A)
          d31 <- data.table::merge.data.table(open_events, senderallevents, by = "receiver", all.x = T, all.y = T, allow.cartesian = T)
          if(collapse::fnrow( d31) == 0){next}
          d31t <- d31[stats::complete.cases(d31)] # complete cases, removing the cases that dont fully merge (i.e., are missing one of the weights)
          if(collapse::fnrow( d31t) == 0){next}
          ### merging the event cases from w(i,j) and w(i,A) to now include w(U, j)
          d31t3 <- data.table::merge.data.table(d31t, targetallevents, by.x = "sender.x", by.y = "sender", all.x = T, all.y = T, allow.cartesian = T)
          if(collapse::fnrow( d31t3) == 0){next}
          d31t4 <- d31t3[stats::complete.cases(d31t3)] # complete cases, removing the cases that dont fully merge (i.e., are missing one of the weights)
          if(collapse::fnrow(d31t4) == 0){next}
          d31t5 <- unique(d31t4) # getting only non-redudant unique cases
          if(collapse::fnrow( d31t5) == 0){next}
          # d31t5[, product := with(d31t5, weightij * weightia * weightbj)] # now computing the product of all three weights, but removing the instances whereby events are recycled from the past 1 million event dataset
          #fourcycle[i] <- sum(d31t5$product)^(1 / 3) # sum all values, take the cubic root, and then finally, add that value to the computed vector ith spot
          fourcycle[i] <- sum(d31t5$weightij * d31t5$weightia * d31t5$weightbj )^(1 / 3) # sum all values, take the cubic root, and then finally, add that value to the computed vector ith spot
          ###### clearing memory for consumption purposes
          rm(list = c(
            "open_events", "senders_with_past_targets", "targets_with_past_senders",
            "targetallevents", "pastsenderloop", "pasttargetloop", "senderallevents",
            "d31", "d31t", "d31t3", "d31t4", "d31t5"
          ))
        }

      }
    }
  }

  ########################################################
  #
  #    If the User Did Want To Use the Sliding Windows Framework
  #
  ########################################################
  if (sliding_windows == TRUE) {
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
    eventSet$EVENT <- processed_seqIDs # the sampled eventnet IDS

    sliding_window$eventENDS <- eventSet[ending_blocks]$EVENT # getting the ending event from the block [the index in which the block ends]
    sliding_window$eventSTARTS <- rep(1, length(starting_blocks)) # creating a dummy variable for the starting blocks for each events (in other words, the starting of relationally relevancy)
    sliding_window$time <- eventSet[starting_blocks]$time # the eventtimes of the starting blocks (of the eventtime of the first event in each block)
    #sliding_window$minimum_times <- minimum_effective_time(sliding_window$time, dyadicweight = dyadic_weight, halflife = halflife, Lerneretal_2013 = Lerneretal_2013) # compute the minimum time weight

    if(dyadic_weight == 0){
      sliding_window$minimum_times <- 0 #no minimum event time since dyadic weight is 0
    }else{
      #compute dyadic weight
      sliding_window$minimum_times <- remExpWeights(sliding_window$time,
                                                         dyadic_weight = dyadic_weight,
                                                         halflife = halflife,
                                                         Lerneretal_2013 = Lerneretal_2013,
                                                         exp.weights = FALSE) # compute the minimum time weight

    }


    ###### Doing a quick test: it should be noted that too large of a halflife parameter in the Lerner et al. 2013
    ###### weighting function results in minimum effective time greater than the eventTime. Therefore, it cannot be
    ###### properly approximated in this case.
    test_vec <- sliding_window$minimum_times - sliding_window$time # getting the difference between the values, these should all be negative
    areanygreater <- collapse::fsum((test_vec) > 0) # checking if we have any positive values,
    if (areanygreater != 0) { # if we do, stop the function, and tell the user this
      base::stop("Error: Unfortunately the combination of the provided halflife parameter and the weighting function (i.e., using the \n
                  Lerner et al. 2013 specification) resulted in a minimum effective time that is greater than the eventTimes, therefore, \n
                  the sliding windows framework cannot be used. Please restart with sliding_windows == FALSE. Please see the documentation \n
                  for this function and the minimum effective time documentation. As always, we hope you know what you're doing.....") # stop computation and tell the user
    }
    for (c in 2:nrow(sliding_window)) { # from 2 to all other windows
      sliding_window$eventSTARTS[c] <- base::min(which(events$time >= sliding_window$minimum_times[c])) # getting the minimum event that provides relevancy for each block
    }

    ######################################################################
    ####### Creating the Vector to Store Event Weights
    ######################################################################
    n_blocks <- nrow(sliding_window) # number of observed events in dataset

    if (priorStats == TRUE) { # if the user provided prior computations

      for (j in 1:n_blocks) { # for all blocks from the sliding windows framework
        eventsj <- events[sliding_window$eventSTARTS[j]:sliding_window$eventENDS[j]] # get the current events information for the jth block
        new_start <- sliding_window$start_block[j] # the starting sampled event for block j
        new_stop <- sliding_window$stop_block[j] # the ending sampled event for block j

        for (i in new_start:new_stop) { # for all sampled events in block j

          #  If outdegree or indegree is equal to 0, then there will be no possible four-cycle values,
          #  therefore we skip this event
          if (eventSet[i]$outdegree == 0 | eventSet[i]$indegree == 0) {
            next
          } else {
            senderi <- eventSet$sender[i] # current sender
            timei <- eventSet$time[i] # current time
            receiveri <- eventSet$receiver[i] # current receiver

            #### getting all level 2 ties!
            targetallevents <- eventsj[eventsj$receiver == receiveri] # past events that the receiver (receiver) is involved in
            if(collapse::fnrow(targetallevents) == 0){next}
            senderallevents <- eventsj[eventsj$sender == senderi] # past events that the sender (sender) is involved in
            if(collapse::fnrow(senderallevents) == 0){next}
            #targetallevents <- targetallevents[time < timei, weightbj := remExpWeights(current = timei, past = time, halflife = halflife, dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013)] # compute the weight of prior events
            senderallevents <- senderallevents[senderallevents$time < timei]
            if(collapse::fnrow( senderallevents) == 0){next}
            senderallevents$weightia <- remExpWeights(current = timei,
                                                      past = senderallevents$time,
                                                      halflife = halflife,
                                                      dyadic_weight = dyadic_weight,
                                                      Lerneretal_2013 = Lerneretal_2013)
            senderallevents <- senderallevents[senderallevents$weightia != 0] # keep if only weight is not equal to 0 (this is important as we are removing events before computation)
            if(collapse::fnrow(senderallevents) == 0){next}
            targetallevents <- targetallevents[targetallevents$time < timei]
            if(collapse::fnrow(targetallevents) == 0){next}
            # targetallevents <- targetallevents[time < timei, weightbj := remExpWeights(current = timei, past = time, halflife = halflife, dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013)] # compute the weight of prior events
            targetallevents$weightbj <- remExpWeights(current = timei,
                                                      past = targetallevents$time,
                                                      halflife = halflife,
                                                      dyadic_weight = dyadic_weight,
                                                      Lerneretal_2013 = Lerneretal_2013)

            targetallevents <- targetallevents[targetallevents$weightbj != 0] # keep if only weight is not equal to 0 (this is important as we are removing events before computation)
            if(collapse::fnrow(targetallevents) == 0){next}
            #### getting only the past events, and then computing the event weight
            ### if event weight is not equal to 0, then keep them! (this lessens the overall computation time by removing potential event orders that would sum to 0)
            #### getting only the past events, and then computing the event weight
            ### if event weight is not equal to 0, then keep them! (this lessens the overall computation time by removing potential event orders that would sum to 0)
            if (collapse::fnrow(senderallevents) == 0 | collapse::fnrow(targetallevents) == 0) {
              next
            } # if there are no events with anyone involved (i.e, the same events)
            ## getting only unique cases
            pastsenderloop <- collapse::funique(targetallevents$sender)
            ## getting only unique cases
            pasttargetloop <- collapse::funique(senderallevents$receiver)
            ##### removing the current sender from the unique list
            inloop <- collapse::whichv(pastsenderloop, senderi)
            if (length(inloop) != 0) {
              pastsenderloop <- pastsenderloop[-(inloop)] # remove sender
            }
            ##### removing the current target from the unique list
            inloop <- collapse::whichv(pasttargetloop, receiveri)
            if (length(inloop) != 0) {
              pasttargetloop <- pasttargetloop[-(inloop)] # remove target
            }
            if (length(pastsenderloop) == 0 | length(pasttargetloop) == 0) {
              next
            } # if there are no events with anyone involved (i.e, the same events)
            # data.table::setkey(events, receiver) # using the data.table set key function [no need to set key anymore]
            targets_with_past_senders <- eventsj[eventsj$receiver %in% pasttargetloop]#DT original code: events[J(pasttargetloop)]
            if(collapse::fnrow( targets_with_past_senders) == 0){next}
            targets_with_past_senders <- targets_with_past_senders[targets_with_past_senders$time < timei] # give me only the prior time events
            if(collapse::fnrow( targets_with_past_senders) == 0){next}
            #targets_with_past_senders[, weight := remExpWeights(current = timei, past = time, halflife = halflife, dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013)] # compute the weight of all events
            targets_with_past_senders$weight <- remExpWeights(current = timei,
                                                              past = targets_with_past_senders$time,
                                                              halflife = halflife,
                                                              dyadic_weight = dyadic_weight,
                                                              Lerneretal_2013 = Lerneretal_2013)
            #targets_with_past_senders[, weight := remExpWeights(current = timei, past = time, halflife = halflife, dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013)] # compute the weight of all events
            targets_with_past_senders <- targets_with_past_senders[targets_with_past_senders$weight != 0] # keep if only weight is not equal to 0 (this is important as we are removing events before computation)
            if(collapse::fnrow( targets_with_past_senders) == 0){next}
            # all events that have the unique targets within for the creation of four cycles
            ### if event weight is not equal to 0, then keep them! (this lessens the overall computation time by removing potential event orders that would sum to 0)
            senders_with_past_targets <- eventsj[eventsj$sender %in% pastsenderloop] # all events that have the unique targets within for the creation of four cycles
            if(collapse::fnrow( senders_with_past_targets ) == 0){next}
            senders_with_past_targets <- senders_with_past_targets[senders_with_past_targets$time < timei] # give me only the prior time events
            if(collapse::fnrow( senders_with_past_targets ) == 0){next}
            #senders_with_past_targets <- #senders_with_past_targets[, weight := remExpWeights(current = timei, past = time, halflife = halflife, dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013)] # compute the weight of all events
            senders_with_past_targets$weight <-  remExpWeights(current = timei,
                                                               past = senders_with_past_targets$time,
                                                               halflife = halflife,
                                                               dyadic_weight = dyadic_weight,
                                                               Lerneretal_2013 = Lerneretal_2013)
            senders_with_past_targets <- senders_with_past_targets[senders_with_past_targets$weight != 0] # keep if only weight is not equal to 0 (this is important as we are removing events before computation)
            if(collapse::fnrow( senders_with_past_targets ) == 0){next}
            ### if event weight is not equal to 0, then keep them! (this lessens the overall computation time by removing potential event orders that would sum to 0)
            # if there are no events with anyone involved (i.e, the same events)
            #### combining the two data sets to give me all prior events with senders that share the same receiver as the current and receivers that share the same sender as past receivers
            open_events <- senders_with_past_targets[targets_with_past_senders,
                                                     on = colnames(senders_with_past_targets)[-which(colnames(senders_with_past_targets) == "weight")],
                                                     nomatch = 0
            ]
            if(collapse::fnrow( open_events ) == 0){next}
            #### getting only the past events, and then computing the event weight
            #open_events <- open_events[time < timei, weightij := remExpWeights(current = timei, past = time, halflife = halflife, dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013)]
            open_events <- open_events[open_events$time < timei]
            if(collapse::fnrow( open_events ) == 0){next}
            open_events$weightij <-  remExpWeights(current = timei,
                                                   past = open_events$time,
                                                   halflife = halflife,
                                                   dyadic_weight = dyadic_weight,
                                                   Lerneretal_2013 = Lerneretal_2013)
            open_events <- open_events[open_events$weightij != 0] # keep if only weight is not equal to 0 (this is important as we are removing events before computation)
            ### if event weight is not equal to 0, then keep them! (this lessens the overall computation time by removing potential event orders that would sum to 0)
            if (collapse::fnrow(open_events) == 0) {
              next
            } # if there are no events in which the current event is being closed (i.e, the same events)
            ### merging the event cases from w(i,j) and w(i,A)
            if(counts == TRUE){
              counti <- 0
              pairs <- unique(open_events,by = c("sender","receiver")) #getting the unique pairs for the summation count
              for(z in 1:nrow(pairs)){
                counti <- counti + base::min(sum(senderallevents$receiver == pairs$receiver[z]), #d(s,r')
                                             sum(targetallevents$sender == pairs$sender[z]),#d(r,s')
                                             sum(open_events$sender == pairs$sender[z] &
                                                   open_events$receiver  == pairs$receiver[z]))#d(s',r')
              }
              fourcycle[i]  <- counti
              ###### clearing memory for consumption purposes
              rm(list = c(
                "open_events", "senders_with_past_targets", "targets_with_past_senders",
                "targetallevents", "pastsenderloop", "pasttargetloop", "senderallevents"
              ))
            }else{

              ### merging the event cases from w(i,j) and w(i,A)
              d31 <- data.table::merge.data.table(open_events, senderallevents, by = "receiver", all.x = T, all.y = T, allow.cartesian = T)
              if(collapse::fnrow( d31) == 0){next}
              d31t <- d31[stats::complete.cases(d31)] # complete cases, removing the cases that dont fully merge (i.e., are missing one of the weights)
              if(collapse::fnrow( d31t) == 0){next}
              ### merging the event cases from w(i,j) and w(i,A) to now include w(U, j)
              d31t3 <- data.table::merge.data.table(d31t, targetallevents, by.x = "sender.x", by.y = "sender", all.x = T, all.y = T, allow.cartesian = T)
              if(collapse::fnrow( d31t3) == 0){next}
              d31t4 <- d31t3[stats::complete.cases(d31t3)] # complete cases, removing the cases that dont fully merge (i.e., are missing one of the weights)
              if(collapse::fnrow(d31t4) == 0){next}
              d31t5 <- unique(d31t4) # getting only non-redudant unique cases
              if(collapse::fnrow( d31t5) == 0){next}
              # d31t5[, product := with(d31t5, weightij * weightia * weightbj)] # now computing the product of all three weights, but removing the instances whereby events are recycled from the past 1 million event dataset
              #fourcycle[i] <- sum(d31t5$product)^(1 / 3) # sum all values, take the cubic root, and then finally, add that value to the computed vector ith spot
              fourcycle[i] <- sum(d31t5$weightij * d31t5$weightia * d31t5$weightbj )^(1 / 3) # sum all values, take the cubic root, and then finally, add that value to the computed vector ith spot
              ###### clearing memory for consumption purposes
              rm(list = c(
                "open_events", "senders_with_past_targets", "targets_with_past_senders",
                "targetallevents", "pastsenderloop", "pasttargetloop", "senderallevents",
                "d31", "d31t", "d31t3", "d31t4", "d31t5"
              ))
            }                ###### clearing memory for consumption purposes

          }











        }
      }
    }
    if (priorStats == FALSE) { # if the user did not provide prior computations

      for (j in 1:n_blocks) { # for all blocks from the sliding windows framework
        eventsj <- events[sliding_window$eventSTARTS[j]:sliding_window$eventENDS[j]] # get the current events information for the jth block
        new_start <- sliding_window$start_block[j] # the starting sampled event for block j
        new_stop <- sliding_window$stop_block[j] # the ending sampled event for block j

        for (i in new_start:new_stop) { # for all sampled events in block j













          senderi <- eventSet$sender[i] # current sender
          timei <- eventSet$time[i] # current time
          receiveri <- eventSet$receiver[i] # current receiver

          #### getting all level 2 ties!
          targetallevents <- eventsj[eventsj$receiver == receiveri] # past events that the receiver (receiver) is involved in
          if(collapse::fnrow(targetallevents) == 0){next}
          senderallevents <- eventsj[eventsj$sender == senderi] # past events that the sender (sender) is involved in
          if(collapse::fnrow(senderallevents) == 0){next}
          #targetallevents <- targetallevents[time < timei, weightbj := remExpWeights(current = timei, past = time, halflife = halflife, dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013)] # compute the weight of prior events
          senderallevents <- senderallevents[senderallevents$time < timei]
          if(collapse::fnrow( senderallevents) == 0){next}
          senderallevents$weightia <- remExpWeights(current = timei,
                                                    past = senderallevents$time,
                                                    halflife = halflife,
                                                    dyadic_weight = dyadic_weight,
                                                    Lerneretal_2013 = Lerneretal_2013)
          senderallevents <- senderallevents[senderallevents$weightia != 0] # keep if only weight is not equal to 0 (this is important as we are removing events before computation)
          if(collapse::fnrow(senderallevents) == 0){next}
          targetallevents <- targetallevents[targetallevents$time < timei]
          if(collapse::fnrow(targetallevents) == 0){next}
          # targetallevents <- targetallevents[time < timei, weightbj := remExpWeights(current = timei, past = time, halflife = halflife, dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013)] # compute the weight of prior events
          targetallevents$weightbj <- remExpWeights(current = timei,
                                                    past = targetallevents$time,
                                                    halflife = halflife,
                                                    dyadic_weight = dyadic_weight,
                                                    Lerneretal_2013 = Lerneretal_2013)

          targetallevents <- targetallevents[targetallevents$weightbj != 0] # keep if only weight is not equal to 0 (this is important as we are removing events before computation)
          if(collapse::fnrow(targetallevents) == 0){next}
          #### getting only the past events, and then computing the event weight
          ### if event weight is not equal to 0, then keep them! (this lessens the overall computation time by removing potential event orders that would sum to 0)
          #### getting only the past events, and then computing the event weight
          ### if event weight is not equal to 0, then keep them! (this lessens the overall computation time by removing potential event orders that would sum to 0)
          if (collapse::fnrow(senderallevents) == 0 | collapse::fnrow(targetallevents) == 0) {
            next
          } # if there are no events with anyone involved (i.e, the same events)
          ## getting only unique cases
          pastsenderloop <- collapse::funique(targetallevents$sender)
          ## getting only unique cases
          pasttargetloop <- collapse::funique(senderallevents$receiver)
          ##### removing the current sender from the unique list
          inloop <- collapse::whichv(pastsenderloop, senderi)
          if (length(inloop) != 0) {
            pastsenderloop <- pastsenderloop[-(inloop)] # remove sender
          }
          ##### removing the current target from the unique list
          inloop <- collapse::whichv(pasttargetloop, receiveri)
          if (length(inloop) != 0) {
            pasttargetloop <- pasttargetloop[-(inloop)] # remove target
          }
          if (length(pastsenderloop) == 0 | length(pasttargetloop) == 0) {
            next
          } # if there are no events with anyone involved (i.e, the same events)
          # data.table::setkey(events, receiver) # using the data.table set key function [no need to set key anymore]
          targets_with_past_senders <- eventsj[eventsj$receiver %in% pasttargetloop]#DT original code: events[J(pasttargetloop)]
          if(collapse::fnrow( targets_with_past_senders) == 0){next}
          targets_with_past_senders <- targets_with_past_senders[targets_with_past_senders$time < timei] # give me only the prior time events
          if(collapse::fnrow( targets_with_past_senders) == 0){next}
          #targets_with_past_senders[, weight := remExpWeights(current = timei, past = time, halflife = halflife, dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013)] # compute the weight of all events
          targets_with_past_senders$weight <- remExpWeights(current = timei,
                                                            past = targets_with_past_senders$time,
                                                            halflife = halflife,
                                                            dyadic_weight = dyadic_weight,
                                                            Lerneretal_2013 = Lerneretal_2013)
          #targets_with_past_senders[, weight := remExpWeights(current = timei, past = time, halflife = halflife, dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013)] # compute the weight of all events
          targets_with_past_senders <- targets_with_past_senders[targets_with_past_senders$weight != 0] # keep if only weight is not equal to 0 (this is important as we are removing events before computation)
          if(collapse::fnrow( targets_with_past_senders) == 0){next}
          # all events that have the unique targets within for the creation of four cycles
          ### if event weight is not equal to 0, then keep them! (this lessens the overall computation time by removing potential event orders that would sum to 0)
          senders_with_past_targets <- eventsj[eventsj$sender %in% pastsenderloop] # all events that have the unique targets within for the creation of four cycles
          if(collapse::fnrow( senders_with_past_targets ) == 0){next}
          senders_with_past_targets <- senders_with_past_targets[senders_with_past_targets$time < timei] # give me only the prior time events
          if(collapse::fnrow( senders_with_past_targets ) == 0){next}
          #senders_with_past_targets <- #senders_with_past_targets[, weight := remExpWeights(current = timei, past = time, halflife = halflife, dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013)] # compute the weight of all events
          senders_with_past_targets$weight <-  remExpWeights(current = timei,
                                                             past = senders_with_past_targets$time,
                                                             halflife = halflife,
                                                             dyadic_weight = dyadic_weight,
                                                             Lerneretal_2013 = Lerneretal_2013)
          senders_with_past_targets <- senders_with_past_targets[senders_with_past_targets$weight != 0] # keep if only weight is not equal to 0 (this is important as we are removing events before computation)
          if(collapse::fnrow( senders_with_past_targets ) == 0){next}
          ### if event weight is not equal to 0, then keep them! (this lessens the overall computation time by removing potential event orders that would sum to 0)
          # if there are no events with anyone involved (i.e, the same events)
          #### combining the two data sets to give me all prior events with senders that share the same receiver as the current and receivers that share the same sender as past receivers
          open_events <- senders_with_past_targets[targets_with_past_senders,
                                                   on = colnames(senders_with_past_targets)[-which(colnames(senders_with_past_targets) == "weight")],
                                                   nomatch = 0
          ]
          if(collapse::fnrow( open_events ) == 0){next}
          #### getting only the past events, and then computing the event weight
          #open_events <- open_events[time < timei, weightij := remExpWeights(current = timei, past = time, halflife = halflife, dyadic_weight = dyadic_weight, Lerneretal_2013 = Lerneretal_2013)]
          open_events <- open_events[open_events$time < timei]
          if(collapse::fnrow( open_events ) == 0){next}
          open_events$weightij <-  remExpWeights(current = timei,
                                                 past = open_events$time,
                                                 halflife = halflife,
                                                 dyadic_weight = dyadic_weight,
                                                 Lerneretal_2013 = Lerneretal_2013)
          open_events <- open_events[open_events$weightij != 0] # keep if only weight is not equal to 0 (this is important as we are removing events before computation)
          ### if event weight is not equal to 0, then keep them! (this lessens the overall computation time by removing potential event orders that would sum to 0)
          if (collapse::fnrow(open_events) == 0) {
            next
          } # if there are no events in which the current event is being closed (i.e, the same events)
          ### merging the event cases from w(i,j) and w(i,A)
          if(counts == TRUE){
            counti <- 0
            pairs <- unique(open_events,by = c("sender","receiver")) #getting the unique pairs for the summation count
            for(z in 1:nrow(pairs)){
              counti <- counti + base::min(sum(senderallevents$receiver == pairs$receiver[z]), #d(s,r')
                                           sum(targetallevents$sender == pairs$sender[z]),#d(r,s')
                                           sum(open_events$sender == pairs$sender[z] &
                                                 open_events$receiver  == pairs$receiver[z]))#d(s',r')
            }
            fourcycle[i]  <- counti
            ###### clearing memory for consumption purposes
            rm(list = c(
              "open_events", "senders_with_past_targets", "targets_with_past_senders",
              "targetallevents", "pastsenderloop", "pasttargetloop", "senderallevents"
            ))
          }else{

            ### merging the event cases from w(i,j) and w(i,A)
            d31 <- data.table::merge.data.table(open_events, senderallevents, by = "receiver", all.x = T, all.y = T, allow.cartesian = T)
            if(collapse::fnrow( d31) == 0){next}
            d31t <- d31[stats::complete.cases(d31)] # complete cases, removing the cases that dont fully merge (i.e., are missing one of the weights)
            if(collapse::fnrow( d31t) == 0){next}
            ### merging the event cases from w(i,j) and w(i,A) to now include w(U, j)
            d31t3 <- data.table::merge.data.table(d31t, targetallevents, by.x = "sender.x", by.y = "sender", all.x = T, all.y = T, allow.cartesian = T)
            if(collapse::fnrow( d31t3) == 0){next}
            d31t4 <- d31t3[stats::complete.cases(d31t3)] # complete cases, removing the cases that dont fully merge (i.e., are missing one of the weights)
            if(collapse::fnrow(d31t4) == 0){next}
            d31t5 <- unique(d31t4) # getting only non-redudant unique cases
            if(collapse::fnrow( d31t5) == 0){next}
            # d31t5[, product := with(d31t5, weightij * weightia * weightbj)] # now computing the product of all three weights, but removing the instances whereby events are recycled from the past 1 million event dataset
            #fourcycle[i] <- sum(d31t5$product)^(1 / 3) # sum all values, take the cubic root, and then finally, add that value to the computed vector ith spot
            fourcycle[i] <- sum(d31t5$weightij * d31t5$weightia * d31t5$weightbj )^(1 / 3) # sum all values, take the cubic root, and then finally, add that value to the computed vector ith spot
            ###### clearing memory for consumption purposes
            rm(list = c(
              "open_events", "senders_with_past_targets", "targets_with_past_senders",
              "targetallevents", "pastsenderloop", "pasttargetloop", "senderallevents",
              "d31", "d31t", "d31t3", "d31t4", "d31t5"
            ))
          }                ###### clearing memory for consumption purposes

        }






















        }
      }






    }

  ########################################################
  #
  #   Returning the values back to the user
  #
  ########################################################
  return(fourcycle)
}
