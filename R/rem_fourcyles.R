## The Creation of Four-Cycle Statistics for Large Relational Events into R
## Purpose: Compute Four Cycle Statistics
## See Butts 2008; Lerner and Lomi 2020; REM R Package (Brandenberger 2018)
## Code written by Kevin Carson (kacarson@arizona.edu) and Deigo Leal (https://www.diegoleal.info/)
## Last Updated: 07-31-24
#' @title Compute the Four-Cycles Network Statistic for Event Dyads in a Relational Event Sequence
#' @name  remstats_fourcycles
#' @param time The vector of event times from the post-processing event sequence.
#' @param sender The vector of event senders from the post-processing event sequence.
#' @param receiver The vector of event receivers from the post-processing event sequence
#' @param observed A vector for the post-processing event sequence where i is equal to 1 if the dyadic event is observed and 0 if not.
#' @param sampled A vector for the post-processing event sequence where i is equal to 1 if the observed dyadic event is sampled and 0 if not.
#' @param counts TRUE/FALSE. TRUE indicates that the counts of past events should be computed (see the details section). FALSE indicates that the temporal exponential weighting function should be used to downweigh past events (see the details section). Set to FALSE by default.
#' @param halflife A numerical value that is the halflife value to be used in the exponential weighting function (see details section). Preset to 2 (should be updated by the user based on substantive context).
#' @param dyadic_weight A numerical value for the dyadic cutoff weight that represents the numerical cutoff value for temporal relevancy based on the exponential weighting function. For example, a numerical value of 0.01, indicates that an exponential weight less than 0.01 will become 0 and that events with such value (or smaller values) will not be included in the sum of the past event weights (see the details section). Set to 0 by default.
#' @param exp_weight_form TRUE/FALSE. TRUE indicates that the Lerner et al. (2013) exponential weighting function will be used (see the details section). FALSE indicates that the Lerner and Lomi (2020) exponential weighting function will be used (see the details section). Set to FALSE by default
#' @import Rcpp
#' @return The vector of four-cycle statistics for the two-mode relational event sequence.
#' @export
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' The function computes the four-cycles network sufficient statistic for a two-mode relational
#' sequence with the exponential weighting function (Lerner and Lomi 2020). In essence, the
#' four-cycles measure captures the tendency for clustering to occur in the network of past
#' events, whereby an event is more likely to occur between a sender node *a* and receiver
#' node *b* given that *a* has interacted with other receivers in past events who have
#' received events from other senders that interacted with *b* (e.g., Duxbury and Haynie 2021, Lerner and Lomi 2020). The function
#' also allows users to use two different weighting functions, return the counts of past events, reduce
#' computational runtime, and specify a dyadic cutoff for relational relevancy.
#'
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
#'\code{\link{remstats_dyadcut}} function.
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
#'EventSet <-create_riskset(type = "two-mode",
#'  time = WikiEvent2018$time, # The Time Variable
#'  eventID = WikiEvent2018$eventID, # The Event Sequence Variable
#'  sender = WikiEvent2018$user, # The Sender Variable
#'  receiver = WikiEvent2018$article, # The Receiver Variable
#'  p_samplingobserved = 0.01, # The Probability of Selection
#'  n_controls = 8, # The Number of Controls to Sample from the Full Risk Set
#'  combine = TRUE,
#'  seed = 9999) # The Seed for Replication
#'
#'#Computing the four-cycles statistics for the relational event sequence with
#'#the exponential weights of past events returned
#'cycle4_weights <- remstats_fourcycles(
#'    time = EventSet$time,
#'    sender = EventSet$sender,
#'    receiver = EventSet$receiver,
#'    sampled = EventSet$sampled,
#'    observed = EventSet$observed,
#'    halflife = 2.592e+09, #halflife parameter
#'    dyadic_weight = 0,
#'    exp_weight_form = FALSE)
#'
#'
#'
#'#Computing the four-cycles statistics for the relational event sequence with
#'#the counts of past events returned
#'cycle4_counts <- remstats_fourcycles(
#'    time = EventSet$time,
#'    sender = EventSet$sender,
#'    receiver = EventSet$receiver,
#'    sampled = EventSet$sampled,
#'    observed = EventSet$observed,
#'    halflife = 2.592e+09, #halflife parameter
#'    dyadic_weight = 0,
#'    counts = TRUE)
#'
#'cbind(cycle4_weights, cycle4_counts)
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

remstats_fourcycles <- function( time,# variable (column) name that contains the time variable
                               sender,# variable (column) name that contains the sender variable
                               receiver,# variable (column) name that contains the target variable
                               observed,# variable (column) name that contains the observed variable
                               sampled,# variable (column) name that contains the sampled variable
                               halflife=2, # the half life value for the weighting function
                               counts = FALSE, #Logical indicating if the raw counts of events should be returned or the exponential weighting function should be used (TRUE = counts; FALSE = exponential weighting)
                               dyadic_weight= 0.00, # dyadic cutoff weight for events that no longer matter
                               exp_weight_form = FALSE
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
  if(halflife < 0){
    base::stop("Error: Halflife values must be positive.") # stop computation and tell the user
  }
  if(sum(observed) == 0){
    base::stop("Error: There are no observed events based upon the 'observed' input. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(sum(sampled) == 0){
    base::stop("Error: There are no sampled events based upon the 'observed' input. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(typeof(time) != "double"){
    base::stop("Error: The 'time' argument is not a numeric value. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  ########################################################
  #
  #   Prepping the data to be sent to c++ for speedy computation
  #
  ########################################################
  appender <- "__NIKOACAR2020__" # a (hopefully) unique joiner for the string!
  dyad.idR <- (base::paste0(sender,appender,receiver)) #this is arguably very inefficent at scale
  Lerneretal_2013 <- exp_weight_form
  weightSchemeR <- ifelse(Lerneretal_2013, 0, 1) #setting this argument up for c++ computation
  countsR <- ifelse(counts, 1, 0) #setting this argument up for c++ computation
  controlR <- 1 - observed #making it such that dummy events have a 1 and real events have a value of 0
  weights <- computefourcyclesrem(time = time,
                                  sampledevent = sampled,
                                  controlevent = controlR,
                                  cutweight = dyadic_weight,
                                  halflife = halflife,
                                  sender = base::paste0(sender),
                                  target = base::paste0(receiver),
                                  dyad_id = dyad.idR,
                                  weightScheme = weightSchemeR,
                                  counts = countsR,
                                  delim = appender)
  ########################################################
  #
  #   Returning the values back to the user
  #
  ########################################################
  return(weights)
}
