## The Creation of Indegree Statistic Scores for Relational Event Models
## See Butts 2008; Lerner and Lomi 2020; REM R Package (Brandenberger 2018)
## Code written by Kevin Carson (kacarson@arizona.edu) and Deigo Leal (https://www.diegoleal.info/)
## Last Updated: 07-15-24.
#' @title Compute Degree Network Statistics for Event Senders and Receivers in a Relational Event Sequence
#' @name remstats_degree
#' @param formation The degree statistic to be computed. "sender-indegree" computes the indegree statistic for the event senders. "receiver-indegree" computes the
#' indegree statistic for the event receivers. "sender-outdegree" computes the outdegree statistic for the event senders. "receiver-outdegree" computes the
#' outdegree statistic for the event receivers.
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
#' @return The vector of degree statistics for the relational event sequence.
#' @export


#' @description
#' `r lifecycle::badge("stable")`
#'
#' The function computes the indegree network sufficient statistic for event senders
#' in a relational event sequence (see Lerner and Lomi 2020; Butts 2008).
#' This measure allows for indegree scores to be only  computed for the sampled
#' events, while creating the weights based on the full event sequence (see
#' Lerner and Lomi 2020; Vu et al. 2015). The function also allows users to use two
#' different weighting functions, return the counts of past events, reduce computational
#' runtime, and specify a dyadic cutoff for relational relevancy.
#'
#'
#'
#'@details The function calculates sender indegree scores for relational event
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
#'
#'**Sender-Indegree Statistic**:
#'
#'The formula for sender indegree for event \eqn{e_i} is:
#'\deqn{sender indegree_{e_{i}} = w(s', s, t) }
#'
#'That is, all past events in which the event receiver is the current sender. Following Butts (2008), if the
#'counts of the past events are requested, the formula for sender indegree for
#'event \eqn{e_i} is:
#'\deqn{sender indegree_{e_{i}} = d(r' = s, t') }
#'Where, \eqn{d()} is the number of past events where the event receiver, *r'*, is the current event sender *s* .
#'
#'
#'**Sender-Outdegree Statistic**:
#'
#'The formula for sender outdegree for event \eqn{e_i} is:
#'\deqn{sender outdegree_{e_{i}} = w(s, r', t) }
#'
#'That is, all past events in which the past sender is the current sender and
#'the event target can be any past user. Following Butts (2008), if the counts
#'of the past events are requested, the formula for sender outdegree for
#'event \eqn{e_i} is:
#'\deqn{sender outdegree_{e_{i}} = d(s = s', t') }
#'Where, \eqn{d()} is the number of past events where the sender *s'* is the current event sender, *s*
#'
#'
#'**Receiver-Outdegree Statistic**:
#'
#'The formula for receiver outdegree for event \eqn{e_i} is:
#'\deqn{receiver outdegree_{e_{i}} = w(r', r, t) }
#'
#'Following Butts (2008), if the counts of the past events are requested, the formula for receiver outdegree for
#'event \eqn{e_i} is:
#'\deqn{receiver outdegree{e_{i}} = d(s' = r, t') }
#'Where, \eqn{d()} is the number of past events where the event sender, *s'*, is the current event receiver, *r'*.
#'
#'
#'**Receiver-Indegree Statistic**:
#'
#'The formula for receiver indegree for event \eqn{e_i} is:
#' \deqn{reciever indegree_{e_{i}} = w(s', r, t) }
#'
#'That is, all past events in which the event receiver is the current receiver.
#'Following Butts (2008), if the counts of the past events are requested, the formula for receiver indegree for
#'event \eqn{e_i} is:
#'\deqn{reciever indegree_{e_{i}} = d(r' = r, t') }
#'where, \eqn{d()} is the number of past events where the past event receiver, *r'*, is the
#'current event receiver (target).
#'
#'
#'Lastly, researchers interested in modeling temporal relevancy (see Quintane,
#'Mood, Dunn, and Falzone 2022; Lerner and Lomi 2020) can specify the dyadic
#'weight cutoff, that is, the minimum value for which the weight is considered
#'relationally relevant. Users who do not know the specific dyadic cutoff value to use, can use the
#'\code{\link{remstats_dyadcut}} function.
#'

#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#' @references
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
#'events <- data.frame(time = 1:18, eventID = 1:18,
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
#'eventSet <- create_riskset(type = "one-mode",
#'                       time = events$time,
#'                       eventID = events$eventID,
#'                       sender = events$sender,
#'                       receiver = events$target,
#'                       p_samplingobserved = 1.00,
#'                       n_controls = 1,
#'                       seed = 9999)
#'
#'#Computing the sender indegree statistic for the relational event sequence
#'eventSet$senderind <- remstats_degree(
#'    formation = "sender-indegree",
#'    time = as.numeric(eventSet$time),
#'    observed = eventSet$observed,
#'    sampled = rep(1,nrow(eventSet)),
#'    sender = eventSet$sender,
#'    receiver = eventSet$receiver,
#'    halflife = 2, #halflife parameter
#'    dyadic_weight = 0,
#'    exp_weight_form = FALSE)
#'
#'#Computing the sender outdegree statistic for the relational event sequence
#'eventSet$senderout <- remstats_degree(
#'    formation = "sender-outdegree",
#'    time = as.numeric(eventSet$time),
#'    observed = eventSet$observed,
#'    sampled = rep(1,nrow(eventSet)),
#'    sender = eventSet$sender,
#'    receiver = eventSet$receiver,
#'    halflife = 2, #halflife parameter
#'    dyadic_weight = 0,
#'    exp_weight_form = FALSE)
#'
#'
#'#Computing the receiver outdegree statistic for the relational event sequence
#'eventSet$recieverout <- remstats_degree(
#'    formation = "receiver-outdegree",
#'    time = as.numeric(eventSet$time),
#'    observed = eventSet$observed,
#'    sampled = rep(1,nrow(eventSet)),
#'    sender = eventSet$sender,
#'    receiver = eventSet$receiver,
#'    halflife = 2, #halflife parameter
#'    dyadic_weight = 0,
#'    exp_weight_form = FALSE)
#'
#'
#'#Computing the receiver indegree statistic for the relational event sequence
#'eventSet$recieverind <- remstats_degree(
#'    formation = "receiver-indegree",
#'    time = as.numeric(eventSet$time),
#'    observed = eventSet$observed,
#'    sampled = rep(1,nrow(eventSet)),
#'    sender = eventSet$sender,
#'    receiver = eventSet$receiver,
#'    halflife = 2, #halflife parameter
#'    dyadic_weight = 0,
#'    exp_weight_form = FALSE)

########################################################################################################
#  Events = the full event sequence
#  eventSet = the sampled event sequence
#  Time = the name for the time variable
#  sender = the name for the sender variable
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

remstats_degree <- function(formation = c("sender-indegree", "receiver-indegree",
                                      "sender-outdegree", "receiver-outdegree"),
                                  time,# variable (column) name that contains the time variable
                                  sender,# variable (column) name that contains the sender variable
                                  receiver,# variable (column) name that contains the target variable
                                  observed,# variable (column) name that contains the observed variable
                                  sampled,# variable (column) name that contains the sampled variable
                                  halflife=2, # the half life value for the weighting function
                                  counts = FALSE, #Logical indicating if the raw counts of events should be returned or the exponential weighting function should be used (TRUE = counts; FALSE = exponential weighting)
                                  dyadic_weight= 0.00, # dyadic cutoff weight for events that no longer matter
                                  exp_weight_form = FALSE # Do we want to use the weighting function of Lerner et al. 2013 (alsoused in the rem R package)?
) {

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
    base::stop("Error: There are no sampled events based upon the 'sampled' input. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(typeof(time) != "double"){
    base::stop("Error: The 'time' argument is not a numeric value. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(length(formation) != 1){
    base::stop("Error: The 'type' argument is not of length 1. Please only input one type at a time! Happy computing!") # stop computation and tell the user
  }
  if(!(formation %in% c("sender-indegree", "receiver-indegree","sender-outdegree", "receiver-outdegree"))){
    base::stop("Error: The 'type' argument was not correctly specific. Please see the arguments section for more details. Happy computing!") # stop computation and tell the user
  }

  ########################################################
  #
  #   Computing the weights in c++
  #
  ########################################################
  Lerneretal_2013 <- exp_weight_form
  if(formation == "sender-indegree"){ #the degree formation for event sender indegree
    #computing the weights
    ########################################################
    #
    #   Prepping the data to be sent to c++ for speedy computation
    #
    ########################################################
    appender <- "__NIKOACAR2020__" # a (hopefully) unique joiner for the string!
    weightSchemeR <- ifelse(Lerneretal_2013, 0, 1) #setting this argument up for c++ computation
    countsR <- ifelse(counts, 1, 0) #setting this argument up for c++ computation
    controlR <- 1 - observed #making it such that dummy events have a 1 and real events have a value of 0
    ########################################################
    #
    #   Computing the weights in c++
    #
    ########################################################
    weights <- computeremweightsv2(time = time,
                                   sampledevent = sampled,
                                   controlevent = controlR,
                                   cutweight = dyadic_weight,
                                   halflife = halflife,
                                   dyad_id = (base::paste0(sender)),
                                   dyad_idOpposite = (base::paste0(receiver)) ,
                                   weightScheme = weightSchemeR,
                                   counts = countsR) #if the value is false then
  } #ending the if statement
  if(formation == "receiver-indegree"){#the degree formation for event target indegree
    ########################################################
    #
    #   Prepping the data to be sent to c++ for speedy computation
    #
    ########################################################
    appender <- "__NIKOACAR2020__" # a (hopefully) unique joiner for the string!
    dyad.idR <- (base::paste0(receiver)) #this is arguably very inefficent at scale
    weightSchemeR <- ifelse(Lerneretal_2013, 0, 1) #setting this argument up for c++ computation
    countsR <- ifelse(counts, 1, 0) #setting this argument up for c++ computation
    controlR <- 1 - observed #making it such that dummy events have a 1 and real events have a value of 0
    ########################################################
    #
    #   Computing the weights in c++
    #
    ########################################################
    weights <- computeREMweightsv1(time = time,
                                   sampledevent = sampled,
                                   controlevent = controlR,
                                   cutweight = dyadic_weight,
                                   halflife = halflife,
                                   dyad_id = dyad.idR,
                                   weightScheme = weightSchemeR,
                                   counts = countsR) #if the value is false then
  } #ending the if statement
  if(formation == "sender-outdegree"){#the degree formation for event sender outdegree
    ########################################################
    #
    #   Prepping the data to be sent to c++ for speedy computation
    #
    ########################################################
    appender <- "__NIKOACAR2020__" # a (hopefully) unique joiner for the string!
    dyad.idR <- (base::paste0(sender)) #this is arguably very inefficent at scale
    weightSchemeR <- ifelse(Lerneretal_2013, 0, 1) #setting this argument up for c++ computation
    countsR <- ifelse(counts, 1, 0) #setting this argument up for c++ computation
    controlR <- 1 - observed #making it such that dummy events have a 1 and real events have a value of 0
    ########################################################
    #
    #   Computing the weights in c++
    #
    ########################################################
    weights <- computeREMweightsv1(time = time,
                                   sampledevent = sampled,
                                   controlevent = controlR,
                                   cutweight = dyadic_weight,
                                   halflife = halflife,
                                   dyad_id = dyad.idR,
                                   weightScheme = weightSchemeR,
                                   counts = countsR) #if the value is false then
  } #ending the if statement
  if(formation == "receiver-outdegree"){#the degree formation for event target outdegree
    ########################################################
    #
    #   Prepping the data to be sent to c++ for speedy computation
    #
    ########################################################
    appender <- "__NIKOACAR2020__" # a (hopefully) unique joiner for the string!
    weightSchemeR <- ifelse(Lerneretal_2013, 0, 1) #setting this argument up for c++ computation
    countsR <- ifelse(counts, 1, 0) #setting this argument up for c++ computation
    controlR <- 1 - observed #making it such that dummy events have a 1 and real events have a value of 0
    ########################################################
    #
    #   Computing the weights in c++
    #
    ########################################################
    weights <- computeremweightsv2(time = time,
                                   sampledevent = sampled,
                                   controlevent = controlR,
                                   cutweight = dyadic_weight,
                                   halflife = halflife,
                                   dyad_id = (base::paste0(receiver)),
                                   dyad_idOpposite = (base::paste0(sender)),
                                   weightScheme = weightSchemeR,
                                   counts = countsR) #if the value is false then
  } #ending the if statement

  ########################################################
  #
  #   Returning the values back to the user
  #
  ########################################################
  return(weights)# return the vector of values
}
