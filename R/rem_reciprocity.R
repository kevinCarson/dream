## The Computation of Reciprocity Network Statistic for Relational Event Sequences
## Code written by Kevin Carson (https://kevincarson.github.io/) and Deigo Leal (https://www.diegoleal.info/)
## Last Updated: 01-09-2025



#' @title Compute the Reciprocity Network Statistic for Event Dyads in a Relational Event Sequence
#' @name remstats_reciprocity
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
#' @return The vector of reciprocity statistics for the relational event sequence.
#' @export
#'
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#'This function calculates the reciprocity network sufficient statistic for a
#'relational event sequence (see Lerner and Lomi 2020; Butts 2008). The reciprocity
#'statistic captures the tendency for a sender a to ‘send a tie’ to (e.g., initiate
#'a communication event with) receiver b given that b sent a tie to a in the
#'past (i.e., an exchange between two medical doctors). This function allows
#'for reciprocity scores to be only computed for the sampled events, while
#'creating the weights based on the full event sequence (see Lerner and
#'Lomi 2020; Vu et al. 2015). The function also allows users to use two
#'different weighting functions, return the counts of past events, reduce computational runtime, and specify
#'a dyadic cutoff for relational relevancy.
#'
#'
#'@details This function calculates reciprocity scores for relational event models
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
#'past event times that meet the weight subset, and \eqn{T_{1/2}} is the halflife parameter.
#'
#'The formula for reciprocity for event \eqn{e_i} is:
#'\deqn{reciprocity_{e_{i}} = w(r, s, t) }
#'
#'That is, all past events in which the past sender is the current receiver and
#'the past receiver is the current sender.
#'
#'Moreover, researchers interested in modeling temporal relevancy (see Quintane,
#'Mood, Dunn, and Falzone 2022; Lerner and Lomi 2020) can specify the dyadic
#'weight cutoff, that is, the minimum value for which the weight is considered
#'relationally relevant. Users who do not know the specific dyadic cutoff value to use, can use the
#'\code{\link{remstats_dyadcut}} function.
#'
#'Following Butts (2008), if the counts of the past events are requested, the formula for reciprocity for
#'event \eqn{e_i} is:
#'\deqn{reciprocity_{e_{i}} = d(r = s', s = r', t') }
#'Where, \eqn{d()} is the number of past events where the event sender, *s'*, is the current event receiver, *r*, and the event
#'receiver (target), *r'*, is the current event sender, *s*. Moreover, the counting equation
#'can be used in tandem with relational relevancy, by specifying the halflife parameter, exponential
#'weighting function, and the dyadic cut off weight values. If the user is not interested in modeling
#'relational relevancy, then those value should be left at their baseline values.


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
#'eventSet <-create_riskset(type = "one-mode",
#'                       time = events$time,
#'                       eventID = events$eventID,
#'                       sender = events$sender,
#'                       receiver = events$target,
#'                       p_samplingobserved = 1.00,
#'                       n_controls = 1,
#'                       seed = 9999)
#'
#'#Computing the reciprocity statistics for the relational event sequence
#'eventSet$recip <- remstats_reciprocity(
#'    time = as.numeric(eventSet$time),
#'    observed = eventSet$observed,
#'    sampled = rep(1,nrow(eventSet)),
#'    sender = eventSet$sender,
#'    receiver = eventSet$receiver,
#'    halflife = 2, #halflife parameter
#'    dyadic_weight = 0,
#'    exp_weight_form = FALSE)
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
########################################################################################################

remstats_reciprocity <- function( time,# variable (column) name that contains the time variable
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
  Lerneretal_2013 <- exp_weight_form
  appender <- "__NIKOACAR2020__" # a (hopefully) unique joiner for the string!
  dyad.idR <- (base::paste0(sender,appender,receiver)) #this is arguably very inefficent at scale
  dyad.idR.oops <- (base::paste0(receiver,appender,sender)) #this is arguably very inefficent at scale
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
                                 dyad_id = dyad.idR.oops,
                                 dyad_idOpposite = dyad.idR,
                                 weightScheme = weightSchemeR,
                                 counts = countsR) #if the value is false then
  ########################################################
  #
  #   Returning the values back to the user
  #
  ########################################################
  return(weights)# return the vector of values
}
