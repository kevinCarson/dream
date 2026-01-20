## Code written by Kevin Carson (https://kevincarson.github.io/) and Deigo Leal (https://www.diegoleal.info/)
## Last Updated: 01-09-2025

#' @title Compute Butts' (2008) Triadic Formation Statistics for Relational Event Sequences
#' @name remstats_triads
#' @param formation The specific triadic formation the statistic will be based on (see details section). "ISP" = incoming shared partners. "OSP" = outgoing shared partners. "OTP" = outgoing two-paths. "ITP" = incoming two-paths.
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
#' @return The vector of triadic formation statistics for the relational event sequence.
#' @export

#' @description
#' `r lifecycle::badge("stable")`
#'
#' The function computes the set of one-mode triadic formation statistics discussed in Butts (2008) for a one-mode relational
#' event sequence (see also Lerner and Lomi 2020). The function can compute the following triadic formations: 1) incoming shared partners (ISP),
#' 2) outgoing shared partners (OSP), 3) incoming two-paths (ITP), and 4) outgoing two-paths (OTP). Importantly, this function allows for the triadic formation
#' statistics to be computed only for the sampled events, while creating the weights based on the full event sequence (see
#' Lerner and Lomi 2020; Vu et al. 2015). The function also allows users to use two different
#' weighting functions, return the counts of past events, reduce computational
#' runtime, and specify a dyadic cutoff for relational relevancy.
#'
#'
#'@details The function calculates the triadic formation statistics discussed in Butts (2008) for relational
#'event sequences based on the exponential weighting function used in either Lerner
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
#'**Outgoing Shared Partners**:
#'
#'The general formula for outgoing shared partners for event \eqn{e_i} is:
#'\deqn{OSP_{e_{i}} = \sqrt{ \sum_h w(s, h, t) \cdot w(r, h, t) }}
#'
#'That is, as discussed in Butts (2008), outgoing shared partners finds all
#'past events where the current sender and target sent a relational tie (i.e.,
#'were a sender in a relational event) to the same *h* node.
#'
#'Following Butts (2008), if the counts of the past events are requested, the formula for outgoing shared partners for
#'event \eqn{e_i} is:
#'\deqn{OSP{e_{i}} = \sum_{i=1}^{|H|} \min\left[d(s,h,t), d(s,h,t)\right]}
#'Where, \eqn{d()} is the number of past events that meet the specific set operations. \eqn{d(s,h,t)} is the number
#'of past events where the current event sender sent a tie to a third actor, *h*, and \eqn{d(r,h,t)} is the number
#'of past events where the current event receiver sent a tie to a third actor, *h*. The sum loops through all
#'unique actors that have formed past outgoing shared partners structures with the current event sender and receiver.
#'Moreover, the counting equation can be used in tandem with relational relevancy, by specifying the halflife parameter, exponential
#'weighting function, and the dyadic cut off weight values. If the user is not interested in modeling
#'relational relevancy, then those value should be left at their defaults.
#'
#'
#'**Outgoing Two-Paths**:
#'
#'The general formula for outgoing two-paths for event \eqn{e_i} is:
#'\deqn{OTP_{e_{i}} = \sqrt{ \sum_h w(s, h, t) \cdot w(h, r, t) }}
#'
#'That is, as discussed in Butts (2008), outgoing two-paths finds all
#'past events where the current sender sends a relational tie to node *h* and
#'the current target receives a relational tie from the same *h* node.
#'
#'Following Butts (2008), if the counts of the past events are requested, the formula for outgoing two paths for
#'event \eqn{e_i} is:
#'\deqn{OTP_{e_{i}} = \sum_{i=1}^{|H|} \min\left[d(s,h,t), d(h,r,t)\right]}
#'Where, \eqn{d()} is the number of past events that meet the specific set operations. \eqn{d(s,h,t)} is the number
#'of past events where the current event sender sent a tie to a third actor, *h*, and \eqn{d(h,r,t)} is the number
#'of past events where the third actor *h* sent a tie to the current event receiver. The sum loops through all
#'unique actors that have formed past outgoing two-path structures with the current event sender and receiver.
#'
#'
#'**Incoming Two-Paths**:
#'
#'The general formula for incoming two-paths for event \eqn{e_i} is:
#'\deqn{ITP_{e_{i}} = \sqrt{ \sum_h w(r, h, t) \cdot w(h, s, t) }}
#'
#'That is, as discussed in Butts (2008), incoming two-paths finds all past events
#'where the current sender was the receiver in a relational event where the sender
#'was a node h and the current target was the sender in a past relational event
#'where the target was the same node h.
#'
#'Following Butts (2008), if the counts of the past events are requested, the formula for incoming two paths for
#'event \eqn{e_i} is:
#'\deqn{ITP_{e_{i}} = \sum_{i=1}^{|H|} \min\left[d(r,h,t), d(h,s,t\right]}
#'Where, \eqn{d()} is the number of past events that meet the specific set operations. \eqn{d(r,h,t)} is the number
#'of past events where the current event receiver sent a tie to a third actor, *h*, and \eqn{d(h,s,t} is the number
#'of past events where the third actor *h* sent a tie to the current event sender. The sum loops through all
#'unique actors that have formed past incoming two-path structures with the current event sender and receiver.
#'
#'
#'**Incoming Shared Partners**:
#'
#'The general formula for incoming shared partners for event \eqn{e_i} is:
#'\deqn{ISP_{e_{i}} = \sqrt{ \sum_h w(h, s, t) \cdot w(h, r, t) }}
#'
#'That is, as discussed in Butts (2008), incoming shared partners finds all
#'past events where the current sender and target were themselves the target
#'in a relational event from the same *h* node.
#'
#'Following Butts (2008), if the counts of the past events are requested, the formula for incoming shared partners for
#'event \eqn{e_i} is:
#'\deqn{ISP_{e_{i}} = \sum_{i=1}^{|H|} \min\left[d(h,s,t), d(h,r,t)\right]}
#'Where, \eqn{d()} is the number of past events that meet the specific set operations, \eqn{d(h,s,t)} is the number
#'of past events where the current event sender received a tie from a third actor, *h*, and \eqn{d(h,r,t)} is the number
#'of past events where the current event receiver received a tie from a third actor, *h*. The sum loops through all
#'unique actors that have formed past incoming shared partners structures with the current event sender and receiver.
#'
#'
#'Lastly, researchers interested in modeling temporal relevancy (see Quintane,
#'Mood, Dunn, and Falzone 2022; Lerner and Lomi 2020) can specify the dyadic
#'weight cutoff, that is, the minimum value for which the weight is considered
#'relationally relevant. Users who do not know the specific dyadic cutoff value to use, can use the
#'\code{\link{remstats_dyadcut}} function.
#'
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
#'eventSet <-create_riskset(type = "one-mode",
#'                       time = events$time,
#'                       eventID = events$eventID,
#'                       sender = events$sender,
#'                       receiver = events$target,
#'                       p_samplingobserved = 1.00,
#'                       n_controls = 1,
#'                       seed = 9999)
#'
#'#compute the triadic statistic for the outgoing shared partners formation
#'eventSet$OSP <- remstats_triads(
#'    formation = "OSP", #outgoing shared partners argument
#'    time = as.numeric(eventSet$time),
#'    observed = eventSet$observed,
#'    sampled = rep(1,nrow(eventSet)),
#'    sender = eventSet$sender,
#'    receiver = eventSet$receiver,
#'    halflife = 2, #halflife parameter
#'    dyadic_weight = 0,
#'    exp_weight_form = FALSE)
#'
#'#compute the triadic statistic for the incoming shared partners formation
#'eventSet$ISP <- remstats_triads(
#'    formation = "ISP", #incoming shared partners argument
#'    time = as.numeric(eventSet$time),
#'    observed = eventSet$observed,
#'    sampled = rep(1,nrow(eventSet)),
#'    sender = eventSet$sender,
#'    receiver = eventSet$receiver,
#'    halflife = 2, #halflife parameter
#'    dyadic_weight = 0,
#'    exp_weight_form = FALSE)
#'
#'#compute the triadic statistic for the outgoing two-paths formation
#'eventSet$OTP <- remstats_triads(
#'    formation = "OTP", #outgoing two-paths argument
#'    time = as.numeric(eventSet$time),
#'    observed = eventSet$observed,
#'    sampled = rep(1,nrow(eventSet)),
#'    sender = eventSet$sender,
#'    receiver = eventSet$receiver,
#'    halflife = 2, #halflife parameter
#'    dyadic_weight = 0,
#'    exp_weight_form = FALSE)
#'
#'#compute the triadic statistic for the incoming two-paths formation
#'eventSet$ITP <- remstats_triads(
#'    formation = "ITP", #incoming two-paths argument
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

remstats_triads <- function(formation = c("ISP", "OSP", "ITP", "OTP"), #the type of traidic formations
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
    base::stop("Error: There are no sampled events based upon the 'observed' input. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(typeof(time) != "double"){
    base::stop("Error: The 'time' argument is not a numeric value. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(length(formation) != 1){
    base::stop("Error: The 'type' argument is not of length 1. Please only input one type at a time! Happy computing!") # stop computation and tell the user
  }
  if(!(formation %in% c("ISP", "OSP", "ITP", "OTP"))){
    base::stop("Error: The 'type' argument was not correctly specific. The input should be of one of four: 'ISP', 'OSP', 'ITP', or 'OTP'. Happy computing!") # stop computation and tell the user
  }

  ########################################################
  #
  #   Prepping the data to be sent to c++ for speedy computation
  #
  ########################################################
  Lerneretal_2013 <- exp_weight_form
  appender <- "__NIKOACAR2020__" # a (hopefully) unique joiner for the string!
  dyad.idR <- (base::paste0(sender,appender,receiver)) #this is arguably very inefficent at scale
  weightSchemeR <- ifelse(Lerneretal_2013, 0, 1) #setting this argument up for c++ computation
  countsR <- ifelse(counts, 1, 0) #setting this argument up for c++ computation
  controleventsR <- 1 - observed #making it such that dummy events have a 1 and real events have a value of 0

  ########################################################
  #
  #   Computing the weights in c++
  #
  ########################################################
  if(formation == "OSP"){ #if outgoing shared partners is requested!
        weights <- computeoutsharedpart(time = time,
                                        sampledevent = sampled,
                                        controlevent = controleventsR,
                                        sender = base::paste0(sender),
                                        target = base::paste0(receiver),
                                        dyad_id = dyad.idR,
                                        weightScheme = weightSchemeR,
                                        counts = countsR,
                                        cutweight = dyadic_weight,
                                        halflife = halflife,
                                        appender = appender)
  }

  if(formation == "ISP"){ #if outgoing shared partners is requested!
      weights <- computeincomingsharedparts(time = time,
                                            sampledevent = sampled,
                                            controlevent = controleventsR,
                                            sender = base::paste0(sender),
                                            target = base::paste0(receiver),
                                            dyad_id = dyad.idR,
                                            weightScheme = weightSchemeR,
                                            counts = countsR,
                                            cutweight = dyadic_weight,
                                            halflife = halflife,
                                            appender = appender)
  }

  if(formation == "ITP"){ #if outgoing shared partners is requested!
      weights <- computeincomingtwopaths(time = time,
                                         sampledevent = sampled,
                                         controlevent = controleventsR,
                                         sender = base::paste0(sender),
                                         target = base::paste0(receiver),
                                         dyad_id = dyad.idR,
                                         weightScheme = weightSchemeR,
                                         counts = countsR,
                                         cutweight = dyadic_weight,
                                         halflife = halflife,
                                         appender = appender)
  }

  if(formation == "OTP"){ #if outgoing shared partners is requested!
    weights <- computeouttwopaths(time = time,
                                  sampledevent = sampled,
                                  controlevent = controleventsR,
                                  sender = base::paste0(sender),
                                  target = base::paste0(receiver),
                                  dyad_id = dyad.idR,
                                  weightScheme = weightSchemeR,
                                  counts = countsR,
                                  cutweight = dyadic_weight,
                                  halflife = halflife,
                                  appender = appender)
  }


  ########################################################
  #
  #   Returning the values back to the user
  #
  ########################################################
  return(weights)# return the vector of values

}
