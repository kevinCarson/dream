# the master file for dream R package functions related to the creation
# of rem-related network statistics

#' @title Compute Butts' (2008) Triadic Formation Statistics for Relational Event Sequences
#' @name dreamstats_triads
#' @param formation The specific triadic formation the statistic will be based on (see details section). "ISP" = incoming shared partners. "OSP" = outgoing shared partners. "OTP" = outgoing two-paths. "ITP" = incoming two-paths.
#' @param data A `dream_sequence` object that contains the processed relational event sequence.
#' @param counts TRUE/FALSE. TRUE indicates that the counts of past events should be computed (see the details section). FALSE indicates that the temporal exponential weighting function should be used to downweigh past events (see the details section). Set to FALSE by default.
#' @param halflife A numerical value that is the halflife value to be used in the exponential weighting function (see details section). Preset to 2 (should be updated by the user based on substantive context).
#' @param dyadic_weight A numerical value for the dyadic cutoff weight that represents the numerical cutoff value for temporal relevancy based on the exponential weighting function. For example, a numerical value of 0.01, indicates that an exponential weight less than 0.01 will become 0 and that events with such value (or smaller values) will not be included in the sum of the past event weights (see the details section). Set to 0 by default.
#' @param exp_weight_form TRUE/FALSE. TRUE indicates that the Lerner et al. (2013) exponential weighting function will be used (see the details section). FALSE indicates that the Lerner and Lomi (2020) exponential weighting function will be used (see the details section). Set to FALSE by default
#' @param return_stats TRUE/FALSE. TRUE indicates that the vector of computed
#' statistics will be returned. FALSE indicates that the vector of computed
#' statistics will be added to the `statistics` element of the `data` argument.
#' @import Rcpp
#' @importFrom methods is
#' @return The vector of triadic formation statistics for the relational event sequence or the updated `data` argument.
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
#'\code{\link{dreamstats_dyadcut}} function.
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
#'eventSet <- create_res(type = "one-mode",
#'                       ordinal = TRUE,
#'                       riskset = "constant_sample",
#'                       time = events$time,
#'                       sender = events$sender,
#'                       receiver = events$target,
#'                       p_samplingobserved = 1.00,
#'                       n_controls = 1,
#'                       seed = 9999)
#'
#'#compute the triadic statistic for the outgoing shared partners formation
#'eventSet <- dreamstats_triads(formation = "OSP",
#'                              data = eventSet,
#'                              halflife = 2, #halflife parameter
#'                              dyadic_weight = 0)
#'#printing the post-processed relational event sequence
#'eventSet
#'#printing the vector of computed values
#'eventSet$statistics$outgoing.shared.partners
#'
#'
#'#Computing theoutgoing shared partners statistic for the relational event sequence
#'#and returning only the vector of computed values
#'osp.stat <- dreamstats_triads(formation = "OSP",
#'                              data = eventSet,
#'                              halflife = 2, #halflife parameter
#'                              dyadic_weight = 0,
#'                              return_stats = TRUE)
#'
#'cor(osp.stat, eventSet$statistics$outgoing.shared.partners)
#'
#'#compute the triadic statistic for the incoming shared partners formation
#'eventSet <- dreamstats_triads(
#'    formation = "ISP",
#'    data = eventSet,
#'    halflife = 2, #halflife parameter
#'    dyadic_weight = 0)
#'
#'#compute the triadic statistic for the outgoing two-paths formation
#'eventSet <- dreamstats_triads(
#'    formation = "OTP",
#'    data = eventSet,
#'    halflife = 2, #halflife parameter
#'    dyadic_weight = 0)
#'
#'#compute the triadic statistic for the incoming two-paths formation
#'eventSet <- dreamstats_triads(
#'    formation = "ITP",
#'    data = eventSet,
#'    halflife = 2, #halflife parameter
#'    dyadic_weight = 0)
#'
#'#extracting the relational event information
#'triad.rems <- as.data.frame(eventSet)
#'triad.rems

dreamstats_triads <- function(formation = c("ISP", "OSP", "ITP", "OTP"), #the type of traidic formations
                              data,
                              halflife=2, # the half life value for the weighting function
                              counts = FALSE, #Logical indicating if the raw counts of events should be returned or the exponential weighting function should be used (TRUE = counts; FALSE = exponential weighting)
                              dyadic_weight= 0.00, # dyadic cutoff weight for events that no longer matter
                              exp_weight_form = FALSE, # Do we want to use the weighting function of Lerner et al. 2013 (alsoused in the rem R package)?
                              return_stats = FALSE
) {

  ########################################################
  #
  #   Checking for Errors in User Inputs
  #
  ########################################################
  if(!inherits(data,"dream_sequence")) base::stop("Error: the object for the `data` argument is not a `dream_sequence` object.")
  time <- data$processed_sequence$time # variable (column) name that contains the time variable
  sender<- data$processed_sequence$sender# variable (column) name that contains the sender variable
  receiver<- data$processed_sequence$receiver# variable (column) name that contains the target variable
  observed<- data$processed_sequence$observed# variable (column) name that contains the observed variable
  sampled<- data$processed_sequence$sampled# variable (column) name that contains the sampled variable

  if(data$type != "one-mode") base::stop("Error: Currently, triadic statistics are only defined for one-mode networks.") # stop computation and tell the user

  if(halflife < 0){
    base::stop("Error: Halflife values must be positive.") # stop computation and tell the user
  }
  if(sum(observed) == 0){
    base::stop("Error: There are no observed events based upon the 'observed' input. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(sum(sampled) == 0){
    base::stop("Error: There are no sampled events based upon the 'observed' input. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(!is.numeric(time)){
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
    weights <- computetriadsrem(time = time,
                                sampledevent = sampled,
                                controlevent = controleventsR,
                                sender = base::paste0(sender),
                                target = base::paste0(receiver),
                                dyad_id = dyad.idR,
                                weightScheme = weightSchemeR,
                                counts = countsR,
                                cutweight = dyadic_weight,
                                halflife = halflife,
                                appender = appender,
                                tri_type = "outgoing.shared.partners")
  }

  if(formation == "ISP"){ #if outgoing shared partners is requested!
    weights <-computetriadsrem(time = time,
                               sampledevent = sampled,
                               controlevent = controleventsR,
                               sender = base::paste0(sender),
                               target = base::paste0(receiver),
                               dyad_id = dyad.idR,
                               weightScheme = weightSchemeR,
                               counts = countsR,
                               cutweight = dyadic_weight,
                               halflife = halflife,
                               appender = appender,
                               tri_type = "incoming.shared.partners")
  }

  if(formation == "ITP"){ #if outgoing shared partners is requested!
    weights <- computetriadsrem(time = time,
                                sampledevent = sampled,
                                controlevent = controleventsR,
                                sender = base::paste0(sender),
                                target = base::paste0(receiver),
                                dyad_id = dyad.idR,
                                weightScheme = weightSchemeR,
                                counts = countsR,
                                cutweight = dyadic_weight,
                                halflife = halflife,
                                appender = appender,
                                tri_type = "incoming.two.paths")
  }

  if(formation == "OTP"){ #if outgoing shared partners is requested!
    weights <- computetriadsrem(time = time,
                                sampledevent = sampled,
                                controlevent = controleventsR,
                                sender = base::paste0(sender),
                                target = base::paste0(receiver),
                                dyad_id = dyad.idR,
                                weightScheme = weightSchemeR,
                                counts = countsR,
                                cutweight = dyadic_weight,
                                halflife = halflife,
                                appender = appender,
                                tri_type = "outgoing.two.paths")
  }


  ########################################################
  #
  #   Returning the values back to the user as an attachment to the class
  #
  ########################################################
  if(return_stats){
    data <- weights #returning just the vector of computed stats
  }else{
    if(formation == "OTP"){  data$statistics$outgoing.two.paths=weights }
    if(formation == "ITP"){  data$statistics$incoming.two.paths=weights}
    if(formation == "ISP"){  data$statistics$incoming.shared.partners=weights}
    if(formation == "OSP"){  data$statistics$outgoing.shared.partners=weights}
  }
  return(data) # output the file to the user!

}






#' @title Add Event-Level Statistics for a Relational Event Sequence
#' @description
#' `r lifecycle::badge("stable")`
#'
#' This function allows users to add event-level statistics that impact the event rates in
#' interval timing relational event model, such as statistics that impact the
#' waiting times between events.
#'
#' @name dreamstats_event
#' @param data A `dream_sequence` object that contains the processed relational event sequence.
#' @param var_name A string that is the name of the variable from the `event_info` argument that represents the actor-level statistic to be added.
#' @param event_info A \eqn{N} x 2 `data.frame` object that contains two columns with the number of observations
#' being the number of observed time points actors. The object should contain
#' two named columns: (1) `time_id`, which is the time associated with each event and (2) `var_name` (from the argument), where
#' the var_name column is the vector of statistics (each time point has one associated value; see details for more information).
#' @param make_factor TRUE/FALSE. TRUE indicates that the vector of values will be made a factor, and FALSE if not.
#' @param return_stats TRUE/FALSE. TRUE indicates that the vector of computed
#' statistics will be returned. FALSE indicates that the vector of computed
#' statistics will be added to the `statistics` element of the `data` argument.
#' @import Rcpp
#' @return The vector of event-level statistics for the relational event sequence or the updated `data` argument.
#' @export
#'
#'@details This function adds user-provided time-varying and time-invariant actor-level
#' statistics to a `dream_sequence` object for relational event models. The `event_info`
#' argument should be a \eqn{N} x 2 `data.frame` object with two named columns. The first
#' column should be named `time_id`, which represents the observed time points
#' based upon the relational event sequence. The second column should be named after the
#' `var_name` argument and represents the event-level statistic for that i time point (i.e., the ith `time_id`).
#'
#'
#'@examples
#'events <- data.frame(time = 1:18, eventID = 1:18,
#'                     sender = c("A", "B", "C",
#'                                "A", "D", "E",
#'                                "F", "B", "A",
#'                                "F", "D", "B",
#'                                "G", "B", "D",
#'                                "H", "A", "D"),
#'                     target = c("B", "C", "D",
#'                                "E", "A", "F",
#'                               "D", "A", "C",
#'                                "G", "B", "C",
#'                                "H", "J", "A",
#'                                "F", "C", "B"))
#'
#'processed <- create_res(type = "one-mode",
#'                       ordinal = TRUE,
#'                       time = events$time,
#'                       riskset = "constant_sample",
#'                       sender = events$sender,
#'                       receiver = events$target,
#'                       p_samplingobserved = 1.00,
#'                       n_controls = 5,
#'                       seed = 9999)
#'
#'#reconstructing the data.frame object to store the time-varying
#'#event-level statistic
#'event_stats <- data.frame(time_id = events$time,
#'                          oliver = rnorm(nrow(events)))
#'
#'#reconstructing the data.frame object to store the time-varying
#'#event-level statistic
#'processed <- dreamstats_event(data = processed,
#'                              var_name = "oliver",
#'                              event_info = event_stats)
#'processed
#'processed$statistics$oliver.event
#'
#'extract.data <- as.data.frame(processed)
#'extract.data
#'
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#'
dreamstats_event <- function(data, #the relational event sequence
                             var_name, #the variable name
                             event_info,
                             make_factor = FALSE, #should the variable be made categorical (a factor in R)?
                             return_stats = FALSE
) {

  ########################################################
  #
  #   Checking for Errors in User Inputs
  #
  ########################################################
  if(!inherits(data,"dream_sequence")) base::stop("Error: the object for the `data` argument is not a `dream_sequence` object.")
  time <- data$processed_sequence$time # variable (column) name that contains the time variable
  sampled<- data$processed_sequence$sampled# variable (column) name that contains the sampled variable
  if(sum(sampled) == 0){
    base::stop("Error: There are no sampled events based upon the 'sampled' input. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(!is.numeric(time)){
    base::stop("Error: The 'time' argument is not a numeric value. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(!inherits(event_info,"data.frame")){
    base::stop("Error: The 'event_info' argument should be a `data.frame` object. Please update this argument and return.") # stop computation and tell the user
  }
  var.names <- colnames(event_info) #extracting the column names of the actor information dataset
  need.names <- c("time_id", var_name)
  #checking the inputs
  if(!all(need.names %in% var.names)){
    base::stop("Error: The `var_name` argument and/or a `time` column are missing in the `event_info` argument. Please update this argument and return.") # stop computation and tell the user
  }
  if(!is.numeric(event_info[,"time_id"])){
    base::stop("Error: The 'time_id' column in the `event_info` argument is not a numeric value. Please update this argument.") # stop computation and tell the user
  }
  #checking to make sure that for all sampled actors, we have information, if not... tell the user
  time <- time[sampled == 1]
  check.information <- all(time %in% event_info[,"time_id"])
  if(!check.information)   base::stop("Error: There are event times in the post-processing event sequence (`data` argument) that do not have observations in the `event_info` argument. Please update this argument and return.") # stop computation and tell the user
  #subsetting to only get those relevant actors information that are needed for the computation
  event_info <- base::subset(event_info, event_info[,"time_id"] %in% time) #only those sampled event times
  #inputting the associated time point values with the match function
  weights <- event_info[,var_name][match(time, event_info[,"time_id"], nomatch = NA)]
  if(make_factor) weights <- base::factor(weights)
  ########################################################
  #
  #   Returning the values back to the user as an attachment to the class
  #
  ########################################################
  if(return_stats){
    data <- weights #returning just the vector of computed stats
  }else{
    new.name <- paste0(var_name,".event", sep="")
    data$statistics[[new.name]] <-weights
  }
  return(data) # output the file to the user!s
}



#' @title Add Actor-Level Fixed Effects for Event Dyads in a Relational Event Sequence
#' @description
#' `r lifecycle::badge("stable")`
#'
#' This function allows users to add event sender and receiver fixed effects
#' for relational event models to a `dream_sequence` object.

#' @name dreamstats_actorfe
#' @param data A `dream_sequence` object that contains the processed relational event sequence.
#' @param sender TRUE/FALSE. TRUE creates event sender fixed effects and FALSE created event receiver fixed effects.
#' @param return_stats TRUE/FALSE. TRUE indicates that the vector of computed
#' statistics will be returned. FALSE indicates that the vector of computed
#' statistics will be added to the `statistics` element of the `data` argument.
#' @import Rcpp
#' @return The vector of actor-level fixed effects for the relational event sequence or the updated `data` argument.
#' @export
#'
#'@details This function adds sender or receiver actor-level fixed effects
#'to a `dream_sequence` object. Internally, the function creates a new variable
#'(`senderFE` for sender fixed effects and `receiverFE` for receiver fixed effects)
#'that is a factor of the event sender/receiver ids.
#'

#'@examples
#'events <- data.frame(time = 1:18, eventID = 1:18,
#'                     sender = c("A", "B", "C",
#'                                "A", "D", "E",
#'                                "F", "B", "A",
#'                                "F", "D", "B",
#'                                "G", "B", "D",
#'                                "H", "A", "D"),
#'                     target = c("B", "C", "D",
#'                                "E", "A", "F",
#'                               "D", "A", "C",
#'                                "G", "B", "C",
#'                                "H", "J", "A",
#'                                "F", "C", "B"))
#'
#'processed <- create_res(type = "one-mode",
#'                       ordinal = TRUE,
#'                       riskset = "constant_sample",
#'                       time = events$time,
#'                       sender = events$sender,
#'                       receiver = events$target,
#'                       p_samplingobserved = 1.00,
#'                       n_controls = 20,
#'                       seed = 9999)
#'
#'#adding the sender fixed effects to the statistics list
#'processed <- dreamstats_actorfe(data=processed,sender=TRUE)
#'
#'#adding the receiver fixed effects to the statistics list
#'processed <- dreamstats_actorfe(data=processed,sender=FALSE)
#'
#'processed #the effects are added
#'
#'#estimating the fixed effects only ordinal timing relational event model
#'model <- estimate_rem(~senderFE + receiverFE, data = processed)
#'
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#'
dreamstats_actorfe <- function(data, #the relational event sequence
                               sender = TRUE,
                               return_stats = FALSE) {

  ########################################################
  #
  #   Checking for Errors in User Inputs
  #
  ########################################################
  if(!inherits(data,"dream_sequence")) base::stop("Error: the object for the `data` argument is not a `dream_sequence` object.")
  sampled<- data$processed_sequence$sampled# variable (column) name that contains the sampled variable
  if(sum(sampled) == 0){
    base::stop("Error: There are no sampled events based upon the 'observed' input. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(sender){
    actor <- data$processed_sequence$sender[sampled == 1]
  }else{
    actor <- data$processed_sequence$receiver[sampled == 1]
  }
  weights <- base::factor(actor)
  ########################################################
  #
  #   Returning the values back to the user as an attachment to the class
  #
  ########################################################
  if(return_stats){
    data <- weights #returning just the vector of computed stats
  }else{
    if(sender) new.name <- "senderFE"
    if(!sender) new.name <- "receiverFE"
    data$statistics[[new.name]] <-weights
  }
  return(data) # output the file to the user!s
}








#' @title Add Dyadic-Level Fixed Effects for Event Dyads in a Relational Event Sequence
#' @description
#' `r lifecycle::badge("stable")`
#'
#' This function allows users to add dyad-level fixed effects
#' for relational event models to a `dream_sequence` object.

#' @name dreamstats_dyadfe
#' @param data A `dream_sequence` object that contains the processed relational event sequence.
#' @param directed TRUE/FALSE. TRUE indicates that the dyadic-level fixed effects
#' will be generated based upon the ordering of the sending and receiving actors (i.e.,
#' AB != BA). FALSE indicates that the dyadic-level fixed effects will be generated
#' based upon the combo of the sending and receiving actors (i.e., AB == BA). Set to
#' TRUE by default. Of course, these will lead to (1) a different number of fixed
#' effects and (2) a different interpretation of the results.
#' @param return_stats TRUE/FALSE. TRUE indicates that the vector of computed
#' statistics will be returned. FALSE indicates that the vector of computed
#' statistics will be added to the `statistics` element of the `data` argument.
#' @import Rcpp
#' @return The vector of actor-level fixed effects for the relational event sequence or the updated `data` argument.
#' @export
#'
#'@details This function adds dyad-level fixed effects to a `dream_sequence`
#' object. Internally, the function creates a new variable named "dyad.fe"
#' which is a factor of the combined event sender/receiver ids.
#'

#'@examples
#'events <- data.frame(time = 1:18, eventID = 1:18,
#'                     sender = c("A", "B", "C",
#'                                "A", "D", "E",
#'                                "F", "B", "A",
#'                                "F", "D", "B",
#'                                "G", "B", "D",
#'                                "H", "A", "D"),
#'                     target = c("B", "C", "D",
#'                                "E", "A", "F",
#'                               "D", "A", "C",
#'                                "G", "B", "C",
#'                                "H", "J", "A",
#'                                "F", "C", "B"))
#'
#'processed <- create_res(type = "one-mode",
#'                       ordinal = TRUE,
#'                       riskset = "complete",
#'                       time = events$time,
#'                       sender = events$sender,
#'                       receiver = events$target,
#'                       seed = 9999)
#'
#'#adding the dyadic fixed effects to the statistics list
#'processed <- dreamstats_dyadfe(data=processed)
#'
#'#estimating the dyadic fixed effects only ordinal timing relational event model
#'model <- estimate_rem(~dyad.fe, data = processed)
#'
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#'
dreamstats_dyadfe <- function(data, #the relational event sequence
                              directed=TRUE,
                              return_stats = FALSE) {

  ########################################################
  #
  #   Checking for Errors in User Inputs
  #
  ########################################################
  if(!inherits(data,"dream_sequence")) base::stop("Error: the object for the `data` argument is not a `dream_sequence` object.")
  sampled<- data$processed_sequence$sampled# variable (column) name that contains the sampled variable
  if(sum(sampled) == 0){
    base::stop("Error: There are no sampled events based upon the 'observed' input. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  sender <- data$processed_sequence$sender[sampled == 1]
  receiver <- data$processed_sequence$receiver[sampled == 1]
  if(directed)  dyad.ids <- paste0(sender, "_+_", receiver)
  if(!directed) dyad.ids <- paste0(pmin(sender, receiver), "_+_", pmax(sender, receiver))
  weights <- base::factor(dyad.ids)
  ########################################################
  #
  #   Returning the values back to the user as an attachment to the class
  #
  ########################################################
  if(return_stats){
    data <- weights #returning just the vector of computed stats
  }else{
    data$statistics[["dyad.fe"]] <-weights
  }
  return(data) # output the file to the user!s
}








#' @title Compute Degree Network Statistics for Event Senders and Receivers in a Post-Processing Relational Event Sequence
#' @name dreamstats_degree
#' @param formation The degree statistic to be computed. "sender-indegree" computes the indegree statistic for the event senders. "receiver-indegree" computes the
#' indegree statistic for the event receivers. "sender-outdegree" computes the outdegree statistic for the event senders. "receiver-outdegree" computes the
#' outdegree statistic for the event receivers.
#' @param data An object of class `dream_sequence` that contains the processed relational event sequence.
#' @param counts TRUE/FALSE. TRUE indicates that the counts of past events should be computed (see the details section). FALSE indicates that
#' the temporal exponential weighting function should be used to downweigh past events (see the details section). Set to FALSE by default.
#' @param halflife A numerical value that is the halflife value to be used in the exponential weighting function (see details section). Preset to
#' 2 (should be updated by the user based on substantive context).
#' @param dyadic_weight A numerical value for the dyadic cutoff weight that represents the numerical cutoff value for temporal relevancy
#' based on the exponential weighting function. For example, a numerical value of 0.01, indicates that an exponential weight less than 0.01
#' will become 0 and that events with such value (or smaller values) will not be included in the sum of the past event
#' weights (see the details section). Set to 0 by default.
#' @param exp_weight_form TRUE/FALSE. TRUE indicates that the Lerner et al. (2013) exponential weighting function will be used (see the details section). FALSE indicates that the Lerner and Lomi (2020) exponential weighting function will be used (see the details section). Set to FALSE by default
#' @param return_stats TRUE/FALSE. TRUE indicates that the vector of computed
#' statistics will be returned. FALSE indicates that the vector of computed
#' statistics will be added to the `statistics` element of the `data` argument.
#' @import Rcpp
#' @return The vector of degree statistics or the updated `data` argument.
#' @export


#' @description
#' `r lifecycle::badge("stable")`
#'
#' The function computes the various degree network sufficient statistic for event senders
#' in a relational event sequence (see Lerner and Lomi 2020; Butts 2008).
#' This measure allows for the degree values to be only computed for the sampled
#' events, while creating the weights based on the full event sequence (see
#' Lerner and Lomi 2020; Vu et al. 2015). The function also allows users to use two
#' different weighting functions, return the counts of past events, reduce computational
#' runtime, and specify a dyadic cutoff for relational relevancy.
#'
#'
#'@details The function calculates degree values for relational event
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
#'\code{\link{dreamstats_dyadcut}} function.
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
#'eventSet <- create_res(type = "one-mode",
#'                       ordinal = TRUE,
#'                       riskset = "constant_sample",
#'                       time = events$time,
#'                       sender = events$sender,
#'                       receiver = events$target,
#'                       p_samplingobserved = 1.00,
#'                       n_controls = 1,
#'                       seed = 9999)
#'
#'#Computing the sender indegree statistic for the relational event sequence
#'eventSet <- dreamstats_degree(
#'    formation = "sender-indegree",
#'    data = eventSet,
#'    halflife = 2, #halflife parameter
#'    dyadic_weight = 0,
#'    exp_weight_form = FALSE)
#'
#'eventSet #printing the post-processed relational event sequence
#'eventSet$statistics$sender.indegree #printing the vector of computed values
#'
#'#Computing the sender indegree statistic for the relational event sequence
#'#and returning only the vector of computed sender indegree values
#'degree.stat <- dreamstats_degree(
#'    formation = "sender-indegree",
#'    data = eventSet,
#'    halflife = 2, #halflife parameter
#'    dyadic_weight = 0,
#'    exp_weight_form = FALSE,
#'    return_stats = TRUE)
#'
#'cor(degree.stat, eventSet$statistics$sender.indegree)
#'
#'#Computing the sender outdegree statistic for the relational event sequence
#'eventSet <- dreamstats_degree(
#'    formation = "sender-outdegree",
#'    data = eventSet,
#'    halflife = 2, #halflife parameter
#'    dyadic_weight = 0,
#'    exp_weight_form = FALSE)
#'
#'
#'#Computing the receiver outdegree statistic for the relational event sequence
#'eventSet <- dreamstats_degree(
#'    formation = "receiver-outdegree",
#'    data = eventSet,
#'    halflife = 2, #halflife parameter
#'    dyadic_weight = 0,
#'    exp_weight_form = FALSE)
#'
#'
#'#Computing the receiver indegree statistic for the relational event sequence
#'eventSet <- dreamstats_degree(
#'    formation = "receiver-indegree",
#'    data = eventSet,
#'    halflife = 2, #halflife parameter
#'    dyadic_weight = 0,
#'    exp_weight_form = FALSE)
#'
#'#printing the post-processed relational event sequence that contains all computed degree statistics
#'degree.info <- as.data.frame(eventSet)
#'degree.info #printing the information to the user
#'

dreamstats_degree <- function(formation = c("sender-indegree",
                                            "receiver-indegree",
                                            "sender-outdegree",
                                            "receiver-outdegree"),
                              data,
                              halflife=2, # the half life value for the weighting function
                              counts = FALSE, #Logical indicating if the raw counts of events should be returned or the exponential weighting function should be used (TRUE = counts; FALSE = exponential weighting)
                              dyadic_weight= 0.00, # dyadic cutoff weight for events that no longer matter
                              exp_weight_form = FALSE, # Do we want to use the weighting function of Lerner et al. 2013 (alsoused in the rem R package)?
                              return_stats = FALSE #returning the statistics to the user
) {

  ########################################################
  #
  #   Checking for Errors in User Inputs
  #
  ########################################################
  if(!inherits(data,"dream_sequence")) base::stop("Error: the object for the `data` argument is not a `dream_sequence` object.")
  time <- data$processed_sequence$time # variable (column) name that contains the time variable
  sender<- data$processed_sequence$sender# variable (column) name that contains the sender variable
  receiver<- data$processed_sequence$receiver# variable (column) name that contains the target variable
  observed<- data$processed_sequence$observed# variable (column) name that contains the observed variable
  sampled<- data$processed_sequence$sampled# variable (column) name that contains the sampled variable

  if(data$type == "two-mode" & formation %in% c("sender-indegree","receiver-outdegree")) base::stop("Error: The sender-indegree and receiver-outdegree statistics are not defined for two-mode event sequences.") # stop computation and tell the user


  if(halflife < 0){
    base::stop("Error: Halflife values must be positive.") # stop computation and tell the user
  }
  if(sum(observed) == 0){
    base::stop("Error: There are no observed events based upon the 'observed' input. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(sum(sampled) == 0){
    base::stop("Error: There are no sampled events based upon the 'sampled' input. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(!is.numeric(time)){
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
  if(return_stats){
    data <- weights
  }else{
    if(formation == "sender-indegree"){data$statistics$sender.indegree=weights }
    if(formation == "receiver-indegree"){data$statistics$receiver.indegree=weights}
    if(formation == "sender-outdegree"){data$statistics$sender.outdegree=weights}
    if(formation == "receiver-outdegree"){data$statistics$receiver.outdegree=weights}
  }
  return(data) # output the file to the user!
}








#' @title Compute Butts' (2008) Recency Network Statistic for Event Dyads in a Relational Event Sequence
#' @description
#' `r lifecycle::badge("stable")`
#'
#'This function computes the recency network sufficient statistic for a relational
#'event sequence (see Butts 2008; Vu et al. 2015; Meijerink-Bosman et al. 2022). The
#'recency statistic captures the tendency for more recent events (i.e., an exchange
#'between two medical doctors) are more likely to re-occur in comparison to events
#'that happened in the more distant past (see Butts 2008 for a discussion). This
#'measure allows for recency scores to be only computed for the sampled events,
#'while computing the statistics based on the full event sequence.
#'
#' @name dreamstats_recency
#' @param data A `dream_sequence` object that contains the processed relational event sequence.
#' @param type A string value that specifies which recency formula will be used to compute the statistics. The options are "raw.diff", "inv.diff.plus1", "rank.ordered.count" (see details section).
#' @param i_neighborhood TRUE/FALSE. TRUE indicates that the recency statistic will be computed in reference to the sender’s past relational history (see details section). FALSE indicates that the recency statistic will be computed in reference to the target’s past relational history (see details section). Set to TRUE by default.
#' @param nopastEvents The numerical value that specifies what value should be given to events in which the sender was not active as a sender in the past (i’s neighborhood when i_neighborhood = TRUE) or was not the recipient of a past event (j’s neighborhood when i_neighborhood = FALSE). Set to NA by default.
#' @param dependency TRUE/FALSE. TRUE indicates that temporal relevancy will be modeled (see details section). FALSE indicates that temporal relevancy will not be modeled, that is, all past events are relevant (see details section). Set to FALSE by default.
#' @param relationalTimeSpan If dependency = TRUE, a numerical value that corresponds to the temporal span for relational relevancy, which must be the same measurement unit as the observed_time and processed_time objects. When dependency = TRUE, the relevant events are events that have occurred between current event time, *t*, and *t - relationalTimeSpan*. For example, if the time measurement is the number of days since the first event and the value for relationalTimeSpan is set to 10, then only those events which occurred in the past 10 days are included in the computation of the statistic.
#' @param return_stats TRUE/FALSE. TRUE indicates that the vector of computed
#' statistics will be returned. FALSE indicates that the vector of computed
#' statistics will be added to the `statistics` element of the `data` argument.
#' @import Rcpp
#' @importFrom methods is
#' @return The vector of recency network statistics for the relational event sequence or the updated `data` argument.
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
#'post.processing <- create_res(type = "one-mode",
#'                            riskset = "constant_sample",
#'                       ordinal = TRUE,
#'                       time = events$time,
#'                       sender = events$sender,
#'                       receiver = events$target,
#'                       p_samplingobserved = 1.00,
#'                       n_controls = 1,
#'                       seed = 9999)
#'
#'#Computing the recency statistics (with raw time difference) for the relational event sequence
#'post.processing <- dreamstats_recency(data = post.processing,
#'                                      type = "raw.diff")
#'
#'#printing the post-processed relational event sequence
#'post.processing
#'
#'#Computing the recency statistics (with raw time difference) for the relational event sequence
#'#and returning only the vector of computed values
#'recency.stat <- dreamstats_recency(data = post.processing,
#'                                   type = "raw.diff",
#'                                   return_stats = TRUE)
#'
#'cor(recency.stat, post.processing$statistics$recency)
#'
#'#Computing the recency statistics (with inverse of time difference) for the
#'#relational event sequence
#'post.processing <- dreamstats_recency(data = post.processing,
#'                                      type = "inv.diff.plus1")
#'
#'#Computing the rank-based recency statistics for the relational event sequence
#'post.processing <- dreamstats_recency(data = post.processing,
#'                                      type = "rank.ordered.count")
#'

dreamstats_recency <-   function( data,
                                  type = c("raw.diff", "inv.diff.plus1", "rank.ordered.count"),
                                  i_neighborhood = TRUE, #should the recency be computed on the i's neighborhood or j's neighborhood
                                  dependency = FALSE, #Boolean for temporal dependency
                                  relationalTimeSpan = NULL, #if dependency == TRUE, this should specific the associated temporal time span
                                  nopastEvents = NA, #the value given to events where there are no past events (defaults to NA)
                                  return_stats = FALSE
) {

  ########################################################
  #
  #   Checking for Errors in User Inputs
  #
  ########################################################
  if(!inherits(data,"dream_sequence")) base::stop("Error: the object for the `data` argument is not a `dream_sequence` object.")
  time <- data$processed_sequence$time # variable (column) name that contains the time variable
  sender<- data$processed_sequence$sender# variable (column) name that contains the sender variable
  receiver<- data$processed_sequence$receiver# variable (column) name that contains the target variable
  observed<- data$processed_sequence$observed# variable (column) name that contains the observed variable
  sampled<- data$processed_sequence$sampled# variable (column) name that contains the sampled variable
  if(dependency == TRUE & is.null(relationalTimeSpan)){
    base::stop("Error: Temporal dependency was requested, however, the relationalTimeSpan is missing. Please add this and restart the function.") # stop computation and tell the user
  }


  if(sum(observed) == 0){
    base::stop("Error: There are no observed events based upon the 'observed' input. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(sum(sampled) == 0){
    base::stop("Error: There are no sampled events based upon the 'sampled' input. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(!is.numeric(time)){
    base::stop("Error: The 'time' argument is not a numeric value. Stopping computation! Please update this argument.") # stop computation and tell the user
  }


  appender <- "__NIKOACAR2020__" # a simple (hopefully) unique string
  controlevents <- 1 - observed #making the control events 1 and the observed events 0
  if(type != "rank.ordered.count"){
    #if raw difference is requested, update the object to be TRUE, (if inverse), then make FALSE
    raw_diff <- ifelse(type == "raw.diff", TRUE, FALSE)
    weights <- computerecencynorank(time=time,
                                    sampledevent=sampled,
                                    controlevent=controlevents,
                                    sender=paste0(sender),
                                    target=paste0(receiver),
                                    dyad_id=paste0(sender, appender, receiver),
                                    raw_diff=raw_diff,
                                    i_neighborhood=i_neighborhood,
                                    appender = appender,
                                    nopastEvents = nopastEvents)
  }else{ #computing rank based recency
    weights <- computerecencyrank(time=time,
                                  sampledevent=sampled,
                                  controlevent=controlevents,
                                  sender=paste0(sender),
                                  target=paste0(receiver),
                                  i_neighborhood=i_neighborhood,
                                  appender = appender,
                                  nopastEvents = nopastEvents)
  }
  ########################################################
  #
  #   Returning the values back to the user as an attachment to the class
  #
  ########################################################
  if(return_stats){
    data <- weights #returning just the vector of computed stats
  }else{
    data$statistics$recency <-weights
  }
  return(data) # output the file to the user!

}





#' @title Compute Butts' (2008) Preferential Attachment Network Statistic for Event Dyads in a Relational Event Sequence
#' @description
#' `r lifecycle::badge("stable")`
#'
#'
#' The function computes the preferential attachment network sufficient statistic for
#' a relational event sequence (see Butts 2008). Preferential attachment measures the tendency towards a
#' positive feedback loop in which actors involved in more past events are more likely to be involved
#' in future events (see Butts 2008 for an empirical example and discussion).This measure allows
#' for preferential attachment scores to be only computed for the sampled events, while creating the statistics based on the full event
#' sequence. Moreover, the function allows users to specify relational relevancy for the resulting statistics.

#' @name dreamstats_prefattachment
#' @param data A `dream_sequence` object that contains the processed relational event sequence.
#' @param dependency TRUE/FALSE. TRUE indicates that temporal relevancy will be modeled (see the details section). FALSE indicates that temporal relevancy will not be modeled, that is, all past events are relevant (see the details section). Set to FALSE by default.
#' @param relationalTimeSpan If dependency = TRUE, a numerical value that corresponds to the temporal span for relational relevancy, which must be the same measurement unit as the observed_time and processed_time objects. When dependency = TRUE, the relevant events are events that have occurred between current event time, *t*, and *t - relationalTimeSpan*. For example, if the time measurement is the number of days since the first event and the value for relationalTimeSpan is set to 10, then only those events which occurred in the past 10 days are included in the computation of the statistic.
#' @param return_stats TRUE/FALSE. TRUE indicates that the vector of computed
#' statistics will be returned. FALSE indicates that the vector of computed
#' statistics will be added to the `statistics` element of the `data` argument.
#' @import Rcpp
#' @importFrom methods is
#' @return The vector of event preferential attachment statistics for the relational event sequence or the updated `data` argument.
#' @export
#'
#'
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#' @references
#' Butts, Carter T. 2008. "A relational event framework for social action." *Sociological Methodology* 38(1): 155-200.
#'
#' Quintane, Eric, Martin Wood, John Dunn, and Lucia Falzon. 2022. “Temporal
#' Brokering: A Measure of Brokerage as a Behavioral Process.” *Organizational Research Methods*
#' 25(3): 459-489.
#'
#'@details The function calculates preferential attachment for a relational event sequence based on Butts (2008).
#'
#'Following Butts (2008), the formula for preferential attachment for event \eqn{e_i} is:
#'\deqn{PA_{e_{i}} = \frac{d^{+}(r(e_{i}), A_t)+d^{-}(r(e_{i}), A_t)}{\sum_{i=1}^{|S|} (d^{+}(i, A_t)+d^{-}(i, A_t))} }
#'
#'where  \eqn{d^{+}(r(e_{i}), A_t)} is the past outdegree of the receiver for \eqn{e_i},  \eqn{d^{-}(r(e_{i}), A_t)} is the past indegree of the receiver for \eqn{e_i},
#'\eqn{\sum_{i=1}^{|S|} (d^{+}(i, A_t)+d^{-}(i, A_t))} is the sum of the past outdegree and indegree for all past event senders in the relational history.
#'
#'Moreover, researchers interested in modeling temporal relevancy (see Quintane,
#'Mood, Dunn, and Falzone 2022) can specify the relational time span, that is, length of time for which events are considered
#'relationally relevant. This should be specified via the option *relationalTimeSpan* with *dependency* set to TRUE.
#'
#'
#'
#'@examples
#'
#'
#'# A Psuedo One-Mode Event Dataset
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
#'post.processing <- create_res(type = "one-mode",
#'                       ordinal = TRUE,
#'                       riskset = "constant_sample",
#'                       time = events$time,
#'                       sender = events$sender,
#'                       receiver = events$target,
#'                       p_samplingobserved = 1.00,
#'                       n_controls = 1,
#'                       seed = 9999)
#'
#'#Computing the preferential attachment statistic for the relational event sequence
#'post.processing <- dreamstats_prefattachment(data = post.processing)
#'
#'#printing the post-processed relational event sequence
#'post.processing
#'
#'
#'#Computing the preferential attachment statistic for the relational event sequence
#'#and returning only the vector of computed values
#'prefattach.stat <- dreamstats_prefattachment(data = post.processing,
#'                                             return_stats = TRUE)
#'
#'cor(prefattach.stat, post.processing$statistics$pref.attachment)



dreamstats_prefattachment <-   function(data,
                                        dependency = FALSE,
                                        relationalTimeSpan =0, # the sizes of the windows that we will use, if NA, we will compute it internally
                                        return_stats = FALSE
) {

  ########################################################
  #
  #   Checking for Errors in User Inputs
  #
  ########################################################
  if(!inherits(data,"dream_sequence")) base::stop("Error: the object for the `data` argument is not a `dream_sequence` object.")
  time <- data$processed_sequence$time # variable (column) name that contains the time variable
  sender<- data$processed_sequence$sender# variable (column) name that contains the sender variable
  receiver<- data$processed_sequence$receiver# variable (column) name that contains the target variable
  observed<- data$processed_sequence$observed# variable (column) name that contains the observed variable
  sampled<- data$processed_sequence$sampled# variable (column) name that contains the sampled variable

  if(dependency == TRUE & is.null(relationalTimeSpan)){
    base::stop("Error: Temporal dependency was requested, however, the relationalTimeSpan is missing. Please add this and restart the function.") # stop computation and tell the user
  }


  if(sum(observed) == 0){
    base::stop("Error: There are no observed events based upon the 'observed' input. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(sum(sampled) == 0){
    base::stop("Error: There are no sampled events based upon the 'sampled' input. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(!is.numeric(time)){
    base::stop("Error: The 'time' argument is not a numeric value. Stopping computation! Please update this argument.") # stop computation and tell the user
  }

  controlevents <- 1 - observed #making the control events 1 and the observed events 0
  if(!dependency){ #if temporal dependency is not requested
    weights <- computeremprefattach(time=time,
                                    sampledevent= sampled,
                                    controlevent= controlevents,
                                    sender= paste0(sender),
                                    target= paste0(receiver))
  }else{ #if temporal dependency is requested
    weights <- prefattachrelspanrem(time=time,
                                    sampledevent= sampled,
                                    controlevent= controlevents,
                                    sender= paste0(sender),
                                    target= paste0(receiver),
                                    reltimespan = relationalTimeSpan)

  }
  ########################################################
  #
  #   Returning the values back to the user as an attachment to the class
  #
  ########################################################
  if(return_stats){
    data <- weights #returning just the vector of computed stats
  }else{
    data$statistics$pref.attachment <-weights
  }
  return(data) # output the file to the user!
}






#' @title Compute Butts' (2008) Persistence Network Statistic for Event Dyads in a Relational Event Sequence
#' @description
#' `r lifecycle::badge("stable")`
#'
#' This function computes the persistence network sufficient statistic for
#' a relational event sequence (see Butts 2008). Persistence measures the proportion of past ties sent from the event sender that went to the current event receiver.
#' Furthermore, this measure allows for persistence scores to be only
#' computed for the sampled events, while creating the weights based on the full event
#' sequence. Moreover, the function allows users to specify relational relevancy for the resulting statistic.

#' @name dreamstats_persistence
#' @param data A `dream_sequence` object that contains the processed relational event sequence.
#' @param ref_sender TRUE/FALSE. TRUE indicates that the persistence statistic will be computed in reference to the sender’s past relational history (see details section). FALSE indicates that the persistence statistic will be computed in reference to the target’s past relational history (see details section). Set to TRUE by default.
#' @param nopastEvents The numerical value that specifies what value should be given to events in which the sender has sent not past ties (i's neighborhood when sender = TRUE) or has not received any past ties (j's neighborhood when sender = FALSE). Set to NA by default.
#' @param dependency TRUE/FALSE. TRUE indicates that temporal relevancy will be modeled (see the details section). FALSE indicates that temporal relevancy will not be modeled, that is, all past events are relevant (see the details section). Set to FALSE by default.
#' @param relationalTimeSpan If dependency = TRUE, a numerical value that corresponds to the temporal span for relational relevancy, which must be the same measurement unit as the observed_time and processed_time objects. When dependency = TRUE, the relevant events are events that have occurred between current event time, *t*, and *t-relationalTimeSpan*. For example, if the time measurement is the number of days since the first event and the value for relationalTimeSpan is set to 10, then only those events which occurred in the past 10 days are included in the computation of the statistic.
#' @param return_stats TRUE/FALSE. TRUE indicates that the vector of computed
#' statistics will be returned. FALSE indicates that the vector of computed
#' statistics will be added to the `statistics` element of the `data` argument.
#' @import Rcpp
#' @importFrom methods is
#' @return The vector of persistence network statistics for the relational event sequence or the updated `data` argument.
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
#'# A Psuedo One-Mode Event Dataset
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
#'post.processing <- create_res(type = "one-mode",
#'                       ordinal = TRUE,
#'                       riskset = "constant_sample",
#'                       time = events$time,
#'                       sender = events$sender,
#'                       receiver = events$target,
#'                       p_samplingobserved = 1.00,
#'                       n_controls = 1,
#'                       seed = 9999)
#'
#'#Computing the persistence statistic for the relational event sequence
#'post.processing <- dreamstats_persistence(data = post.processing)
#'
#'#printing the post-processed relational event sequence
#'post.processing
#'
#'
#'#Computing the persistence statistic for the relational event sequence
#'#and returning only the vector of computed values
#'persistence.stat <- dreamstats_persistence(data = post.processing,
#'                                           return_stats = TRUE)
#'
#'cor(persistence.stat, post.processing$statistics$persistence)



dreamstats_persistence <-   function(data,
                                     ref_sender = TRUE,
                                     nopastEvents = NA, #the value given to events where there are no past events (defaults to NA)
                                     dependency = FALSE,
                                     relationalTimeSpan = 0,
                                     return_stats = FALSE
) {

  ########################################################
  #
  #   Checking for Errors in User Inputs
  #
  ########################################################
  if(!inherits(data,"dream_sequence")) base::stop("Error: the object for the `data` argument is not a `dream_sequence` object.")
  time <- data$processed_sequence$time # variable (column) name that contains the time variable
  sender<- data$processed_sequence$sender# variable (column) name that contains the sender variable
  target<- data$processed_sequence$receiver# variable (column) name that contains the target variable
  observed<- data$processed_sequence$observed# variable (column) name that contains the observed variable
  sampled<- data$processed_sequence$sampled# variable (column) name that contains the sampled variable

  if(dependency == TRUE & is.null(relationalTimeSpan)){
    base::stop("Error: Temporal dependency was requested, however, the relationalTimeSpan is missing. Please add this and restart the function.") # stop computation and tell the user
  }
  if(sum(observed) == 0){
    base::stop("Error: There are no observed events based upon the 'observed' input. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(sum(sampled) == 0){
    base::stop("Error: There are no sampled events based upon the 'sampled' input. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(!is.numeric(time)){
    base::stop("Error: The 'time' argument is not a numeric value. Stopping computation! Please update this argument.") # stop computation and tell the user
  }


  ########################################################
  #
  #   Prepping for c++ computation
  #
  ########################################################
  dyad_sep <- "__NIKOACAR2020__"
  if(ref_sender){ #if the sender is the reference
    actor <- sender
    dyad_id <- paste0(sender, dyad_sep, target)
  }else{#if the sender is the not reference
    actor <- target
    dyad_id <- paste0(sender, dyad_sep, target)
  }

  ########################################################
  #
  #   Computing in c++ for speedy computation
  #
  ########################################################
  persistence <- persistencerem(time =time,
                                sampledevent=sampled,
                                controlevent=observed,
                                dyad_id=dyad_id,
                                actor=actor,
                                timedependency=dependency,
                                cuttime=relationalTimeSpan,
                                nopastEvents=nopastEvents)

  ########################################################
  #
  #   Returning the values back to the user as an attachment to the class
  #
  ########################################################
  if(return_stats){
    data <- persistence #returning just the vector of computed stats
  }else{
    data$statistics$persistence <-persistence
  }
  return(data) # output the file to the user!

}








#' @title Compute Butts' (2008) Repetition Network Statistic for Event Dyads in a Relational Event Sequence
#' @description
#' `r lifecycle::badge("stable")`
#'
#' This function computes the repetition network sufficient statistic for a relational event
#' sequence (see Lerner and Lomi 2020; Butts 2008). Repetition measures the increased tendency
#' for events between S and R to occur given that S and R have interacted
#' in the past. Furthermore, this function allows for repetition scores to be only
#' computed for the sampled events, while creating the weights based on the full event
#' sequence (see Lerner and Lomi 2020; Vu et al. 2015). The function also allows
#' users to use two different weighting functions, return the counts of past events,
#' reduce computational runtime, and specify a dyadic cutoff for relational relevancy.

#' @name dreamstats_repetition
#' @param data A `dream_sequence` object that contains the processed relational event sequence.
#' @param counts TRUE/FALSE. TRUE indicates that the counts of past events should be computed (see the details section). FALSE indicates that the temporal exponential weighting function should be used to downweigh past events (see the details section). Set to FALSE by default.
#' @param halflife A numerical value that is the halflife value to be used in the exponential weighting function (see details section). Preset to 2 (should be updated by the user based on substantive context).
#' @param dyadic_weight A numerical value for the dyadic cutoff weight that represents the numerical cutoff value for temporal relevancy based on the exponential weighting function. For example, a numerical value of 0.01, indicates that an exponential weight less than 0.01 will become 0 and that events with such value (or smaller values) will not be included in the sum of the past event weights (see the details section). Set to 0 by default.
#' @param exp_weight_form TRUE/FALSE. TRUE indicates that the Lerner et al. (2013) exponential weighting function will be used (see the details section). FALSE indicates that the Lerner and Lomi (2020) exponential weighting function will be used (see the details section). Set to FALSE by default
#' @param return_stats TRUE/FALSE. TRUE indicates that the vector of computed
#' statistics will be returned. FALSE indicates that the vector of computed
#' statistics will be added to the `statistics` element of the `data` argument.
#' @import Rcpp
#' @importFrom methods is
#' @return The vector of repetition statistics for the relational event sequence or the updated `data` argument.
#' @export
#'
#'
#'
#'@details This function calculates the repetition scores for relational event models
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
#'The formula for repetition for event \eqn{e_i} is:
#'\deqn{repetition_{e_{i}} = w(s, r, t) }
#'
#'Moreover, researchers interested in modeling temporal relevancy (see Quintane,
#'Mood, Dunn, and Falzone 2022; Lerner and Lomi 2020) can specify the dyadic
#'weight cutoff, that is, the minimum value for which the weight is considered
#'relationally relevant. Users who do not know the specific dyadic cutoff value to use, can use the
#'\code{\link{dreamstats_dyadcut}} function.
#'
#'Following Butts (2008), if the counts of the past events are requested, the formula for repetition for
#'event \eqn{e_i} is:
#'\deqn{repetition_{e_{i}} = d(s = s', r = r', t') }
#'Where, \eqn{d()} is the number of past events where the event sender, *s'*, is the current event sender, *s*, the event
#'receiver (target), *r'*, is the current event receiver, *r*. Moreover, the counting equation
#'can be used in tandem with relational relevancy, by specifying the halflife parameter, exponential
#'weighting function, and the dyadic cut off weight values. If the user is not interested in modeling
#'relational relevancy, then those value should be left at their baseline values.
#'
#'@examples
#'data("WikiEvent2018.first100k", package = "dream")
#'WikiEvent2018 <- WikiEvent2018.first100k[1:1000,] #the first one thousand events
#'WikiEvent2018$time <- as.numeric(WikiEvent2018$time) #making the variable numeric
#'### Creating the EventSet By Employing Case-Control Sampling With M = 5 and
#'### Sampling from the Observed Event Sequence with P = 0.01
#'post.processing <- create_res(type = "two-mode",
#'  ordinal = TRUE,
#'  riskset = "constant_sample",
#'  time = WikiEvent2018$time, # The Time Variable
#'  sender = WikiEvent2018$user, # The Sender Variable
#'  receiver = WikiEvent2018$article, # The Receiver Variable
#'  p_samplingobserved = 0.01, # The Probability of Selection
#'  n_controls = 8, # The Number of Controls to Sample from the Full Risk Set
#'  seed = 9999)
#'
#'
#'#Computing the repetition statistics for the relational event sequence with
#'#the exponential weights of past events returned
#'post.processing <- dreamstats_repetition(data = post.processing,
#'                                         halflife = 2.592e+09)
#'
#'#printing the post-processed relational event sequence
#'post.processing
#'
#'
#'#Computing the repetition statistic for the relational event sequence
#'#and returning only the vector of computed values
#'repetition.stat <-  dreamstats_repetition(data = post.processing,
#'                                          halflife = 2.592e+09,
#'                                          return_stats = TRUE)
#'
#'cor(repetition.stat, post.processing$statistics$repetition)
#'
#'#Computing the repetition statistics for the relational event sequence with
#'#the counts of past events returned
#'post.processing <- dreamstats_repetition(data = post.processing,
#'                                         halflife = 2.592e+09,
#'                                         counts = TRUE)
#'
#'cbind(post.processing$statistics$repetition, repetition.stat)
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
#'
#'
#
#'
#'
#'
dreamstats_repetition <-    function(data,
                                     halflife=2, # the half life value for the weighting function
                                     counts = FALSE, #Logical indicating if the raw counts of events should be returned or the exponential weighting function should be used (TRUE = counts; FALSE = exponential weighting)
                                     dyadic_weight= 0.00, # dyadic cutoff weight for events that no longer matter
                                     exp_weight_form = FALSE, # Do we want to use the weighting function of Lerner et al. 2013 (alsoused in the rem R package)?
                                     return_stats = FALSE
) {

  ########################################################
  #
  #   Checking for Errors in User Inputs
  #
  ########################################################
  if(!inherits(data,"dream_sequence")) base::stop("Error: the object for the `data` argument is not a `dream_sequence` object.")
  time <- data$processed_sequence$time # variable (column) name that contains the time variable
  sender<- data$processed_sequence$sender# variable (column) name that contains the sender variable
  receiver<- data$processed_sequence$receiver# variable (column) name that contains the target variable
  observed<- data$processed_sequence$observed# variable (column) name that contains the observed variable
  sampled<- data$processed_sequence$sampled# variable (column) name that contains the sampled variable
  if(halflife < 0){
    base::stop("Error: Halflife values must be positive.") # stop computation and tell the user
  }
  if(sum(observed) == 0){
    base::stop("Error: There are no observed events based upon the 'observed' input. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(sum(sampled) == 0){
    base::stop("Error: There are no sampled events based upon the 'observed' input. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(!is.numeric(time)){
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
  weightSchemeR <- ifelse(Lerneretal_2013, 0, 1) #setting this argument up for c++ computation
  countsR <- ifelse(counts, 1, 0) #setting this argument up for c++ computation
  controlR <- 1- observed #making it such that dummy events have a 1 and real events have a value of 0
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

  ########################################################
  #
  #   Returning the values back to the user as an attachment to the class
  #
  ########################################################
  if(return_stats){
    data <- weights #returning just the vector of computed stats
  }else{
    data$statistics$repetition <-weights
  }
  return(data) # output the file to the user!s
}






#' @title Compute the Reciprocity Network Statistic for Event Dyads in a Relational Event Sequence
#' @name dreamstats_reciprocity
#' @param data An object of class `dream_sequence` that contains the processed relational event sequence.
#' @param counts TRUE/FALSE. TRUE indicates that the counts of past events should be computed (see the details section). FALSE indicates that the temporal exponential weighting function should be used to downweigh past events (see the details section). Set to FALSE by default.
#' @param halflife A numerical value that is the halflife value to be used in the exponential weighting function (see details section). Preset to 2 (should be updated by the user based on substantive context).
#' @param dyadic_weight A numerical value for the dyadic cutoff weight that represents the numerical cutoff value for temporal relevancy based on the exponential weighting function. For example, a numerical value of 0.01, indicates that an exponential weight less than 0.01 will become 0 and that events with such value (or smaller values) will not be included in the sum of the past event weights (see the details section). Set to 0 by default.
#' @param exp_weight_form TRUE/FALSE. TRUE indicates that the Lerner et al. (2013) exponential weighting function will be used (see the details section). FALSE indicates that the Lerner and Lomi (2020) exponential weighting function will be used (see the details section). Set to FALSE by default
#' @param return_stats TRUE/FALSE. TRUE indicates that the vector of computed
#' statistics will be returned. FALSE indicates that the vector of computed
#' statistics will be added to the `statistics` element of the `data` argument.
#' @import Rcpp
#' @return The vector of reciprocity statistics for the relational event sequence or the updated `data` argument.
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
#'\code{\link{dreamstats_dyadcut}} function.
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
#'# A Psuedo One-Mode Event Dataset
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
#'post.processing <- create_res(type = "one-mode",
#'                       ordinal = TRUE,
#'                       riskset = "constant_sample",
#'                       time = events$time,
#'                       sender = events$sender,
#'                       receiver = events$target,
#'                       p_samplingobserved = 1.00,
#'                       n_controls = 1,
#'                       seed = 9999)
#'
#'#Computing the reciprocity statistic for the relational event sequence
#'post.processing <- dreamstats_reciprocity(data = post.processing,
#'                                          halflife = 2)
#'
#'#printing the post-processed relational event sequence
#'post.processing
#'
#'
#'#Computing the reciprocity statistic for the relational event sequence
#'#and returning only the vector of computed values
#'reciprocity.stat <-  dreamstats_reciprocity(data = post.processing,
#'                                          halflife = 2,
#'                                          return_stats = TRUE)
#'
#'cor(reciprocity.stat, post.processing$statistics$reciprocity)

dreamstats_reciprocity <- function(data, #a dream class object
                                   halflife=2, # the half life value for the weighting function
                                   counts = FALSE, #Logical indicating if the raw counts of events should be returned or the exponential weighting function should be used (TRUE = counts; FALSE = exponential weighting)
                                   dyadic_weight= 0.00, # dyadic cutoff weight for events that no longer matter
                                   exp_weight_form = FALSE, # Do we want to use the weighting function of Lerner et al. 2013 (alsoused in the rem R package)?
                                   return_stats = FALSE) {


  ########################################################
  #
  #   Checking for Errors in User Inputs
  #
  ########################################################
  if(data$type == "two-mode") base::stop("Error: The reciprocity statistic is not defined for two-mode event sequences.") # stop computation and tell the user
  if(!inherits(data,"dream_sequence")) base::stop("Error: the object for the `data` argument is not a `dream_sequence` object.")
  time <- data$processed_sequence$time # variable (column) name that contains the time variable
  sender<- data$processed_sequence$sender# variable (column) name that contains the sender variable
  receiver<- data$processed_sequence$receiver# variable (column) name that contains the target variable
  observed<- data$processed_sequence$observed# variable (column) name that contains the observed variable
  sampled<- data$processed_sequence$sampled# variable (column) name that contains the sampled variable

  if(halflife < 0){
    base::stop("Error: Halflife values must be positive.") # stop computation and tell the user
  }
  if(sum(observed) == 0){
    base::stop("Error: There are no observed events based upon the 'observed' input. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(sum(sampled) == 0){
    base::stop("Error: There are no sampled events based upon the 'observed' input. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(!is.numeric(time)){
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
  #   Returning the values back to the user as an attachment to the class
  #
  ########################################################
  if(return_stats){
    data <- weights #returning just the vector of computed stats
  }else{
    data$statistics$reciprocity <-weights
  }
  return(data) # output the file to the user!s
}




#' @title Compute the Four-Cycles Network Statistic for Event Dyads in a Relational Event Sequence
#' @name  dreamstats_fourcycles
#' @param data An object of class `dream_sequence` that contains the processed relational event sequence.
#' @param counts TRUE/FALSE. TRUE indicates that the counts of past events should be computed (see the details section). FALSE indicates that the temporal exponential weighting function should be used to downweigh past events (see the details section). Set to FALSE by default.
#' @param halflife A numerical value that is the halflife value to be used in the exponential weighting function (see details section). Preset to 2 (should be updated by the user based on substantive context).
#' @param dyadic_weight A numerical value for the dyadic cutoff weight that represents the numerical cutoff value for temporal relevancy based on the exponential weighting function. For example, a numerical value of 0.01, indicates that an exponential weight less than 0.01 will become 0 and that events with such value (or smaller values) will not be included in the sum of the past event weights (see the details section). Set to 0 by default.
#' @param exp_weight_form TRUE/FALSE. TRUE indicates that the Lerner et al. (2013) exponential weighting function will be used (see the details section). FALSE indicates that the Lerner and Lomi (2020) exponential weighting function will be used (see the details section). Set to FALSE by default
#' @param return_stats TRUE/FALSE. TRUE indicates that the vector of computed
#' statistics will be returned. FALSE indicates that the vector of computed
#' statistics will be added to the `statistics` element of the `data` argument.
#' @import Rcpp
#' @importFrom methods is
#' @return The vector of four-cycle statistics for the two-mode relational event sequence or the updated `data` argument.
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
#'\code{\link{dreamstats_dyadcut}} function.
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
#'data("WikiEvent2018.first100k", package = "dream")
#'WikiEvent2018 <- WikiEvent2018.first100k[1:1000,] #the first one thousand events
#'WikiEvent2018$time <- as.numeric(WikiEvent2018$time) #making the variable numeric
#'### Creating the EventSet By Employing Case-Control Sampling With M = 5 and
#'### Sampling from the Observed Event Sequence with P = 0.01
#'post.processing <- create_res(type = "two-mode",
#'  ordinal = TRUE,
#'  riskset = "constant_sample",
#'  time = WikiEvent2018$time, # The Time Variable
#'  sender = WikiEvent2018$user, # The Sender Variable
#'  receiver = WikiEvent2018$article, # The Receiver Variable
#'  p_samplingobserved = 0.01, # The Probability of Selection
#'  n_controls = 8, # The Number of Controls to Sample from the Full Risk Set
#'  seed = 9999)
#'
#'
#'#Computing the four-cycles statistics for the relational event sequence with
#'#the exponential weights of past events returned
#'post.processing <- dreamstats_fourcycles(data = post.processing,
#'                                         halflife = 2.592e+09)
#'
#'#printing the post-processed relational event sequence
#'post.processing
#'
#'
#'#Computing the four-cycles statistic for the relational event sequence
#'#and returning only the vector of computed values
#'cycle4.stat <-  dreamstats_fourcycles(data = post.processing,
#'                                      halflife = 2.592e+09,
#'                                      return_stats = TRUE)
#'
#'cor(cycle4.stat, post.processing$statistics$four.cycles)
#'
#'#Computing the four-cycles statistics for the relational event sequence with
#'#the counts of past events returned
#'post.processing <- dreamstats_fourcycles(data = post.processing,
#'                                         halflife = 2.592e+09,
#'                                         counts = TRUE)
#'
#'cbind(post.processing$statistics$four.cycles, cycle4.stat)
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
dreamstats_fourcycles <- function(data,
                                  halflife=2, # the half life value for the weighting function
                                  counts = FALSE, #Logical indicating if the raw counts of events should be returned or the exponential weighting function should be used (TRUE = counts; FALSE = exponential weighting)
                                  dyadic_weight= 0.00, # dyadic cutoff weight for events that no longer matter
                                  exp_weight_form = FALSE,
                                  return_stats = FALSE
) {
  if(!inherits(data,"dream_sequence")) base::stop("Error: the object for the `data` argument is not a `dream_sequence` object.")
  time <- data$processed_sequence$time # variable (column) name that contains the time variable
  sender<- data$processed_sequence$sender# variable (column) name that contains the sender variable
  receiver<- data$processed_sequence$receiver# variable (column) name that contains the target variable
  observed<- data$processed_sequence$observed# variable (column) name that contains the observed variable
  sampled<- data$processed_sequence$sampled# variable (column) name that contains the sampled variable

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
  if(!is.numeric(time)){
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
  #   Returning the values back to the user as an attachment to the class
  #
  ########################################################
  if(return_stats){
    data <- weights #returning just the vector of computed stats
  }else{
    data$statistics$four.cycles <-weights
  }
  return(data) # output the file to the user!s
}







#' @title Add Dyadic-Level Statistics for Event Dyads in a Relational Event Sequence
#' @description
#' `r lifecycle::badge("stable")`
#'
#' This function allows users to add time-varying and time-invariant
#' dyadic-level statistics that impact the dyadic event rates in
#' relational event models.

#' @name dreamstats_dyadic
#' @param data A `dream_sequence` object that contains the processed relational event sequence.
#' @param var_name A string that is the name of the variable from the `statistics` element in the `dream_sequence` argument for which the statistic should be created.
#' @param transformation The type of transformation for how the sender and receiver values will be compared (see details). The following arguments are provided: "same", "abs.diff", and "inv.abs.diff".
#' @param return_stats TRUE/FALSE. TRUE indicates that the vector of computed
#' statistics will be returned. FALSE indicates that the vector of computed
#' statistics will be added to the `statistics` element of the `data` argument.
#' @import Rcpp
#' @return The vector of dyadic-level statistics for the relational event sequence or the updated `data` argument.
#' @export
#'
#'@details This function adds user-provided time-varying and time-invariant dyadic-level
#' statistics to a `dream_sequence` object for relational event models. For this
#' function to work, in the object provided to the `data` argument, the `statistics`
#' element must contain the following previously computed statistics: `var_name.sender` and
#' `var_name.receiver`, where `var_name` is the user-provided argument. For instance,
#' if the `var_name` is `male`, then the two variables need to be included in the
#' `statistics` list: `male.sender` and `male.receiver`.
#'
#' The function allows for three types of transformations to compare the values
#' for the event senders and event receivers. When the `transformation` argument
#' is `same`, then the values are 1 if the event sender and receivers are the
#' same, and 0 if not. When the `transformation` argument
#' is `abs.diff`, then the values are \eqn{|sender - receiver|}. Finally, when
#' the `transformation` argument is `inv.abs.diff`, then the values are \eqn{1/|sender - receiver|} (if the difference is 0, the value is set to 1).
#'
#'@examples
#'events <- data.frame(time = 1:18, eventID = 1:18,
#'                     sender = c("A", "B", "C",
#'                                "A", "D", "E",
#'                                "F", "B", "A",
#'                                "F", "D", "B",
#'                                "G", "B", "D",
#'                                "H", "A", "D"),
#'                     target = c("B", "C", "D",
#'                                "E", "A", "F",
#'                               "D", "A", "C",
#'                                "G", "B", "C",
#'                                "H", "J", "A",
#'                                "F", "C", "B"))
#'
#'processed <- create_res(type = "one-mode",
#'                       ordinal = TRUE,
#'                       riskset = "constant_sample",
#'                       time = events$time,
#'                       sender = events$sender,
#'                       receiver = events$target,
#'                       p_samplingobserved = 1.00,
#'                       n_controls = 1,
#'                       seed = 9999)
#'
#'#reconstructing the data.frame object to store the time-varying
#'#actor-level statistic for receivers
#'ids <- unique(c(events$sender,events$target))
#'actor_stats <- data.frame(actor_id = ids,
#'                          time_start = 0,
#'                          time_end = 18,
#'                          male = sample(0:1, length(ids), TRUE))
#'
#'#adding the value to the post-processing relational event sequence where
#'#the new variable is named "male"
#'processed <- dreamstats_actor(data = processed,
#'                              var_name = "male",
#'                              sender_ref = TRUE,
#'                              actor_info = actor_stats)
#'processed
#'processed$statistics$male.sender
#'#adding the value for the event receivers
#'processed <- dreamstats_actor(data = processed,
#'                              var_name = "male",
#'                              sender_ref = FALSE,
#'                              actor_info = actor_stats)
#'processed
#'processed$statistics$male.receiver
#'
#'
#'#adding the dyadic same value to the post-processing relational event sequence where
#'#the new variable is 1 if the event sender and receiver have the same value
#'# and 0 if not
#'processed <- dreamstats_dyadic(data = processed,
#'                              var_name = "male",
#'                              transformation = "same")
#'processed
#'processed$statistics$male.same
#'
#'extract.data <- as.data.frame(processed)
#'extract.data #checking the values are the same!
#'
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#'
dreamstats_dyadic <- function(data, #the relational event sequence
                              var_name, #the variable name
                              transformation = c("same", "abs.diff", "inv.abs.diff"),
                              return_stats = FALSE
) {

  ########################################################
  #
  #   Checking for Errors in User Inputs
  #
  ########################################################
  if(!inherits(data,"dream_sequence")) base::stop("Error: the object for the `data` argument is not a `dream_sequence` object.")

  pre.stats <- names(data$statistics)
  #getting the event sender values
  new.var1 <- paste0(var_name, ".sender")
  if(!(new.var1 %in% pre.stats)){
    base::stop("Error: There is not a previously computed statistic for the `var_name` argument for event senders. The variable should be named `var_name.sender`. Please update this argument.")
  }
  #getting the event receiver values
  new.var2 <- paste0(var_name, ".receiver")
  if(!(new.var2 %in% pre.stats)){
    base::stop("Error: There is not a previously computed statistic for the `var_name` argument for event receivers. The variable should be named `var_name.receiver`. Please update this argument.")
  }
  sender.values <- data$statistics[[new.var1]]
  receiver.values <- data$statistics[[new.var2]]

  if(length(transformation) != 1){
    base::stop("Error: There should only be one argument specified for the `transformation` argument. Please update this argument.")
  }
  if(!(transformation %in% c("same", "abs.diff", "inv.abs.diff"))){
    base::stop("Error: The `transformation` argument should be of one of the following values: `same`, `abs.diff`, `inv.abs.diff`. Please update this argument.")
  }

  if(transformation == "same"){ #the same transformation should be applied: 1 = the values are the same, 0 if not
    weights <- ifelse(sender.values == receiver.values, 1 , 0)
    new.name <- paste0(var_name,".same", sep="")
  }

  if(transformation == "abs.diff"){ #the absdiff transformation should be applied: |sender_value - receiver_value|
    weights <- abs(sender.values-receiver.values)
    new.name <- paste0(var_name,".abs.diff", sep="")
  }

  if(transformation == "inv.abs.diff"){ #the absdiff transformation should be applied: 1/|sender_value - receiver_value|
    weights <- abs(sender.values-receiver.values)
    weights[weights == 0] <- 1 #a check to not get division by 0
    weights <- 1/weights
    new.name <- paste0(var_name,".inv.abs.diff", sep="")

  }

  ########################################################
  #
  #   Returning the values back to the user as an attachment to the class
  #
  ########################################################
  if(return_stats){
    data <- weights #returning just the vector of computed stats
  }else{
    data$statistics[[new.name]] <-weights
  }
  return(data) # output the file to the user!s
}






#' @title Add Actor-Level Statistics for Event Dyads in a Relational Event Sequence
#' @description
#' `r lifecycle::badge("stable")`
#'
#' This function allows users to add time-varying and time-invariant
#' actor-level statistics to be used in the estimation of ordinal and interval timing relational event models.
#'
#' @name dreamstats_actor
#' @param data A `dream_sequence` object that contains the processed relational event sequence.
#' @param var_name A string that is the name of the variable from the `actor_info` argument that represents the actor-level statistic to be added.
#' @param sender_ref TRUE/FALSE. TRUE indicates that the variable should be added with respect to the sender (i.e., the variable is associated with the sender). FALSE indicates that the variable is representative of the event receivers.
#' @param actor_info A \eqn{N} x 4 `data.frame` object that contains four columns with the number of observations
#' being the number of referenced actors (i.e., if sender = TRUE, then the senders in the relational event sequence). The object should contain
#' four named columns: (1) `actor_id`, which is the id that is associated with the actors, (2) `var_name` (from the argument), where
#' the var_name column is the vector of statistics, and for time-varying statistics, (3) `time_start`, a column that contains the
#' time that represents the time at which the actor has that statistic, and (4) `time_end`, a column that contains the
#' time for which the actor stops possessing that statistic value (see details for more information).
#' @param make_factor TRUE/FALSE. TRUE indicates that the vector of values will be made a factor, and FALSE if not.
#' @param return_stats TRUE/FALSE. TRUE indicates that the vector of computed
#' statistics will be returned. FALSE indicates that the vector of computed
#' statistics will be added to the `statistics` element of the `data` argument.
#' @import Rcpp
#' @return The vector of actor-level statistics for the relational event sequence or the updated `data` argument.
#' @export
#'
#'@details This function adds user-provided time-varying and time-invariant actor-level
#' statistics to a `dream_sequence` object for relational event models. The `actor_info`
#' argument should be a \eqn{N} x 4 `data.frame` object with four named columns. The first
#' column should be named `actor_id`, which represents the sender/receiver unique ids
#' based upon the observed relational event sequence. The second column should be named
#' `time_start` and represents the time in which the specific actor adopts the specific
#' value. The third column should be named `time_end` and represents the time in which the specific actor adopts the specific
#' value. For time-invariant actor-level statistics, the `time_start` value should be
#' 0 and the `time_end` value should be the time that marks the end of the relational
#' event sequence (e.g., the time of the last observed event). The last column should be named after the
#' `var_name` argument and represents the actor-level statistic for the actor i during the
#' times between `time_start` and `time_end`.
#'

#'@examples
#'events <- data.frame(time = 1:18, eventID = 1:18,
#'                     sender = c("A", "B", "C",
#'                                "A", "D", "E",
#'                                "F", "B", "A",
#'                                "F", "D", "B",
#'                                "G", "B", "D",
#'                                "H", "A", "D"),
#'                     target = c("B", "C", "D",
#'                                "E", "A", "F",
#'                               "D", "A", "C",
#'                                "G", "B", "C",
#'                                "H", "J", "A",
#'                                "F", "C", "B"))
#'
#'processed <- create_res(type = "one-mode",
#'                       ordinal = TRUE,
#'                       riskset = "constant_sample",
#'                       time = events$time,
#'                       sender = events$sender,
#'                       receiver = events$target,
#'                       p_samplingobserved = 1.00,
#'                       n_controls = 1,
#'                       seed = 9999)
#'
#'#reconstructing the data.frame object to store the time-varying
#'#actor-level statistic for sender
#'ids <- unique(c(events$sender,events$target))
#'actor_stats <- data.frame(actor_id = ids,
#'                          time_start = 0,
#'                          time_end = 18,
#'                          rv = sample(0:1, length(ids), TRUE))
#'
#'#adding the value to the post-processing relational event sequence where
#'#the new variable is named "rv"
#'processed <- dreamstats_actor(data = processed,
#'                              var_name = "rv",
#'                              sender_ref = TRUE,
#'                              actor_info = actor_stats)
#'processed
#'processed$statistics$rv.sender
#'
#'#constructing the data.frame object to store the time-invariant
#'#actor-level statistic for receivers
#'ids <- unique(c(events$sender,events$target))
#'actor_stats <- data.frame(actor_id = rep(ids,2),
#'                          time_start = c(rep(0, 9), rep(10,9)),
#'                          time_end = c(rep(9,9),rep(18,9)),
#'                          trv = rnorm(length(ids)*2))
#'
#'#adding the value to the post-processing relational event sequence where
#'#the new variable is named "trv" and the values is for the event receivers
#'processed <- dreamstats_actor(data = processed,
#'                              var_name = "trv",
#'                              sender_ref = FALSE,
#'                              actor_info = actor_stats)
#'processed
#'processed$statistics$trv.receiver
#'
#'extract.data <- as.data.frame(processed)
#'extract.data
#'
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#'
dreamstats_actor <- function(data, #the relational event sequence
                             var_name, #the variable name
                             sender_ref = TRUE,
                             actor_info,
                             make_factor = FALSE, #should the variable be made categorical (a factor in R)?
                             return_stats = FALSE
) {

  ########################################################
  #
  #   Checking for Errors in User Inputs
  #
  ########################################################
  if(!inherits(data,"dream_sequence")) base::stop("Error: the object for the `data` argument is not a `dream_sequence` object.")
  time <- data$processed_sequence$time # variable (column) name that contains the time variable
  sender<- data$processed_sequence$sender# variable (column) name that contains the sender variable
  receiver<- data$processed_sequence$receiver# variable (column) name that contains the target variable
  sampled<- data$processed_sequence$sampled# variable (column) name that contains the sampled variable
  if(sum(sampled) == 0){
    base::stop("Error: There are no sampled events based upon the 'observed' input. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(!is.numeric(time)){
    base::stop("Error: The 'time' argument is not a numeric value. Stopping computation! Please update this argument.") # stop computation and tell the user
  }
  if(!inherits(actor_info,"data.frame")){
    base::stop("Error: The 'actor_info' argument should be a `data.frame` object. Please update this argument and return.") # stop computation and tell the user
  }
  var.names <- colnames(actor_info) #extracting the column names of the actor information dataset
  need.names <- c("time_start", "time_end", var_name, "actor_id")
  #checking the inputs
  if(!all(need.names %in% var.names)){
    base::stop("Error: The `var_name` argument and/or a `time` column are missing in the `actor_info` argument. Please update this argument and return.") # stop computation and tell the user
  }
  if(!is.numeric(actor_info[,"time_start"])){
    base::stop("Error: The 'time_start' column in the `actor_info` argument is not a numeric value. Please update this argument.") # stop computation and tell the user
  }
  if(!is.numeric(actor_info[,"time_end"])){
    base::stop("Error: The 'time_end' column in the `actor_info` argument is not a numeric value. Please update this argument.") # stop computation and tell the user
  }
  info.actors <- as.character(actor_info[,"actor_id"])
  if(sender_ref){ #we should be checking the sender argument
    observed.actors <- sender[sampled==1]
  }else{ #we should be checking the target argument
    observed.actors <- receiver[sampled==1]
  }
  #checking to make sure that for all sampled actors, we have information, if not... tell the user
  time <- time[sampled == 1]
  check.information <- all(observed.actors %in% info.actors)
  if(!check.information)   base::stop("Error: There are actors in the post-processing event sequence (`data` argument) that do not have observations in the `actor_info` argument. Please update this argument and return.") # stop computation and tell the user
  #subsetting to only get those relevant actors information that are needed for the computation
  actor_info <- base::subset(actor_info, actor_info[,"actor_id"] %in% observed.actors)
  weights <- dream_nodestats(time=time,
                             time_start = actor_info[,"time_start"],
                             time_end = actor_info[,"time_end"],
                             observed_actor = observed.actors,
                             actor_id = actor_info[,"actor_id"],
                             value = actor_info[,var_name])
  if(make_factor) weights <- base::factor(weights)
  ########################################################
  #
  #   Returning the values back to the user as an attachment to the class
  #
  ########################################################
  if(return_stats){
    data <- weights #returning just the vector of computed stats
  }else{
    if(sender_ref) new.name <- paste0(var_name,".sender", sep="")
    if(!sender_ref) new.name <- paste0(var_name,".receiver", sep="")
    data$statistics[[new.name]] <-weights
  }
  return(data) # output the file to the user!s
}








#' @title A Helper Function to Assist Researchers in Finding Dyadic Weight Cutoff Values
#' @name dreamstats_dyadcut
#' @param halflife A numerical value that is the halflife value to be used in the exponential weighting function (see details section). Preset to 2 (should be updated by the user based on substantive context).
#' @param relationalWidth The numerical value that corresponds to the time range for which the user specifies for temporal relevancy.
#' @param exp_weight_form TRUE/FALSE. TRUE indicates that the Lerner et al. (2013) exponential weighting function will be used (see the details section). FALSE indicates that the Lerner and Lomi (2020) exponential weighting function will be used (see the details section). Set to FALSE by default
#' @return The dyadic weight cutoff based on user specified values.
#' @export
#'
#
#' @description
#' `r lifecycle::badge("stable")`
#'
#' A user-helper function to assist researchers in finding the dyadic
#' cutoff value to compute sufficient statistics for relational event models based upon temporal dependency.
#'
#'@details
#' This function is specifically designed as a user-helper function to assist
#' researchers in finding the dyadic cutoff value for creating sufficient statistics
#' based upon temporal dependency. In other words, this function estimates a dyadic
#' cutoff value for relational relevance, that is, the minimum dyadic weight for past
#' events to be potentially relevant (i.e., to possibly have an impact) on the current
#' event. All non-relevant events (i.e., events too distant in the past from the
#' current event to be considered relevant, that is, those below the cutoff value)
#' will have a weight of 0. This cutoff value is based upon two user-specified
#' values: the events' halflife (i..e, `halflife`) and the relationally relevant event
#' or time span (i.e., `relationalWidth`). Ideally, both the values for `halflife` and
#' `relationalWidth` would be based on the researcher’s command of the relevant
#' substantive literature. Importantly, `halflife` and `relationalWidth` must be in
#' the same units of measurement (e.g., days). If not, the function will not return
#' the correct answer.
#'
#' For example, let’s say that the user defines the `halflife` to be 15
#' days (i.e., two weeks) and the relationally relevant event or time
#' span (i.e., `relationalWidth`) to be 30 days (i.e., events that occurred
#' more than 1 month in the past are not considered relationally relevant
#' for the current event). The user would then specify `halflife` = 15 and `relationalWidth` = 30.

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
#'The task of this function is to find the weight, \eqn{ w(s, r, t)}, that corresponds to the
#'time difference provided by the user.
#'
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#' @references
#'Lerner, Jürgen and Alessandro Lomi. 2020. “Reliability of relational event
#'model estimates under sampling: How to fit a relational event model to 360
#'million dyadic events.” *Network Science* 8(1): 97-135.
#'
#'Lerner, Jürgen, Margit Bussman, Tom A.B. Snijders, and Ulrik Brandes. 2013. "
#'Modeling Frequency and Type of Interaction in Event Networks."
#'*The Corvinus Journal of Sociology and Social Policy* 4(1): 3-32.

#'@examples
#' #To replicate the example in the details section:
#' # with the Lerner et al. 2013 weighting function
#' dreamstats_dyadcut(halflife = 15,
#'                  relationalWidth = 30,
#'                  exp_weight_form = TRUE)
#'
#' # without the Lerner et al. 2013 weighting function
#' dreamstats_dyadcut(halflife = 15,
#'                  relationalWidth = 30,
#'                  exp_weight_form = FALSE)
#'
#'# A result to test the function (should come out to 0.50)
#' dreamstats_dyadcut(halflife = 30,
#'                  relationalWidth = 30,
#'                  exp_weight_form = FALSE)
#'
#'
#'# Replicating Lerner and Lomi (2020):
#'#"We set T1/2 to 30 days so that an event counts as (close to) one in the very next instant of time,
#'#it counts as 1/2 one month later, it counts as 1/4 two months after the event, and so on. To reduce
#'#the memory consumption needed to store the network of past events, we set a dyadic weight to
#'#zero if its value drops below 0.01. If a single event occurred in some dyad this would happen after
#'#6.64×T1/2, that is after more than half a year." (Lerner and Lomi 2020: 104).
#'
#'# Based upon Lerner and Lomi (2020: 104), the result should be around 0.01. Since the
#'# time values in Lerner and Lomi (2020) are in milliseconds, we have to change
#'# all measurements into milliseconds
#'dreamstats_dyadcut(halflife = (30*24*60*60*1000), #30 days in milliseconds
#'                 relationalWidth = (6.64*30*24*60*60*1000), #Based upon the paper
#'                 #using the Lerner and Lomi (2020) weighting function
#'                 exp_weight_form = FALSE)
#'
#'
dreamstats_dyadcut <- function(halflife = 2,  #the user specificed halflife
                               relationalWidth, #this is a value that measures how long the time span is
                               exp_weight_form = FALSE){ #should the Lerner et al. 2013 weighting function be used
  message("You are employing this function to find the corresponding dyadic cutoff value
          for temporal relevancy. The eventTime, relationalWidth, and halflife parameters must
          all be in the same measurement unit (e.g., hours, days).

          We hope you are providing the correct values...")
  Lerneretal_2013 <- exp_weight_form
  if(Lerneretal_2013 == FALSE){ #originally, 1 was an arugment, however, the function is not dependent upon this
    X <- 1 - relationalWidth #a corresponding time difference
    dyadCutoff <- exp((-(1 - X) * log(2)/(halflife))) #the weighting function
  }
  if(Lerneretal_2013 == TRUE){
    X <- 1 - relationalWidth #a corresponding time difference
    dyadCutoff <- exp((-(1 - X) * log(2)/(halflife))) * log(2)/(halflife)#the weighting function
  }
  return(dyadCutoff) #return the specificed value
}






