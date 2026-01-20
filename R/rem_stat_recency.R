## Code written by Kevin Carson (kacarson@arizona.edu) and Deigo Leal (https://www.diegoleal.info/)
## Last Updated: 07-15-24

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
#' @name remstats_recency
#' @param time The vector of event times from the post-processing event sequence.
#' @param sender The vector of event senders from the post-processing event sequence.
#' @param receiver The vector of event receivers from the post-processing event sequence
#' @param observed A vector for the post-processing event sequence where i is equal to 1 if the dyadic event is observed and 0 if not.
#' @param sampled A vector for the post-processing event sequence where i is equal to 1 if the observed dyadic event is sampled and 0 if not.
#' @param type A string value that specifies which recency formula will be used to compute the statistics. The options are "raw.diff", "inv.diff.plus1", "rank.ordered.count" (see details section).
#' @param i_neighborhood TRUE/FALSE. TRUE indicates that the recency statistic will be computed in reference to the sender’s past relational history (see details section). FALSE indicates that the recency statistic will be computed in reference to the target’s past relational history (see details section). Set to TRUE by default.
#' @param nopastEvents The numerical value that specifies what value should be given to events in which the sender was not active as a sender in the past (i’s neighborhood when i_neighborhood = TRUE) or was not the recipient of a past event (j’s neighborhood when i_neighborhood = FALSE). Set to NA by default.
#' @param dependency TRUE/FALSE. TRUE indicates that temporal relevancy will be modeled (see details section). FALSE indicates that temporal relevancy will not be modeled, that is, all past events are relevant (see details section). Set to FALSE by default.
#' @param relationalTimeSpan If dependency = TRUE, a numerical value that corresponds to the temporal span for relational relevancy, which must be the same measurement unit as the observed_time and processed_time objects. When dependency = TRUE, the relevant events are events that have occurred between current event time, *t*, and *t - relationalTimeSpan*. For example, if the time measurement is the number of days since the first event and the value for relationalTimeSpan is set to 10, then only those events which occurred in the past 10 days are included in the computation of the statistic.
#' @import Rcpp
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
#'eventSet <- create_riskset(type = "one-mode",
#'                           time = events$time,
#'                           eventID = events$eventID,
#'                           sender = events$sender,
#'                           receiver = events$target,
#'                           p_samplingobserved = 1.00,
#'                          n_controls = 6,
#'                          seed = 9999)
#'
#'#Computing the recency statistics (with raw time difference) for the relational event sequence
#'eventSet$recency_rawdiff <- remstats_recency(
#'    time = as.numeric(eventSet$time),
#'    observed = eventSet$observed,
#'    sampled = rep(1,nrow(eventSet)),
#'    sender = eventSet$sender,
#'    receiver = eventSet$receiver,
#'    type = "raw.diff")
#'
#'#Computing the recency statistics (with inverse of time difference) for the
#'#relational event sequence
#'eventSet$recency_rawdiff <- remstats_recency(
#'    time = as.numeric(eventSet$time),
#'    observed = eventSet$observed,
#'    sampled = rep(1,nrow(eventSet)),
#'    sender = eventSet$sender,
#'    receiver = eventSet$receiver,
#'    type = "inv.diff.plus1")
#'
#'#Computing the rank-based recency statistics for the relational event sequence
#'eventSet$recency_rawdiff <- remstats_recency(
#'    time = as.numeric(eventSet$time),
#'    observed = eventSet$observed,
#'    sampled = rep(1,nrow(eventSet)),
#'    sender = eventSet$sender,
#'    receiver = eventSet$receiver,
#'    type = "rank.ordered.count")
#'

remstats_recency <-   function(   time,
                                sender,
                                receiver,
                                sampled,
                                observed,
                                type = c("raw.diff", "inv.diff.plus1", "rank.ordered.count"),
                                i_neighborhood = TRUE, #should the recency be computed on the i's neighborhood or j's neighborhood
                                dependency = FALSE, #Boolean for temporal dependency
                                relationalTimeSpan = NULL, #if dependency == TRUE, this should specific the associated temporal time span
                                nopastEvents = NA #the value given to events where there are no past events (defaults to NA)
 ) {

  ########################################################
  #
  #   Checking for Errors in User Inputs
  #
  ########################################################
  if(dependency == TRUE & is.null(relationalTimeSpan)){
    base::stop("Error: Temporal dependency was requested, however, the relationalTimeSpan is missing. Please add this and restart the function.") # stop computation and tell the user
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
  return(weights)# return the vector of values
}

