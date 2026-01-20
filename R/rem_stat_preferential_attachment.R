## Code written by Kevin Carson (kacarson@arizona.edu) and Deigo Leal (https://www.diegoleal.info/)
## Last Updated: 07-15-24

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

#' @name remstats_prefattachment
#' @param time The vector of event times from the post-processing event sequence.
#' @param sender The vector of event senders from the post-processing event sequence.
#' @param receiver The vector of event receivers from the post-processing event sequence
#' @param observed A vector for the post-processing event sequence where i is equal to 1 if the dyadic event is observed and 0 if not.
#' @param sampled A vector for the post-processing event sequence where i is equal to 1 if the observed dyadic event is sampled and 0 if not.
#' @param dependency TRUE/FALSE. TRUE indicates that temporal relevancy will be modeled (see the details section). FALSE indicates that temporal relevancy will not be modeled, that is, all past events are relevant (see the details section). Set to FALSE by default.
#' @param relationalTimeSpan If dependency = TRUE, a numerical value that corresponds to the temporal span for relational relevancy, which must be the same measurement unit as the observed_time and processed_time objects. When dependency = TRUE, the relevant events are events that have occurred between current event time, *t*, and *t - relationalTimeSpan*. For example, if the time measurement is the number of days since the first event and the value for relationalTimeSpan is set to 10, then only those events which occurred in the past 10 days are included in the computation of the statistic.
#' @import Rcpp
#' @return The vector of event preferential attachment statistics for the relational event sequence.
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
#'eventSet <- create_riskset( type = "one-mode",
#'                           time = events$time,
#'                           eventID = events$eventID,
#'                           sender = events$sender,
#'                           receiver = events$target,
#'                           p_samplingobserved = 1.00,
#'                           n_controls = 6,
#'                           seed = 9999)
#'
#'#Computing the preferential attachment statistic for the relational event sequence
#'eventSet$pref <- remstats_prefattachment(
#'    time = as.numeric(eventSet$time),
#'    observed = eventSet$observed,
#'    sampled = rep(1,nrow(eventSet)),
#'    sender = eventSet$sender,
#'    receiver = eventSet$receiver)
#'


remstats_prefattachment <-   function(time,
                                sampled,
                                observed,
                                sender,
                                receiver,
                                dependency = FALSE,
                                relationalTimeSpan =0 # the sizes of the windows that we will use, if NA, we will compute it internally
) {

  ########################################################
  #
  #   Checking for Errors in User Inputs
  #
  ########################################################
  if(dependency == TRUE & is.null(relationalTimeSpan)){
    base::stop("Error: Temporal dependency was requested, however, the relationalTimeSpan is missing. Please add this and restart the function.") # stop computation and tell the user
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
  return(weights) #returning the computed scores to the user
}


