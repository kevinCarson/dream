## The Creation of Repetition Scores for Relational Event Models (Can be Used in Both one-mode and two-mode networks)
## See Butts 2008; Lerner and Lomi 2020; REM R Package (Brandenberger 2018)
## Code written by Kevin Carson (kacarson@arizona.edu) and Deigo Leal (https://www.diegoleal.info/)
## Last Updated: 07-15-24

#' @title Compute Butts' (2008) Persistence Network Statistic for Event Dyads in a Relational Event Sequence
#' @description
#' `r lifecycle::badge("stable")`
#'
#' This function computes the persistence network sufficient statistic for
#' a relational event sequence (see Butts 2008). Persistence measures the proportion of past ties sent from the event sender that went to the current event receiver.
#' Furthermore, this measure allows for persistence scores to be only
#' computed for the sampled events, while creating the weights based on the full event
#' sequence. Moreover, the function allows users to specify relational relevancy for the resulting statistic.

#' @name remstats_persistence
#' @param time The vector of event times from the post-processing event sequence.
#' @param sender The vector of event senders from the post-processing event sequence.
#' @param target The vector of event targets from the post-processing event sequence
#' @param observed A vector for the post-processing event sequence where i is equal to 1 if the dyadic event is observed and 0 if not.
#' @param sampled A vector for the post-processing event sequence where i is equal to 1 if the observed dyadic event is sampled and 0 if not.
#' @param ref_sender TRUE/FALSE. TRUE indicates that the persistence statistic will be computed in reference to the sender’s past relational history (see details section). FALSE indicates that the persistence statistic will be computed in reference to the target’s past relational history (see details section). Set to TRUE by default.
#' @param nopastEvents The numerical value that specifies what value should be given to events in which the sender has sent not past ties (i's neighborhood when sender = TRUE) or has not received any past ties (j's neighborhood when sender = FALSE). Set to NA by default.
#' @param dependency TRUE/FALSE. TRUE indicates that temporal relevancy will be modeled (see the details section). FALSE indicates that temporal relevancy will not be modeled, that is, all past events are relevant (see the details section). Set to FALSE by default.
#' @param relationalTimeSpan If dependency = TRUE, a numerical value that corresponds to the temporal span for relational relevancy, which must be the same measurement unit as the observed_time and processed_time objects. When dependency = TRUE, the relevant events are events that have occurred between current event time, *t*, and *t-relationalTimeSpan*. For example, if the time measurement is the number of days since the first event and the value for relationalTimeSpan is set to 10, then only those events which occurred in the past 10 days are included in the computation of the statistic.
#' @import Rcpp
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
#'eventSet <- create_riskset(type = "one-mode",
#'                           time = events$time,
#'                           eventID = events$eventID,
#'                           sender = events$sender,
#'                           receiver = events$target,
#'                           p_samplingobserved = 1.00,
#'                           n_controls = 6,
#'                           seed = 9999)
#'
#'#Computing the persistence statistic for the relational event sequence
#'eventSet$remstats_persistence <- remstats_persistence(
#'    time = as.numeric(eventSet$time),
#'    observed = eventSet$observed,
#'    sampled = rep(1,nrow(eventSet)),
#'    sender = eventSet$sender,
#'    target = eventSet$receiver,
#'    ref_sender = TRUE)
#'


remstats_persistence <-   function(time, # variable (column) name that contains the time variable
                                 sender,
                                 target,
                                 sampled, # variable (column) name that contains the sender variable
                                 observed, # variable (column) name that contains the receiver variable
                                 ref_sender = TRUE,
                                 nopastEvents = NA, #the value given to events where there are no past events (defaults to NA)
                                 dependency = FALSE,
                                 relationalTimeSpan = 0
) {

  #base::cat("Checking Data Structure and User Inputs.......") # outputting status to user

  ########################################################
  #
  #   Checking for Errors in User Inputs
  #
  ########################################################

  if(dependency == TRUE & is.null(relationalTimeSpan)){
    base::stop("Error: Temporal dependency was requested, however, the relationalTimeSpan is missing. Please add this and restart the function.") # stop computation and tell the user
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

  return(persistence)# return the vector of values

}
