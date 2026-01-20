## The Creation of Dynamic Risk Sets
## Code written by Kevin Carson (https://kevincarson.github.io/) and Deigo Leal (https://www.diegoleal.info/)
## Last Updated: 01-09-2025

#' @title Process and Create Risk Sets for a One- and Two-Mode Relational Event Sequences
#' @name create_riskset
#' @param type "two-mode" indicates that this is a two-mode event sequence. "one-mode" indicates that the event sequence is one-mode.
#' @param time The vector of event time values from the observed event sequence.
#' @param sender The vector of event senders from the observed event sequence.
#' @param receiver The vector of event receivers from the observed event sequence.
#' @param eventID The vector of event IDs from the observed event sequence (typically a numerical event sequence that goes from 1 to *n*).
#' @param p_samplingobserved The numerical value for the probability of selection for sampling from the observed event sequence. Set to 1 by default indicating that all observed events from the event sequence will be included in the post-processing event sequence.
#' @param n_controls The numerical value for the number of null event controls for each (sampled) observed event.
#' @param combine TRUE/FALSE. TRUE indicates that the post-sampling (processing) event sequence should be merged with the pre-processing dataset. FALSE only returns the post-processing event sequence (that is, only the sampled events).
#' @param seed The random number seed for user replication.
#' @import Rcpp
#' @importFrom data.table rbindlist
#' @importFrom data.table data.table
#' @return A post-processing data.table object with the following columns:
#' \itemize{
#'   \item \code{time} - The event time for the sampled and observed events.
#'   \item \code{eventID} - The numerical event sequence ID for the sampled and observed events.
#'   \item \code{sender} - The event senders of the sampled and observed events.
#'   \item \code{receiver} - The event targets (receivers) of the sampled and observed events.
#'   \item \code{observed} - Boolean indicating if the event is an observed or control event. (1 = observed; 0 = control)
#'   \item \code{sampled} - Boolean indicating if the event is sampled or not sampled. (1 = sampled; 0 = not sampled)
#' }
#' @export
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#'
#' This function creates one- and two-mode post-sampling eventset with options for case-control
#' sampling (Vu et al. 2015) and sampling from the observed event sequence (Lerner and Lomi 2020). Case-control
#' sampling samples an arbitrary *m* number of controls from the risk set for any event
#' (Vu et al. 2015). Lerner and Lomi (2020) proposed sampling from the observed event sequence
#' where observed events are sampled with probability *p*. Importantly, this function generates risk sets
#' that assume that the risk set for each event is fixed across all time points, that is, all actors active
#' at any time point across the event sequence are in the set of potential events. Users interested in
#' generating time-/event-varying risks sets should consult the \code{\link[dream]{processOMEventSeq}} function
#' for one-mode event sequences and the \code{\link[dream]{processTMEventSeq}} function for two-mode event
#' sequences. Future versions of the `dream` package will incorporate this option into this function in a
#' principled manner.
#'
#' @details This function processes observed events from the set \eqn{E}, where each event \eqn{e_i} is
#' defined as:
#' \deqn{e_{i} \in E = (s_i, r_i, t_i, G[E;t])}
#' where:
#' \itemize{
#'   \item \eqn{s_i} is the sender of the event.
#'   \item \eqn{r_i} is the receiver of the event.
#'   \item \eqn{t_i} represents the time of the event.
#'   \item \eqn{G[E;t] = \{e_1, e_2, \ldots, e_{t'} \mid t' < t\}} is the network of past events, that is, all events that occurred prior to the current event, \eqn{e_i}.
#' }
#'
#' Following Butts (2008) and Butts and Marcum (2017), for one-mode event sequences, the risk (support)
#' set is defined as all possible  events at time \eqn{t}, \eqn{A_t}, as the full Cartesian
#' product of prior senders and receivers in the set \eqn{G[E;t]} that could have
#' occurred at time \eqn{t}. Formally:
#' \deqn{A_t = \{ (s, r) \mid s \in S \times r \in R\}}
#' where \eqn{S} is the set of potential event senders and \eqn{R} is the set of potential event receivers. In this function,
#' the full risk set is considered fixed across all time points.
#'
#' For two-mode event sequences, the risk (support) set is defined as all possible
#' events at time \eqn{t}, \eqn{A_t}, as the cross product of two disjoint sets, namely, prior senders and receivers,
#' in the set \eqn{G[E;t]} that could have occurred at time \eqn{t}. Formally:
#' \deqn{A_t = \{ (s, r) \mid s \in S \times r \in R\}}
#' where \eqn{S} is the set of potential event senders and \eqn{R} is the set of potential event receivers. In this function,
#' the full risk set is considered fixed across all time points.
#'
#' Case-control sampling maintains the full set of observed events, that is, all events in \eqn{E}, and
#' samples an arbitrary number \eqn{m} of non-events from the support set \eqn{A_t} (Vu et al. 2015; Lerner
#' and Lomi 2020). This process generates a new support set, \eqn{SA_t}, for any relational event
#' \eqn{e_i} contained in \eqn{E} given a network of past events \eqn{G[E;t]}. \eqn{SA_t} is formally defined as:
#' \deqn{SA_t \subseteq \{ (s, r) \mid s \in S \times r \in R \}}
#' and in the process of sampling from the observed events, \eqn{n} number of observed events are
#' sampled from the set \eqn{E} with known probability \eqn{0 < p \le 1}. More formally, sampling from
#' the observed set generates a new set \eqn{SE \subseteq E}.
#'
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#' @references
#' Butts, Carter T. 2008. "A Relational Event Framework for Social Action." *Sociological Methodology* 38(1): 155-200.
#'
#' Butts, Carter T. and Christopher Steven Marcum. 2017. "A Relational Event Approach to Modeling Behavioral Dynamics." In A.
#' Pilny & M. S. Poole (Eds.), *Group processes: Data-driven computational approaches*. Springer International Publishing.
#'
#' Lerner, Jürgen and Alessandro Lomi. 2020. "Reliability of relational event model estimates under sampling: How to
#' fit a relational event model to 360 million dyadic events." *Network Science* 8(1): 97–135.
#'
#' Vu, Duy, Philippa Pattison, and Garry Robins. 2015. "Relational event models for social learning in MOOCs." *Social Networks* 43: 121-135.
#' @examples
#'
#' data("WikiEvent2018.first100k")
#' WikiEvent2018.first100k$time <- as.numeric(WikiEvent2018.first100k$time)
#' ### Creating the EventSet By Employing Case-Control Sampling With M = 10 and
#' ### Sampling from the Observed Event Sequence with P = 0.01
#' EventSet <- create_riskset(
#'   type = "two-mode",
#'   time = WikiEvent2018.first100k$time, # The Time Variable
#'   eventID = WikiEvent2018.first100k$eventID, # The Event Sequence Variable
#'   sender = WikiEvent2018.first100k$user, # The Sender Variable
#'   receiver = WikiEvent2018.first100k$article, # The Receiver Variable
#'   p_samplingobserved = 0.01, # The Probability of Selection
#'   n_controls = 10, # The Number of Controls to Sample from the Full Risk Set
#'   seed = 9999) # The Seed for Replication
#'
#'
#' ### Creating A New EventSet with more observed events and less control events
#' ### Sampling from the Observed Event Sequence with P = 0.02
#' ### Employing Case-Control Sampling With M = 2
#' EventSet1 <- create_riskset(
#'   type = "two-mode",
#'   time = WikiEvent2018.first100k$time, # The Time Variable
#'   eventID = WikiEvent2018.first100k$eventID, # The Event Sequence Variable
#'   sender = WikiEvent2018.first100k$user, # The Sender Variable
#'   receiver = WikiEvent2018.first100k$article, # The Receiver Variable
#'   p_samplingobserved = 0.02, # The Probability of Selection
#'   n_controls = 2, # The Number of Controls to Sample from the Full Risk Set
#'   seed = 9999) # The Seed for Replication
#'

create_riskset <-  function(type = c("two-mode", "one-mode"), #the type of risk set to be created
                           time, # variable (column) name that contains the time variable
                           eventID, # variable (column) name that contains event Sequence ID
                           sender, # variable (column) name that contains the sender variable
                           receiver, # variable (column) name that contains the receiver variable
                           p_samplingobserved = 1, # probability of selection for case control sampling
                           n_controls, # number of controls for each selected event
                           combine = TRUE,
                           seed = 9999) { # seed for replication (user can change this value)

  #base::cat("Checking Data Structure and User Inputs.......") # outputting status to user

  ########################################################
  #
  #   Checking for Errors in User Inputs
  #
  #######################################################
  if (p_samplingobserved > 1 | p_samplingobserved < 0) { # if the probability is not a probability (i.e., bounded by 0 and 1)
    base::stop("Error: Probabilty of Selection Must be within the interval: 0 < p <= 1. We hope you know what you're doing. Happy computing!") # stop computation and tell the user
  }
  if (!(n_controls > 0)) { # if number of controls is equal to 0, that is, no null events
    base::stop("Error: Number of Controls Must Be At Least 1. We hope you know what you're doing. Happy computing!") # stop computation and tell the user
  }
  if(!(type %in% c("two-mode", "one-mode"))){
    base::stop("Error: The type argument is not valid. Please see the help page and retry! Happy computing!") # stop computation and tell the user
  }
  sender <- as.character(sender)
  receiver <- as.character(receiver)
  ########################################################
  #
  #   Creating the dataset in c++
  #
  #######################################################
  if(type == "two-mode"){ #if it is a two-mode event sequence
      #then create a two-mode event sequence!
      bigeventlist <- processREMseqTM(time = time,
                                      seqid = eventID,
                                      sender = (sender),
                                      target = (receiver),
                                      pobserved = p_samplingobserved,
                                      ncontrols = n_controls,
                                      rseed = seed)
  }else{#then create a one-mode event sequence!
      bigeventlist <- processREMseqOM(time = time,
                                      seqid = eventID,
                                      sender = (sender),
                                      target = (receiver),
                                      pobserved = p_samplingobserved,
                                      ncontrols = n_controls,
                                      rseed = seed)
  }
  #cleaning the post-processing dataset!
  eventSeq <- data.table::rbindlist(bigeventlist) # merging everything into a nice dataframe to be exported to the user!
  eventSeq$sampled <- 1
  colnames(eventSeq) <- c("time", "eventID", "sender", "receiver", "observed", "sampled")

  if(combine){ #if the user wants the complete dataframe returned!
    #creating the data.table object to be all non-sampled observed events!
    data <- data.table::data.table(time = time[!(eventID %in% eventSeq$eventID)], #the non-sampled event times
                                   eventID = eventID[!(eventID %in% eventSeq$eventID)], #the non-sampled event senders
                                   sender = sender[!(eventID %in% eventSeq$eventID)],#the non-sampled event senders
                                   receiver = receiver[!(eventID %in% eventSeq$eventID)],#the non-sampled event targets
                                   observed = 1,#they are observed
                                   sampled = 0)#they are not sampled
    eventSeq <- data.table::rbindlist(list(eventSeq,data)) #combining the objects together!
    eventSeq <- eventSeq[order(eventSeq$time)] #temporal (ordering) sorting the post-processing event sequence by time!
  }
  return(eventSeq) # output the file to the user!
}

