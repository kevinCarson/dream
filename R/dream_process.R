#' @title Process One- and Two-Mode Relational Event Sequences and Create Post-Processing Relational Event Sequences
#' @name create_res
#' @param time The vector of event time values from the observed event sequence,
#' where the jth entry is the relative time at which the jth event occurred. The event times
#' should be relative to the onset (start) of the relational event sequence.
#' @param sender The vector of event senders from the observed event sequence, where
#' the jth entry is the event sender for the jth observed/realized event.
#' @param receiver The vector of event receivers from the observed event sequence where
#' the jth entry is the event receiver for the jth observed/realized event.
#' @param type "two-mode" indicates that this is a two-mode event sequence (i.e., observed actors
#' can only be either event senders or event receivers). The option "one-mode"
#' indicates that the observed event sequence is one-mode (i.e., observed actors
#' can be event senders and receivers)
#' @param ordinal TRUE/FALSE. TRUE indicates that observed timing of the events is ordinal (and the ordinal timing likelihood function will be used). FALSE denotes
#' that the observed timing is observed, relative to the start of the event sequence (and the interval timing likelihood function will be used). The interval
#' timing option adds the right-censored events to the post-processed relational event sequence (i.e., the set of (sampled) controls events for the
#' time point t, that marks the end of the relational event sequence.) Please see the references for more information.
#' @param t If ordinal is set to FALSE, the time that marks the end of the relational event sequence, relative to the start of the event
#' sequence. If t is NULL and ordinal is set to FALSE, then the right-censored events are not added to the post-processing event sequence and
#' the \code{\link{estimate_rem}} function will not add the right-censoring term to the log-likelihood. The log-likelihood thus
#' reduces to that used in the MLE interval timing relational event model in the \code{\link[remstimate]{remstimate}} function. The default
#' value is NULL.
#' @param riskset The argument should be one of the following strings: "complete", "constant_sample", "dynamic_sample",
#' "actor_varying", "actor_varying_sample". "complete" will create a post-processing relational event
#' sequence where the set of null events for each (sampled) realized/observed event is the full
#' set of actors that were active at anytime in the event sequence (for one-mode sequences, this is the
#' full Cartesian plot of actors, and for two-mode sequences, this is the full cross-product of the
#' event sender and receiver sets). "constant_sample" will create a post-processing relational event
#' sequence where the set of null events for each (sampled) realized/observed event is a random sample from the full
#' set of actors that were active at anytime in the event sequence (for one-mode sequences, this is the
#' full Cartesian plot of actors, and for two-mode sequences, this is the full cross-product of the
#' event sender and receiver sets), where the number of sampled events is dependent upon the `n_controls`
#' argument. "dynamic_sample" will create a post-processing relational event
#' sequence where the set of null events for each (sampled) realized/observed event at time t is a random sample from the full
#' set of actors that have been active up to and including t. "actor_varying" will create a post-processing relational event
#' sequence where the set of null events for each (sampled) realized/observed event at time t is the full
#' set of actors that are considered relationally active at time t dependent upon the `active_times`
#' argument, whereas the "actor_varying_sample" option will create a post-processing relational event
#' sequence where the set of null events for each (sampled) realized/observed event is a random
#' sample from the set of relationally active actors.
#' @param p_samplingobserved The numerical value for the probability of selection for sampling
#' from the observed event sequence. Set to 1 by default indicating that all observed
#' events from the event sequence will be included in the post-processing event sequence.
#' @param n_controls The numerical value for the number of null event controls for
#' each (sampled) observed event. This argument should be specified when one of the
#' following `riskset` types is provided: "constant_sample", "dynamic_sample",
#' and "actor_varying_sample".
#' @param active_times This argument is either a `data.frame` object for "one-mode" relational event
#' sequences or a `list` that contains 2 `data.frame` objects when `type` is set to "two-mode". For
#' "one-mode" event sequence types, the `data.frame` object must have three named
#' columns: "actor_id", "time_start", "time_end", and the number of rows (observations)
#' is the number of actor active spells, where an actor spell is defined as the
#' time in which the ith actor becomes relationally active (the ith entry for the "time_start"
#' column) and then becomes relationally inactive (the ith entry for the "time_end"
#' column). The "time_start" and "time_end" vectors should be relative
#' to the start of the relational event sequence (commonly set to start at time 0). Importantly, a single
#' actor (denoted by the "actor_id" column) can have  multiple active spell rows, as actors may exit and re-enter the relational event
#' sequence. In comparison, actors who do not leave the relational event sequence and always have the
#' capacity to become relationally active should have one row (entry) where the
#' "time_start" value is set to 0 and the "time_end" value should be set to the time point
#' that marks the end of the relational event sequence (t for interval event sequences
#' and the time of the last observed event for ordinal event sequences). Finally,
#' for "two-mode" event sequences, the `active_times` argument should be a list of
#' 2 `data.frame` objects where the list elements are named "senders" and "receivers", where
#' "senders" is a `data.frame` object that follows the same design as the "one-mode"
#' arguments and defines how event senders become relationally active and inactive. Similarly,
#' the second element should be named "receivers" and is a `data.frame` object that follows the above
#' stated design and defines how event receivers become relationally active and inactive.
#' @param seed The random number seed for user replication. This argument is set to
#' NULL be default.
#' @import Rcpp
#' @importFrom data.table rbindlist
#' @importFrom data.table data.table
#' @return An object of class `dream_sequence` that contains a list of the following elements:
#' \itemize{
#'   \item \code{processed_sequence} - A `data.table` object that contains the post-processing relational event sequence with
#'   the following columns: "time", "eventID", "sender", "receiver", "sampled", and "observed". The "time" column
#'   is the vector of event times for the realized and control events. The "eventID" column represents the
#'   order that the event occurred in the relational event sequence. The "sender" and "receiver" columns
#'   are the specific dyad for that row. The "observed" vector takes a value of 1 if the
#'   dyadic pair is the observed dyad at the specific time. The "sampled" vector takes a value
#'   of 1 if the dyad was sampled at that event time and 0 if not (relevant for case-control sampling and sampling from the observed event sequence).
#'   \item \code{ordinal} - Based upon the user's input. (`ordinal` and `interval`)
#'   \item \code{t} - Based upon the user's input.
#'   \item \code{riskset} - Based upon the user's input.
#'   \item \code{p} - The probability of sampling from the observed event sequence. Based upon the user's input.
#'   \item \code{m} - The number of null event controls for each (sampled) observed event. Based upon the user's input.
#'   \item \code{type} - Based upon the user's input. (`two-mode` and `one-mode`)
#'   \item \code{n} - The number of observed events.
#'   \item \code{sampled_events} - The number of sampled observed events.
#'   \item \code{null} - The number of sampled null (control) events.
#'   \item \code{statistics} - An empty list to store the future computed relational event network statistics.
#'   \item \code{interevent.times} - The vector of interevent times (the time difference between observed events).
#' }
#' @export
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' This function creates one- and two-mode post-processing (sampled) relational event
#' sequences with options for case-control sampling (Vu et al. 2015; Butts 2008), sampling
#' from the observed event sequence (Lerner and Lomi 2020), dynamic time-varying supports
#' sets, actor-varying supports sets, and the full Cartesian product for one-mode
#' sequences and the full cross-product for two-mode event sequences (Butts 2008). The created post-processing
#' relational event sequences are designed to be modeled by relational event models (Butts 2008).
#' Case-control sampling samples an arbitrary *m* number of controls from the risk set for any event
#' (Vu et al. 2015; Butts 2008). Lerner and Lomi (2020) proposed sampling from the observed event sequence
#' where observed events are sampled with probability *p*. Importantly, this function allows
#' users to generate post-processing relational event sequences for *ordinal* and *interval* relational
#' event model likelihoods. Lastly, the post-processing relational event sequence
#' is a `dream_sequence` object that is the required object for this package's
#' functions to compute exogenous and endogenous network statistics, alongside
#' the function to estimate Maximum Likelihood relational event models.
#'
#'
#'
#' @details This function processes observed events from the set \eqn{A_t}, where each event \eqn{a_i} is
#' defined as:
#' \deqn{a_{i} \in A_t = (s_i, r_i, \tau_i, G[A_t; \tau_i])}
#' where:
#' \itemize{
#'   \item \eqn{s_i} is the sender of the event.
#'   \item \eqn{r_i} is the receiver of the event.
#'   \item \eqn{\tau_i} represents the time of the event.
#'   \item \eqn{G[A_t; \tau_i] = \{a_1, a_2, \ldots, a_{t'} \mid t' < \tau_i\}} is the network of past events, that is, all events that occurred prior to the current event, \eqn{a_i}.
#' }
#'
#'
#' For the post-processing event sequences where the `ordinal` argument is set to FALSE, the
#' last set of processed events, marked with the time point \eqn{t}, represent the
#' right-censoring events. The function generates post-processing relational event
#' sequences across three axises: (1) the inclusion of sampling from the observed/realized
#' relational event sequence, \eqn{A_t}, (2) one-mode vs. two-mode event types,
#' where the relevant actors can be either senders or receivers (in the case of two-mode)
#' sequences, or where the relevant actors can be both (in the case of one-mode
#' sequences), and (3) how the processed support set for each event should be
#' constructed. The third axis is based upon the `riskset` argument, which is
#' one of the following: "complete", "constant_sample", "dynamic_sample",
#' "actor_varying", "actor_varying_sample". "complete" will create a post-processing relational event
#' sequence where the set of null events for each (sampled) realized/observed event is the full
#' set of actors that were active at anytime in the event sequence (for one-mode sequences, this is the
#' full Cartesian plot of actors, and for two-mode sequences, this is the full cross-product of the
#' event sender and receiver sets). "constant_sample" will create a post-processing relational event
#' sequence where the set of null events for each (sampled) realized/observed event is a random sample from the full
#' set of actors that were active at anytime in the event sequence (for one-mode sequences, this is the
#' full Cartesian plot of actors, and for two-mode sequences, this is the full cross-product of the
#' event sender and receiver sets), where the number of sampled events is dependent upon the `n_controls`
#' argument. "dynamic_sample" will create a post-processing relational event
#' sequence where the set of null events for each (sampled) realized/observed event at time t is a random sample from the full
#' set of actors that have been active up to and including t. "actor_varying" will create a post-processing relational event
#' sequence where the set of null events for each (sampled) realized/observed event at time t is the full
#' set of actors that are considered relationally active at time t dependent upon the `active_times`
#' argument, whereas the "actor_varying_sample" option will create a post-processing relational event
#' sequence where the set of null events for each (sampled) realized/observed event is a random
#' sample from the set of relationally active actors.
#'
#'
#' **Complete Risk Sets**:
#'
#' Following Butts (2008) and Butts and Marcum (2017), for one-mode
#' event sequences (`type` = "one-mode"), the risk (support) set is defined as all possible
#' events at time \eqn{t}, \eqn{M_t}, as the full Cartesian
#' product of actors active in the relational event sequence and is the defined as the
#' set \eqn{N}.  Formally:
#' \deqn{M_t = \{ (s, r) \mid s \in N \times r \in N\}}
#' where \eqn{N} is the set of all possible actors in the sequence. In this function,
#' the full risk set is considered fixed (constant) across all time points.
#'
#'
#' For two-mode event sequences (`type` = "two-mode"), the risk (support) set is defined as all possible
#' events at time \eqn{t}, \eqn{M_t}, as the cross product of two disjoint sets, namely, event senders (i.e., \eqn{S}) and
#' event receivers (i.e., \eqn{R}). Formally:
#' \deqn{M_t = \{ (s, r) \mid s \in S \times r \in R\}}
#' where \eqn{S} is the set of potential event senders and \eqn{R} is the set of potential event receivers. In this function,
#' the full risk set is considered fixed across all time points.
#'
#' **Constant Sample Risk Sets**:
#'
#' Following Butts (2008), Vu et al. (2015), and Lerner and Lomi (2020), case-control sampling
#' samples an arbitrary number \eqn{m} of non-events from the above risk/support set definitions \eqn{M_t}. This
#' process generates a new support set, \eqn{\tilde{M}_t}, for any relational event
#' \eqn{a_i} contained in \eqn{A_t}. \eqn{\tilde{M}_t}, for one-mode
#' relational event sequences, is formally defined as:
#' \deqn{\tilde{M}_t \subseteq \{ (s, r) \mid s \in N \times r \in N \}}
#'
#' **Dynamic Sample Risk Sets**:
#'
#' For dynamic risk sets, for one-mode event sequences (`type` = "one-mode"), the risk (support)
#' set at time \eqn{t}, that is, \eqn{M_t}, is defined as a sample of \eqn{m} dyads from the full Cartesian product
#' of all past actors who have been involved in a relational event at or before time \eqn{t}.
#' Formally:
#' \deqn{\tilde{M}_t \subseteq \{ (s, r) \mid s \in N_t \times r \in N_t\}}
#' where \eqn{N_t} is the set of potential event senders and targets at and time \eqn{t}. The
#' definition follows the same as above for two-mode event sequences, where the sets
#' are now defined as \eqn{S_t} and \eqn{R_t}.
#'
#'
#' **Actor-varying Risk Sets**:
#'
#' Actor-varying support sets allows for actors to enter, exit, and re-enter
#' the relational event sequence as time progresses. For one-mode event sequences,
#' the node set \eqn{Y_t} is defined as those actors who are considered
#' relationally active at time \eqn{t}. For two-mode event sequences,
#' the node sets \eqn{S_t} and \eqn{R_t} are defined, respectively, as the
#' senders who can be active at time \eqn{t} and the receivers who can be
#' active at time \eqn{t}. For one-mode sequences, the formal definition is:
#' \deqn{V_t = \{ (s, r) \mid s \in Y_t \times r \in Y_t\}}
#'
#'
#' **Actor-varying Sampling Risk Sets**:
#'
#' Based upon the formal definition above, sampling from the full actor-varying
#' risk set generates a new support set definition:
#' \deqn{\tilde{V_t} \subseteq \{ (s, r) \mid s \in Y_t \times r \in Y_t\}}
#' where \eqn{\tilde{V_t}} represents the \eqn{m} number of dyads sampled
#' (with equal probability) from the set \eqn{V_t}.
#'
#'
#' @author Kevin A. Carson <kacarson@arizona.edu> and Diego F. Leal <dflc@arizona.edu>
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
#'
#' @examples
#'
#'# Generating a psuedo one-mode relational event sequence
#'set.seed(9999)
#'events <- data.frame(time = sort(rexp(18)),
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
#'# Creating the full one-mode relational risk set with p = 1.00 (all true events)
#'# based upon the ordinal timing relational event framework
#'full.process <- create_res(ordinal = TRUE,
#'                       type = "one-mode",
#'                       riskset = "complete",
#'                       time = events$time,
#'                       sender = events$sender,
#'                       receiver = events$target,
#'                       p_samplingobserved = 1.00,
#'                       seed = 9999)
#'
#'# Creating a fixed one-mode relational risk set with p = 1.00 (all true events)
#'# and 5 controls based upon the ordinal timing relational event framework
#'sample.process <- create_res(ordinal = TRUE,
#'                             type = "one-mode",
#'                             riskset = "constant_sample",
#'                             time = events$time,
#'                             sender = events$sender,
#'                             receiver = events$target,
#'                             p_samplingobserved = 1.00,
#'                             n_controls = 10,
#'                             seed = 9999)
#'
#'# Creating a dynamic one-mode relational risk set with p = 1.00 (all true events)
#'# and 5 controls based upon the interval timing relational event framework
#'dynamic.process <- create_res(ordinal = FALSE,
#'                       t = max(events$time) + rexp(1),
#'                       type = "one-mode",
#'                       riskset = "dynamic_sample",
#'                       time = events$time,
#'                       sender = events$sender,
#'                       receiver = events$target,
#'                       p_samplingobserved = 1.00,
#'                       n_controls = 5,
#'                       seed = 9999)
#'
#'
#'# Creating a actor-varying one-mode relational event sequence where actors
#'# enter, exit, and re-enter the event sequence dependent upon the user-specified
#'# active times. Each row contains the actor id, the time for which, in a specific relevant
#'#spell, become active (enter the sequence) and become inactive (exit the sequence)
#'#for actors who are always relevant, the active time ranges from 0 to the
#'#end of the sequence
#'t <- max(events$time) + rexp(1)
#'active.times <- data.frame(actor_id = c("A", "B", "B", "C", "D", "E", "F", "G",
#'                                        "H" ,"J"),
#'                          time_start = c(0,0,0.65,0,0,0.45,0.60,0,0,0),
#'                          time_end = c(t,0.45,t,t,t,1.00,t,t,t,t))
#'
#'actor.varying.process <- create_res(ordinal = TRUE,
#'                                    type= "one-mode",
#'                                    riskset = "actor_varying",
#'                                    time = events$time,
#'                                    sender = events$sender,
#'                                    receiver = events$target,
#'                                    active_times = active.times)
#'
#'
create_res <-  function(type = c("two-mode", "one-mode"), #the type of risk set to be created
                        ordinal = TRUE, #the ordinal timing rem (in the interval rem, if the ending time is known, create an additional risk set point with only controls)
                        t = NULL, #the time for the end of the relational sequence (if different from the last event)
                        time, # variable (column) name that contains the time variable
                        sender, # variable (column) name that contains the sender variable
                        receiver, # variable (column) name that contains the receiver variable
                        riskset = c("complete", "constant_sample", "dynamic_sample", "actor_varying", "actor_varying_sample"),
                        p_samplingobserved = 1, # probability of selection for case control sampling
                        n_controls=NULL, # number of controls for each selected event
                        active_times = NULL,
                        seed = NULL) { # seed for replication (user can change this value)

  ########################################################
  #
  #   Checking for Errors in User Inputs
  #
  #######################################################
  if (p_samplingobserved > 1 | p_samplingobserved < 0) { # if the probability is not a probability (i.e., bounded by 0 and 1)
    base::stop("Error: Probabilty of Selection Must be within the interval: 0 < p <= 1. We hope you know what you're doing. Happy dreaming!") # stop computation and tell the user
  }
  if(!is.null(n_controls)){
        if (!(n_controls > 0)) { # if number of controls is equal to 0, that is, no null events
          base::stop("Error: Number of Controls Must Be At Least 1. We hope you know what you're doing. Happy dreaming!") # stop computation and tell the user
        }
    }
  if(!(type %in% c("two-mode", "one-mode"))){
    base::stop("Error: The type argument is not valid. Please see the help page and retry! Happy dreaming!") # stop computation and tell the user
  }
  if(!(riskset %in% c("complete", "constant_sample", "dynamic_sample", "actor_varying", "actor_varying_sample"))){
    base::stop("Error: The `riskset` argument is not valid. Please see the help page and retry! Happy dreaming!") # stop computation and tell the user
  }
  if(c(riskset %in% c("dynamic_sample", "constant_sample", "actor_varying_sample")) & is.null(n_controls)){
    base::stop("Error: The `n_controls` argument was not specified and one of the case-control sampling riskset formulations was specified. Please see the help page and retry! Happy dreaming!") # stop computation and tell the user
  }

  if(riskset %in% c("actor_varying", "actor_varying_sample")){
    if(type=="one-mode"){
      if(!inherits(active_times, "data.frame")){
        base::stop("Error: A one-mode actor-varying riskset was requested, however, the `active_times` argument is not a data.frame. Please see the help page and retry! Happy dreaming!") # stop computation and tell the user
      }
      x <- colnames(active_times)
      need <- c("actor_id", "time_start", "time_end")
      if(!all(need %in% x)){
        base::stop("Error: An actor-varying riskset was requested The following column names must be included: actor_id, time_start, and time_end. Please see the help page and retry! Happy dreaming!") # stop computation and tell the user
      }
      actors <- active_times[,"actor_id"] #the actor ids for relevant actors
      if(!inherits(actors, "character")) base::stop("Error: The `actor_id` column in the `active_times` argument must be a vector of character values. Please see the help page and retry! Happy dreaming!")
      timestart <- active_times[,"time_start"] #the times for which actors enter
      if(!inherits(timestart, "numeric")) base::stop("Error: The `time_start` column in the `active_times` argument must be a vector of numeric values. Please see the help page and retry! Happy dreaming!")
      timeend <- active_times[,"time_end"] #the times for which actors exit
      if(!inherits(timeend, "numeric")) base::stop("Error: The `time_end` column in the `active_times` argument must be a vector of numeric values. Please see the help page and retry! Happy dreaming!")
    }
    if(type=="two-mode"){
      if(!inherits(active_times, "list")){
        base::stop("Error: A two-mode actor-varying riskset was requested, however, the `active_times` argument is not a list. Please see the help page and retry! Happy dreaming!") # stop computation and tell the user
      }
      x1 <- active_times[["senders"]] #the senders
      x2 <- active_times[["receivers"]] #the receivers
      if(!inherits(x1, "data.frame")){
        base::stop("Error: A two-mode actor-varying riskset was requested, however, the `senders` element in the `active_times` list is not a data.frame object. Please see the help page and retry! Happy dreaming!") # stop computation and tell the user
      }
      if(!inherits(x2, "data.frame")){
        base::stop("Error: A two-mode actor-varying riskset was requested, however, the `receivers` element in the `active_times` list is not a data.frame object. Please see the help page and retry! Happy dreaming!") # stop computation and tell the user
      }
      need <- c("actor_id", "time_start", "time_end")
      if(!all(need %in% colnames(x1))){
        base::stop("Error:  two-mode actor-varying riskset was requested, however, the `senders` element must be a data.frame object with the following column names: actor_id, time_start, and time_end. Please see the help page and retry! Happy dreaming!") # stop computation and tell the user
      }
      if(!all(need %in% colnames(x2))){
        base::stop("Error:  two-mode actor-varying riskset was requested, however, the `receivers` element must be a data.frame object with the following column names: actor_id, time_start, and time_end. Please see the help page and retry! Happy dreaming!") # stop computation and tell the user
      }

      sendingactors <- x1[,"actor_id"] #the actor ids for relevant actors
      if(!inherits(sendingactors, "character")) base::stop("Error: The `actor_id` column in the `senders` data.frame (taken from the `active_times` list) must be a vector of character values. Please see the help page and retry! Happy dreaming!")
      sendingtimestart <- x1[,"time_start"] #the times for which actors enter
      if(!inherits(sendingtimestart, "numeric")) base::stop("Error: The `time_start` column in the `senders` data.frame (taken from the `active_times` list) must be a vector of numeric values. Please see the help page and retry! Happy dreaming!")
      sendingtimeend <- x1[,"time_end"] #the times for which actors exit
      if(!inherits(sendingtimeend, "numeric")) base::stop("Error: The `time_end` column in the `senders` data.frame (taken from the `active_times` list) must be a vector of numeric values. Please see the help page and retry! Happy dreaming!")

      receivingactors <- x1[,"actor_id"] #the actor ids for relevant actors
      if(!inherits(receivingactors, "character")) base::stop("Error: The `actor_id` column in the `receivers` data.frame (taken from the `active_times` list) must be a vector of character values. Please see the help page and retry! Happy dreaming!")
      receivingtimestart <- x1[,"time_start"] #the times for which actors enter
      if(!inherits(receivingtimestart, "numeric")) base::stop("Error: The `time_start` column in the `receivers` data.frame (taken from the `active_times` list) must be a vector of numeric values. Please see the help page and retry! Happy dreaming!")
      receivingtimeend <- x1[,"time_end"] #the times for which actors exit
      if(!inherits(receivingtimeend, "numeric")) base::stop("Error: The `time_end` column in the `receivers` data.frame (taken from the `active_times` list) must be a vector of numeric values. Please see the help page and retry! Happy dreaming!")

    }

  }

  stopifnot(is.numeric(time))
  eventID <- 1:length(time)
  sender <- as.character(sender)
  receiver <- as.character(receiver)
  if(is.null(t)){
    t <- -1 #presetting it to a negative value
    interevent.times <- c(time[1],diff(time))
  }else{
    interevent.times <- c(time[1],diff(time),(t - max(time)))
  }

  if(!is.null(seed)) set.seed(seed) #based upon the user provided inputs, set the seed for reproducibility

  ########################################################
  #
  #   Creating the dataset in c++
  #
  #######################################################


  if(type=="one-mode"){ #for one-mode relational event sequences (the cartesian)
    bigeventlist <- switch(riskset,
                           "complete"=fullrisksetom(time = time,
                                                    seqid = eventID,
                                                    sender = (sender),
                                                    target = (receiver),
                                                    pobserved = p_samplingobserved,
                                                    interval = !ordinal,
                                                    t=t),
                           "constant_sample"=processREMseqOM(time = time,
                                                             seqid = eventID,
                                                             sender = (sender),
                                                             target = (receiver),
                                                             pobserved = p_samplingobserved,
                                                             ncontrols = n_controls,
                                                             interval = !ordinal,
                                                             t=t),
                           "dynamic_sample"=processREMseqOM_varying(time = time,
                                                                    seqid = eventID,
                                                                    sender = (sender),
                                                                    target = (receiver),
                                                                    pobserved = p_samplingobserved,
                                                                    ncontrols = n_controls,
                                                                    interval = !ordinal,
                                                                    t=t),
                           "actor_varying"=timevaryingomriskset(time=time,
                                                                seqid=eventID,
                                                                sender=sender,
                                                                target=receiver,
                                                                timestart=timestart,
                                                                timeend = timeend,
                                                                actors=actors,
                                                                pobserved = 1,
                                                                interval = TRUE,
                                                                t=t),
                           "actor_varying_sample"=timevaryingomrisksetwithsample(time=time,
                                                                                 seqid=eventID,
                                                                                 sender=sender,
                                                                                 target=receiver,
                                                                                 timestart=timestart,
                                                                                 timeend = timeend,
                                                                                 actors=actors,
                                                                                 pobserved = 1,
                                                                                 ncontrols = n_controls,
                                                                                 interval = TRUE,
                                                                                 t=t))
  }
  if(type=="two-mode"){#for one-mode relational event sequences (the cross-product)
    bigeventlist <- switch(riskset,
                           "complete"=fullrisksettm(time = time,
                                                    seqid = eventID,
                                                    sender = (sender),
                                                    target = (receiver),
                                                    pobserved = p_samplingobserved,
                                                    interval = !ordinal,
                                                    t=t),
                           "constant_sample"=processREMseqTM(time = time,
                                                             seqid = eventID,
                                                             sender = (sender),
                                                             target = (receiver),
                                                             pobserved = p_samplingobserved,
                                                             ncontrols = n_controls,
                                                             interval = !ordinal,
                                                             t=t),
                           "dynamic_sample"=processREMseqTM_varying(time = time,
                                                                    seqid = eventID,
                                                                    sender = (sender),
                                                                    target = (receiver),
                                                                    pobserved = p_samplingobserved,
                                                                    ncontrols = n_controls,
                                                                    interval = !ordinal,
                                                                    t=t),
                           "actor_varying"=timevaryingtmriskset(time=time,
                                                                seqid=eventID,
                                                                sender=sender,
                                                                target=receiver,
                                                                timestartsender=sendingtimestart,
                                                                timeendsenders=sendingtimeend,
                                                                actorsenders=sendingactors,
                                                                timestarttarget=receivingtimestart,
                                                                timeendtargets=receivingtimeend,
                                                                actortargets=receivingactors,
                                                                pobserved = 1,
                                                                interval = TRUE,
                                                                t=t),
                           "actor_varying_sample"=timevaryingtmrisksetwithsample(time=time,
                                                                                 seqid=eventID,
                                                                                 sender=sender,
                                                                                 target=receiver,
                                                                                 timestartsender=sendingtimestart,
                                                                                 timeendsenders=sendingtimeend,
                                                                                 actorsenders=sendingactors,
                                                                                 timestarttarget=receivingtimestart,
                                                                                 timeendtargets=receivingtimeend,
                                                                                 actortargets=receivingactors,
                                                                                 pobserved = 1,
                                                                                 interval = TRUE,
                                                                                 t=t))
  }

  #cleaning the post-processing dataset!
  eventSeq <- data.table::rbindlist(bigeventlist) # merging everything into a nice dataframe to be exported to the user!
  eventSeq$sampled <- 1
  colnames(eventSeq) <- c("time", "eventID", "sender", "receiver", "observed", "sampled")
  ### getting the interevent event times
  if(ordinal) interevent.times <- 1
  if(!ordinal) interevent.times <- interevent.times[unique(eventSeq$eventID)]
  #creating the data.table object to be all non-sampled observed events!
  data <- data.table::data.table(time = time[!(eventID %in% eventSeq$eventID)], #the non-sampled event times
                                 eventID = eventID[!(eventID %in% eventSeq$eventID)], #the non-sampled event senders
                                 sender = sender[!(eventID %in% eventSeq$eventID)],#the non-sampled event senders
                                 receiver = receiver[!(eventID %in% eventSeq$eventID)],#the non-sampled event targets
                                 observed = 1,#they are observed
                                 sampled = 0)#they are not sampled
  x <- data.table::rbindlist(list(eventSeq,data)) #combining the objects together!
  x <- x[order(x$time)] #temporal (ordering) sorting the post-processing event sequence by time!
  n <- sum(x$observed)
  sampled_events <- sum(x$sampled * x$observed)
  m<-n_controls
  n_controls <- sum(x$sampled * (1  - x$observed))
  if(type=="one-mode"){
    n_actors <- length(unique(c(sender, receiver))) #the total number of unique actors
    n_senders <- NA
    n_receivers <- NA
  }
  if(type=="two-mode"){
    n_senders <- length(unique(sender))
    n_receivers <- length(unique(receiver))
    n_actors <- n_senders+n_receivers
  }

  #creating the dreamsequence object
  post <- new_dream_sequence(x=x,
                             statistics=list(),
                             ordinal=ordinal,
                             type=type,
                             t=t,
                             p=p_samplingobserved,
                             m=m,
                             n=n,
                             null=n_controls,
                             riskset=riskset,
                             sampled_events=sampled_events,
                             interevent_times=interevent.times,
                             n_actors=n_actors,
                             n_senders=n_senders,
                             n_receivers=n_receivers)
  #validate the object
  post<-validate_dream_sequence(post)
  #returning the created object to the user
  return(post) # output the file to the user!

}




















#' @title Helper Function to Create `dream_sequence` Objects
#' @description
#'
#' The function `dream_sequence()` is a user helper function that transforms
#' user-created processed event sequences into `dream_sequence` objects to be
#' used in the `dream` functions to compute sufficient network statistics
#' and estimate relational event models. This function may also be helpful
#' to user who need to computed network statistics for the estimation of
#' relational outcome models (see [Lerner and Hâncean (2023)](https://www.cambridge.org/core/journals/network-science/article/microlevel-network-dynamics-of-scientific-collaboration-and-impact-relational-hyperevent-models-for-the-analysis-of-coauthor-networks/375932B5B86D2033A0A290DE8198BB32)).
#'
#' If `ordinal` is FALSE, that is, if the relational event sequence is to use the
#' interval timing likelihood, then the events for the last observation time point (the
#' set of realized and null events at the last time point) should all be control events, as
#' they represent the set of right-censoring non-realized events. The `t` should
#' specify the time point that marks the end of the relational event sequence \eqn{A_t}. If `t`
#' is not known, then the value should be left as `NULL`.
#'
#' @name dream_sequence
#' @param type "two-mode" indicates that this is a two-mode event sequence. "one-mode" indicates that the event sequence is one-mode.
#' @param ordinal TRUE/FALSE. TRUE indicates that observed timing of the events is ordinal (and the ordinal timing likelihood function will be used). FALSE denotes
#' that the observed timing is observed, relative to the start of the event sequence (and the interval timing likelihood function will be used). Please see the references for more information.
#' @param t If ordinal is set to FALSE, the time that marks the end of the relational event sequence, relative
#' to the start of the event sequence.
#' @param time A numeric vector that contains the timing of the events in the relational event sequence.
#' @param sender A character vector that contains the sender of the events in the relational event sequence.
#' @param receiver A character vector that contains the receiver/target of the events in the relational event sequence.
#' @param observed A numeric vector that is 1 if the observation is an observed event in the relational event sequence, or 0 if the observation is a control event in the relational event sequence (see \code{\link{create_res}} for more information).
#' @param sampled A numeric vector that is 1 if the observation is a sampled event in the relational event sequence, or 0 if the observation is a non-sampled event in the relational event sequence (see \code{\link{create_res}} for more information).
#' @param statistics A `data.frame` object that contains the previously computed statistics
#' for the processed events. If processing has occurred (i.e., the `sampled` argument is
#' specified; control events have been created and/or if sampling from the observed event sequence),
#' the the `statistics` argument must have a number of observations equal to the sum of the
#' `sampled` argument. If processing has not occurred, then the `statistics` argument must have
#' a number of observations equal to the length of the `time` argument.
#' @param ... Additional arguments (currently unused).
#' @importFrom data.table data.table
#' @return A `dream_sequence` object that contains the user-provided information.
#' @export
#'
#' @examples
#' #a pseudo event sequence
#' events <- data.frame(time = 1:18, eventID = 1:18,
#'                      sender = c("A", "B", "C",
#'                                 "A", "D", "E",
#'                                 "F", "B", "A",
#'                                 "F", "D", "B",
#'                                "G", "B", "D",
#'                                "H", "A", "D"),
#'                    target = c("B", "C", "D",
#'                                "E", "A", "F",
#'                                "D", "A", "C",
#'                                 "G", "B", "C",
#'                                 "H", "J", "A",
#'                                 "F", "C", "B"))
#' #making a post-processing event sequence
#' eventSet <- create_res(type = "one-mode",
#'                        ordinal = TRUE,
#'                        riskset = "constant_sample",
#'                        time = events$time,
#'                        sender = events$sender,
#'                        receiver = events$target,
#'                        p_samplingobserved = 1.00,
#'                        n_controls = 1,
#'                        seed = 9999)
#'
#' #computing the sender indegree statistics
#' eventSet <- dreamstats_degree(formation = "sender-indegree", data = eventSet)
#' #making the dream_sequence object a data.frame object
#' new.data <- as.data.frame(eventSet) #making the event set a data.frame object
#' stats <- new.data["sender.indegree"] #the computed statistics
#'
#' #reconverting the object to a dream_sequence object
#' psuedo.data <- dream_sequence(ordinal = TRUE,
#'                               time = new.data$time,
#'                               sender = new.data$sender,
#'                               receiver = new.data$receiver,
#'                               observed = new.data$observed,
#'                               sampled = new.data$sampled,
#'                               type= "one-mode",
#'                               statistics = stats)
#' psuedo.data #printing the object
#'
#' #reconverting we original event sequence to a dream_sequence object
#' #(this is helpful for the estimation of relational outcome models!)
#' psuedo.data1 <-  dream_sequence(ordinal = TRUE,
#'                                 time = events$time,
#'                                 sender = events$sender,
#'                                 receiver = events$target,
#'                                 type= "one-mode")
#'
#' psuedo.data1 #printing the object
#'
#' #computing a statistic on the data
#' psuedo.data1 <- dreamstats_degree(formation = "sender-outdegree",
#'                                   data = psuedo.data1,
#'                                   counts = TRUE)
#'
#' psuedo.data1 #printing the object with the computed statistics
#' psuedo.data1$statistics #printing the object with the computed statistics
dream_sequence <- function(ordinal = TRUE,
                           t = NULL,
                           time,
                           sender,
                           receiver,
                           observed = NULL,
                           sampled = NULL,
                           type = "one-mode",
                           statistics = NULL,
                           ...){

  #checking for errors in user inputs
  if(length(time) != length(sender) ||
     length(time) != length(receiver)){
    base::stop("The `time`, `sender`, and `receiver` arguments must all be the same length.")
  }

  #checking for errors in user inputs
  if(!is.null(observed) & length(time) != length(observed)){
    base::stop("The `observed` argument must be the same length as the other entries.")
  }
  #checking user inputs
  if(is.null(observed)) base::message("An entry was not specified for the `observed` argument, so we are assuming that all entries are realized (observed) events.")
  if(is.null(observed)) observed <- rep(1,length(time)) #making the vector for the argument

  #checking for errors in user inputs
  if(!is.null(sampled) & length(sampled) != length(observed)){
    base::stop("The `sampled` argument must be the same length as the other entries.")
  }

  #checking for errors in user inputs
  if(!is.null(statistics)){
    if(!inherits(statistics, "data.frame")){
      base::stop("The `statistics` argument must be either a `data.frame` or a `data.table` object.", call. = FALSE)
    }
  }

  #checking user inputs
  if(is.null(sampled)) base::message("An entry was not specified for the `sampled` argument, so we are assuming that the included events represent the entire relational event sequence.")
  if(is.null(sampled)) sampled <- rep(1,length(time)) #making the vector for the argument

  if(!is.null(statistics)){
    if(nrow(statistics) != sum(sampled)) base::stop("The number of observations for the `statistics` argument must correspond to the number of `sampled` events. If the `sampled` argument was not specified, then the number of observations should correspond to the length of the `time` argument.")
  }

  #converting objects for later use
  if(!is.null(statistics)) statistics <- as.list(statistics) #making the object a list
  if(is.null(statistics)) statistics <- list() #making the object a list

  #making the object
  if(is.null(t)) t <- -1

  interevent.times1 <- c(time[1],diff(time))
  interevent.times <- c(interevent.times1[sampled==1],0) #only those sampled events that will be used in the REM later
  interevent.times <- unique(interevent.times) #the unique time points

  events <- data.table::data.table(time,sender,receiver,observed,sampled)
  n <- sum(observed)
  sampled_events <- sum(observed*sampled)
  n_controls <- sum(1-observed)
  #these are just attributes that do not matter for actual variables
  p <- NA
  m <- NA
  seed <- NA
  riskset <- "user-provided"
  if(type=="one-mode"){
    n_actors <- length(unique(c(sender, receiver))) #the total number of unique actors
    n_senders <- NA
    n_receivers <- NA
  }
  if(type=="two-mode"){
    n_senders <- length(unique(sender))
    n_receivers <- length(unique(receiver))
    n_actors <- n_senders+n_receivers
  }

  #creating the dreamsequence object

  #creating the dreamsequence object
  post <- new_dream_sequence(x=events,
                             statistics=statistics,
                             ordinal=ordinal,
                             type=type,
                             t=t,
                             p=p,
                             m=m,
                             n=n,
                             null=n_controls,
                             riskset=riskset,
                             sampled_events=sampled_events,
                             interevent_times=interevent.times,
                             n_actors=n_actors,
                             n_senders=n_senders,
                             n_receivers=n_receivers)
  #validate the object
  post<-validate_dream_sequence(post)
  #returning the created object to the user
  return(post)
}





#' Constructor for `dream_sequence` objects
#' @param x `data.frame` for the processed event sequence
#' @param statistics a `list` of the computed network statistics.
#' @param ordinal `logical` indicating if the sequence type is interval or ordinal timing.
#' @param type either of "one-mode" or "two-mode" characterizing the type of relational event sequence.
#' @param t If interval timing event sequence, the time point that marks the end of the sequence.
#' @param p The probability of sampling from the observed event sequence.
#' @param m The number of controls for each event.
#' @param n The number of realized (observed) events.
#' @param null The number of non-realized (control) events.
#' @param riskset `logical` indicating the definition of the risk set.
#' @param sampled_events The number of sampled observed events.
#' @param n_actors The number of actors in the relational event sequence.
#' @param n_receivers The number of event receivers in the relational event sequence.
#' @param n_senders The number of event senders in the relational event sequence.
#' @param interevent_times The numeric vector of interevent times (the timing between events).
#' @keywords internal
new_dream_sequence <- function(x,
                               statistics,
                               ordinal,
                               type,
                               riskset,
                               t,
                               p,
                               n,
                               sampled_events,
                               null,
                               n_actors,
                               n_senders,
                               n_receivers,
                               m,
                               interevent_times
){
  #checking the class of the processed event sequence (x is built elsewher)
  stopifnot(inherits(x,"data.frame") | inherits(x,"data.table"))
  #checking the argument of the ordinal value
  stopifnot(inherits(ordinal,"logical"))
  #checking the class of the statistics argument
  stopifnot(inherits(statistics,"list"))
  #checking the argument of the type value
  stopifnot(type %in% c("one-mode", "two-mode"))
  #checking the argument of the interevent type value
  stopifnot(is.numeric(interevent_times))
  #creating the object
  x <- list(processed_sequence = x, #the post processing event sequence
            statistics = statistics, #the list of computed statistics
            ordinal = ordinal, #ordinal vs. interval timing
            type = type, #the type of event sequence
            n_actor=n_actors,
            n_senders=n_senders,
            n_receivers=n_receivers,
            t = t, #the ending time point in interval timing log likelihood
            p = p, #the probability of sampling from the observed event sequence
            m = m, #the number of sampled null event controls
            n = n, #the number of observed events
            null = null, #the number of null events
            sampled_events=sampled_events, #the number of sampled observed events
            riskset = riskset, #the type of risk set employed
            interevent_times = interevent_times) #the interevent times between events
  post <- structure(x, class = "dream_sequence")
  return(post) #returning the object
}






#' Validator for `dream_sequence` objects
#' @param x a created `dream_sequence` object to be checked.
#' @keywords internal
validate_dream_sequence <- function(x){

  elements <- unclass(x)

  if(!inherits(elements$processed_sequence,"data.frame") |
     !inherits(elements$processed_sequence,"data.table")){
    stop("The `processed_sequence` argument must be an object of class `data.frame` or `data.table`.",call. = FALSE)
  }

  d <- colnames(elements$processed_sequence)
  q <- c("time", "sender", "receiver", "sampled", "observed")
  if(!all(q %in% d)){
    stop("The `processed_sequence` argument must be have the following columms: time, sender, receiver, sampled, and observed.",call. = FALSE)
  }

  #checking for errors in user inputs
  if(!is.numeric(elements$processed_sequence$time)) base::stop("The `time` vector must be a vector of numeric value.",call. = FALSE)
  if(!is.character(elements$processed_sequence$sender)) base::stop("The `sender` vector must be a vector of character values.",call. = FALSE)
  if(!is.character(elements$processed_sequence$receiver)) base::stop("The `receiver` vector must be a vector of character values.",call. = FALSE)

  if(!is.numeric(elements$processed_sequence$observed)) base::stop("The `observed` argument must be a vector of numeric value.",call. = FALSE)

  the.values <- unique(elements$processed_sequence$observed) #the unique values
  if(length(the.values) == 1){
    if(!(the.values %in% c(0,1))) base::stop("The entries in the `observed` argument must be only be 0 (null/control event) and 1 (realized event).",call. = FALSE)
  }else{
    if(length(the.values) != 2)  base::stop("The entries in the `observed` argument must be only be 0 (null/control event) and 1 (realized event).",call. = FALSE)
    if(!(0 %in% the.values))  base::stop("The entries in the `observed` argument must be only be 0 (null/control event) and 1 (realized event).",call. = FALSE)
    if(!(1 %in% the.values))  base::stop("The entries in the `observed` argument must be only be 0 (null/control event) and 1 (realized event).",call. = FALSE)
  }

  if(!is.numeric(elements$processed_sequence$sampled))base::stop("The `sampled` argument must be a vector of numeric value.",call. = FALSE)
  the.values <- unique(elements$processed_sequence$sampled) #the unique values
  if(length(the.values) == 1){
    if(!(the.values %in% c(0,1))) base::stop("The entries in the `sampled` argument must be only be 0 (null/control event) and 1 (realized event).",call. = FALSE)
  }else{
    if(length(the.values) != 2)  base::stop("The entries in the `sampled` argument must be only be 0 (non-sampled event) and 1 (sampled event).",call. = FALSE)
    if(!(0 %in% the.values))  base::stop("The entries in the `sampled` argument must be only be 0 (non-sampled event) and 1 (sampled event).",call. = FALSE)
    if(!(1 %in% the.values))  base::stop("The entries in the `sampled` argument must be only be 0 (non-sampled event) and 1 (sampled event).",call. = FALSE)
  }

  if(sum(elements$processed_sequence$sampled[elements$processed_sequence$observed == 0] == 0) != 0) base::stop("You supplied entries for the `sampled` and `observed` arguments. There are values of 0 for the `sampled` argument for events (entries) that are `control` events (i.e., have a value of 0 in the observed vector). Please fix this, since by default, in this package, it is assumed that control events are sampled (i.e., a value of 1 on the sampled vector).",call. = FALSE)

  if(!inherits(elements$statistics,"list")){
    stop("The `statistics` argument must be a list (i.e., object of class `list`).",call. = FALSE)
  }

  if(!inherits(elements$ordinal, "logical")){
    stop("The `ordinal` argument must be a TRUE (ordinal timing) / FALSE (interval timing).",call. = FALSE)
  }

  if(!(elements$type %in% c("one-mode", "two-mode"))){
    stop("The `type` argument must be either `one-mode` or `two-mode`.",call. = FALSE)
  }
  #checking for errors in user inputs
  if(length(elements$type) != 1) base::stop("The `type` argument must be of length 1 (i.e., either one-mode or two-mode).",call. = FALSE)

  if(!is.numeric(elements$t)){
    stop("The `t` argument must be a numeric value.",call. = FALSE)
  }

  if(!is.numeric(elements$interevent_times)){
    stop("The `interevent_times` argument must be a vector of numeric values.",call. = FALSE)
  }

  if(any(elements$interevent_times < 0)){
    stop("All of the `interevent_times` entries must be non-negative values.",call. = FALSE)
  }
  return(x)
}













#' Coerce a `dream_sequence` Object into a `data.frame` Object
#'
#' This function will create a `data.frame` object from a `dream_sequence` object, where
#' the generated data.frame includes the processed event sequence and the
#' computed statistics. If `all.events` is set to FALSE, then the created
#' object will only contains the sampled events, whereas, if it is TRUE, the
#' returned object will contain all events, and if any computed statistics are attached,
#' the non-sampled event entries will be NAs.
#'
#'
#'
#' @param x An object of class `dream_sequence`.
#' @param row.names The `row.names` argument from the \code{\link[base]{as.data.frame}} function. Please
#' see that function for more details on this argument.
#' @param optional The `optional` argument from the \code{\link[base]{as.data.frame}} function. Please
#' see that function for more details on this argument.
#' @param all.events TRUE/FALSE. If sampling from the observed event sequence has occurred, TRUE
#' returns all of the sampled and non-sampled observed events, whereas FALSE returns only those
#' sampled observed events. FALSE by default.
#' @param ... Additional arguments for other methods..
#' @return Returns a `data.frame` object that contains the processed relational event statistics and
#' the associated sufficient network statistics for the event sequence.
#' @export
#'
#' @examples
#'
#' #a pseudo event sequence
#' events <- data.frame(time = 1:18, eventID = 1:18,
#'                      sender = c("A", "B", "C",
#'                                 "A", "D", "E",
#'                                 "F", "B", "A",
#'                                 "F", "D", "B",
#'                                "G", "B", "D",
#'                                "H", "A", "D"),
#'                    target = c("B", "C", "D",
#'                                "E", "A", "F",
#'                                "D", "A", "C",
#'                                 "G", "B", "C",
#'                                 "H", "J", "A",
#'                                 "F", "C", "B"))
#' #making a post-processing event sequence
#' processed <- create_res(type = "one-mode",
#'                        ordinal = TRUE,
#'                        riskset = "constant_sample",
#'                        time = events$time,
#'                        sender = events$sender,
#'                        receiver = events$target,
#'                        p_samplingobserved = 1.00,
#'                        n_controls = 1,
#'                        seed = 9999)
#'
#' #computing the sender indegree statistics
#' processed <- dreamstats_degree(formation = "sender-indegree", data = processed)
#'  #computing the outgoing two paths statistics
#' processed <- dreamstats_triads(formation = "OTP", data = processed)
#' #computing the repetition/inertia statistics
#' processed <- dreamstats_repetition(data = processed)
#' #making the dream_sequence object a data.frame object
#' new.data <- as.data.frame(processed) #making the event set a data.frame object
#' new.data #the processed event sequence returned with the computed statistics

as.data.frame.dream_sequence <- function(x, row.names = NULL, optional = FALSE, all.events=FALSE, ...) {

  if(isTRUE(all.events)) data <- x$processed_sequence
  if(!isTRUE(all.events)) data <- x$processed_sequence[x$processed_sequence$sampled == 1,]
  if(length(x$statistics) > 0){
    if(!isTRUE(all.events)){
      data <- cbind(data, data.frame(x$statistics))
    }else{
      data1 <- data[data$sampled==1,]
      data1 <- data.table::data.table(cbind(data1,data.frame(x$statistics)))
      data2 <- data[data$sampled==0,]
      data3 <- data.table::rbindlist(list(data1,data2),fill = TRUE)
      data <- data3[order(data3$eventID),]
    }
  }
  x <- as.data.frame(data, row.names = row.names, optional = optional)
  return(x)
}










#' Print Method for 'dream'  object
#'
#' @param x An object of class 'dream_support' .
#' @param digits The number of digits to print after the decimal point.
#' @param ... Additional arguments (currently unused).
#' @return No return value. Prints out the main results of a 'dream' object.
#' @export

print.dream_sequence <- function(x,digits=4,...) {
  cat("Processed Relational Event Sequence: \n")
  cat(" -> Relational event sequence type:", x$type,"\n")
  cat(" -> Relational event sequence timing:", ifelse(x$ordinal, "ordinal", "interval"),"\n")
  if(!x$ordinal & x$t > 0){
    cat(" -> The end of the relational event sequence:", x$t,"\n")
  }
  if(x$type == "one-mode")   cat(" -> Number of actors:", x$n_actor,"\n")
  if(x$type == "two-mode")   cat(" -> Number of senders:", x$n_senders,"\n -> Number of receivers:", x$n_receivers,"\n")
  cat(" -> Number of realized events:", x$n,"\n")
  n <- x$n
  cat(" -> Sampling from the realized event sequence?:", ifelse(x$p == 1, "no", "yes"),"\n")
  if(x$p != 1 & !is.na(x$p)){
    cat(" -> The probabilty of sampling from the realized event sequence:",x$p,"\n")
    cat(" -> Number of sampled realized events:",x$sampled_events," \n")
    n <- x$sampled_events
  }
  cat(" -> The risk/support set defintion:", x$riskset,"\n")
  cat(" -> Case-control sampling of control events?:", ifelse(x$riskset %in% c("constant_sample", "dynamic_sample",
                                                                               "actor_varying_sample"), "yes", "no"),"\n")

  if(x$riskset %in% c("constant_sample", "dynamic_sample","actor_varying_sample"))  cat(" -> Case-control sampling m:", x$m,"\n")
  cat(" -> Number of non-realized/control events:",x$null," \n")
  cat(" -> The total number of processed realized and non-realized events:",x$n + x$null," \n")
  if(length(x$statistics) > 0){
    cat("\nComputed Relational Event Statistics:\n")
    for(i in 1:length(x$statistics)){
      cat(" -> ",names(x$statistics)[i],"\n")
    }
  }
  invisible(x)
}

#' Summary Method for dream_sequence Objects
#'
#' Summarizes the main components of a processed relational event sequence.
#'
#' @param object An object of class "dream_sequence".
#' @param digits The number of digits to print after the decimal point.
#' @param ... Additional arguments (currently unused).
#' @return A list of descripitive statistics for a processed relational event sequence.
#' @export
#'
summary.dream_sequence <- function(object,digits=4,...) {
  res <- object
  class(res) <- "summary.dream_sequence"
  return(res)
}

#' Print Method for dream Model
#'
#' @param x An object of class "dream_sequence".
#' @param digits The number of digits to print after the decimal point.
#' @param ... Additional arguments (currently unused).
#' @return No return value. Prints out the main results of a 'dream_sequence' summary object.
#' @export
print.summary.dream_sequence<- function(x,digits=3,...) {
  cat("Processed Relational Event Sequence: \n")
  cat(" -> Relational event sequence type:", x$type,"\n")
  cat(" -> Relational event sequence timing:", ifelse(x$ordinal, "ordinal", "interval"),"\n")
  if(!x$ordinal & x$t > 0){
    cat(" -> The end of the relational event sequence:", x$t,"\n")
  }
  if(x$type == "one-mode")   cat(" -> Number of actors:", x$n_actor,"\n")
  if(x$type == "two-mode")   cat(" -> Number of senders:", x$n_senders,"\n -> Number of receivers:", x$n_receivers,"\n")
  cat(" -> Number of realized events:", x$n,"\n")
  n <- x$n
  cat(" -> Sampling from the realized event sequence?:", ifelse(x$p == 1, "no", "yes"),"\n")
  if(x$p != 1 & !is.na(x$p)){
    cat(" -> The probabilty of sampling from the realized event sequence:",x$p,"\n")
    cat(" -> Number of sampled realized events:",x$sampled_events," \n")
    n <- x$sampled_events
  }
  cat(" -> The risk/support set defintion:", x$riskset,"\n")
  cat(" -> Case-control sampling of control events?:", ifelse(x$riskset %in% c("constant_sample", "dynamic_sample",
                                                                               "actor_varying_sample"), "yes", "no"),"\n")

  if(x$riskset %in% c("constant_sample", "dynamic_sample","actor_varying_sample"))  cat(" -> Case-control sampling m:", x$m,"\n")
  cat(" -> Number of non-realized/control events:",x$null," \n")
  cat(" -> The total number of processed realized and non-realized events:",x$n + x$null," \n")

  if(length(x$statistics) > 0){
    cat("\nComputed Relational Event Statistics:\n")
    for(i in 1:length(x$statistics)){
      cat("   -> ",names(x$statistics)[i],"\n")
      print(round(summary(x$statistics[[i]]), digits = digits))
      cat("\n")
    }
  }
}





