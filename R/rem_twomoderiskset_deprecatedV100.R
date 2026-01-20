## The Creation of  Dynamic Two-Mode Risk Sets
## See Butts 2008; Lerner and Lomi 2020; REM R Package (Brandenberger 2018)
## Code written by Kevin Carson (kacarson@arizona.edu) and Deigo Leal (https://www.diegoleal.info/)
## Last Updated: 07-15-24
#' @title Process and Create Risk Sets for a Two-Mode Relational Event Sequence
#' @name processTMEventSeq
#' @param data The full relational event sequence dataset.
#' @param time The vector of event time values from the observed event sequence.
#' @param sender The vector of event senders from the observed event sequence.
#' @param receiver The vector of event receivers from the observed event sequence.
#' @param eventID The vector of event IDs from the observed event sequence (typically a numerical event sequence that goes from 1 to *n*).
#' @param p_samplingobserved The numerical value for the probability of selection for sampling from the observed event sequence. Set to 1 by default indicating that all observed events from the event sequence will be included in the post-processing event sequence.
#' @param n_controls The numerical value for the number of null event controls for each (sampled) observed event.
#' @param time_dependent TRUE/FALSE. TRUE indicates that a time- or event-dependent dynamic risk set will be created in which only actors involved in a user-specified relationally relevant (time or event) span (i.e., the ‘stretch’ of relational relevancy, such as one month for a time-dependent risk set or 100 events for an event-dependent risk set) are included in the potential risk set. FALSE indicates the complete set of actors involved in past events will be included in the risk set (see the details section). Set to FALSE by default.
#' @param timeDV If time_dependent = TRUE, the vector of event time values that corresponds to the creation of the time- *or* event-dependent dynamic risk set (see the details section). *This may or may not be the same vector provided to the time argument*. The *timeDV* vector can be the same vector provided to the *time* argument, in which the relational time span will be based on the event timing within the dataset. In contrast, the *timeDV* vector can also be the vector of numerical event IDs which correspond to the number sequence of events. Moreover, the *timeDV* can also be another measurement that is not the *time* argument or a numerical event ID sequence, such as the number of days, months, years, etc. since the first event.
#' @param timeDif If time_dependent = TRUE, the numerical value that represents the time or event span for the creation of the risk set (see the details section). This argument must be in the same measurement unit as the `timeDV` argument. For instance, in an event-dependent dynamic risk set, if `timeDV` is the number of events since the first event (i.e., a numerical event ID sequence) and only those actors involved in the past, say, 100 events, are considered relationally relevant for the creation of the null events for the current observed event, then `timeDIF` should be set to 100. In the time-dependent dynamic risk set case, let’s say that only those actors involved in events that occurred in the past month are considered relationally relevant for the risk set. Let’s also assume that the `timeDV` vector is measured in the number of days since the first event. Then `timeDif` should be set to 30 in this particular case.
#' @param seed The random number seed for user replication.
#' @import data.table
#' @importFrom collapse whichv
#' @importFrom collapse fduplicated
#' @importFrom collapse funique
#' @importFrom collapse `%==%`
#' @importFrom dqrng dqsample
#' @importFrom fastmatch `%fin%`
#' @return A post-processing data table with the following columns:
#' \itemize{
#'   \item \code{sender} - The event senders of the sampled and observed events.
#'   \item \code{receiver} - The event targets (receivers) of the sampled and observed events.
#'   \item \code{time} - The event time for the sampled and observed events.
#'   \item \code{sequenceID} - The numerical event sequence ID for the sampled and observed events.
#'   \item \code{observed} - Boolean indicating if the event is a sampled event or observed event. (1 = observed; 0 = sampled)
#' }
#' @export
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `processTMEventSeq()` has been deprecated starting on version 1.0.0 of the `dream` package. Please use the `create_riskset()` function and see the `NEWS.md` file for more details.
#'
#'
#' This function creates a two-mode post-sampling eventset with options for case-control
#' sampling (Vu et al. 2015), sampling from the observed event sequence (Lerner and Lomi 2020), and time- or event-dependent
#' risk sets. Case-control sampling samples an arbitrary *m* number of controls from the risk set for any event
#' (Vu et al. 2015). Lerner and Lomi (2020) proposed sampling from the observed event sequence
#' where observed events are sampled with probability *p*. The time- and event-dependent risk sets generate risk sets where the
#' potential null events are based upon a specified past relational time window, such as events that have occurred in the past month.
#' Users interested in generating risk sets that assume all actors active at any time point within the event
#' sequence are in the risk set at every time point should consult the \code{\link[rem]{createRemDataset}}
#' and \code{\link[remify]{remify}} functions. Future versions of this package will
#' incorporate this option into the function.
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
#' Following Butts (2008) and Butts and Marcum (2017), we define the risk (support) set of all possible
#' events at time \eqn{t}, \eqn{A_t}, as the cross product of two disjoint sets, namely, prior senders and receivers,
#' in the set \eqn{G[E;t]} that could have occurred at time \eqn{t}. Formally:
#' \deqn{A_t = \{ (s, r) \mid s \in G[E;t] \text{ X } r \in G[E;t] \}}
#' where \eqn{G[E;t]} is the set of events up to time \eqn{t}.
#'
#' Case-control sampling maintains the full set of observed events, that is, all events in \eqn{E}, and
#' samples an arbitrary number \eqn{m} of non-events from the support set \eqn{A_t} (Vu et al. 2015; Lerner
#' and Lomi 2020). This process generates a new support set, \eqn{SA_t}, for any relational event
#' \eqn{e_i} contained in \eqn{E} given a network of past events \eqn{G[E;t]}. \eqn{SA_t} is formally defined as:
#' \deqn{SA_t \subseteq \{ (s, r) \mid s \in G[E;t] \text{ X } r \in G[E;t] \}}
#' and in the process of sampling from the observed events, \eqn{n} number of observed events are
#' sampled from the set \eqn{E} with known probability \eqn{0 < p \le 1}. More formally, sampling from
#' the observed set generates a new set \eqn{SE \subseteq E}.
#'
#' A time *or* event-dependent dynamic risk set can be created where the set of potential events,
#' that is, all events in the risk set, At, is based only on the set of actors active in a
#' specified event or time span from the current event (e.g., such as within the past month
#' or within the past 100 events). In other words, the specified event or time span can be
#' based on either: a) a specified time span based upon the actual timing of the past events
#' (e.g., years, months, days or even milliseconds as in the case of Lerner and Lomi 2020),
#' or b) a specified number of events based on the ordering of the past events (e.g., such
#' as all actors involved in the past 100 events). Thus, if time- or event-dependent dynamic
#' risk sets are desired, the user should set time_dependent to TRUE, and then specify the
#' accompanying time vector, `timeDV`, defined as the number of time units (e.g., days) or the
#' number of events since the first event. Moreover, the user should also specify the cutoff
#' threshold with the `timeDif` value that corresponds directly to the measurement unit of
#' `timeDV` (e.g., days). For example, let’s say you wanted to create a time-dependent dynamic
#' risk set that only includes actors active within the past month, then you should create a
#' vector of values `timeDV`, which for each event represents the number of days since the first
#' event, and then specify `timeDif` to 30. Similarly, let’s say you wanted to create an event-dependent
#' dynamic risk set that only includes actors involved in the past 100 events, then you should create
#' a vector of values `timeDV`, that is, the counts of events since the first event (e.g., 1:n), and
#' then specify `timeDif` to 100.
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
#' EventSet <- processTMEventSeq(
#'   data = WikiEvent2018.first100k, # The Event Dataset
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
#' EventSet1 <- processTMEventSeq(
#'   data = WikiEvent2018.first100k, # The Event Dataset
#'   time = WikiEvent2018.first100k$time, # The Time Variable
#'   eventID = WikiEvent2018.first100k$eventID, # The Event Sequence Variable
#'   sender = WikiEvent2018.first100k$user, # The Sender Variable
#'   receiver = WikiEvent2018.first100k$article, # The Receiver Variable
#'   p_samplingobserved = 0.02, # The Probability of Selection
#'   n_controls = 2, # The Number of Controls to Sample from the Full Risk Set
#'   seed = 9999) # The Seed for Replication
#'
#' ### Creating An Event-Dependent EventSet with P = 0.001 and m = 5 with
#' ### where only actors involved in the past 20 events are involved in the
#' ### creation of the risk set.
#'event_dependent <- processTMEventSeq(
#'  data = WikiEvent2018.first100k,
#'  time = WikiEvent2018.first100k$time,
#'  sender = WikiEvent2018.first100k$user,
#'  receiver = WikiEvent2018.first100k$article,
#'  eventID = WikiEvent2018.first100k$eventID,
#'  p_samplingobserved = 0.001,
#'  n_controls = 5,
#'  time_dependent = TRUE,
#'  timeDV = 1:nrow(WikiEvent2018.first100k),
#'  timeDif = 20, #20 past events
#'  seed = 9999)

#' ### Creating An Time-Dependent EventSet with P = 0.001 and m = 5 with
#' ### where only actors involved in the past 30 days are involved in the
#' ### creation of the risk set.
# Generate the time difference: one month in milliseconds
# 30 days * 24 hours * 60 minutes * 60 seconds * 1000
#'timeSinceStart <- WikiEvent2018.first100k$time-WikiEvent2018.first100k$time[1]
#'timeDifMonth <- 30*24*60*60*1000
#'timedependent <- processTMEventSeq(
#'  data = WikiEvent2018.first100k,
#'  time = WikiEvent2018.first100k$time,
#'  sender = WikiEvent2018.first100k$user,
#'  receiver = WikiEvent2018.first100k$article,
#'  eventID = WikiEvent2018.first100k$eventID,
#'  p_samplingobserved = 0.001,
#'  n_controls = 5,
#'  time_dependent = TRUE,
#'  timeDV = timeSinceStart,
#'  timeDif = timeDifMonth,
#'  seed = 9999)



########################################################################################################
#  Data = the full event sequence
#  Time = the name for the time variable
#  eventID = the event sequence that holds the event number (or ID) for each event (must be sequentially ordered)
#  sender = the name for the sender variable
#  receiver = the name for the receiver variable
#  p_samplingobserved = probability of selection for sampling from the observed event sequence (see Lerner and Lomi 2020)
#  n_controls = number of sampled controls from the null events per each event (see Lerner and Lomi 2020; Vu et al. 2015)
#  seed = random seed for replication
########################################################################################################
processTMEventSeq <- function(    data, # full event data sequence
                           time, # variable (column) name that contains the time variable
                           eventID, # variable (column) name that contains event Sequence ID
                           sender, # variable (column) name that contains the sender variable
                           receiver, # variable (column) name that contains the receiver variable
                           p_samplingobserved = 1, # probability of selection for case control sampling
                           n_controls, # number of controls for each selected event
                           time_dependent = FALSE, #should a dyanmic time-varying risk set be created?
                           timeDV = NULL, #the dynamic time-varying vector? A vector of standarized time values, could be an increasing minutes, days, months etc
                           timeDif = NULL, #the value of time difference based on the timeDV vector of values
                           seed = 9999) { # seed for replication (user can change this value)

  #base::cat("Checking Data Structure and User Inputs.......") # outputting status to user
  lifecycle::deprecate_warn("1.0.0", " processTMEventSeq()", "create_riskset()")
  ########################################################
  #
  #   Checking for Errors in User Inputs
  #
  ########################################################

  if (!data.table::is.data.table(data)) { # if data is not a data.table, (this speeds up computation)
    data <- data.table::data.table(data) # then, make it a data.table
  }
  n_forrealevents <- base::nrow(data) #the number of real events provided by user

  if (length(time) != n_forrealevents) { # if the length of the time vector is not the same as the full dataset
    base::stop("Error: The length of the provided events time is not the same length as the events dataset") # stop computation and tell the user
  }
  if (length( sender) != n_forrealevents) { # if the length of the time vector is not the same as the full dataset
    base::stop("Error: The length of the provided events sender is not the same length as the events dataset") # stop computation and tell the user
  }
  if (length( receiver) != n_forrealevents) { # if the length of the time vector is not the same as the full dataset
    base::stop("Error: The length of the provided events receiver is not the same length as the events dataset") # stop computation and tell the user
  }
  if (length(eventID) != n_forrealevents) { # if the length of the time vector is not the same as the full dataset
    base::stop("Error: The length of the provided events id is not the same length as the events dataset") # stop computation and tell the user
  }
  p <- p_samplingobserved # copying value for sample selection
  if (p > 1 | p < 0) { # if the probability is not a probability (i.e., bounded by 0 and 1)
    base::stop("Error: Probabilty of Selection Must be within the interval: 0 < p <= 1") # stop computation and tell the user
  }
  m <- n_controls # copying value for number of controls
  if (n_controls <= 0) { # if number of controls is equal to 0, that is, no null events
    base::stop("Error: Number of Controls Must Be At Least 0") # stop computation and tell the user
  }

  ########################################################
  # Renaming the columns to match the user inputs
  ########################################################
  data$sender <- sender # renaming the sender column
  data$receiver <- receiver # renaming the receiver column
  data$time <- time # renaming the time column
  data$EVENT <- eventID # renaming the time column

  ########################################################
  #### Clearing User Inputs for Memory
  ########################################################
  rm(list = c("sender", "time", "eventID", "receiver"))
  ########################################################
  #
  #   Sampling from Observed Events Based on User Probability Entry
  #
  ########################################################

  #base::cat("Sampling from Observed Event Sequence.......") # outputting status to user
  n <- base::nrow(data) # getting number of data events

  data$sampled <- 0 # creating a new column to indicate the event is sampled or not (1 = yes; 0 = no)
  set.seed(seed) # setting random seed
  how.many <- base::round(n * p) # how many observed cases do we want?
  eventsequence <- 1:n # the event sequence, 1:n all events (essentially an indicator)
  observed <- dqrng::dqsample(eventsequence, size = how.many, replace = F) # sampling from values
  observed <- base::sort(observed) # sorting the events from largest to smallest
  data$sampled[observed] <- 1 # make theses values be equal to 1 for sampled cases
  events <- data[observed] # only getting the observed sampled events
  ########################################################
  # Creating a list to store the id of observed events and the sampled event information
  ########################################################
  sampled_events <- list(
    events = events,
    observed = observed
  )
  rm(list = c("events", "observed")) # removing objects that are no longer required to clean memory

  ################################################################
  #
  #   Sampling from Null Events Based on Case-Control Sampling
  #
  ################################################################

  #base::cat("Sampling from the Risk Set.......") # outputting status to user

  events_sampling <- list() # an empty list to store the sampled events for each observed event
  observed <- sampled_events$observed # the ids of the sampled observed events
  events <- data.table::as.data.table(sampled_events$events) # making the events dataframe a data.table
  ##### only getting the variables that we need in order to create the null events
  events <- events[, which(colnames(events) %in% c(
    "sender", "receiver", "sampled", "EVENT",
    "OBSERVED", "time"
  )), with = FALSE]
  events$OBSERVED <- 1 # making the variable 1 for all true events; 1 = real event, 0 = control event
  sender <- data$sender # making the full event sequence sender variable a vector
  receiver <- data$receiver # making the full event sequence receiver variable a vector
  time <- data$time # making the full event sequence time variable a vector

  ################################################################################
  #### Case Control Sampling from Unobserved Events with Known M
  ################################################################################

  if(time_dependent == FALSE){ #creating a non-dynamic risk set
  #### Step 3: Randomly Sample Control Events (decide an n for selection)
  ## we should store all these values into a list
  #### Now lets sample m possible controls per event
  i <- 1
  #####
  ## Below we use the dqrng R package which is faster for random sampling of integers, however, it is
  ## exactly the same as saying: data$user[sample(1:observed, 5, replace = T)]
  #####

  past_senders <- collapse::funique(sender[1:(observed[i]-1)]) # starting a vector to store all the unique senders in the past event sequence
  past_reciever <- collapse::funique(receiver[1:(observed[i]-1)]) # starting a vector to store all the  unique receivers in the past event sequence

  if (observed[i] == 1) { # Following prior REM software, if the first selected event is the first event, then
    # the observed event cannot have any controls!
    events_sampling[[i]] <- events[i]
    enough_actors <- FALSE ### a logical vector to make sure we have enough actors to create

  } else {
    ################################################################################
    #### Null Events Dataframe
    ################################################################################
    # Note: the riskset is the full Cartesian plot at time i, however, due to the computational
    # time that comes with fully creating and updating the risk set (for instance, once we reach 10000 unique
    # senders and 10000 unique receivers, the riskset set is 100 million events (which becomes unmanageable)).
    # To circumvent this, we sample from the x (sender) and y (receiver) axis.
    nullEventsi <- data.table::data.table( # combine this frame with the prior sampled ones
      sender = dqrng::dqsample(past_senders, m, replace = T), # sample from the past senders (i.e., only senders that have appears in the past)
      receiver = dqrng::dqsample(past_reciever, m, replace = T), # sample from the past receivers (i.e., only receivers that have appears in the past)
      EVENT = rep(observed[i], m), # current observed value should be the event they are stored with
      OBSERVED = rep(0, m), # OBSERVED is 0, meaning that they are not are real event
      time = rep(time[observed[i]], m), # current time should be from the event that are being sorted based on
      sampled = rep(1, m)
    ) # add a sampled value equal to one for these cases!

    ################################################################################
    #### Let's Check for Duplicated Events
    ################################################################################
    Duplicated <- collapse::fduplicated(nullEventsi) # the indicators of duplicated control events
    nDuplicated <- sum(Duplicated) # the number of duplicated cases

    enough_actors <- FALSE ### a logical vector to make sure we have enough actors to create
    ### a unique (n controls) x (n controls) (sender x receiver) null set,
    ### if there is not enough, then skip the null events

    # if we do have enough cases to fill the full null event: that is there is more actors in general than controls
    # For instance, if n_controls = 5, and we have at least 5 past senders or recievers, then we have enough actors to fill
    # that null event set
    if ((length(past_senders) >= n_controls) | (length(past_reciever) >= n_controls)) {
      enough_actors <- TRUE
    }
    # if we do not have at least n_control number of past senders or receivers, then let's check if we can fill the
    # full combinations for the n_controls requested. Again, if n_controls = 10, but if we have 4 past senders
    # and 3 past receivers, then we have at 12 different combinations.
    if (enough_actors == FALSE) {
      #### Provide the full Cartesian plot of past senders x past receivers
      possible_combinations <- base::expand.grid(past_senders, past_reciever)
      if (nrow(possible_combinations) >= n_controls) { # if the number of combinations >= n_controls
        enough_actors <- TRUE # we have enough actors
      }
    }

    if (enough_actors == FALSE) { # Now, if after those two checks, we still do not have enough past actors
      #### Provide the full Cartesian plot of past senders x past receivers
      possible_combinations <- data.table::data.table(base::expand.grid(
        sender = past_senders,
        receiver = past_reciever
      ))

      nullEventsi <- data.table::data.table( # combine this frame with the prior sampled ones
        sender = base::as.character(possible_combinations$sender), # sample from the past senders (i.e., only senders that have appears in the past)
        receiver = base::as.character(possible_combinations$receiver), # sample from the past receivers (i.e., only receivers that have appears in the past)
        EVENT = rep(observed[i], nrow(possible_combinations)), # current observed value should be the event they are stored with
        OBSERVED = rep(0, nrow(possible_combinations)), # OBSERVED is 0, meaning that they are not are real event
        time = rep(time[observed[i]], nrow(possible_combinations)), # current time should be from the event that are being sorted based on
        sampled = rep(1, nrow(possible_combinations))
      ) # add a sampled value equal to one for these cases!
    }

    if (nDuplicated != 0 & enough_actors == TRUE) { # Now, if we do have enough actors, AND we have duplicated controls (that is,
      # out of our n controls, we do not have n unique controls)
      while (nDuplicated != 0) { # while we still have duplicated cases
        whichduplicated <- collapse::whichv(Duplicated, TRUE) # which values are duplicated?

        ##### Create a new data set, called redoNull, that contains the controls to replace the duplicated ones!
        redoNull <- data.table::data.table( # combine this frame with the prior sampled ones
          sender = dqrng::dqsample(past_senders, nDuplicated, replace = T), # sample from the past senders (i.e., only senders that have appears in the past)
          receiver = dqrng::dqsample(past_reciever, nDuplicated, replace = T), # sample from the past receivers (i.e., only receivers that have appears in the past)
          EVENT = rep(observed[i], nDuplicated), # current observed value should be the event they are stored with
          OBSERVED = rep(0, nDuplicated), # OBSERVED is 0, meaning that they are not are real event
          time = rep(time[observed[i]], nDuplicated), # current time should be from the event that are being sorted based on
          sampled = rep(1, nDuplicated)
        ) # add a sampled value equal to one for these cases!

        nullEventsi[whichduplicated] <- redoNull # replace the duplicated rows with the new "non-duplicated" rows
        Duplicated <- collapse::fduplicated(nullEventsi) # now, recheck for duplicated rows
        nDuplicated <- sum(Duplicated) # how many duplicated values do we have?
      }
    }



    ################################################################################################
    #### One Last Check for the Selected Control: Is the observed event included in the null events
    ################################################################################################
    fulleventsi <- data.table::rbindlist(
      list(
        events[i], # current event i
        nullEventsi
      ), # null control events
      use.names = TRUE
    ) # combine the dataframe together
    #### Now, we have to make sure that the observed event is not included in the null events
    Duplicated <- base::duplicated(fulleventsi, by = c("sender", "receiver")) # are any of these values duplicates of the observed event?
    nDuplicated <- sum(Duplicated) # the number of duplicated events (this should be equal to 1 if there is a duplicate; 0 else)
    if (nDuplicated != 0) { # if we do a duplicated case: that is, a control is the same event as the observed event

      if (enough_actors == FALSE) { # if we do not have enough cases to populate the full control risk set
        # remove the duplicated case from the control set
        fulleventsi <- fulleventsi[-collapse::whichv(Duplicated, TRUE)] # which values are duplicated?
      }

      if (enough_actors == TRUE) { # if we do have enough cases to populate the full control risk set
        duplicatedValue <- 1 # set this value to 1, this is considered a switch until we not longer have any duplicates
        while (duplicatedValue == 1) { # while we still have duplicated events
          duplicatedID <- collapse::whichv(Duplicated, TRUE)
          fulleventsi$sender[duplicatedID] <- dqrng::dqsample(past_senders, nDuplicated, replace = T)
          fulleventsi$receiver[duplicatedID] <- dqrng::dqsample(past_reciever, nDuplicated, replace = T)
          Duplicated <- base::duplicated(fulleventsi, by = c("sender", "receiver")) # are any of these values duplicates of the observed event?
          nDuplicated <- sum(Duplicated) # the number of duplicated events (this should be equal to 1 if there is a duplicate; 0 else)
          if (nDuplicated == 0) {
            duplicatedValue <- 0 # there are no longer duplicated events!
          }
        }
      }
    }


    ################################################################################
    #### Merging the final cases!
    ################################################################################
    events_sampling[[i]] <- fulleventsi
  }
  #################################################################################
  #### Now, for 2 to the rest of the observed events! (i.e., 2, 3, ....., n - 1, n)
  #################################################################################



  for (i in 2:length(observed)) { # for all events from 2 to the last event
    past_sendersi <- collapse::funique(sender[observed[i - 1]:(observed[i]-1)]) # getting the unique senders from the last observed event to the current
    different <-  collapse::`%==%`(fastmatch::`%fin%`(past_sendersi, past_senders),FALSE) # getting all the senders that have not previously occurred
    if (length(different) != 0) { # if there are different cases, then add them to the past senders vector
      past_senders[(length(past_senders) + 1):(length(past_senders) + length(past_sendersi[different]))] <- past_sendersi[different]
    }
    past_recieveri <- collapse::funique(receiver[observed[i - 1]:(observed[i]-1)]) # getting the unique receivers from the last observed event to the current
    different <-  collapse::`%==%`(fastmatch::`%fin%`(past_recieveri, past_reciever), FALSE) # getting all the receivers that have not previously occurred
    if (length(different) != 0) { # if there are different cases, then add them to the past receivers vector
      past_reciever[(length(past_reciever) + 1):(length(past_reciever) + length(past_recieveri[different]))] <- past_recieveri[different]
    }

    nullEventsi <- data.table::data.table( # combine this frame with the prior sampled ones
      sender = dqrng::dqsample(past_senders, m, replace = T), # sample from the past senders (i.e., only senders that have appears in the past)
      receiver = dqrng::dqsample(past_reciever, m, replace = T), # sample from the past receivers (i.e., only receivers that have appears in the past)
      EVENT = rep(observed[i], m), # current observed value should be the event they are stored with
      OBSERVED = rep(0, m), # OBSERVED is 0, meaning that they are not are real event
      time = rep(time[observed[i]], m), # current time should be from the event that are being sorted based on
      sampled = rep(1, m)
    ) # add a sampled value equal to one for these cases!

    Duplicated <- collapse::fduplicated(nullEventsi) # the indicators of duplicated control events
    nDuplicated <- sum(Duplicated) # the number of duplicated cases
    # if we do not have at least n_control number of past senders or receivers, then let's check if we can fill the
    # full combinations for the n_controls requested. Again, if n_controls = 10, but if we have 4 past senders
    # and 3 past receivers, then we have at 12 different combinations.
    if (enough_actors == FALSE) { # if
      if ((length(past_senders) * length(past_reciever)) > (n_controls - 1)) { # f we have enough actors
        enough_actors <- TRUE # change value to TRUE
      }

      if (enough_actors == FALSE & nDuplicated != 0) { # if value is still flase and we have duplicated events
        #### Provide the full Cartesian plot of past senders x past receivers
        possible_combinations <- data.table::data.table(base::expand.grid(
          sender = past_senders,
          receiver = past_reciever
        ))

        nullEventsi <- data.table::data.table( # combine this frame with the prior sampled ones
          sender = base::as.character(possible_combinations$sender), # sample from the past senders (i.e., only senders that have appears in the past)
          receiver = base::as.character(possible_combinations$receiver), # sample from the past receivers (i.e., only receivers that have appears in the past)
          EVENT = rep(observed[i], nrow(possible_combinations)), # current observed value should be the event they are stored with
          OBSERVED = rep(0, nrow(possible_combinations)), # OBSERVED is 0, meaning that they are not are real event
          time = rep(time[observed[i]], nrow(possible_combinations)), # current time should be from the event that are being sorted based on
          sampled = rep(1, nrow(possible_combinations))
        ) # add a sampled value equal to one for these cases!
      }
    }

    if (nDuplicated != 0 & enough_actors == TRUE) { # Now, if we do have enough actors, AND we have duplicated controls (that is,
      # out of our n controls, we do not have n unique controls)
      while (nDuplicated != 0) { # while we still have duplicated cases
        whichduplicated <- collapse::whichv(Duplicated, TRUE) # which values are duplicated?

        redoNull <- data.table::data.table( # combine this frame with the prior sampled ones
          sender = dqrng::dqsample(past_senders, nDuplicated, replace = T), # sample from the past senders (i.e., only senders that have appears in the past)
          receiver = dqrng::dqsample(past_reciever, nDuplicated, replace = T), # sample from the past receivers (i.e., only receivers that have appears in the past)
          EVENT = rep(observed[i], nDuplicated), # current observed value should be the event they are stored with
          OBSERVED = rep(0, nDuplicated), # OBSERVED is 0, meaning that they are not are real event
          time = rep(time[observed[i]], nDuplicated), # current time should be from the event that are being sorted based on
          sampled = rep(1, nDuplicated)
        ) # add a sampled value equal to one for these cases!

        nullEventsi[whichduplicated] <- redoNull # replace the duplicated rows with the new "non-duplicated" rows
        Duplicated <- collapse::fduplicated(nullEventsi) # now, recheck for duplicated rows
        nDuplicated <- sum(Duplicated) # how many duplicated values do we have?
      }
    }

    ################################################################################################
    #### One Last Check for the Selected Control: Is the observed event included in the null events
    ################################################################################################
    fulleventsi <- data.table::rbindlist(
      list(
        events[i], # current event i
        nullEventsi
      ), # null control events
      use.names = TRUE
    ) # combine the dataframe together
    #### Now, we have to make sure that the observed event is not included in the null events
    Duplicated <- base::duplicated(fulleventsi, by = c("sender", "receiver")) # are any of these values duplicates of the observed event?
    nDuplicated <- sum(Duplicated) # the number of duplicated events (this should be equal to 1 if there is a duplicate; 0 else)
    if (nDuplicated != 0) { # if we do a duplicated case: that is, a control is the same event as the observed event

      if (enough_actors == FALSE) { # if we do not have enough cases to populate the full control risk set
        # remove the duplicated case from the control set
        fulleventsi <- fulleventsi[-collapse::whichv(Duplicated, TRUE)] # which values are duplicated?
      }

      if (enough_actors == TRUE) { # if we do have enough cases to populate the full control risk set
        duplicatedValue <- 1 # set this value to 1, this is considered a switch until we not longer have any duplicates
        while (duplicatedValue == 1) { # while we still have duplicated events
          duplicatedID <- collapse::whichv(Duplicated, TRUE)
          fulleventsi$sender[duplicatedID] <- dqrng::dqsample(past_senders, nDuplicated, replace = T)
          fulleventsi$receiver[duplicatedID] <- dqrng::dqsample(past_reciever, nDuplicated, replace = T)
          Duplicated <- base::duplicated(fulleventsi, by = c("sender", "receiver")) # are any of these values duplicates of the observed event?
          nDuplicated <- sum(Duplicated) # the number of duplicated events (this should be equal to 1 if there is a duplicate; 0 else)
          if (nDuplicated == 0) {
            duplicatedValue <- 0 # there are no longer duplicated events!
          }
        }
      }
    }


    ################################################################################
    #### Merging the final cases!
    ################################################################################
    events_sampling[[i]] <- fulleventsi
    rm(list = c("past_sendersi", "different", "past_recieveri")) # clearing memory for future computations
  }
  }

  #creating a dynamic risk set based on the user time inputs and difference values
  if(time_dependent == TRUE){

    if(is.null(timeDV) | is.null(timeDif)){#checking the dimensions of the user inputs
      base::stop("Error: You requested that a dynamic time risk set to be created. However
                 you did not include a time difference vector or did not specify
                 a timeDif value")
    }
    if(length(timeDif) != 1){#checking the dimensions of the user inputs
      base::stop("Error: The time difference value should be a scalar.")
    }
    if(length(timeDV) != length(sender)){ #checking the dimensions of the user inputs
      base::stop("Error: The length of the provided timeDV vector is not the same length as the events dataset") # stop computation and tell the user
    }

    #### Step 3: Randomly Sample Control Events (decide an n for selection)
    ## we should store all these values into a list
    #### Now lets sample m possible controls per event
    i <- 1
    #####
    ## Below we use the dqrng R package which is faster for random sampling of integers, however, it is
    ## exactly the same as saying: data$user[sample(1:observed, 5, replace = T)]
    #####


    #past_senders <- collapse::funique(sender[1:(observed[i]-1)]) # starting a vector to store all the unique senders in the past event sequence
    #past_reciever <- collapse::funique(receiver[1:(observed[i]-1)]) # starting a vector to store all the  unique receivers in the past event sequence


    if (observed[i] == 1) { # Following prior REM software, if the first selected event is the first event, then
      # the observed event cannot have any controls!
      events_sampling[[i]] <- events[i]
      enough_actors <- FALSE ### a logical vector to make sure we have enough actors to create

    } else {
      #------------------------------------#
      #updated earlier for reference with respect to event i not event i - 1
      #------------------------------------#
      time_Update <- timeDV[observed[i]] - timeDif #updated this
      if(time_Update < 0){
        time_Update <- 0
      }
      timesGood <- (time_Update <= timeDV & timeDV < timeDV[observed[i]])
      past_senders <- collapse::funique(sender[timesGood]) # starting a vector to store all the unique senders in the past event sequence
      past_reciever <- collapse::funique( receiver[timesGood]) # starting a vector to store all the unique senders in the past event sequence


      ################################################################################
      #### Null Events Dataframe
      ################################################################################
      # Note: the riskset is the full Cartesian plot at time i, however, due to the computational
      # time that comes with fully creating and updating the risk set (for instance, once we reach 10000 unique
      # senders and 10000 unique receivers, the riskset set is 100 million events (which becomes unmanageable)).
      # To circumvent this, we sample from the x (sender) and y (receiver) axis.
      nullEventsi <- data.table::data.table( # combine this frame with the prior sampled ones
        sender = dqrng::dqsample(past_senders, m, replace = T), # sample from the past senders (i.e., only senders that have appears in the past)
        receiver = dqrng::dqsample(past_reciever, m, replace = T), # sample from the past receivers (i.e., only receivers that have appears in the past)
        EVENT = rep(observed[i], m), # current observed value should be the event they are stored with
        OBSERVED = rep(0, m), # OBSERVED is 0, meaning that they are not are real event
        time = rep(time[observed[i]], m), # current time should be from the event that are being sorted based on
        sampled = rep(1, m)
      ) # add a sampled value equal to one for these cases!

      ################################################################################
      #### Let's Check for Duplicated Events
      ################################################################################
      Duplicated <- collapse::fduplicated(nullEventsi) # the indicators of duplicated control events
      nDuplicated <- sum(Duplicated) # the number of duplicated cases

      enough_actors <- FALSE ### a logical vector to make sure we have enough actors to create
      ### a unique (n controls) x (n controls) (sender x receiver) null set,
      ### if there is not enough, then skip the null events

      # if we do have enough cases to fill the full null event: that is there is more actors in general than controls
      # For instance, if n_controls = 5, and we have at least 5 past senders or recievers, then we have enough actors to fill
      # that null event set
      if ((length(past_senders) >= n_controls) | (length(past_reciever) >= n_controls)) {
        enough_actors <- TRUE
      }
      # if we do not have at least n_control number of past senders or receivers, then let's check if we can fill the
      # full combinations for the n_controls requested. Again, if n_controls = 10, but if we have 4 past senders
      # and 3 past receivers, then we have at 12 different combinations.
      if (enough_actors == FALSE) {
        #### Provide the full Cartesian plot of past senders x past receivers
        possible_combinations <- base::expand.grid(past_senders, past_reciever)
        if (nrow(possible_combinations) >= n_controls) { # if the number of combinations >= n_controls
          enough_actors <- TRUE # we have enough actors
        }
      }

      if (enough_actors == FALSE) { # Now, if after those two checks, we still do not have enough past actors
        #### Provide the full Cartesian plot of past senders x past receivers
        possible_combinations <- data.table::data.table(base::expand.grid(
          sender = past_senders,
          receiver = past_reciever
        ))

        nullEventsi <- data.table::data.table( # combine this frame with the prior sampled ones
          sender = base::as.character(possible_combinations$sender), # sample from the past senders (i.e., only senders that have appears in the past)
          receiver = base::as.character(possible_combinations$receiver), # sample from the past receivers (i.e., only receivers that have appears in the past)
          EVENT = rep(observed[i], nrow(possible_combinations)), # current observed value should be the event they are stored with
          OBSERVED = rep(0, nrow(possible_combinations)), # OBSERVED is 0, meaning that they are not are real event
          time = rep(time[observed[i]], nrow(possible_combinations)), # current time should be from the event that are being sorted based on
          sampled = rep(1, nrow(possible_combinations))
        ) # add a sampled value equal to one for these cases!
      }

      if (nDuplicated != 0 & enough_actors == TRUE) { # Now, if we do have enough actors, AND we have duplicated controls (that is,
        # out of our n controls, we do not have n unique controls)
        while (nDuplicated != 0) { # while we still have duplicated cases
          whichduplicated <- collapse::whichv(Duplicated, TRUE) # which values are duplicated?

          ##### Create a new data set, called redoNull, that contains the controls to replace the duplicated ones!
          redoNull <- data.table::data.table( # combine this frame with the prior sampled ones
            sender = dqrng::dqsample(past_senders, nDuplicated, replace = T), # sample from the past senders (i.e., only senders that have appears in the past)
            receiver = dqrng::dqsample(past_reciever, nDuplicated, replace = T), # sample from the past receivers (i.e., only receivers that have appears in the past)
            EVENT = rep(observed[i], nDuplicated), # current observed value should be the event they are stored with
            OBSERVED = rep(0, nDuplicated), # OBSERVED is 0, meaning that they are not are real event
            time = rep(time[observed[i]], nDuplicated), # current time should be from the event that are being sorted based on
            sampled = rep(1, nDuplicated)
          ) # add a sampled value equal to one for these cases!

          nullEventsi[whichduplicated] <- redoNull # replace the duplicated rows with the new "non-duplicated" rows
          Duplicated <- collapse::fduplicated(nullEventsi) # now, recheck for duplicated rows
          nDuplicated <- sum(Duplicated) # how many duplicated values do we have?
        }
      }



      ################################################################################################
      #### One Last Check for the Selected Control: Is the observed event included in the null events
      ################################################################################################
      fulleventsi <- data.table::rbindlist(
        list(
          events[i], # current event i
          nullEventsi
        ), # null control events
        use.names = TRUE
      ) # combine the dataframe together
      #### Now, we have to make sure that the observed event is not included in the null events
      Duplicated <- base::duplicated(fulleventsi, by = c("sender", "receiver")) # are any of these values duplicates of the observed event?
      nDuplicated <- sum(Duplicated) # the number of duplicated events (this should be equal to 1 if there is a duplicate; 0 else)
      if (nDuplicated != 0) { # if we do a duplicated case: that is, a control is the same event as the observed event

        if (enough_actors == FALSE) { # if we do not have enough cases to populate the full control risk set
          # remove the duplicated case from the control set
          fulleventsi <- fulleventsi[-collapse::whichv(Duplicated, TRUE)] # which values are duplicated?
        }

        if (enough_actors == TRUE) { # if we do have enough cases to populate the full control risk set
          duplicatedValue <- 1 # set this value to 1, this is considered a switch until we not longer have any duplicates
          while (duplicatedValue == 1) { # while we still have duplicated events
            duplicatedID <- collapse::whichv(Duplicated, TRUE)
            fulleventsi$sender[duplicatedID] <- dqrng::dqsample(past_senders, nDuplicated, replace = T)
            fulleventsi$receiver[duplicatedID] <- dqrng::dqsample(past_reciever, nDuplicated, replace = T)
            Duplicated <- base::duplicated(fulleventsi, by = c("sender", "receiver")) # are any of these values duplicates of the observed event?
            nDuplicated <- sum(Duplicated) # the number of duplicated events (this should be equal to 1 if there is a duplicate; 0 else)
            if (nDuplicated == 0) {
              duplicatedValue <- 0 # there are no longer duplicated events!
            }
          }
        }
      }


      ################################################################################
      #### Merging the final cases!
      ################################################################################
      events_sampling[[i]] <- fulleventsi
    }
    #################################################################################
    #### Now, for 2 to the rest of the observed events! (i.e., 2, 3, ....., n - 1, n)
    #################################################################################



    for (i in 2:length(observed)) { # for all events from 2 to the last event
      enough_actors<-FALSE
      #------------------------------------#
      #updated earlier for reference with respect to event i not event i - 1
      #------------------------------------#
      time_Update <- timeDV[observed[i]] - timeDif #updated this
      if (time_Update < 0) {
        time_Update <- 0
      }
      timesGood <- (time_Update <= timeDV & timeDV < timeDV[observed[i]])
      past_senders <- collapse::funique(sender[timesGood])
      past_reciever <- collapse::funique(receiver[timesGood])
      if(length(past_senders) == 0 | length(   past_reciever) == 0){
        events_sampling[[i]] <- events[i]
        enough_actors <- FALSE
        next
      }
      nullEventsi <- data.table::data.table( # combine this frame with the prior sampled ones
        sender = dqrng::dqsample(past_senders, m, replace = T), # sample from the past senders (i.e., only senders that have appears in the past)
        receiver = dqrng::dqsample(past_reciever, m, replace = T), # sample from the past receivers (i.e., only receivers that have appears in the past)
        EVENT = rep(observed[i], m), # current observed value should be the event they are stored with
        OBSERVED = rep(0, m), # OBSERVED is 0, meaning that they are not are real event
        time = rep(time[observed[i]], m), # current time should be from the event that are being sorted based on
        sampled = rep(1, m)
      ) # add a sampled value equal to one for these cases!

      Duplicated <- collapse::fduplicated(nullEventsi) # the indicators of duplicated control events
      nDuplicated <- sum(Duplicated) # the number of duplicated cases
      # if we do not have at least n_control number of past senders or receivers, then let's check if we can fill the
      # full combinations for the n_controls requested. Again, if n_controls = 10, but if we have 4 past senders
      # and 3 past receivers, then we have at 12 different combinations.
      if (enough_actors == FALSE) { # if
        if ((length(past_senders) * length(past_reciever)) > (n_controls - 1)) { # f we have enough actors
          enough_actors <- TRUE # change value to TRUE
        }

        if (enough_actors == FALSE & nDuplicated != 0) { # if value is still flase and we have duplicated events
          #### Provide the full Cartesian plot of past senders x past receivers
          possible_combinations <- data.table::data.table(base::expand.grid(
            sender = past_senders,
            receiver = past_reciever
          ))

          nullEventsi <- data.table::data.table( # combine this frame with the prior sampled ones
            sender = base::as.character(possible_combinations$sender), # sample from the past senders (i.e., only senders that have appears in the past)
            receiver = base::as.character(possible_combinations$receiver), # sample from the past receivers (i.e., only receivers that have appears in the past)
            EVENT = rep(observed[i], nrow(possible_combinations)), # current observed value should be the event they are stored with
            OBSERVED = rep(0, nrow(possible_combinations)), # OBSERVED is 0, meaning that they are not are real event
            time = rep(time[observed[i]], nrow(possible_combinations)), # current time should be from the event that are being sorted based on
            sampled = rep(1, nrow(possible_combinations))
          ) # add a sampled value equal to one for these cases!
        }
      }

      if (nDuplicated != 0 & enough_actors == TRUE) { # Now, if we do have enough actors, AND we have duplicated controls (that is,
        # out of our n controls, we do not have n unique controls)
        while (nDuplicated != 0) { # while we still have duplicated cases
          whichduplicated <- collapse::whichv(Duplicated, TRUE) # which values are duplicated?

          redoNull <- data.table::data.table( # combine this frame with the prior sampled ones
            sender = dqrng::dqsample(past_senders, nDuplicated, replace = T), # sample from the past senders (i.e., only senders that have appears in the past)
            receiver = dqrng::dqsample(past_reciever, nDuplicated, replace = T), # sample from the past receivers (i.e., only receivers that have appears in the past)
            EVENT = rep(observed[i], nDuplicated), # current observed value should be the event they are stored with
            OBSERVED = rep(0, nDuplicated), # OBSERVED is 0, meaning that they are not are real event
            time = rep(time[observed[i]], nDuplicated), # current time should be from the event that are being sorted based on
            sampled = rep(1, nDuplicated)
          ) # add a sampled value equal to one for these cases!

          nullEventsi[whichduplicated] <- redoNull # replace the duplicated rows with the new "non-duplicated" rows
          Duplicated <- collapse::fduplicated(nullEventsi) # now, recheck for duplicated rows
          nDuplicated <- sum(Duplicated) # how many duplicated values do we have?
        }
      }

      ################################################################################################
      #### One Last Check for the Selected Control: Is the observed event included in the null events
      ################################################################################################
      fulleventsi <- data.table::rbindlist(
        list(
          events[i], # current event i
          nullEventsi
        ), # null control events
        use.names = TRUE
      ) # combine the dataframe together
      #### Now, we have to make sure that the observed event is not included in the null events
      Duplicated <- base::duplicated(fulleventsi, by = c("sender", "receiver")) # are any of these values duplicates of the observed event?
      nDuplicated <- sum(Duplicated) # the number of duplicated events (this should be equal to 1 if there is a duplicate; 0 else)
      if (nDuplicated != 0) { # if we do a duplicated case: that is, a control is the same event as the observed event

        if (enough_actors == FALSE) { # if we do not have enough cases to populate the full control risk set
          # remove the duplicated case from the control set
          fulleventsi <- fulleventsi[-collapse::whichv(Duplicated, TRUE)] # which values are duplicated?
        }

        if (enough_actors == TRUE) { # if we do have enough cases to populate the full control risk set
          duplicatedValue <- 1 # set this value to 1, this is considered a switch until we not longer have any duplicates
          while (duplicatedValue == 1) { # while we still have duplicated events
            duplicatedID <- collapse::whichv(Duplicated, TRUE)
            fulleventsi$sender[duplicatedID] <- dqrng::dqsample(past_senders, nDuplicated, replace = T)
            fulleventsi$receiver[duplicatedID] <- dqrng::dqsample(past_reciever, nDuplicated, replace = T)
            Duplicated <- base::duplicated(fulleventsi, by = c("sender", "receiver")) # are any of these values duplicates of the observed event?
            nDuplicated <- sum(Duplicated) # the number of duplicated events (this should be equal to 1 if there is a duplicate; 0 else)
            if (nDuplicated == 0) {
              duplicatedValue <- 0 # there are no longer duplicated events!
            }
          }
        }
      }


      ################################################################################
      #### Merging the final cases!
      ################################################################################
      events_sampling[[i]] <- fulleventsi
      rm(list = c("past_senders",  "past_reciever")) # clearing memory for future computations
    }


  }

  ################################################################################
  #### Merging all sampled events to make the final outputted file!
  ################################################################################
  #base::cat("Combining all sampled files together.......") # outputting status to user
  eventSet <- data.table::rbindlist(events_sampling) # combine and stack all the sampled events
  eventSet <- eventSet[, c("sender", "receiver", "time", "EVENT", "OBSERVED")]
  colnames(eventSet) <- c("sender", "receiver", "time", "sequenceID", "observed")
  return(eventSet) # output the file to the user!
}
