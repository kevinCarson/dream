pkgname <- "redr"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('redr')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("om.constraint")
### * om.constraint

flush(stderr()); flush(stdout())

### Name: om.constraint
### Title: Compute Burt's (1992) Constraint for Ego Networks
### Aliases: om.constraint

### ** Examples


# For this example, we recreate the ego network provided in Burt (1992: 56):
BurtEgoNet <- matrix(c(
  0,1,0,0,1,1,1,
 1,0,0,1,0,0,1,
 0,0,0,0,0,0,1,
 0,1,0,0,0,0,1,
 1,0,0,0,0,0,1,
 1,0,0,0,0,0,1,
 1,1,1,1,1,1,0),
 nrow = 7, ncol = 7)
colnames(BurtEgoNet) <- rownames(BurtEgoNet) <- c("A", "B", "C", "D", "E",
                                                 "F", "ego")
#the constraint value for the ego replicates that provided in Burt (1992: 56)
om.constraint(BurtEgoNet)





cleanEx()
nameEx("om.effective")
### * om.effective

flush(stderr()); flush(stdout())

### Name: om.effective
### Title: Compute Burt's (1992) Effective Size for Ego Networks from a
###   Sociomatrix
### Aliases: om.effective

### ** Examples

# For this example, we recreate the ego network provided in Borgatti (1997):
BorgattiEgoNet <- matrix(
 c(0,1,0,0,0,0,0,0,1,
   1,0,0,0,0,0,0,0,1,
   0,0,0,1,0,0,0,0,1,
   0,0,1,0,0,0,0,0,1,
   0,0,0,0,0,1,0,0,1,
  0,0,0,0,1,0,0,0,1,
  0,0,0,0,0,0,0,1,1,
   0,0,0,0,0,0,1,0,1,
   1,1,1,1,1,1,1,1,0),
 nrow = 9, ncol = 9, byrow = TRUE)
colnames(BorgattiEgoNet) <- rownames(BorgattiEgoNet) <- c("A", "B", "C",
                                                         "D", "E", "F",
                                                        "G", "H", "ego")
#the effective size value for the ego replicates that provided in Borgatti (1997)
om.effective(BorgattiEgoNet)

# For this example, we recreate the ego network provided in Burt (1992: 56):
BurtEgoNet <- matrix(c(
  0,1,0,0,1,1,1,
 1,0,0,1,0,0,1,
 0,0,0,0,0,0,1,
 0,1,0,0,0,0,1,
 1,0,0,0,0,0,1,
 1,0,0,0,0,0,1,
 1,1,1,1,1,1,0),
 nrow = 7, ncol = 7)
colnames(BurtEgoNet) <- rownames(BurtEgoNet) <- c("A", "B", "C", "D", "E",
                                                 "F", "ego")
#the effective size value for the ego replicates that provided in Burt (1992: 56)
om.effective(BurtEgoNet)



cleanEx()
nameEx("om.npaths")
### * om.npaths

flush(stderr()); flush(stdout())

### Name: om.npaths
### Title: Compute the Number of Paths of Length K
### Aliases: om.npaths

### ** Examples


# For this example, we generate a random one-mode graph with the sna package.
#creating the random network with 10 actors
set.seed(9999)
rnet <- matrix(sample(c(0,1), 10*10, replace = TRUE, prob = c(0.8,0.2)),
               nrow = 10, ncol = 10, byrow = TRUE)
diag(rnet) <- 0 #setting self ties to 0
#counting the paths of length 2
om.npaths(rnet, k = 2)
#counting the paths of length 5
om.npaths(rnet, k = 5)



cleanEx()
nameEx("om.pi.brokerage")
### * om.pi.brokerage

flush(stderr()); flush(stdout())

### Name: om.pi.brokerage
### Title: Compute Interpersonal Cultural Brokerage Based on Leal (2025)
### Aliases: om.pi.brokerage

### ** Examples


# For this example, we recreate Figure 3 in Leal (2025)
LealNet <- matrix( c(
 0,1,0,0,0,0,0,
 1,0,1,1,0,0,0,
 0,1,0,0,1,1,0,
 0,1,0,0,1,0,0,
 0,0,1,1,0,0,0,
 0,0,1,0,0,0,1,
 0,0,0,0,0,1,0),
 nrow = 7, ncol = 7, byrow = TRUE)

colnames(LealNet) <- rownames(LealNet) <- c("A", "B", "C","D",
                                           "E", "F", "G")
categorical_variable <- c(0,0,1,0,0,0,0)
#These values are exactly the same as reported by Leal (2025)
om.pi.brokerage(LealNet,
   symmetric = TRUE,
   g.mem = categorical_variable)






cleanEx()
nameEx("rem.riskset.om")
### * rem.riskset.om

flush(stderr()); flush(stdout())

### Name: rem.riskset.om
### Title: Create Dynamic One-Mode Risk Sets for Relational Event Sequences
### Aliases: rem.riskset.om

### ** Examples

# A random one-mode relational event sequence
events <- data.table::data.table(time = 1:18,
                                eventID = 1:18,
                                sender = c("A", "B", "C",
                                           "A", "D", "E",
                                           "F", "B", "A",
                                           "F", "D", "B",
                                           "G", "B", "D",
                                          "H", "A", "D"),
                               target = c("B", "C", "D",
                                          "E", "A", "F",
                                          "D", "A", "C",
                                          "G", "B", "C",
                                          "H", "J", "A",
                                          "F", "C", "B"))

# Creating a one-mode relational risk set with p = 1.00 (all true events) and 5 controls
eventSet <- rem.riskset.om(data = events,
                      time = events$time,
                      eventID = events$eventID,
                      sender = events$sender,
                      receiver = events$target,
                      p_samplingobserved = 1.00,
                      n_controls = 5,
                      seed = 9999)

# Creating a event-dependent one-mode relational risk set with p = 1.00 (all
# true events) and 3 controls based upon the past 5 events prior to the current event.
events$timeseq <- 1:nrow(events)
eventSetT <- rem.riskset.om(data = events,
                       time = events$time,
                       eventID = events$eventID,
                       sender = events$sender,
                       receiver = events$target,
                       p_samplingobserved = 1.00,
                       time_dependent = TRUE,
                       timeDV = events$timeseq,
                       timeDif = 5,
                       n_controls = 3,
                       seed = 9999)



cleanEx()
nameEx("rem.riskset.tm")
### * rem.riskset.tm

flush(stderr()); flush(stdout())

### Name: rem.riskset.tm
### Title: Create Dynamic Two-Mode Risk Sets for Relational Event Sequences
### Aliases: rem.riskset.tm

### ** Examples


data("WikiEvent2018.first100k")
WikiEvent2018.first100k$time <- as.numeric(WikiEvent2018.first100k$time)
### Creating the EventSet By Employing Case-Control Sampling With M = 10 and
### Sampling from the Observed Event Sequence with P = 0.01
EventSet <- rem.riskset.tm(
  data = WikiEvent2018.first100k, # The Event Dataset
  time = WikiEvent2018.first100k$time, # The Time Variable
  eventID = WikiEvent2018.first100k$eventID, # The Event Sequence Variable
  sender = WikiEvent2018.first100k$user, # The Sender Variable
  receiver = WikiEvent2018.first100k$article, # The Receiver Variable
  p_samplingobserved = 0.01, # The Probability of Selection
  n_controls = 10, # The Number of Controls to Sample from the Full Risk Set
  seed = 9999) # The Seed for Replication


### Creating A New EventSet with more observed events and less control events
### Sampling from the Observed Event Sequence with P = 0.10
### Employing Case-Control Sampling With M = 2
EventSet1 <- rem.riskset.tm(
  data = WikiEvent2018.first100k, # The Event Dataset
  time = WikiEvent2018.first100k$time, # The Time Variable
  eventID = WikiEvent2018.first100k$eventID, # The Event Sequence Variable
  sender = WikiEvent2018.first100k$user, # The Sender Variable
  receiver = WikiEvent2018.first100k$article, # The Receiver Variable
  p_samplingobserved = 0.02, # The Probability of Selection
  n_controls = 2, # The Number of Controls to Sample from the Full Risk Set
  seed = 9999) # The Seed for Replication





cleanEx()
nameEx("rem.simulate")
### * rem.simulate

flush(stderr()); flush(stdout())

### Name: rem.simulate
### Title: Simulate A One-Mode Random Relational Event Sequence
### Aliases: rem.simulate

### ** Examples

#Creating a random relational sequence with 5 actors and 25 events
rem1<- rem.simulate(n_actors = 25,
                     n_events = 1000,
                     inertia = TRUE,
                     inertia_p = 0.12,
                     recip = TRUE,
                     recip_p = 0.08,
                     sender_outdegree = TRUE,
                     sender_outdegree_p = 0.09,
                     target_indegree = TRUE,
                     target_indegree_p = 0.05,
                     assort = TRUE,
                     assort_p = -0.01,
                     trans_trips = TRUE,
                     trans_trips_p = 0.09,
                     three_cycles = TRUE,
                     three_cycles_p = 0.04,
                     starting_events = NULL,
                     returnStats = TRUE)
rem1

#Creating a random relational sequence with 100 actors and 1000 events with
#only inertia and reciprocity
rem2 <- rem.simulate(n_actors = 100,
                     n_events = 1000,
                     inertia = TRUE,
                     inertia_p = 0.12,
                     recip = TRUE,
                     recip_p = 0.08,
                     returnStats = TRUE)
rem2

#Creating a random relational sequence based on the starting sequence with
#only inertia and reciprocity
rem3 <- rem.simulate(n_actors = 100, #does not matter can be any value, this is
                                    #overridden by the starting event sequence
                    n_events = 100,
                    inertia = TRUE,
                    inertia_p = 0.12,
                    recip = TRUE,
                    recip_p = 0.08,
                    #a random starting event sequence
                    starting_events = matrix(c(1:10, 10:1),
                    nrow = 10, ncol = 2,  byrow = FALSE),
                    returnStats = TRUE)
rem3



cleanEx()
nameEx("rem.stat.ISP")
### * rem.stat.ISP

flush(stderr()); flush(stdout())

### Name: rem.stat.ISP
### Title: Compute Incoming Shared Partner Statistic for Relational Event
###   Sequences
### Aliases: rem.stat.ISP

### ** Examples

events <- data.table::data.table(time = 1:18,
                                eventID = 1:18,
                                sender = c("A", "B", "C",
                                           "A", "D", "E",
                                           "F", "B", "A",
                                           "F", "D", "B",
                                           "G", "B", "D",
                                          "H", "A", "D"),
                               target = c("B", "C", "D",
                                          "E", "A", "F",
                                          "D", "A", "C",
                                          "G", "B", "C",
                                          "H", "J", "A",
                                          "F", "C", "B"))

eventSet <- rem.riskset.om(data = events,
                      time = events$time,
                      eventID = events$eventID,
                      sender = events$sender,
                      receiver = events$target,
                      p_samplingobserved = 1.00,
                      n_controls = 1,
                      seed = 9999)

# Computing Incoming Shared Partners Statistic without the sliding windows framework
eventSet$ISP <- rem.stat.ISP(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   dyadic_weight = 0,
   Lerneretal_2013 = FALSE)

# Computing Incoming Shared Partners Statistic with the sliding windows framework
eventSet$ISP_SW <- rem.stat.ISP(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   processed_seqIDs = eventSet$sequenceID,
   dyadic_weight = 0,
   sliding_window = TRUE,
   Lerneretal_2013 = FALSE)

#The results with and without the sliding windows are the same (see correlation below).
#Using the sliding windows method is recommended when the data are 'big' so
#that memory allotment is more efficient.
cor(eventSet$ISP , eventSet$ISP_SW)

# Computing Reciprocity Statistics with the counts of events being returned
eventSet$ISPC <- rem.stat.ISP(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   dyadic_weight = 0,
   sliding_window = FALSE,
   counts = TRUE,
   Lerneretal_2013 = FALSE)

cbind(eventSet$ISP,
     eventSet$ISP_SW,
     eventSet$ISPC)



cleanEx()
nameEx("rem.stat.ITP")
### * rem.stat.ITP

flush(stderr()); flush(stdout())

### Name: rem.stat.ITP
### Title: Compute Incoming Two Path Statistic for Relational Event
###   Sequences
### Aliases: rem.stat.ITP

### ** Examples

events <- data.table::data.table(time = 1:18,
                                eventID = 1:18,
                                sender = c("A", "B", "C",
                                           "A", "D", "E",
                                           "F", "B", "A",
                                           "F", "D", "B",
                                           "G", "B", "D",
                                          "H", "A", "D"),
                               target = c("B", "C", "D",
                                          "E", "A", "F",
                                          "D", "A", "C",
                                          "G", "B", "C",
                                          "H", "J", "A",
                                          "F", "C", "B"))

eventSet <- rem.riskset.om(data = events,
                      time = events$time,
                      eventID = events$eventID,
                      sender = events$sender,
                      receiver = events$target,
                      p_samplingobserved = 1.00,
                      n_controls = 1,
                      seed = 9999)

# Computing Incoming Two Paths Statistics without the sliding windows framework
eventSet$ITP <- rem.stat.ITP(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   dyadic_weight = 0,
   Lerneretal_2013 = FALSE)

# Computing Incoming Two Paths Statistics with the sliding windows framework
eventSet$ITP_SW <- rem.stat.ITP(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   processed_seqIDs = eventSet$sequenceID,
   dyadic_weight = 0,
   sliding_window = TRUE,
   Lerneretal_2013 = FALSE)

#The results with and without the sliding windows are the same (see correlation below).
#Using the sliding windows method is recommended when the data are 'big' so
#that memory allotment is more efficient.
cor(eventSet$ITP, eventSet$ITP_SW)

# Computing Reciprocity Statistics with the counts of events being returned
eventSet$ITPC <- rem.stat.ITP(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   dyadic_weight = 0,
   sliding_window = FALSE,
   counts = TRUE,
   Lerneretal_2013 = FALSE)

cbind(eventSet$ITP,
     eventSet$ITP_SW,
     eventSet$ITPC)



cleanEx()
nameEx("rem.stat.OSP")
### * rem.stat.OSP

flush(stderr()); flush(stdout())

### Name: rem.stat.OSP
### Title: Compute Outgoing Shared Partner Statistic for Relational Event
###   Sequences
### Aliases: rem.stat.OSP

### ** Examples

events <- data.table::data.table(time = 1:18,
                                eventID = 1:18,
                                sender = c("A", "B", "C",
                                           "A", "D", "E",
                                           "F", "B", "A",
                                           "F", "D", "B",
                                           "G", "B", "D",
                                          "H", "A", "D"),
                               target = c("B", "C", "D",
                                          "E", "A", "F",
                                          "D", "A", "C",
                                          "G", "B", "C",
                                          "H", "J", "A",
                                          "F", "C", "B"))

eventSet <- rem.riskset.om(data = events,
                      time = events$time,
                      eventID = events$eventID,
                      sender = events$sender,
                      receiver = events$target,
                      p_samplingobserved = 1.00,
                      n_controls = 1,
                      seed = 9999)

# Computing Outgoing Shared Partners Statistics without the sliding windows framework
eventSet$OSP <- rem.stat.OSP(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   dyadic_weight = 0,
   Lerneretal_2013 = FALSE)

# Computing Outgoing Shared Partners Statistics with the sliding windows framework
eventSet$OSP_SW <- rem.stat.OSP(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   processed_seqIDs = eventSet$sequenceID,
   dyadic_weight = 0,
   sliding_window = TRUE,
   Lerneretal_2013 = FALSE)

#The results with and without the sliding windows are the same (see correlation below).
#Using the sliding windows method is recommended when the data are 'big' so
#that memory allotment is more efficient.
cor(eventSet$OSP , eventSet$OSP_SW)

# Computing  Outgoing Shared Partners Statistics with the counts of events being returned
eventSet$OSP_C <- rem.stat.OSP(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   dyadic_weight = 0,
   sliding_window = FALSE,
   counts = TRUE,
   Lerneretal_2013 = FALSE)

cbind(eventSet$OSP,
     eventSet$OSP_SW,
     eventSet$OSP_C)



cleanEx()
nameEx("rem.stat.OTP")
### * rem.stat.OTP

flush(stderr()); flush(stdout())

### Name: rem.stat.OTP
### Title: Compute Outgoing Two Path Statistic for Relational Event
###   Sequences
### Aliases: rem.stat.OTP

### ** Examples

events <- data.table::data.table(time = 1:18,
                                eventID = 1:18,
                                sender = c("A", "B", "C",
                                           "A", "D", "E",
                                           "F", "B", "A",
                                           "F", "D", "B",
                                           "G", "B", "D",
                                          "H", "A", "D"),
                               target = c("B", "C", "D",
                                          "E", "A", "F",
                                          "D", "A", "C",
                                          "G", "B", "C",
                                          "H", "J", "A",
                                          "F", "C", "B"))

eventSet <- rem.riskset.om(data = events,
                      time = events$time,
                      eventID = events$eventID,
                      sender = events$sender,
                      receiver = events$target,
                      p_samplingobserved = 1.00,
                      n_controls = 1,
                      seed = 9999)

# Computing Outgoing Two Paths Statistics without the sliding windows framework
eventSet$OTP <- rem.stat.OTP(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   dyadic_weight = 0,
   Lerneretal_2013 = FALSE)

# Computing Outgoing Two Paths Statistics with the sliding windows framework
eventSet$OTP_SW <- rem.stat.OTP(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   processed_seqIDs = eventSet$sequenceID,
   dyadic_weight = 0,
   sliding_window = TRUE,
   Lerneretal_2013 = FALSE)

#The results with and without the sliding windows are the same (see correlation below).
#Using the sliding windows method is recommended when the data are 'big' so
#that memory allotment is more efficient.
cor(eventSet$OTP , eventSet$OTP_SW)

# Computing Reciprocity Statistics with the counts of events being returned
eventSet$OTPC <- rem.stat.OTP(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   dyadic_weight = 0,
   sliding_window = FALSE,
   counts = TRUE,
   Lerneretal_2013 = FALSE)

cbind(eventSet$OTP,
     eventSet$OTP_SW,
     eventSet$OTPC)



cleanEx()
nameEx("rem.stat.cycle4")
### * rem.stat.cycle4

flush(stderr()); flush(stdout())

### Name: rem.stat.cycle4
### Title: Compute The Four-Cycles Statistic for Two-Mode Relational Event
###   Sequences
### Aliases: rem.stat.cycle4

### ** Examples

data("WikiEvent2018.first100k")
WikiEvent2018 <- WikiEvent2018.first100k[1:10000,] #the first ten thousand events
WikiEvent2018$time <- as.numeric(WikiEvent2018$time) #making the variable numeric
### Creating the EventSet By Employing Case-Control Sampling With M = 5 and
### Sampling from the Observed Event Sequence with P = 0.01
EventSet <- rem.riskset.tm(
 data = WikiEvent2018, # The Event Dataset
 time = WikiEvent2018$time, # The Time Variable
 eventID = WikiEvent2018$eventID, # The Event Sequence Variable
 sender = WikiEvent2018$user, # The Sender Variable
 receiver = WikiEvent2018$article, # The Receiver Variable
 p_samplingobserved = 0.01, # The Probability of Selection
 n_controls = 5, # The Number of Controls to Sample from the Full Risk Set
 seed = 9999) # The Seed for Replication

#### Estimating the Four-Cycle Statistic Without the Sliding Windows Framework
EventSet$fourcycle <- rem.stat.cycle4(
   observed_time = WikiEvent2018$time,
   observed_sender = WikiEvent2018$user,
   observed_receiver = WikiEvent2018$article,
   processed_time = EventSet$time,
   processed_sender = EventSet$sender,
   processed_receiver = EventSet$receiver,
   halflife = 2.592e+09, #halflife parameter
   dyadic_weight = 0,
   Lerneretal_2013 = FALSE)

#### Estimating the Four-Cycle Statistic With the Sliding Windows Framework
EventSet$cycle4SW <- rem.stat.cycle4(
   observed_time = WikiEvent2018$time,
   observed_sender = WikiEvent2018$user,
   observed_receiver = WikiEvent2018$article,
   processed_time = EventSet$time,
   processed_sender = EventSet$sender,
   processed_receiver = EventSet$receiver,
   processed_seqIDs = EventSet$sequenceID,
   halflife = 2.592e+09, #halflife parameter
   dyadic_weight = 0,
   sliding_window = TRUE,
   Lerneretal_2013 = FALSE)

#The results with and without the sliding windows are the same (see correlation below).
#Using the sliding windows method is recommended when the data are 'big' so
#that memory allotment is more efficient.
cor(EventSet$fourcycle, EventSet$cycle4SW)

#### Estimating the Four-Cycle Statistic  with the Counts of Events Returned
EventSet$cycle4C <- rem.stat.cycle4(
   observed_time = WikiEvent2018$time,
   observed_sender = WikiEvent2018$user,
   observed_receiver = WikiEvent2018$article,
   processed_time = EventSet$time,
   processed_sender = EventSet$sender,
   processed_receiver = EventSet$receiver,
   processed_seqIDs = EventSet$sequenceID,
   halflife = 2.592e+09, #halflife parameter
   dyadic_weight = 0,
   sliding_window = FALSE,
   counts = TRUE,
   Lerneretal_2013 = FALSE)

cbind(EventSet$fourcycle,
     EventSet$cycle4SW,
     EventSet$cycle4C)




cleanEx()
nameEx("rem.stat.dyadCut")
### * rem.stat.dyadCut

flush(stderr()); flush(stdout())

### Name: rem.stat.dyadCut
### Title: Helper Function to Assist Researchers Finding Dyadic Weight
###   Cutoff Values
### Aliases: rem.stat.dyadCut

### ** Examples

#To replicate the example in the details section:
# with the Lerner et al. 2013 weighting function
rem.stat.dyadCut(halflife = 30,
                 timeValue = 100,
                 timeWidth = 60,
                 Lerneretal_2013 = TRUE)

# without the Lerner et al. 2013 weighting function
rem.stat.dyadCut(halflife = 30,
                 timeValue = 100,
                 timeWidth = 60,
                 Lerneretal_2013 = FALSE)

# A result to test the function (should come out to 0.50)
rem.stat.dyadCut(halflife = 30,
                 timeValue = 100,
                 timeWidth = 30,
                 Lerneretal_2013 = FALSE)


# Replicating Lerner and Lomi (2020):
#"We set T1/2 to 30 days so that an event counts as (close to) one in the very next instant of time,
#it counts as 1/2 one month later, it counts as 1/4 two months after the event, and so on. To reduce
#the memory consumption needed to store the network of past events, we set a dyadic weight to
#zero if its value drops below 0.01. If a single event occurred in some dyad this would happen after
#6.64×T1/2, that is after more than half a year." (Lerner and Lomi 2020: 104).

# Based upon Lerner and Lomi (2020: 104), the result should be around 0.01. Since the
# time values are in milliseconds, we have to change all measurements into milliseconds
rem.stat.dyadCut(halflife = (30*24*60*60*1000), #30 days in milliseconds
                #the first value in the Lerner and Lomi (2020) WikiEvent 2018 dataset
                timeValue = 979686793000,
                timeWidth = (6.64*30*24*60*60*1000), #Based upon the paper
                #using the Lerner and Lomi (2020) weighting function
                Lerneretal_2013 = FALSE)





cleanEx()
nameEx("rem.stat.indegreeR")
### * rem.stat.indegreeR

flush(stderr()); flush(stdout())

### Name: rem.stat.indegreeR
### Title: Compute Indegree Statistic for Event Receivers in Relational
###   Event Sequences
### Aliases: rem.stat.indegreeR

### ** Examples

events <- data.table::data.table(time = 1:18,
                                eventID = 1:18,
                                sender = c("A", "B", "C",
                                           "A", "D", "E",
                                           "F", "B", "A",
                                           "F", "D", "B",
                                           "G", "B", "D",
                                          "H", "A", "D"),
                               target = c("B", "C", "D",
                                          "E", "A", "F",
                                          "D", "A", "C",
                                          "G", "B", "C",
                                          "H", "J", "A",
                                          "F", "C", "B"))

eventSet <- rem.riskset.om(data = events,
                      time = events$time,
                      eventID = events$eventID,
                      sender = events$sender,
                      receiver = events$target,
                      p_samplingobserved = 1.00,
                      n_controls = 1,
                      seed = 9999)

# Computing Target Indegree Statistics without the sliding windows framework
eventSet$target_indegree <- rem.stat.indegreeR(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   dyadic_weight = 0,
   Lerneretal_2013 = FALSE)

# Computing Target Indegree Statistics with the sliding windows framework
eventSet$target_indegreeSW <- rem.stat.indegreeR(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   processed_seqIDs = eventSet$sequenceID,
   dyadic_weight = 0,
   sliding_window = TRUE,
   Lerneretal_2013 = FALSE)

#The results with and without the sliding windows are the same (see correlation below).
#Using the sliding windows method is recommended when the data are 'big' so
#that memory allotment is more efficient.
cor(eventSet$target_indegree , eventSet$target_indegreeSW )

# Computing Target Indegree Statistics with the counts of events being returned
eventSet$target_indegreeC <- rem.stat.indegreeR(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   processed_seqIDs = eventSet$sequenceID,
   dyadic_weight = 0,
   sliding_window = TRUE,
   Lerneretal_2013 = FALSE,
   counts = TRUE)

cbind(eventSet$target_indegree,
     eventSet$target_indegreeSW,
     eventSet$target_indegreeC)



cleanEx()
nameEx("rem.stat.indegreeS")
### * rem.stat.indegreeS

flush(stderr()); flush(stdout())

### Name: rem.stat.indegreeS
### Title: Compute Indegree Statistic for Event Senders in Relational Event
###   Sequences
### Aliases: rem.stat.indegreeS

### ** Examples

events <- data.table::data.table(time = 1:18,
                                eventID = 1:18,
                                sender = c("A", "B", "C",
                                           "A", "D", "E",
                                           "F", "B", "A",
                                           "F", "D", "B",
                                           "G", "B", "D",
                                          "H", "A", "D"),
                               target = c("B", "C", "D",
                                          "E", "A", "F",
                                          "D", "A", "C",
                                          "G", "B", "C",
                                          "H", "J", "A",
                                          "F", "C", "B"))

eventSet <- rem.riskset.om(data = events,
                      time = events$time,
                      eventID = events$eventID,
                      sender = events$sender,
                      receiver = events$target,
                      p_samplingobserved = 1.00,
                      n_controls = 1,
                      seed = 9999)

# Computing Sender Indegree Statistics without the sliding windows framework
eventSet$sender.indegree <- rem.stat.indegreeS(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   dyadic_weight = 0,
   Lerneretal_2013 = FALSE)

# Computing Sender Indegree Statistics with the sliding windows framework
eventSet$sender.indegree.SW <- rem.stat.indegreeS(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   processed_seqIDs = eventSet$sequenceID,
   dyadic_weight = 0,
   sliding_window = TRUE,
   Lerneretal_2013 = FALSE)

#The results with and without the sliding windows are the same (see correlation below).
#Using the sliding windows method is recommended when the data are 'big' so
#that memory allotment is more efficient.
cor(eventSet$sender.indegree.SW,eventSet$sender.indegree)

# Computing Sender Indegree Statistics with the counts of events being returned
eventSet$sender.indegreeC <- rem.stat.indegreeS(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   dyadic_weight = 0,
   Lerneretal_2013 = FALSE,
   counts = TRUE)

cbind(eventSet$sender.indegree.SW,
     eventSet$sender.indegree,
     eventSet$sender.indegreeC)



cleanEx()
nameEx("rem.stat.outdegreeR")
### * rem.stat.outdegreeR

flush(stderr()); flush(stdout())

### Name: rem.stat.outdegreeR
### Title: Compute Outdegree Statistic for Event Receivers in Relational
###   Event Sequences
### Aliases: rem.stat.outdegreeR

### ** Examples

events <- data.table::data.table(time = 1:18,
                                eventID = 1:18,
                                sender = c("A", "B", "C",
                                           "A", "D", "E",
                                           "F", "B", "A",
                                           "F", "D", "B",
                                           "G", "B", "D",
                                          "H", "A", "D"),
                               target = c("B", "C", "D",
                                          "E", "A", "F",
                                          "D", "A", "C",
                                          "G", "B", "C",
                                          "H", "J", "A",
                                          "F", "C", "B"))

eventSet <- rem.riskset.om(data = events,
                      time = events$time,
                      eventID = events$eventID,
                      sender = events$sender,
                      receiver = events$target,
                      p_samplingobserved = 1.00,
                      n_controls = 1,
                      seed = 9999)

# Computing Target Outdegree Statistics without the sliding windows framework
eventSet$target_outdegree <- rem.stat.outdegreeR(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   dyadic_weight = 0,
   Lerneretal_2013 = FALSE)

# Computing Target Outdegree Statistics with the sliding windows framework
eventSet$target_outdegreeSW <- rem.stat.outdegreeR(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   processed_seqIDs = eventSet$sequenceID,
   dyadic_weight = 0,
   sliding_window = TRUE,
   Lerneretal_2013 = FALSE)

#The results with and without the sliding windows are the same (see correlation below).
#Using the sliding windows method is recommended when the data are 'big' so
#that memory allotment is more efficient.
cor(eventSet$target_outdegreeSW , eventSet$target_outdegree)

# Computing  Target Outdegree Statistic with the counts of events being returned
eventSet$target_outdegreeC <- rem.stat.outdegreeR(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   dyadic_weight = 0,
   sliding_window = FALSE,
   counts = TRUE,
   Lerneretal_2013 = FALSE)

cbind(eventSet$target_outdegree,
     eventSet$target_outdegreeSW,
     eventSet$target_outdegreeC)



cleanEx()
nameEx("rem.stat.outdegreeS")
### * rem.stat.outdegreeS

flush(stderr()); flush(stdout())

### Name: rem.stat.outdegreeS
### Title: Compute Outdegree Statistic for Event Senders in Relational
###   Event Sequences
### Aliases: rem.stat.outdegreeS

### ** Examples

events <- data.table::data.table(time = 1:18,
                                eventID = 1:18,
                                sender = c("A", "B", "C",
                                           "A", "D", "E",
                                           "F", "B", "A",
                                           "F", "D", "B",
                                           "G", "B", "D",
                                          "H", "A", "D"),
                               target = c("B", "C", "D",
                                          "E", "A", "F",
                                          "D", "A", "C",
                                          "G", "B", "C",
                                          "H", "J", "A",
                                          "F", "C", "B"))

eventSet <- rem.riskset.om(data = events,
                      time = events$time,
                      eventID = events$eventID,
                      sender = events$sender,
                      receiver = events$target,
                      p_samplingobserved = 1.00,
                      n_controls = 1,
                      seed = 9999)

# Computing Sender Outdegree Statistics without the sliding windows framework
eventSet$sender_outdegree <- rem.stat.outdegreeS(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   halflife = 2, #halflife parameter
   dyadic_weight = 0,
   Lerneretal_2013 = FALSE)

# Computing Sender Outdegree Statistics with the sliding windows framework
eventSet$sender_outdegreeSW <- rem.stat.outdegreeS(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   halflife = 2, #halflife parameter
   processed_seqIDs = eventSet$sequenceID,
   dyadic_weight = 0,
   sliding_window = TRUE,
   Lerneretal_2013 = FALSE)

#The results with and without the sliding windows are the same (see correlation below).
#Using the sliding windows method is recommended when the data are 'big' so
#that memory allotment is more efficient.
cor(eventSet$sender_outdegreeSW , eventSet$sender_outdegree)

# Computing  Sender Outdegree Statistic with the counts of events being returned
eventSet$sender_outdegreeC <- rem.stat.outdegreeS(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   halflife = 2, #halflife parameter
   dyadic_weight = 0,
   sliding_window = FALSE,
   counts = TRUE,
   Lerneretal_2013 = FALSE)

cbind(eventSet$sender_outdegree,
     eventSet$sender_outdegreeSW,
     eventSet$sender_outdegreeC)



cleanEx()
nameEx("rem.stat.pref.attach")
### * rem.stat.pref.attach

flush(stderr()); flush(stdout())

### Name: rem.stat.pref.attach
### Title: Compute the Preferential Attachment Statistic for Relational
###   Event Sequences
### Aliases: rem.stat.pref.attach

### ** Examples



# A Dummy One-Mode Event Dataset
events <- data.table::data.table(time = 1:18,
                                eventID = 1:18,
                                sender = c("A", "B", "C",
                                           "A", "D", "E",
                                           "F", "B", "A",
                                           "F", "D", "B",
                                           "G", "B", "D",
                                           "H", "A", "D"),
                                target = c("B", "C", "D",
                                           "E", "A", "F",
                                           "D", "A", "C",
                                           "G", "B", "C",
                                           "H", "J", "A",
                                           "F", "C", "B"))

# Creating the Post-Processing Event Dataset with Null Events
eventSet <- rem.riskset.om(data = events,
                          time = events$time,
                          eventID = events$eventID,
                          sender = events$sender,
                          receiver = events$target,
                          p_samplingobserved = 1.00,
                         n_controls = 6,
                         seed = 9999)

# Compute Preferential Attachment Statistic without Sliding Windows Framework and
# No Temporal Dependency
eventSet$pref <- rem.stat.pref.attach(observed_time = events$time,
                                     observed_receiver = events$target,
                                     observed_sender = events$sender,
                                     processed_time = eventSet$time,
                                     processed_receiver = eventSet$receiver,
                                     processed_sender = eventSet$sender,
                                    dependency = FALSE)

# Compute Preferential Attachment Statistic with Sliding Windows Framework and
# No Temporal Dependency
eventSet$prefSW <- rem.stat.pref.attach(observed_time = events$time,
                                       observed_receiver = events$target,
                                       observed_sender = events$sender,
                                       processed_time = eventSet$time,
                                       processed_receiver = eventSet$receiver,
                                       processed_sender = eventSet$sender,
                                       dependency = FALSE,
                                       sliding_windows = TRUE,
                                       processed_seqIDs = eventSet$sequenceID)

#The results with and without the sliding windows are the same (see correlation below).
#Using the sliding windows method is recommended when the data are 'big' so
#that memory allotment is more efficient.
cor(eventSet$pref,eventSet$prefSW) #the correlation of the values


# Compute Preferential Attachment Statistic without Sliding Windows Framework and
# Temporal Dependency
eventSet$prefdep <- rem.stat.pref.attach(observed_time = events$time,
                                        observed_receiver = events$target,
                                        observed_sender = events$sender,
                                        processed_time = eventSet$time,
                                        processed_receiver = eventSet$receiver,
                                        processed_sender = eventSet$sender,
                                        dependency = TRUE,
                                        relationalTimeSpan = 10)

# Compute Preferential Attachment Statistic with Sliding Windows Framework and
# Temporal Dependency
eventSet$pref1dep <- rem.stat.pref.attach(observed_time = events$time,
                                         observed_receiver = events$target,
                                         observed_sender = events$sender,
                                         processed_time = eventSet$time,
                                         processed_receiver = eventSet$receiver,
                                         processed_sender = eventSet$sender,
                                         dependency = TRUE,
                                         relationalTimeSpan = 10,
                                        sliding_windows = TRUE,
                                         processed_seqIDs = eventSet$sequenceID)

#The results with and without the sliding windows are the same (see correlation below).
#Using the sliding windows method is recommended when the data are 'big' so
#that memory allotment is more efficient.
cor(eventSet$prefdep,eventSet$pref1dep) #the correlation of the values




cleanEx()
nameEx("rem.stat.presistence")
### * rem.stat.presistence

flush(stderr()); flush(stdout())

### Name: rem.stat.presistence
### Title: Compute Persistence Statistic for Relational Event Sequences
### Aliases: rem.stat.presistence rem.stat.persistence

### ** Examples



# A Dummy One-Mode Event Dataset
events <- data.table::data.table(time = 1:18,
                                eventID = 1:18,
                                sender = c("A", "B", "C",
                                           "A", "D", "E",
                                           "F", "B", "A",
                                           "F", "D", "B",
                                           "G", "B", "D",
                                           "H", "A", "D"),
                                target = c("B", "C", "D",
                                           "E", "A", "F",
                                           "D", "A", "C",
                                           "G", "B", "C",
                                           "H", "J", "A",
                                           "F", "C", "B"))

# Creating the Post-Processing Event Dataset with Null Events
eventSet <- rem.riskset.om(data = events,
                          time = events$time,
                          eventID = events$eventID,
                          sender = events$sender,
                          receiver = events$target,
                          p_samplingobserved = 1.00,
                          n_controls = 6,
                          seed = 9999)

#Compute Persistence with respect to the sender's past relational history without
#the sliding windows framework and no temporal dependency
eventSet$persist <- rem.stat.persistence(observed_time = events$time,
                                        observed_receiver = events$target,
                                        observed_sender = events$sender,
                                        processed_time = eventSet$time,
                                        processed_receiver = eventSet$receiver,
                                        processed_sender = eventSet$sender,
                                        sender = TRUE,
                                        nopastEvents = 0)

#Compute Persistence with respect to the sender's past relational history with
#the sliding windows framework and no temporal dependency
eventSet$persistSW <- rem.stat.persistence(observed_time = events$time,
                                        observed_receiver = events$target,
                                        observed_sender = events$sender,
                                        processed_time = eventSet$time,
                                        processed_receiver = eventSet$receiver,
                                        processed_sender = eventSet$sender,
                                        sender = TRUE,
                                        sliding_windows = TRUE,
                                        processed_seqIDs = eventSet$sequenceID,
                                        nopastEvents = 0)

#The results with and without the sliding windows are the same (see correlation below).
#Using the sliding windows method is recommended when the data are 'big' so
#that memory allotment is more efficient.
cor(eventSet$persist,eventSet$persistSW)


#Compute Persistence with respect to the sender's past relational history without
#the sliding windows framework and temporal dependency
eventSet$persistDep <- rem.stat.persistence(observed_time = events$time,
                                        observed_receiver = events$target,
                                        observed_sender = events$sender,
                                        processed_time = eventSet$time,
                                        processed_receiver = eventSet$receiver,
                                        processed_sender = eventSet$sender,
                                        sender = TRUE,
                                        dependency = TRUE,
                                        relationalTimeSpan = 5, #the past 5 events
                                        nopastEvents = 0)

#Compute Persistence with respect to the receiver's past relational history without
#the sliding windows framework and no temporal dependency
eventSet$persistT <- rem.stat.persistence(observed_time = events$time,
                                        observed_receiver = events$target,
                                        observed_sender = events$sender,
                                        processed_time = eventSet$time,
                                        processed_receiver = eventSet$receiver,
                                        processed_sender = eventSet$sender,
                                        sender = FALSE,
                                        nopastEvents = 0)

#Compute Persistence with respect to the receiver's past relational history with
#the sliding windows framework and no temporal dependency
eventSet$persistSWT <- rem.stat.persistence(observed_time = events$time,
                                        observed_receiver = events$target,
                                        observed_sender = events$sender,
                                        processed_time = eventSet$time,
                                        processed_receiver = eventSet$receiver,
                                        processed_sender = eventSet$sender,
                                        sender = FALSE,
                                        sliding_windows = TRUE,
                                        processed_seqIDs = eventSet$sequenceID,
                                        nopastEvents = 0)

#The results with and without the sliding windows are the same (see correlation below).
#Using the sliding windows method is recommended when the data are 'big' so
#that memory allotment is more efficient.
cor(eventSet$persistT,eventSet$persistSWT)


#Compute Persistence with respect to the receiver's past relational history without
#the sliding windows framework and temporal dependency
eventSet$persistDepT <- rem.stat.persistence(observed_time = events$time,
                                        observed_receiver = events$target,
                                        observed_sender = events$sender,
                                        processed_time = eventSet$time,
                                        processed_receiver = eventSet$receiver,
                                        processed_sender = eventSet$sender,
                                        sender = FALSE,
                                        dependency = TRUE,
                                        relationalTimeSpan = 5, #the past 5 events
                                        nopastEvents = 0)




cleanEx()
nameEx("rem.stat.recency")
### * rem.stat.recency

flush(stderr()); flush(stdout())

### Name: rem.stat.recency
### Title: Compute the Recency Statistic for Relational Event Sequences
### Aliases: rem.stat.recency

### ** Examples



# A Dummy One-Mode Event Dataset
events <- data.table::data.table(time = 1:18,
                                eventID = 1:18,
                                sender = c("A", "B", "C",
                                           "A", "D", "E",
                                           "F", "B", "A",
                                           "F", "D", "B",
                                           "G", "B", "D",
                                           "H", "A", "D"),
                                target = c("B", "C", "D",
                                           "E", "A", "F",
                                           "D", "A", "C",
                                           "G", "B", "C",
                                           "H", "J", "A",
                                           "F", "C", "B"))

# Creating the Post-Processing Event Dataset with Null Events
eventSet <- rem.riskset.om(data = events,
                          time = events$time,
                          eventID = events$eventID,
                          sender = events$sender,
                          receiver = events$target,
                          p_samplingobserved = 1.00,
                         n_controls = 6,
                         seed = 9999)

# Compute Recency Statistic without Sliding Windows Framework and
# No Temporal Dependency
eventSet$recency_rawdiff <- rem.stat.recency(
 observed_time = events$time, # variable (column) name that contains the time variable
 observed_receiver = events$target, # variable (column) name that contains the receiver variable
 observed_sender = events$sender,
 processed_time = eventSet$time,
 processed_receiver = eventSet$receiver,
 processed_sender = eventSet$sender,
 type = "raw.diff",
 dependency = FALSE,
 i_neighborhood = TRUE,
 nopastEvents = 0)

# Compute Recency Statistic without Sliding Windows Framework and
# No Temporal Dependency
eventSet$recency_inv <- rem.stat.recency(
 observed_time = events$time, # variable (column) name that contains the time variable
 observed_receiver = events$target, # variable (column) name that contains the receiver variable
 observed_sender = events$sender,
 processed_time = eventSet$time,
 processed_receiver = eventSet$receiver,
 processed_sender = eventSet$sender,
 type = "inv.diff.plus1",
 dependency = FALSE,
 i_neighborhood = TRUE,
 nopastEvents = 0)


# Compute Recency Statistic without Sliding Windows Framework and
# No Temporal Dependency
eventSet$recency_rank <- rem.stat.recency(
 observed_time = events$time,
 observed_receiver = events$target,
 observed_sender = events$sender,
 processed_time = eventSet$time,
 processed_receiver = eventSet$receiver,
 processed_sender = eventSet$sender,
 type = "rank.ordered.count",
 dependency = FALSE,
 i_neighborhood = TRUE,
 nopastEvents = 0)

# Compute Recency Statistic with Sliding Windows Framework and No Temporal Dependency
eventSet$recency_rawdiffSW <- rem.stat.recency(
 observed_time = events$time,
 observed_receiver = events$target,
 observed_sender = events$sender,
 processed_time = eventSet$time,
 processed_receiver = eventSet$receiver,
 processed_sender = eventSet$sender,
 type = "raw.diff",
 dependency = FALSE,
 i_neighborhood = TRUE,
 sliding_windows = TRUE,
 processed_seqIDs = eventSet$sequenceID,
 nopastEvents = 0)


# Compute Recency Statistic with Sliding Windows Framework and No Temporal Dependency
eventSet$recency_invSW <- rem.stat.recency(
 observed_time = events$time,
 observed_receiver = events$target,
 observed_sender = events$sender,
 processed_time = eventSet$time,
 processed_receiver = eventSet$receiver,
 processed_sender = eventSet$sender,
 type = "inv.diff.plus1",
 dependency = FALSE,
 i_neighborhood = TRUE,
 sliding_windows = TRUE,
 processed_seqIDs = eventSet$sequenceID,
 nopastEvents = 0)


# Compute Recency Statistic with Sliding Windows Framework and No Temporal Dependency
eventSet$recency_rankSW <- rem.stat.recency(
 observed_time = events$time,
 observed_receiver = events$target,
 observed_sender = events$sender,
 processed_time = eventSet$time,
 processed_receiver = eventSet$receiver,
 processed_sender = eventSet$sender,
 type = "rank.ordered.count",
 dependency = FALSE,
 i_neighborhood = TRUE,
 sliding_windows = TRUE,
 processed_seqIDs = eventSet$sequenceID,
 nopastEvents = 0)




cleanEx()
nameEx("rem.stat.recip")
### * rem.stat.recip

flush(stderr()); flush(stdout())

### Name: rem.stat.recip
### Title: Compute Reciprocity Statistic for Relational Event Sequences
### Aliases: rem.stat.recip

### ** Examples

events <- data.table::data.table(time = 1:18,
                                eventID = 1:18,
                                sender = c("A", "B", "C",
                                           "A", "D", "E",
                                           "F", "B", "A",
                                           "F", "D", "B",
                                           "G", "B", "D",
                                          "H", "A", "D"),
                               target = c("B", "C", "D",
                                          "E", "A", "F",
                                          "D", "A", "C",
                                          "G", "B", "C",
                                          "H", "J", "A",
                                          "F", "C", "B"))

eventSet <- rem.riskset.om(data = events,
                      time = events$time,
                      eventID = events$eventID,
                      sender = events$sender,
                      receiver = events$target,
                      p_samplingobserved = 1.00,
                      n_controls = 1,
                      seed = 9999)

# Computing Reciprocity Statistics without the sliding windows framework
eventSet$recip <- rem.stat.recip(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   dyadic_weight = 0,
   Lerneretal_2013 = FALSE)

# Computing Reciprocity Statistics with the sliding windows framework
eventSet$recipSW <- rem.stat.recip(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   processed_seqIDs = eventSet$sequenceID,
   dyadic_weight = 0,
   sliding_window = TRUE,
   Lerneretal_2013 = FALSE)

#The results with and without the sliding windows are the same (see correlation below).
#Using the sliding windows method is recommended when the data are 'big' so
#that memory allotment is more efficient.
cor(eventSet$recipSW , eventSet$recip)

# Computing Reciprocity Statistics with the counts of events being returned
eventSet$recipC <- rem.stat.recip(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   dyadic_weight = 0,
   sliding_window = FALSE,
   counts = TRUE,
   Lerneretal_2013 = FALSE)

cbind(eventSet$recip,
     eventSet$recipSW,
     eventSet$recipC)



cleanEx()
nameEx("rem.stat.repetition")
### * rem.stat.repetition

flush(stderr()); flush(stdout())

### Name: rem.stat.repetition
### Title: Compute Repetition Statistic for Relational Event Sequences
### Aliases: rem.stat.repetition

### ** Examples

data("WikiEvent2018.first100k")
WikiEvent2018 <- WikiEvent2018.first100k[1:10000,] #the first ten thousand events
WikiEvent2018$time <- as.numeric(WikiEvent2018$time) #making the variable numeric
### Creating the EventSet By Employing Case-Control Sampling With M = 5 and
### Sampling from the Observed Event Sequence with P = 0.01
EventSet <- rem.riskset.tm(
 data = WikiEvent2018, # The Event Dataset
 time = WikiEvent2018$time, # The Time Variable
 eventID = WikiEvent2018$eventID, # The Event Sequence Variable
 sender = WikiEvent2018$user, # The Sender Variable
 receiver = WikiEvent2018$article, # The Receiver Variable
 p_samplingobserved = 0.01, # The Probability of Selection
 n_controls = 5, # The Number of Controls to Sample from the Full Risk Set
 seed = 9999) # The Seed for Replication
#### Estimating Repetition Scores Without the Sliding Windows Framework
EventSet$rep <- rem.stat.repetition(
   observed_time = WikiEvent2018$time,
   observed_sender = WikiEvent2018$user,
   observed_receiver = WikiEvent2018$article,
   processed_time = EventSet$time,
   processed_sender = EventSet$sender,
   processed_receiver = EventSet$receiver,
   halflife = 2.592e+09, #halflife parameter
   dyadic_weight = 0,
   Lerneretal_2013 = FALSE)

EventSet$sw_rep <- rem.stat.repetition(
   observed_time = WikiEvent2018$time,
   observed_sender = WikiEvent2018$user,
   observed_receiver = WikiEvent2018$article,
   processed_time = EventSet$time,
   processed_sender = EventSet$sender,
   processed_receiver = EventSet$receiver,
   processed_seqIDs = EventSet$sequenceID,
   halflife = 2.592e+09, #halflife parameter
   dyadic_weight = 0,
   sliding_window = TRUE,
   Lerneretal_2013 = FALSE)

#The results with and without the sliding windows are the same (see correlation below).
#Using the sliding windows method is recommended when the data are 'big' so
#that memory allotment is more efficient.
cor(EventSet$sw_rep, EventSet$rep)

#### Estimating Repetition Scores with the Counts of Events Returned
EventSet$repC <- rem.stat.repetition(
   observed_time = WikiEvent2018$time,
   observed_sender = WikiEvent2018$user,
   observed_receiver = WikiEvent2018$article,
   processed_time = EventSet$time,
   processed_sender = EventSet$sender,
   processed_receiver = EventSet$receiver,
   halflife = 2.592e+09, #halflife parameter
   dyadic_weight = 0,
   Lerneretal_2013 = FALSE,
   counts = TRUE)

cbind(EventSet$rep,
     EventSet$sw_rep,
     EventSet$repC)





cleanEx()
nameEx("rem.stat.triad")
### * rem.stat.triad

flush(stderr()); flush(stdout())

### Name: rem.stat.triad
### Title: Compute Triad Closure Statistics for Relational Event Sequences
### Aliases: rem.stat.triad

### ** Examples

events <- data.table::data.table(time = 1:18,
                                eventID = 1:18,
                                sender = c("A", "B", "C",
                                           "A", "D", "E",
                                           "F", "B", "A",
                                           "F", "D", "B",
                                           "G", "B", "D",
                                          "H", "A", "D"),
                               target = c("B", "C", "D",
                                          "E", "A", "F",
                                          "D", "A", "C",
                                          "G", "B", "C",
                                          "H", "J", "A",
                                          "F", "C", "B"))

eventSet <- rem.riskset.om(data = events,
                      time = events$time,
                      eventID = events$eventID,
                      sender = events$sender,
                      receiver = events$target,
                      p_samplingobserved = 1.00,
                      n_controls = 1,
                      seed = 9999)

# Computing Triadic Statistics without the sliding windows framework
eventSet$triadic <- rem.stat.triad(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   dyadic_weight = 0,
   Lerneretal_2013 = FALSE)

# Computing Triadic Statistics with the sliding windows framework
eventSet$triadicSW <- rem.stat.triad(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   processed_seqIDs = eventSet$sequenceID,
   dyadic_weight = 0,
   sliding_window = TRUE,
   Lerneretal_2013 = FALSE)

#The results with and without the sliding windows are the same (see correlation below).
#Using the sliding windows method is recommended when the data are 'big' so
#that memory allotment is more efficient.
cor(eventSet$triadic , eventSet$triadicSW)

# Computing Triadic Statistics with the counts of events being returned
eventSet$triadicC <- rem.stat.triad(
   observed_time = events$time, # variable (column) name that contains the time variable
   observed_sender = events$sender, # variable (column) name that contains the sender variable
   observed_receiver = events$target, # variable (column) name that contains the receiver variable
   processed_time = eventSet$time,
   processed_sender = eventSet$sender,
   processed_receiver = eventSet$receiver,
   halflife = 2, #halflife parameter
   dyadic_weight = 0,
   sliding_window = FALSE,
   counts = TRUE,
   Lerneretal_2013 = FALSE)

cbind(eventSet$triadic,
     eventSet$triadicSW,
     eventSet$triadicC)



cleanEx()
nameEx("tm.constraint")
### * tm.constraint

flush(stderr()); flush(stdout())

### Name: tm.constraint
### Title: Compute Burchard and Cornwell (2018) Two-Mode Constraint
### Aliases: tm.constraint

### ** Examples


# For this example, we recreate Figure 2 in Burchard and Cornwell (2018: 13)
BCNet <- matrix(
 c(1,1,0,0,
   1,0,1,0,
   1,0,0,1,
   0,1,1,1),
 nrow = 4, ncol = 4, byrow = TRUE)
colnames(BCNet) <- c("1", "2", "3", "4")
rownames(BCNet) <- c("i", "j", "k", "m")
#library(sna) #To plot the two mode network, we use the sna R package
#gplot(BCNet, usearrows = FALSE,
#      gmode = "twomode", displaylabels = TRUE)
tm.constraint(BCNet)

#For this example, we recreate Figure 9 in Burchard and Cornwell (2018:18) for
#weighted two mode networks.
BCweighted <- matrix(c(1,2,1, 1,0,0,
                      0,2,1,0,0,1),
                    nrow = 4, ncol = 3,
                    byrow = TRUE)
rownames(BCweighted) <- c("i", "j", "k", "l")
tm.constraint(BCweighted, weighted = TRUE)







cleanEx()
nameEx("tm.degree")
### * tm.degree

flush(stderr()); flush(stdout())

### Name: tm.degree
### Title: Compute Degree Centrality Values for Two-Mode Networks
### Aliases: tm.degree

### ** Examples

#Replicating the biparitate graph presented in Knoke and Yang (2020: 109)
knoke_yang_PC <- matrix(c(1,1,0,0, 1,1,0,0,
                          1,1,1,0, 0,0,1,1,
                          0,0,1,1), byrow = TRUE,
                          nrow = 5, ncol = 4)
colnames(knoke_yang_PC) <- c("Rubio-R","McConnell-R", "Reid-D", "Sanders-D")
rownames(knoke_yang_PC) <- c("UPS", "MS", "HD", "SEU", "ANA")
tm.degree(knoke_yang_PC, level1 = TRUE) #note: this value matches that of the book
tm.degree(knoke_yang_PC, level1 = FALSE) #note: this value matches that of the book



cleanEx()
nameEx("tm.density")
### * tm.density

flush(stderr()); flush(stdout())

### Name: tm.density
### Title: Compute Density for Two-Mode Networks
### Aliases: tm.density

### ** Examples

#Replicating the biparitate graph presented in Knoke and Yang (2020: 109)
knoke_yang_PC <- matrix(c(1,1,0,0, 1,1,0,0,
                          1,1,1,0, 0,0,1,1,
                          0,0,1,1), byrow = TRUE,
                          nrow = 5, ncol = 4)
colnames(knoke_yang_PC) <- c("Rubio-R","McConnell-R", "Reid-D", "Sanders-D")
rownames(knoke_yang_PC) <- c("UPS", "MS", "HD", "SEU", "ANA")
tm.density(knoke_yang_PC, level1 = TRUE) #note: this value does not match that of the book,
                                  #but does match that of Wasserman and Faust (1995: 317)
                                  #for the ceo dataset.
tm.density(knoke_yang_PC, level1 = FALSE) #note: this value matches that of the book




cleanEx()
nameEx("tm.effective")
### * tm.effective

flush(stderr()); flush(stdout())

### Name: tm.effective
### Title: Compute Burchard and Cornwell (2018) Two-Mode Effective Size
### Aliases: tm.effective

### ** Examples


# For this example, we recreate Figure 2 in Burchard and Cornwell (2018: 13)
BCNet <- matrix(
 c(1,1,0,0,
   1,0,1,0,
   1,0,0,1,
   0,1,1,1),
 nrow = 4, ncol = 4, byrow = TRUE)
colnames(BCNet) <- c("1", "2", "3", "4")
rownames(BCNet) <- c("i", "j", "k", "m")
#library(sna) #To plot the two mode network, we use the sna R package
#gplot(BCNet, usearrows = FALSE,
#      gmode = "twomode", displaylabels = TRUE)
tm.effective(BCNet)

#For this example, we recreate Figure 9 in Burchard and Cornwell (2018:18) f
#or weighted two mode networks.
BCweighted <- matrix(c(1,2,1, 1,0,0,
                      0,2,1,0,0,1),
                      nrow = 4, ncol = 3,
                      byrow = TRUE)
rownames(BCweighted) <- c("i", "j", "k", "l")
tm.effective(BCweighted, weighted = TRUE)




cleanEx()
nameEx("tm.homo4cycles")
### * tm.homo4cycles

flush(stderr()); flush(stdout())

### Name: tm.homo4cycles
### Title: Compute Homophlious Four Cycles in Two-Mode Networks
### Aliases: tm.homo4cycles

### ** Examples


# For this example, we use the Davis Southern Women's Dataset.
data("southern.women")
#creating a random binary membership vector
set.seed(9999)
membership <- sample(0:1, nrow(southern.women), replace = TRUE)
#the homophilous four-cycle values
tm.homo4cycles(southern.women, mem = membership)



cleanEx()
nameEx("tm.homoDis")
### * tm.homoDis

flush(stderr()); flush(stdout())

### Name: tm.homoDis
### Title: Compute Ego Homophily Distance in Two-Mode Networks
### Aliases: tm.homoDis

### ** Examples


# For this example, we use the Davis Southern Women's Dataset.
data("southern.women")
#creating a random binary membership vector
set.seed(9999)
membership <- sample(0:1, nrow(southern.women), replace = TRUE)
#the ego 2 mode distance non-standardized
tm.homoDis(southern.women, mem = membership)
#the ego 2 mode distance standardized
tm.homoDis(southern.women, mem = membership, standardize = TRUE)




cleanEx()
nameEx("tm.redundancy")
### * tm.redundancy

flush(stderr()); flush(stdout())

### Name: tm.redundancy
### Title: Compute Burchard and Cornwell (2018) Two-Mode Redundancy
### Aliases: tm.redundancy

### ** Examples


# For this example, we recreate Figure 2 in Burchard and Cornwell (2018: 13)
BCNet <- matrix(
 c(1,1,0,0,
   1,0,1,0,
   1,0,0,1,
   0,1,1,1),
 nrow = 4, ncol = 4, byrow = TRUE)
colnames(BCNet) <- c("1", "2", "3", "4")
rownames(BCNet) <- c("i", "j", "k", "m")
#library(sna) #To plot the two mode network, we use the sna R package
#gplot(BCNet, usearrows = FALSE,
#      gmode = "twomode", displaylabels = TRUE)
tm.redundancy(BCNet) #this values replicate those reported by Burchard and Cornwell (2018: 14)


#For this example, we recreate Figure 9 in Burchard and Cornwell (2018:18)
#for weighted two mode networks.
BCweighted <- matrix(c(1,2,1, 1,0,0,
                      0,2,1,0,0,1),
                      nrow = 4, ncol = 3,
                      byrow = TRUE)
rownames(BCweighted) <- c("i", "j", "k", "l")
tm.redundancy(BCweighted, weighted = TRUE)





### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
