## The Computation of a Dyadic Weight Cutoff
## See Butts 2008; Lerner and Lomi 2020; REM R Package (Brandenberger 2018)
## Code written by Kevin Carson (kacarson@arizona.edu) and Deigo Leal (https://www.diegoleal.info/)
## Last Updated: 07-15-24
#' @title A Helper Function to Assist Researchers in Finding Dyadic Weight Cutoff Values
#' @name remstats_dyadcut
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
#' remstats_dyadcut(halflife = 15,
#'                  relationalWidth = 30,
#'                  exp_weight_form = TRUE)
#'
#' # without the Lerner et al. 2013 weighting function
#' remstats_dyadcut(halflife = 15,
#'                  relationalWidth = 30,
#'                  exp_weight_form = FALSE)
#'
#'# A result to test the function (should come out to 0.50)
#' remstats_dyadcut(halflife = 30,
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
#'remstats_dyadcut(halflife = (30*24*60*60*1000), #30 days in milliseconds
#'                 relationalWidth = (6.64*30*24*60*60*1000), #Based upon the paper
#'                 #using the Lerner and Lomi (2020) weighting function
#'                 exp_weight_form = FALSE)
#'
#'
remstats_dyadcut <- function(halflife = 2,  #the user specificed halflife
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


