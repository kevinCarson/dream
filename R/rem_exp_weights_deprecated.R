## Helper functions for computing statistics for relational event models
## This script includes two functions: minimum effective time and
## the exponential weighting function
## Written by KC (06-13-2025)

# --------------------------------------------------------------------------------------------------
#' @title Helper Function to Compute Minimum Effective Time and Exponential Weights for REM Statistics
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `remExpWeights()` has been deprecated starting on version 1.0.0 of the `dream` package. Please see the `NEWS.md` file for more details.
#'
#' A helper function for computing exponential decay weights and the corresponding minimum effective time used
#' to calculate network statistics in relational event models within the \pkg{dream} package.
#' This implementation follows the formulations of Lerner et al. (2013) and Lerner & Lomi (2020).
#' Although primarily designed for internal use (e.g., within \code{\link{computeReciprocity}}),
#' it may also be of interest to users working directly with REM statistics (e.g., creating new statistics).
#' @details
#' - **Exponential Weighting Function**:
#'   - Lerner & Lomi (2020): \eqn{w(u,a,t) = \sum \exp(- (t - t') * (\log(2)/T_{1/2}))}
#'   - Lerner et al. (2013): \eqn{w(u,a,t) = \sum \exp(- (t - t') * (\log(2)/T_{1/2})) * (\log(2)/T_{1/2})}
#'
#' - **Minimum Effective Time (MEF)**:
#'   - Lerner & Lomi (2020): \eqn{MEF = t + \log(w) / (\log(2)/T_{1/2})}
#'   - Lerner et al. (2013): \eqn{MEF = t + [T_{1/2} * \log((w * T_{1/2}) / \log(2))] / \log(2)}
#'
#' @param current The current relational event time.
#' @param past The numeric vector of past event times (for exponential weighting only).
#' @param halflife The halflife parameter for exponential weighting.
#' @param dyadic_weight The dyadic (event) weight cutoff for relational relevancy.
#' @param Lerneretal_2013 TRUE/FALSE. If TRUE, the function uses the Lerner et al. (2013) exponential weighting function. If FALSE, the function uses the Lerner and Lomi (2020) exponential weighting function.
#' @param exp.weights TRUE/FALSE. If TRUE, the function computes the exponential weights for past relational events. If FALSE, the function computes the minimum effective time for a relational event (that is, the minimum past time that would result in a 0 value for an exponential weight).
#' @return When exp.weights = TRUE, the numeric vector of exponential decay weights. When exp.weights = FALSE, the scalar for the minimum event cut-off time.
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#' @references
#' Lerner, Jürgen and Alessandro Lomi. 2020. “Reliability of relational event
#' model estimates under sampling: How to fit a relational event model to 360
#' million dyadic events.” *Network Science* 8(1): 97-135.
#'
#' Lerner, Jürgen, Margit Bussman, Tom A.B. Snijders, and Ulrik Brandes. 2013. "
#' Modeling Frequency and Type of Interaction in Event Networks."
#' *The Corvinus Journal of Sociology and Social Policy* 4(1): 3-32.
# --------------------------------------------------------------------------------------------------

remExpWeights <- function(current, #if finding weights: observed event time; if minimum: current time
                        past = NULL, #argument only based for finding exponential weights
                        halflife, #halflife from user
                        dyadic_weight,  #dyadic weight from user
                        Lerneretal_2013 = FALSE, #weighting function
                        exp.weights = TRUE) { #should weights be found (if false: minimum time is found)


  if(exp.weights == TRUE){#finding exponential weights for relational event statistics

    if (Lerneretal_2013 == FALSE) { # If the user wants to use the lerner and lomi 2020 weight
      #
      #                           [ (-(t - t')) * (ln(2) / T1/2)        ]
      #  w(u,a,t) = ∑allevents exp
      #    where, t is the current event time, t' is the past event time, and T1/2 represents the given halflife
      #    ln(2) is a constant in the calculation
      w_uat1 <- exp((-(current - past) * log(2) / (halflife)))
      w_uat1[  w_uat1 < dyadic_weight] <- 0 # if the value is less than the dyadic weight
      # value becomes 0
    }

    if (Lerneretal_2013 == TRUE) { # If the user wants to use the lerner et. al 2013 weight
      #
      #                           [ (-(t - t')) * (ln(2) / T1/2)        ]
      #  w(u,a,t) = ∑allevents exp                                        * (ln(2) / T1/2)
      #    where, t is the current event time, t' is the past event time, and T1/2 represents the given halflife
      #    ln(2) is a constant in the calculation
      w_uat1 <- exp((-(current - past) * log(2) / (halflife)))  * log(2) / (halflife) # exponentiate the weight and multiply by (ln(2) / T1/2)
      w_uat1[  w_uat1 < dyadic_weight] <- 0 # if the value is less than the dyadic weight
      # value becomes 0
    }
    return(w_uat1) # return the calculated weight

  }else{ #finding minimum time for relational relevancy for relational event statistics

    if (Lerneretal_2013 == FALSE) { # if the user does not want to use the Lerner et al. 2013 weight
      #  Lerner et al. 2020: w(s,r,t) = ∑allevents exp[ (-(t - t')) * (ln(2) / T1/2)        ]
      #  mef = currenttime + ln(dyadicweight) / (ln(2)/T1/2)
      minimumtime <- current + (log(dyadic_weight) / (log(2) / halflife))
      minimumtime <- ifelse(minimumtime < 0, 0, minimumtime) # if the minimum time is negative, set cut off to 0, essentially, all past time is relational relevant
    }

    if (Lerneretal_2013 == TRUE) { # if the user does  want to use the Lerner et al. 2013 weight
      #  Lerner et al. 2013: w(s,r,t) = ∑allevents exp[ (-(t - t')) * (ln(2) / T1/2) ] * (ln(2) / T1/2)
      #  mef = currenttime + [ T1/2 * ln((dyadicweight*T1/2) / ln(2))] / ln(2)
      minimumtime <- current + (halflife * log((dyadic_weight * halflife) / log(2))) / log(2)
      minimumtime <- ifelse(minimumtime < 0, 0, minimumtime) # if the minimum time is negative, set cut off to 0, essentially, all past time is relational relevant
    }

    return(round(minimumtime)) #return the effective time
  }

}




