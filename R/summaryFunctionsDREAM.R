# Code Written By Kevin Carson (06-27-25)

#' Summary Method for dream Objects
#'
#' Summarizes the results of an ordinal timing relational event model.
#'
#' @param object An object of class "dream".
#' @param digits The number of digits to print after the decimal point.
#' @param ... Additional arguments (currently unused).
#' @return A list of summary statistics for the relational event model
#' including parameter estimates, (null) likelihoods, and tests of significance
#' for likelihood ratios and estimated parameters.
#' @export
#'
summary.dream <- function(object,digits=4,...) {
  coef_table <- cbind(
    `Estimate` = round(object$parameters,digits),
    `Std. Error` = round(object$se.parameter,digits),
    `z value` = round(object$z.values,digits),
    `Pr(>|z|)` = round(object$p.values,digits)
  )
  res <- list(
    call = object$call,
    coefficients = coef_table,
    logLik = round(object$loglikelihood.full,digits),
    nullLogLik = round(object$loglikelihood.null,digits),
    chiStat = round(object$chi.stat,digits),
    df= object$df.full,
    plrtest = round(object$loglikelihood.test,digits),
    nevents = object$n.events,digits,
    nullevents = object$null.events,
    AIC = round(object$AIC,digits),
    BIC = round(object$BIC,digits),
    optimization.method = object$optimization.method,
    iterations = object$newton.iterations
  )
  class(res) <- "summary.dream"
  return(res)
}

#' Print Method for dream Model
#'
#' @param x An object of class "dream".
#' @param digits The number of digits to print after the decimal point.
#' @param ... Additional arguments (currently unused).
#' @return No return value. Prints out the main results of a 'dream' summary object.
#' @export
print.summary.dream <- function(x,digits=4,...) {
  cat("Ordinal Timing Relational Event Model\n\n")

  cat("Call:\n")
  print(x$call)
  cat("\n n events:",x$nevents, "null events:", x$nullevents,"\n")

  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  cat("\nNull Likelihood:",round(x$nullLogLik,digits),"Model Likelihood:",
      x$logLik,"\n")
  cat("\nLikelihood Ratio Test:",round(x$chiStat,digits)," with df:",x$df,"p-value:",round(x$plrtest,digits) , "\n")
  cat("\nAIC",round(x$AIC,digits),"BIC",round(x$BIC,digits),"\n")
  if(x$optimization.method!="optim"){
    cat("\nNumber of Newton Iterations:",x$iterations,"\n")
  }
}

#' Print Method for Summary of dream Model
#'
#' @param x An object of class "summary.dream".
#' @param digits The number of digits to print after the decimal point.
#' @param ... Additional arguments (currently unused).
#' @return No return value. Prints out the main results of a 'dream' object.
#' @export

print.dream <- function(x,digits=4,...) {
  cat("Ordinal Timing Relational Event Model\n\n")

  cat("Call:\n")
  print(x$call)
  cat("\n n events:",x$n.events, "null events:",x$null.events,"\n")

  cat("\nParameter Estimates:\n")
  coefs <- round(x$parameters,digits)
  print(coefs)

  cat("\nNull Likelihood:",round(x$loglikelihood.null,digits),"Model Likelihood:",
      round(x$loglikelihood.full,digits),"\n")
  cat("\nLikelihood Ratio Test:",round(x$chi.stat,digits)," with df:",x$df.full,"p-value:",round(x$loglikelihood.test,digits) , "\n")
  cat("\nAIC",round(x$AIC,digits),"BIC",round(x$BIC,digits),"\n")
  if(x$optimization.method!="optim"){
    cat("\nNumber of Newton Iterations:",x$iterations,"\n")
  }
  invisible(x)
}
