# Code Written By Kevin Carson (06-27-25)

#' @title Fit a Relational Event Model (REM) to Event Sequence Data
#' @name estimateREM
#' @param formula A formula object with the dependent variable on the left hand side of
#' ~ and the covariates on the right hand side. This is the same argument found in \code{\link[stats]{lm}} and \code{\link[stats]{glm}}.
#' @param event.cluster An integer or factor vector that groups each observed event with its corresponding control (null) events.
#' This vector defines the strata in the event sequence, ensuring that each stratum contains one observed event and its associated
#' null alternatives. It is used to structure the likelihood by stratifying events based on their occurrence in time.
#' @param ordinal TRUE/FALSE. Currently, this function supports only the estimation of ordinal timing relational event models
#' (see Butts 2008). Future versions of the package will include estimation options for interval timing
#' relational event models. At this time, this argument is preset to TRUE and should not be modified by the user.
#' @param multiple.events TRUE/FALSE. Currently, this function assumes that only one event occurs per
#' event cluster (i.e., time point). Future versions of the package will include estimation options for
#' multiple events per time point, commonly referred to as tied events, via the Breslow approximation
#' technique (see Box-Steffensmeier and Jones 2004). At this moment, this argument is preset to FALSE
#' and should not be modified by the user.
#' @param data The data.frame that contains the variable included in the formula argument.
#' @param newton.rhapson TRUE/FALSE. TRUE indicates an internal Newton-Rhapson iteration procedure with line searching is used to
#' find the set of maximum likelihood estimates. FALSE indicates that the log likelihood function will be optimized via the
#' \code{\link{optim}} function. The function defaults to TRUE.
#' @param optim.method If newton.rhapson is FALSE, what optim method should be used in conjunction with the \code{\link{optim}} function. Defaults
#' to "BFGS". See the \code{\link{optim}} function for the set of options.
#' @param optim.control If newton.rhapson is FALSE, a list of control to be used in the \code{\link{optim}} function. See the \code{\link{optim}} function for the set of controls.
#' @param tolerance If newton.rhapson is TRUE, the stopping criterion for the absolute difference in the log likelihoods for each Newton-Rhapson iteration.
#' The optimization procedure stops when the absolute change in the log likelihoods is less than `tolerance` (see Greene 2003).
#' @param maxit If newton.rhapson is TRUE, the maximum number of iterations for the Newton-Rhapson optimization procedure (see Greene 2003).
#' @param starting.beta A numeric vector that represents the starting parameter estimates for the Newton-Rhapson optimization procedure. This may be a beneficial argument
#' if the optimization procedure fails, since the Newton-Rhapson optimization procedure is sensitive to starting values. Preset to NULL.
#' @param ... Additional arguments.
#' @import stats
#' @import Rcpp
#' @return An object of class "dream" as a list containing the following components:
#' \describe{
#'   \item{optimization.method}{The optimzation method used to find the parameters..}
#'   \item{converged}{TRUE/FALSE. TRUE indicates that the REM converged.}
#'   \item{loglikelihood.null}{The log likelihood of the null model (i.e., the model where the parameters are assumed to be 0).}
#'   \item{loglikelihood.full}{The log likelihood of the estimated model.}
#'   \item{chi.stat}{The chi-statistic of the likelihood ratio test.}
#'   \item{loglikelihood.test}{The p-value of the likelihood ratio test.}
#'   \item{df.null}{The degrees of freedom of the null model.}
#'   \item{df.full}{The degrees of freedom of the full model.}
#'   \item{parameters}{The MLE parameter estimates.}
#'   \item{hessian}{The estimated hessian matrix.}
#'   \item{gradient}{The estimated gradient vector.}
#'   \item{se.parameter}{The standard errors of the MLE parameter estimates.}
#'   \item{covariance.mat}{The estimated variance-covariance matrix.}
#'   \item{z.values}{The z-scores for the MLE parameter estimates.}
#'   \item{p.values}{The p-values for the MLE parameter estimates.}
#'   \item{AIC}{The AIC of the estimated REM.}
#'   \item{BIC}{The BIC of the estimated REM.}
#'   \item{n.events}{The number of observed events in the relational event sequence.}
#'   \item{null.events}{The number of control events in the relational event sequence.}
#'   \item{newton.iterations}{The number of Newton-Rhapson iterations.}
#'   \item{search.algo}{A data.frame object that contains the Newton-Rhapson searching algorithm results.}
#'}
#' @export


#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `estimateREM()` has been deprecated starting on version 1.0.0 of the `dream` package. Please use the `estimate_rem_logit()` function and see the `NEWS.md` file for more details.
#'
#'
#'
#' This function estimates the ordinal timing relational event model by maximizing the
#' likelihood function given by Butts (2008) via maximum likelihood estimation. A nice outcome
#' is that the ordinal timing relational event model is equivalent to the conditional logistic
#' regression (see Greene 2003; for R functions, see \code{\link[survival]{clogit}}). In
#' addition, based on this outcome and the structure of the data, this function can estimate
#' the Cox proportional hazards model (see Box-Steffensmeier and Jones 2004; for R functions, see \code{\link[survival]{coxph}})
#' given that the likelihood functions are equivalent. An important assumption this model
#' makes is that only one event occurs at each time point. If this is unfeasible for
#' the user's specific dataset, we encourage the user to see the \code{\link[survival]{clogit}}
#' function for the Breslow approximation technique (Box-Steffensmeier and Jones 2004). Future
#' versions of the package will include options for interval timing relational event
#' models and tied event data (e.g., multiple events at one time point).
#'
#' @details
#' This function maximizes the ordinal timing relational event model likelihood
#' function provided in the seminal REM paper by Butts (2008). The likelihood function
#' is:
#' \deqn{L(E|\beta) = \prod_{i=1}^{|E|} \frac{\lambda_{e_i}}{\sum_{e' \in RS_{e_i}} \lambda_{e'}}}
#' where, following Butts (2008) and Duxbury (2020), \eqn{E} is the relational event sequence,
#' \eqn{\lambda_{e_i}} is the hazard rate for event *i*, which is formulated to be equal to
#' \eqn{exp(\beta^{T}z(x,Y))}, that is, the linear combination of user-specific covariates, \eqn{z(x,Y)}, and associated
#' REM parameters, \eqn{\beta}. Following Duxbury (2020), \eqn{z(x,Y)} is a mapping
#' function that represents the endogenous network statistics computed on the network
#' of past events,\eqn{x}, and exogenous covariates, \eqn{Y}. The user provides these
#' covariates via the `formula` argument.
#'
#' This function provides two numerical optimization techniques to find the maximum
#' likelihood estimates for the associated parameters. First, this function allows
#' the user to use the \code{\link{optim}} function to find the associated parameters
#' based on the above likelihood function. Secondly, and by default, this function
#' employs a Newton-Rhapson iteration algorithm with line-searching to find
#' the unknown parameters (see Greene 2003 for a discussion of this algorithm). If desired, the user can
#' provide the initial searching values for both algorithms with the `starting.beta` argument.
#'
#' It's important to note that the modeling concerns of the conditional logistic regression apply to the
#' ordinal timing relational event model, such as no within-sequence fixed effects, that is,
#' a variable that does not vary within event cluster (i.e., a variable that is the
#' same for both the null and observed events). The function internally checks for
#' this and provides the user with a warning if any requested effects has no total within-event
#' variance. Moreover, any observed events that have no associated control events are
#' removed from the analysis as they provide no information to the log likelihood (see Greene 2003). The function
#' removes these events from the sequence prior to estimation.
#'
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#'
#' @references
#' Box-Steffensmeier, Janet and Bradford S. Jones. 2004. *Event History Modeling: A Guide for Social Scientists*. Cambridge University Press.
#'
#' Butts, Carter T. 2008. "A Relational Event Framework for Social Action." *Sociological Methodology* 38(1): 155-200.
#'
#' Duxbury, Scott. 2020. *Longitudinal Network Models*. Sage University Press. Quantitative Applications in
#' the Social Sciences: 192.
#'
#' Greene, William H. 2003. *Econometric Analysis*. Fifth Edition. Prentice Hall Press.
#'
#'
#'
#'@examples
#'#Creating a psuedo one-mode relational event sequence with ordinal timing
#'relational.seq <- simulateRESeq(n_actors = 8,
#'                                n_events = 50,
#'                                inertia = TRUE,
#'                                inertia_p = 0.10,
#'                                sender_outdegree = TRUE,
#'                                sender_outdegree_p = 0.05)
#'
#'#Creating a post-processing event sequence for the above relational sequence
#'post.processing <- processOMEventSeq(data = relational.seq,
#'                                     time = relational.seq$eventID,
#'                                     eventID = relational.seq$eventID,
#'                                     sender = relational.seq$sender,
#'                                     receiver = relational.seq$target,
#'                                     n_controls = 5)
#'
#'#Computing the sender-outdegree statistic for the above post-processing
#'#one-mode relational event sequence
#'post.processing$sender.outdegree <- computeSenderOutdegree(
#'                                    observed_time = relational.seq$eventID,
#'                                    observed_sender = relational.seq$sender,
#'                                    processed_time = post.processing$time,
#'                                    processed_sender = post.processing$sender,
#'                                    counts = TRUE)
#'
#'#Computing the inertia/repetition statistic for the above post-processing
#'#one-mode relational event sequence
#'post.processing$inertia <- computeRepetition(
#'                           observed_time = relational.seq$eventID,
#'                           observed_sender = relational.seq$sender,
#'                           observed_receiver = relational.seq$target,
#'                           processed_time = post.processing$time,
#'                           processed_sender = post.processing$sender,
#'                           processed_receiver = post.processing$receiver,
#'                           counts = TRUE)
#'
#'#Fitting a (ordinal) relational event model to the above one-mode relational
#'#event sequence
#'rem <- estimateREM(observed~sender.outdegree+inertia,
#'                   event.cluster = post.processing$time,
#'                   data=post.processing)
#'summary(rem) #summary of the relational event model
#'
#'#Fitting a (ordinal) relational event model to the above one-mode relational
#'#event sequence via the optim function
#'rem1 <- estimateREM(observed~sender.outdegree+inertia,
#'                   event.cluster = post.processing$time,
#'                   data=post.processing,
#'                   newton.rhapson=FALSE) #use the optim function
#'summary(rem1) #summary of the relational event model


estimateREM <- function(formula,
                           event.cluster,
                           data,
                           ordinal = TRUE,
                           multiple.events = FALSE, #are there multiple observed events per time
                           newton.rhapson = TRUE, #if the newton rhapson iteration sohuld be used
                           optim.method = "BFGS", #the optim method if desired
                           optim.control = list(), # a list of controls for the optim function
                           tolerance=1e-09, #the stopping absolute tolerance for the netwon rhaspon updated
                           maxit=20,#the maximum number of netwon rhaspon updates
                           starting.beta = NULL,
                           ...){ # a starting search vector
  lifecycle::deprecate_warn("1.0.0", "estimateREM()", "estimate_rem_logit()")
  base::message("Extracting user-provided data.") #send a message to the user
  #--------------------------------------------------------------------#
  # Step 1: The Function Starts By Extracting the Variables based on
  # the formula and data arguments.
  #--------------------------------------------------------------------#
  data.stats <- model.frame(formula, data = data) #extracting the variables
  outcome <- model.extract(data.stats,"response") #the outcome variable (1 = event; 0 = null)
  n.events <- sum(outcome) #the number of true observed events
  null.events <- length(outcome)-n.events #the number of null events
  net_stats <- model.matrix(formula, data = data) #extracting the model matrix
  effects <- colnames(net_stats)
  net_stats <- as.matrix(net_stats[,which(colnames(net_stats) != "(Intercept)")]) #removing the intercept
  effects <- effects[-which(effects== "(Intercept)")]
  #--------------------------------------------------------------------#
  # Step 2: Creating the event-based matrices
  #--------------------------------------------------------------------#
  eventIDS <- base::unique(event.cluster) #the unique event IDs per event
  #--------------------------------------------------------------------#
  # Doing an internal check to remove observed events with no null events
  #--------------------------------------------------------------------#
  event.ID.Null <- tapply(X = 1-outcome, INDEX = event.cluster, FUN = sum)
  extract.non.use <- names(event.ID.Null)[event.ID.Null==0] #removing observed events that have no controls
  if(length(extract.non.use)!=0){ #updating all model statistics
    outcome <- outcome[-which(event.cluster %in% extract.non.use)]
    net_stats <- as.matrix(net_stats[-which(event.cluster %in% extract.non.use),]) #maintaining a matrix structure
    event.cluster <- event.cluster[-which(event.cluster %in% extract.non.use)]
    eventIDS <- base::unique(event.cluster) #the unique event IDs per event
    n.events <- sum(outcome) #the number of true observed events
    null.events <- length(outcome)-n.events #the number of null events
  }
  #--------------------------------------------------------------------#
  # Creating the event-based matrices (kinda slow need to update here)
  #--------------------------------------------------------------------#
  base::message("Prepping data for numerical optimization.") #send a message to the user

  stats.by.event <- extractEventData(stats=net_stats,#the network statistics (and other user-provided covariates)
                                       outcome = outcome,#the event outcome indicator
                                       event_cluster = event.cluster, #the event clustering information
                                       names = c(effects,"dummy","id")) #the names of the data

  #--------------------------------------------------------------------#
  # Making a check to see if the user interested a within-event fixed effect
  #--------------------------------------------------------------------#
  varCheck <- checkVarianceData(stats.by.event,K = length(effects)) #checking the variance of the effects
  withinVar <- ifelse(any(varCheck==0),TRUE,FALSE) #if any of the variables have a total variance of 0
  if(withinVar==TRUE){
    base::warning("Located a variable that has no total within-variance across the event sequence.") #send a message to the user
  }
  whichCheck <- which(varCheck==0)
  #--------------------------------------------------------------------#
  # Step 3: Starting the Optimization Procedure (need to find a better way to start the process)
  #--------------------------------------------------------------------#
  if(is.null(starting.beta)){beta.1 <- rep(0,ncol(net_stats))
  }else{ beta.1 <- starting.beta }

  #--------------------------------------------------------------------#
  # If Newton.Rhapson = TRUE, do it by hand via newton's iteration with line search alog
  #--------------------------------------------------------------------#
  base::message("Starting optimzation for parameters.") #send a message to the user
  if(newton.rhapson == TRUE){
    optim.used <- FALSE
    converged <- TRUE #check model convergence after
    loglike.1 <- remloglike(beta.1,stats.by.event) #The starting estimates
    #--------------------------------------------------------------------#
    # Step 4a: Set up line search algorithm (see Greene 2003: 943)
    #--------------------------------------------------------------------#
    #these lambas are more than likely not optimal (Note: 1 is the traditional step
    #size for the Newton Methods for Numerical Optimization)
    line.search.lamba <- c(0.25,0.50,0.75,1.00,1.25,1.50,1.75,2.00)
    n.lambdas <- length(line.search.lamba)
    # B(i+1) = B(i) - (Hessian)^1*Gradient*lamba(k)
    search.algo <- data.frame(loglike = rep(NA,maxit),
                              loglikeDIF = rep(NA,maxit),
                              line.step = rep(NA,maxit))
    i<-1
    search.algo[i,]<-c(loglike.1,NA,NA)
    for(i in 2:maxit){ #for up to the maximum iterations
      gradient.betai <- remGradient(beta.1,stats.by.event)
      hessian.betai <-  remHessian(beta.1,stats.by.event)
      # Now to find the best update via line search algo
      results <- data.frame(loglike = rep(0,length(line.search.lamba)),
                            linelamba = line.search.lamba) # a dataframe to store the loglikes
      step.increment <- solve(hessian.betai,gradient.betai)
      for(j in 1:n.lambdas){
        update <- beta.1 - step.increment*line.search.lamba[j]
        results$loglike[j] <- remloglike(beta = update,stats.by.event)
      }
      #take the first value if multiple (https://statacumen.com/teach/SC1/SC1_11_LogisticRegression.pdf )
      best.line <- results$linelamba[which(results$loglike == max(results$loglike))[1]]
      beta.2 <- beta.1 - step.increment*best.line
      loglike.2 <-  remloglike(beta.2,stats.by.event) #The starting estimates#the best loglikelihood
      diff.log.like <- loglike.2 - loglike.1 #differences in log likelihoods
      if(abs(diff.log.like) < tolerance ){break} #if the differnece is less than the tolerance, stop (climbing is done)
      beta.1<-beta.2
      loglike.1 <- loglike.2
      search.algo[i,]<-c(loglike.1,diff.log.like,best.line)
    }
    beta.1<-beta.2
    loglike.1 <- loglike.2
    search.algo[i,]<-c(loglike.1,diff.log.like,best.line)
    search.algo <- search.algo[!is.na(search.algo$loglike),]
    if(i == maxit & abs(diff.log.like) > tolerance){converged<-FALSE;
    base::warning("Netwon Rhapson iteration failed to converge!!!")}
    df.null <- 0 #all parameters are equal to 0
    loglikelihood.null <- remloglike(rep(0, ncol(net_stats)),stats.by.event) #The starting estimates
    loglikelihood.full <- loglike.1 #the updated estimates log likelihood for step i + 1
    parameters <- beta.1 #the estimated parameters
    hessian <- remHessian(parameters,stats.by.event) #the current hessian value
    gradient <- remGradient(parameters,stats.by.event) #the current hessian value
    covariance.mat <- solve(-hessian) #the information matrix
    se.parameter <- sqrt(diag(covariance.mat)) #the standard errors
    z.values <- parameters/se.parameter #z-values of parameters
    p.values <- 2*pnorm(abs(z.values),lower.tail = FALSE)#p-values of parameters
    df.full <- length(parameters) #the degrees of freedom
    AIC <- 2*length(parameters) - 2*loglikelihood.full
    BIC <- length(parameters)*log(n.events)-2*loglikelihood.full
    chi.stat<-2*(loglikelihood.full-loglikelihood.null)
    loglikelihood.test <- 1-pchisq(chi.stat,df.full-df.null)
    base::message("Optimzation via Netwon's Method is complete.") #send a message to the user

  }else{
    optim.used <- TRUE
    # optimization via optim
    optim_results <- stats::optim(par = beta.1, #the starting parameters
                                  fn = remloglikeOPTIM, #the log likelihood function
                                  method = optim.method, #the optimization method
                                  event_stats = stats.by.event, #the network stats (and user-provided stats)
                                  hessian = TRUE, #return the hessian matrix
                                  control = optim.control) #the optim control list (user-provided)
    converged <- ifelse(optim_results$convergence==0,TRUE,FALSE) #check model convergence after
    loglikelihood.null <- remloglike(rep(0,length(optim_results$par)),stats.by.event)#likelihood with all effects being 0
    df.null <- 0 #all parameters are equal to 0
    parameters <- optim_results$par #the estimated parameters
    loglikelihood.full <- -optim_results$value    #Since optim minimzes, flip the sign for maximization
    hessian <- -optim_results$hessian #Since optim minimzes, flip the sign for maximization
    gradient <- NULL #optim does not compute the gradient so making it null
    covariance.mat <- solve(-hessian) #the covariance matrix
    se.parameter <- sqrt(diag(covariance.mat))
    z.values <- parameters/se.parameter
    p.values <- 2*pnorm(abs(z.values),lower.tail = FALSE)
    df.full <- length(parameters)
    AIC <- 2*length(parameters) - 2*loglikelihood.full
    BIC <- length(parameters)*log(n.events)-2*loglikelihood.full
    chi.stat<-2*(loglikelihood.full-loglikelihood.null)
    loglikelihood.test <- 1-pchisq(chi.stat,df.full-df.null)
    search.algo <- NULL #no searching algorithm for the optim method
  }
  #--------------------------------------------------------------------#
  # Sending a warning message if the model failed to converge
  #--------------------------------------------------------------------#
  if(converged==FALSE){
    base::warning("The relational event model failed to converge after the maximum iterations. \n
                You may try to restart the function with a greater value of maxit (only when \n
                Newton.Rhapson is set to true) or provide a starting parameter vector \n
                via starting.beta You can also tune the absolute difference tolerance \n
                argument (only when Newton.Rhapson is set to true). Cheers!") }
  parameters<-as.vector(parameters) #the parameter estimates
  names(parameters)<-effects #the names of the parameters
  z.values<-as.vector(z.values) #the z-scores for the parameters
  names(z.values)<-effects#the names of the z-scores
  p.values<-as.vector(p.values) #the p-values
  names(p.values)<-effects #the names of the p-values
  if(withinVar==TRUE){ #cleaning the results based on bad variable input
    parameters[whichCheck]<-0.00 #making bad variables to be NA
    z.values[whichCheck]<-0.00 #making bad variables to be NA
    p.values[whichCheck]<-0.00 #making bad variables to be NA
    se.parameter[whichCheck]<-0.00 #making bad variables to be NA
    covariance.mat[whichCheck,]<-0.00 #making bad variables to be NA
    covariance.mat[,whichCheck]<-0.00 #making bad variables to be NA
  } #and finished!
  results <- list( #combining all results to be outputted
    optimization.method = ifelse(optim.used==TRUE,"optim","Newton-Rhapson with line searching"), #optimization method
    converged = converged,#did the model converge?
    loglikelihood.null=loglikelihood.null, #the null log likelihood
    loglikelihood.full=loglikelihood.full, #the full log likelihood
    chi.stat=chi.stat, #the lr test chi stat
    loglikelihood.test=loglikelihood.test, #the likelihood test pvalue
    df.null=df.null, #the df for the null model
    df.full=df.full, #the df for the full model
    parameters=parameters, #the MLE parameter estimates
    hessian=hessian, #the estimated hessian
    gradient=ifelse(optim.used == TRUE, NA, as.vector(gradient)), #the gradient
    se.parameter=se.parameter, #the standard errors of the MLE esitmates
    covariance.mat=covariance.mat, #the covariance matrix
    z.values=z.values, #the z-values
    p.values=p.values, #the p-values
    AIC=AIC, #the AIC for the model
    BIC=BIC, #the BIC for the model
    call = match.call(), #the formula
    n.events = n.events, #the number of observed events
    null.events=null.events, #the number of null events
    newton.iterations = ifelse(optim.used == TRUE, NA, i), #the number of newton iterations
    search.algo=search.algo) #the search algorithm results
  class(results) <- "dream" #changing the object class
  return(results) #returning the REM
}
