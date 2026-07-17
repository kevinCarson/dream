#' @title Fit a Maximum Likelihood Relational Event Model (REM) to A Processed Relational Event Sequence
#' @name estimate_rem
#' @param formula A formula object where the covariates are on the right hand
#' side of ~. The names of the covariates must follow the names in the
#' `statistics` element of the `data` argument, since the right hand side covariates
#' are taken from this list. The dependent variable does not need to be defined, since it is internally
#' defined based upon the `data` argument. This is the same argument found in \code{\link[stats]{lm}} and \code{\link[stats]{glm}}.
#' @param data An object of class `dream_sequence` that contains the processed relational event sequence and
#' statistics that are included in the `formula` argument.
#' @param multiple.events TRUE/FALSE. Currently, this function assumes that only one event occurs per
#' event cluster (i.e., time point). Future versions of the package will include estimation options for
#' multiple events per time point, commonly referred to as tied events, via the Breslow approximation
#' technique (see Box-Steffensmeier and Jones 2004). At this moment, this argument is preset to FALSE
#' and should not be modified by the user.
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
#' @return An object of class "dream_rem" as a list containing the following components:
#' \itemize{
#'   \item \code{optimization.method} - The optimzation method used to find the parameters.
#'   \item \code{converged} - TRUE/FALSE. TRUE indicates that the REM converged.
#'   \item \code{loglikelihood.null} - The log likelihood of the null model (i.e., the model where the parameters are assumed to be 0).
#'   \item \code{loglikelihood.full} - The log likelihood of the estimated model.
#'   \item \code{chi.stat} - The chi-statistic of the likelihood ratio test.
#'   \item \code{loglikelihood.test} - The p-value of the likelihood ratio test.
#'   \item \code{df.null} - The degrees of freedom of the null model.
#'   \item \code{df.full} - The degrees of freedom of the full model.
#'   \item \code{parameters} - The MLE parameter estimates.
#'   \item \code{hessian} - The estimated hessian matrix.
#'   \item \code{gradient} - The estimated gradient vector.
#'   \item \code{se.parameter} - The standard errors of the MLE parameter estimates.
#'   \item \code{covariance.mat} - The estimated variance-covariance matrix.
#'   \item \code{z.values} - The z-scores for the MLE parameter estimates.
#'   \item \code{p.values} - The p-values for the MLE parameter estimates.
#'   \item \code{AIC} - The AIC of the estimated REM.
#'   \item \code{BIC} - The BIC of the estimated REM.
#'   \item \code{n.events} - The number of observed events in the relational event sequence.
#'   \item \code{null.events} - The number of control events in the relational event sequence.
#'   \item \code{newton.iterations} - The number of Newton-Rhapson iterations.
#'   \item \code{search.algo} - A data.frame object that contains the Newton-Rhapson searching algorithm results.
#'}
#' @export
#' @seealso
#' \code{\link{predict.dream_rem}}, \code{\link{vcov.dream_rem}}, \code{\link{logLik.dream_rem}},
#' \code{\link{AIC}}, \code{\link{gof_rem}}, \code{\link{residuals.dream_rem}},\code{\link{plot.dream_rem}},\code{\link{coef.dream_rem}}.
#'



#' @description
#' `r lifecycle::badge("stable")`
#'
#' This function estimates the ordinal and interval timing relational event model by maximizing the
#' likelihood function given by Butts (2008) via maximum likelihood estimation. A nice outcome
#' is that the ordinal timing relational event model is equivalent to the conditional logistic
#' regression (see Greene 2003; for R functions, see \code{\link[survival]{clogit}}). In
#' addition, based on this outcome and the structure of the data, this function can estimate
#' the Cox proportional hazards model (see Box-Steffensmeier and Jones 2004; for R functions, see \code{\link[survival]{coxph}})
#' given that the likelihood functions are equivalent. An important assumption this model
#' makes is that only one event occurs at each time point. If this is unfeasible for
#' the user's specific dataset, we encourage the user to see the \code{\link[survival]{clogit}}
#' function for the Breslow approximation technique (Box-Steffensmeier and Jones 2004). Future
#' versions of the package will include options for tied event data (e.g., multiple events at one time point).
#'
#' @details
#'
#' This function estimates the ordinal and interval timing relational event model by maximizing the
#' likelihood function provided in the seminal REM paper by Butts (2008) via maximum likelihood estimation. The
#' ordinal timing likelihood function is:
#' \deqn{L(A_t|\beta) = \prod_{i=1}^{|A_t|} \frac{\lambda_{a_i}}{\sum_{a' \in M_t \lambda_{a'}}}}
#' where, following Butts (2008) and Duxbury (2020), \eqn{A_t} is the relational event sequence,
#' \eqn{\lambda_{a_i}} is the hazard rate for event *i*, which is formulated to be equal to
#' \eqn{exp(\beta^{T}z(x,Y))}, that is, the linear combination of user-specific covariates, \eqn{z(x,Y)}, and associated
#' REM parameters, \eqn{\beta}. The user provides these covariates via the `formula` argument. \eqn{M_t} is the support set for event \eqn{a_i \in A_t}. The likelihood function for the interval timing relational event
#' model is:
#' \deqn{L(A_t | \beta) = [ \prod_{i=1}^{|A_t|} \lambda_i \prod_{j \in M_{\tau(i)}} \exp(-\lambda_j \{\tau(i) - \tau(i-1) \}) ] \times [\prod_{j \in M_t} \exp(-\lambda_j \{t - \tau(M) \})]}
#' where \eqn{\tau(i)} is the time of the observed (realized) event *i* and \eqn{t} is the time that marks the end of the relational event sequence.
#' Following Duxbury (2020), \eqn{z(x,Y)} is a mapping function that represents the endogenous network statistics computed on the network
#' of past events,\eqn{x}, and exogenous covariates, \eqn{Y}. In comparison to the ordinal
#' timing relational event formulation, the hazard rate for event *i*, \eqn{\lambda_{a_i}}, includes
#' the baseline hazard rate (the intercept), \eqn{exp(\beta_0 + \beta^{T}z(x,Y))}. If \eqn{t} is not known
#' by the user, then the interval timing likelihood is:
#' \deqn{L(A_t | \beta) =  \prod_{i=1}^{|A_t|} \lambda_i \prod_{j \in M_{\tau(i)}} \exp(-\lambda_j \{\tau(i) - \tau(i-1) \})  }
#' In this case, the likelihood function is the same as employed in the \code{\link[remstimate]{remstimate}}
#' for interval timing relational event models. The values for \eqn{t} are taken from the
#' `data` object.
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
#' @author Kevin A. Carson <kacarson@arizona.edu> and Diego F. Leal <dflc@arizona.edu>
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
#'relational.seq <- simulate_rem_seq(n_actors = 8,
#'                                   n_events = 50,
#'                                   inertia = TRUE,
#'                                   inertia_p = 0.10,
#'                                   sender_outdegree = TRUE,
#'                                   sender_outdegree_p = 0.05)
#'
#'#Creating a post-processing event sequence for the above relational sequence
#'post.processing <-  create_res(type = "one-mode",
#'                               ordinal = TRUE,
#'                               riskset = "fixed",
#'                               time = relational.seq$eventID,
#'                               sender = as.character(relational.seq$sender),
#'                               receiver = as.character(relational.seq$target))
#'
#'#Computing the sender-outdegree statistic for the above post-processing
#'#one-mode relational event sequence
#'post.processing <- dreamstats_degree(formation = "sender-outdegree",
#'                                   data = post.processing,
#'                                   halflife = 2)
#'
#'#Computing the inertia/repetition statistic for the above post-processing
#'#one-mode relational event sequence
#'post.processing <- dreamstats_repetition(data = post.processing,
#'                                       halflife = 2)
#'
#'#Fitting an ordinal timing relational event model to the above one-mode relational
#'#event sequence
#'rem <- estimate_rem(~ sender.outdegree + repetition,
#'                          data=post.processing)
#'summary(rem) #summary of the relational event model
#'
#'vcov(rem) #printing the variance-covariance matrix
#'logLik(rem) #printing the model log-likelihood
#'AIC(rem) #printing the model AIC
#'rates <- predict(rem) #extracting the predicted event rates
#'
#'
#'#Fitting a (ordinal) relational event model to the above one-mode relational
#'#event sequence via the optim function
#'rem1 <- estimate_rem(~ sender.outdegree + repetition,
#'                    data=post.processing,
#'                    newton.rhapson=FALSE)
#'summary(rem1) #summary of the relational event model
#'
#'
#'# a psuedo relational event sequence
#'events <- data.frame(time = 1:18, eventID = 1:18,
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
#'# Creating a dynamic one-mode relational risk set with p = 1.00 (all true events)
#'# and 5 controls based upon the interval timing relational event framework
#'eventSet <- create_res(ordinal = FALSE,
#'                       t = max(events$time) + rexp(1),
#'                       riskset = "cumulative",
#'                       type = "one-mode",
#'                       time = events$time,
#'                       sender = events$sender,
#'                       receiver = events$target,
#'                       p_samplingobserved = 1.00,
#'                       case_control=TRUE,
#'                       n_controls = 5,
#'                       seed = 9999)
#'
#'#Computing the sender indegree statistic for the relational event sequence
#'eventSet <- dreamstats_degree(formation = "sender-indegree",
#'                              data = eventSet,
#'                              halflife = 2)
#'
#'#Computing the outgoing two path statistic for the relational event sequence
#'eventSet <- dreamstats_triads(formation = "OTP",
#'                              data = eventSet,
#'                              halflife = 2)
#'
#'#Fitting an interval timing relational event model to the above one-mode relational
#'#event sequence
#'rem.interval <- estimate_rem(~ sender.indegree + outgoing.two.paths,
#'                             data=eventSet)
#'summary(rem.interval) #summary of the relational event model
#'
#'rem.interval.optim <- estimate_rem(~ sender.indegree + outgoing.two.paths,
#'                                   data=eventSet,
#'                                   newton.rhapson=FALSE)
#'summary(rem.interval.optim) #summary of the relational event model

estimate_rem <- function(formula,
                         data,
                         newton.rhapson = TRUE, #if the newton rhapson iteration sohuld be used
                         optim.method = "BFGS", #the optim method if desired
                         optim.control = list(), # a list of controls for the optim function
                         tolerance=1e-09, #the stopping absolute tolerance for the netwon rhaspon updated
                         maxit=100,#the maximum number of netwon rhaspon updates
                         starting.beta = NULL,
                         multiple.events = FALSE,
                         ...){ # a starting search vector

  #--------------------------------------------------------------------#
  # Step 1: The Function Starts By Extracting the Variables based on
  # the formula and data arguments.
  #--------------------------------------------------------------------#
  if(!inherits(data,"dream_sequence")) base::stop("The `data` argument is not a `dream_sequence` object.")
  ordinal <- data$ordinal
  interevent <- data$interevent_times #the timing between the events
  rem.data <-data
  tknown <- TRUE #presetting the values
  if(data$t < 0) tknown <- FALSE #presetting the values based upon variable inputs
  data <- as.data.frame(data, all.events = FALSE)
  event.cluster <- data$time #the sampled realized event times
  specification <- formula
  formula = as.formula(paste("observed ~" , as.character(formula)[2]))
  data.stats <- model.frame(formula, data = data) #extracting the variables
  outcome <- model.extract(data.stats,"response") #the outcome variable (1 = event; 0 = null)
  out.mat <- data.frame(time=data$time, outcome=outcome)
  n.events <- sum(outcome) #the number of true observed events
  null.events <- length(outcome)-n.events #the number of null events
  net_stats <- model.matrix(formula, data = data) #extracting the model matrix
  effects <- colnames(net_stats)
  if(ordinal){
    net_stats <- as.matrix(net_stats[,which(colnames(net_stats) != "(Intercept)")]) #removing the intercept
    effects <- effects[-which(effects== "(Intercept)")]
  }
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

  stats.by.event <-extractEventData(stats=net_stats,#the network statistics (and other user-provided covariates)
                                             outcome = outcome,#the event outcome indicator
                                             event_cluster = event.cluster, #the event clustering information
                                             names = c(effects,"dummy","id")) #the names of the data

  #--------------------------------------------------------------------#
  # Making a check to see if the user interested a within-event fixed effect
  #--------------------------------------------------------------------#
  varCheck <- checkVarianceData(stats.by.event,K = length(effects)) #checking the variance of the effects
  withinVar <- ifelse(any(varCheck==0),TRUE,FALSE) #if any of the variables have a total variance of 0
  if(withinVar==TRUE & ordinal){
    base::warning("Located a variable that has no total within-variance across the event sequence.") #send a message to the user
  }
  whichCheck <- which(varCheck==0)
  #--------------------------------------------------------------------#
  # Step 3: Starting the Optimization Procedure (need to find a better way to start the process)
  #--------------------------------------------------------------------#
  if(is.null(starting.beta)){beta.1 <- rep(0,ncol(net_stats))
  }else{ beta.1 <- starting.beta }
  rem.type <- ifelse(ordinal,"ordinal", "interval")
  #finding the indices where the event exists
  indices <- lapply(stats.by.event, function(x){as.numeric(which(x[,"dummy"] == 1)) -1})
  indices <- unlist(indices) #unlisting the values
  #extracting only the relevent endogenous and exogenous network statistics
  stats.by.event <- lapply(stats.by.event, function(x){
    x[,which(!(colnames(x) %in% c("dummy","id")))] })
  #--------------------------------------------------------------------#
  # If Newton.Rhapson = TRUE, do it by hand via newton's iteration with line search alog
  #--------------------------------------------------------------------#
  if(newton.rhapson == TRUE){
    optim.used <- FALSE
    converged <- TRUE #check model convergence after
    loglike.1 <- switch(rem.type,
                        ordinal = -remloglik_ordinal(beta.1,stats.by.event,indices),
                        interval = -remloglik_interval(beta.1,stats.by.event,interevent,indices,tknown))
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

      gradient.betai <- switch(rem.type,
                               ordinal = remgrad_ordinal(beta.1,stats.by.event,indices),
                               interval = remgrad_interval(beta.1,stats.by.event,interevent,indices,tknown))
      hessian.betai <- switch(rem.type,
                              ordinal = remhess_ordinal(beta.1,stats.by.event),
                              interval = remhess_interval(beta.1,stats.by.event,interevent,tknown))
      # Now to find the best update via line search algo
      results <- data.frame(loglike = rep(0,length(line.search.lamba)),
                            linelamba = line.search.lamba) # a dataframe to store the loglikes
      step.increment <- solve(hessian.betai,gradient.betai)
      for(j in 1:n.lambdas){
        update <- beta.1 - step.increment*line.search.lamba[j]
        results$loglike[j] <- switch(rem.type,
                                     ordinal = -remloglik_ordinal(update,stats.by.event,indices),
                                     interval = -remloglik_interval(update,stats.by.event,interevent,indices,tknown))

      }
      #take the first value if multiple (https://statacumen.com/teach/SC1/SC1_11_LogisticRegression.pdf )
      best.line <- results$linelamba[which(results$loglike == max(results$loglike))[1]]
      beta.2 <- beta.1 - step.increment*best.line
      loglike.2 <-  switch(rem.type,
                           ordinal = -remloglik_ordinal(beta.2,stats.by.event,indices),
                           interval = -remloglik_interval(beta.2,stats.by.event,interevent,indices,tknown))
      #The starting estimates#the best loglikelihood
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
    null <- rep(0, ncol(net_stats))
    loglikelihood.null <-  switch(rem.type,
                                  ordinal = -remloglik_ordinal(null,stats.by.event,indices),
                                  interval = -remloglik_interval(null,stats.by.event,interevent,indices,tknown))
    loglikelihood.full <- loglike.1 #the updated estimates log likelihood for step i + 1
    parameters <- beta.1 #the estimated parameters
    hessian <- switch(rem.type,
                      ordinal = remhess_ordinal(parameters,stats.by.event),
                      interval = remhess_interval(parameters,stats.by.event,interevent,tknown))
    gradient <- switch(rem.type,
                       ordinal = remgrad_ordinal(parameters,stats.by.event,indices),
                       interval = remgrad_interval(parameters,stats.by.event,interevent,indices,tknown))
    covariance.mat <- solve(-hessian) #the information matrix
    se.parameter <- sqrt(diag(covariance.mat)) #the standard errors
    z.values <- parameters/se.parameter #z-values of parameters
    p.values <- 2*pnorm(abs(z.values),lower.tail = FALSE)#p-values of parameters
    df.full <- length(parameters) #the degrees of freedom
    AIC <- 2*length(parameters) - 2*loglikelihood.full
    BIC <- length(parameters)*log(n.events)-2*loglikelihood.full
    chi.stat<-2*(loglikelihood.full-loglikelihood.null)
    loglikelihood.test <- 1-pchisq(chi.stat,df.full-df.null)

  }else{
    optim.used <- TRUE
    # optimization via optim

    if(ordinal){
      optim_results <- stats::optim(par = beta.1, #the starting parameters
                                    fn = remloglik_ordinal, #the log likelihood function
                                    method = optim.method, #the optimization method
                                    event_stats = stats.by.event, #the network stats (and user-provided stats)
                                    obsindex = indices,
                                    hessian = TRUE, #return the hessian matrix
                                    control = optim.control) #the optim control list (user-provided)
    }else{
      optim_results <- stats::optim(par = beta.1, #the starting parameters
                                    fn = remloglik_interval, #the log likelihood function
                                    method = optim.method, #the optimization method
                                    event_stats = stats.by.event, #the network stats (and user-provided stats)
                                    obsindex = indices,
                                    interevent = interevent,
                                    tknown=tknown,
                                    hessian = TRUE, #return the hessian matrix
                                    control = optim.control) #the optim control list (user-provided)
    }
    converged <- ifelse(optim_results$convergence==0,TRUE,FALSE) #check model convergence after
    null <- rep(0, length(optim_results$par))
    loglikelihood.null <-  switch(rem.type,
                                  ordinal = -remloglik_ordinal(null,stats.by.event,indices),
                                  interval = -remloglik_interval(null,stats.by.event,interevent,indices,tknown))
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
  if(withinVar==TRUE & ordinal){ #cleaning the results based on bad variable input
    parameters[whichCheck]<-0.00 #making bad variables to be NA
    z.values[whichCheck]<-0.00 #making bad variables to be NA
    p.values[whichCheck]<-0.00 #making bad variables to be NA
    se.parameter[whichCheck]<-0.00 #making bad variables to be NA
    covariance.mat[whichCheck,]<-0.00 #making bad variables to be NA
    covariance.mat[,whichCheck]<-0.00 #making bad variables to be NA
  } #and finished!
  results <- list( #combining all results to be outputted
    outcome=outcome,
    rem.data = rem.data, #the relational event sequence dataset
    rem.type = rem.type,
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
    specification = specification,
    n.events = n.events, #the number of observed events
    null.events=null.events, #the number of null events
    newton.iterations = ifelse(optim.used == TRUE, as.numeric(optim_results$counts[1]), i), #the number of newton iterations
    search.algo=search.algo,
    statistics = net_stats) #the search algorithm results
  class(results) <- "dream_rem" #changing the object class
  return(results) #returning the REM
}
















#' Summary Method for dreamrem Objects
#'
#' Summarizes the results of an ordinal timing relational event model.
#'
#' @param object An object of class "dream_rem".
#' @param digits The number of digits to print after the decimal point.
#' @param ... Additional arguments (currently unused).
#' @return A list of summary statistics for the relational event model
#' including parameter estimates, (null) likelihoods, and tests of significance
#' for likelihood ratios and estimated parameters.
#' @export
#'
summary.dream_rem <- function(object,digits=6,...) {
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
    iterations = object$newton.iterations,
    type = object$rem.type
  )
  class(res) <- "summary.dream_rem"
  return(res)
}

#' Print Method for dreamrem Model
#'
#' @param x An object of class "dream_rem".
#' @param digits The number of digits to print after the decimal point.
#' @param ... Additional arguments (currently unused).
#' @return No return value. Prints out the main results of a 'dream' summary object.
#' @export
print.summary.dream_rem  <- function(x,digits=6,...) {
  if(x$type == "ordinal") cat("Ordinal Timing Relational Event Model\n\n")
  if(x$type == "interval") cat("Interval Timing Relational Event Model\n\n")

  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat("Relational Event Sequence Information:\n")

  cat("The number of realized events =",x$nevents,"\nThe number of control events =", x$nullevents,"\n")

  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)

  cat("\nREM Fit Information:\n")
  cat("Null Model Likelihood =",round(x$nullLogLik,digits),"\n")
  cat("Full Model Likelihood =",round(x$logLik,digits),"\n")
  cat("Likelihood Ratio Test: ",round(x$chiStat,digits)," (df=",x$df,"; p-value: ",round(x$plrtest,digits) ,")\n",sep = "")
  cat("AIC =",round(x$AIC,digits),"\n")
  cat("BIC =",round(x$BIC,digits),"\n")
  cat("Number of Newton Iterations =",x$iterations,"\n")

}

#' Print Method for dreamrem Model
#'
#' @param x An object of class "dream_rem".
#' @param digits The number of digits to print after the decimal point.
#' @param ... Additional arguments (currently unused).
#' @return No return value. Prints out the main results of a 'dream' object.
#' @export

print.dream_rem  <- function(x,digits=6,...) {
  if(x$rem.type == "ordinal") cat("Ordinal Timing Relational Event Model\n\n")
  if(x$rem.type == "interval") cat("Interval Timing Relational Event Model\n\n")

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
    cat("\nNumber of Newton Iterations:",x$newton.iterations,"\n")
  }
  invisible(x)
}















#' Extract the model log-likelihood from Relational Event Model Fits
#'
#' This function extracts the model loglikelhood from estimated
#' relational event model fits.
#'
#' @param object An object of class "dream_rem".
#' @param ... Additional arguments for other methods.
#' @param REML From the generic `logLik` function. Set to FALSE and does not
#' need to changed by the user.
#' @export
#' @examples
#'#Creating a psuedo one-mode relational event sequence with ordinal timing
#'relational.seq <- simulate_rem_seq(n_actors = 8,
#'                                   n_events = 50,
#'                                   inertia = TRUE,
#'                                   inertia_p = 0.10,
#'                                   sender_outdegree = TRUE,
#'                                   sender_outdegree_p = 0.05)
#'
#'#Creating a post-processing event sequence for the above relational sequence
#'post.processing <-  create_res(type = "one-mode",
#'                               ordinal = TRUE,
#'                               riskset = "fixed",
#'                               time = relational.seq$eventID,
#'                               sender = as.character(relational.seq$sender),
#'                               receiver = as.character(relational.seq$target),
#'                               case_control=TRUE,
#'                               n_controls = 5)
#'
#'#Computing the sender-outdegree statistic for the above post-processing
#'#one-mode relational event sequence
#'post.processing <- dreamstats_degree(formation = "sender-outdegree",
#'                                   data = post.processing,
#'                                   halflife = 2)
#'
#'#Computing the inertia/repetition statistic for the above post-processing
#'#one-mode relational event sequence
#'post.processing <- dreamstats_repetition(data = post.processing,
#'                                         halflife = 2)
#'
#'#Fitting an ordinal timing relational event model to the above one-mode relational
#'#event sequence
#'rem <- estimate_rem(~ sender.outdegree + repetition,
#'                          data=post.processing)
#'logLik(rem)
#'

logLik.dream_rem <- function(object,...,REML = FALSE){
  val <- object$loglikelihood.full
  attr(val, "df") <- object$df.full
  attr(val, "nobs") <- object$n.events
  class(val) <- "logLik"
  val
}

#' Extract the ML parameter estimates from Relational Event Model Fits
#'
#' This function extracts the Maximum Likelihood (ML) parameter estimates from estimated
#' relational event model fits.
#'
#' @param object An object of class "dream_rem".
#' @param ... Additional arguments for other methods.
#' @export
#' @examples
#'#Creating a psuedo one-mode relational event sequence with ordinal timing
#'relational.seq <- simulate_rem_seq(n_actors = 8,
#'                                   n_events = 50,
#'                                   inertia = TRUE,
#'                                   inertia_p = 0.10,
#'                                   sender_outdegree = TRUE,
#'                                   sender_outdegree_p = 0.05)
#'
#'#Creating a post-processing event sequence for the above relational sequence
#'post.processing <-  create_res(type = "one-mode",
#'                               ordinal = TRUE,
#'                               riskset = "fixed",
#'                               time = relational.seq$eventID,
#'                               sender = as.character(relational.seq$sender),
#'                               receiver = as.character(relational.seq$target),
#'                               case_control=TRUE,
#'                               n_controls = 5)
#'
#'#Computing the sender-outdegree statistic for the above post-processing
#'#one-mode relational event sequence
#'post.processing <- dreamstats_degree(formation = "sender-outdegree",
#'                                   data = post.processing,
#'                                   halflife = 2)
#'
#'#Computing the inertia/repetition statistic for the above post-processing
#'#one-mode relational event sequence
#'post.processing <- dreamstats_repetition(data = post.processing,
#'                                         halflife = 2)
#'
#'#Fitting an ordinal timing relational event model to the above one-mode relational
#'#event sequence
#'rem <- estimate_rem(~ sender.outdegree + repetition,
#'                          data=post.processing)
#'coef(rem) #summary of the relational event model
#'

coef.dream_rem  <- function(object,...){
  object$parameters
}


#' Extract variance-covariance matrix from Relational Event Model Fits
#'
#' This function extracts the variance-covariance matrix from estimated
#' relational event model fits.
#'
#' @param object An object of class "dream_rem".
#' @param ... Additional arguments for other methods.
#' @export
#' @examples
#'#Creating a psuedo one-mode relational event sequence with ordinal timing
#'relational.seq <- simulate_rem_seq(n_actors = 8,
#'                                   n_events = 50,
#'                                   inertia = TRUE,
#'                                   inertia_p = 0.10,
#'                                   sender_outdegree = TRUE,
#'                                   sender_outdegree_p = 0.05)
#'
#'#Creating a post-processing event sequence for the above relational sequence
#'post.processing <-  create_res(type = "one-mode",
#'                               riskset = "fixed",
#'                               ordinal = TRUE,
#'                               time = relational.seq$eventID,
#'                               sender = as.character(relational.seq$sender),
#'                               receiver = as.character(relational.seq$target),
#'                               case_control=TRUE,
#'                               n_controls = 5)
#'
#'#Computing the sender-outdegree statistic for the above post-processing
#'#one-mode relational event sequence
#'post.processing <- dreamstats_degree(formation = "sender-outdegree",
#'                                   data = post.processing,
#'                                   halflife = 2)
#'
#'#Computing the inertia/repetition statistic for the above post-processing
#'#one-mode relational event sequence
#'post.processing <- dreamstats_repetition(data = post.processing,
#'                                         halflife = 2)
#'
#'#Fitting an ordinal timing relational event model to the above one-mode relational
#'#event sequence
#'rem <- estimate_rem(~ sender.outdegree + repetition,
#'                          data=post.processing)
#'vcov(rem)
#'

vcov.dream_rem  <- function(object,...){
  vcov <- object$covariance.mat
  colnames(vcov) <- rownames(vcov) <- names(object$parameters)
  vcov
}

#' Predict method for Relational Event Model Fits
#'
#' Predicted event hazard rates based on `dream_rem` relational event model objects.
#'
#' @param object An object of class "dream_rem".
#' @param newdata If requested, either an object of `data.frame` or `dream_sequence` that contains the transformed
#' variables.
#' @param se.fit TRUE/FALSE. If TRUE, the standard errors of the predicted event rates will be returned.
#' @param ... Additional arguments for other methods.
#' @details
#' Following Butts (2008: 166), the rate for an event \eqn{a_i} at time \eqn{t} is formulated as:
#' \deqn{\lambda(a_i) = \exp[\lambda_0 + \theta^T z(s(a_i), r(a_i), x_i, A_t)]}
#' where \eqn{z()} is a mapping function that represents a set of sufficient
#' statistics for event \eqn{a_i}. \eqn{z()} represents the covariates included
#' in the estimated relational event model (i.e., the `object` argument) and
#' \eqn{\theta} are the estimated parameters/coefficients. \eqn{\lambda_0} is the
#' baseline rate across the relational event sequence and the intercept in the
#' interval relational event models and is arbitrary in the ordinal relational event
#' models. In addition, the standard error for each predicted value can be returned
#' when `se.fit` is set to TRUE, and is based upon the delta method, where the
#' standard error of the predicted values are:
#' \deqn{SE(\hat{\lambda}) = \sqrt(\lambda^2 diag(X\Sigma X^T))}
#' where \eqn{\hat{\lambda}} is the vector of predicted event rates, \eqn{X} is the
#' matrix of event covariates, and \eqn{\Sigma} is the covariance matrix of the
#' relational event model.
#'
#' @references
#' Butts, Carter T. 2008. "A Relational Event Framework for Social Action." *Sociological Methodology* 38(1): 155-200.
#' @examples
#'#Creating a psuedo one-mode relational event sequence with ordinal timing
#'relational.seq <- simulate_rem_seq(n_actors = 8,
#'                                   n_events = 50,
#'                                   inertia = TRUE,
#'                                   inertia_p = 0.10,
#'                                   sender_outdegree = TRUE,
#'                                   sender_outdegree_p = 0.05)
#'
#'#Creating a post-processing event sequence for the above relational sequence
#'post.processing <-  create_res(type = "one-mode",
#'                               ordinal = TRUE,
#'                                riskset = "fixed",
#'                               time = relational.seq$eventID,
#'                               sender = as.character(relational.seq$sender),
#'                               receiver = as.character(relational.seq$target),
#'                               case_control=TRUE,
#'                               n_controls = 5)
#'
#'#Computing the sender-outdegree statistic for the above post-processing
#'#one-mode relational event sequence
#'post.processing <- dreamstats_degree(formation = "sender-outdegree",
#'                                   data = post.processing,
#'                                   halflife = 2)
#'
#'#Computing the inertia/repetition statistic for the above post-processing
#'#one-mode relational event sequence
#'post.processing <- dreamstats_repetition(data = post.processing,
#'                                         halflife = 2)
#'
#'#Fitting an ordinal timing relational event model to the above one-mode relational
#'#event sequence
#'rem <- estimate_rem(~ sender.outdegree + repetition,
#'                          data=post.processing)
#'summary(rem) #summary of the relational event model
#'
#'#the predicted event rates
#'rates <- predict(rem)
#'hist(rates)
#'
#'
#'#adding a new model matrix specification (this is purely an example) based
#'#upon a `data.frame` object
#'changed.data <- data.frame(sender.outdegree = rnorm(100),
#'                           repetition = sample(0:1, 100, TRUE))
#'new.rates <- predict(rem, newdata = changed.data)
#'hist(new.rates)
#'
#'#adding a new model matrix specification (this is purely an example) based
#'#upon the transformed `dream_sequence` object
#'new.sequence <- dreamstats_degree(formation = "sender-outdegree",
#'                                  data = post.processing,
#'                                  halflife = 4) #updating the halflife value
#'new.rates <- predict(rem, newdata = new.sequence)
#'hist(new.rates)
#'
#'@export
predict.dream_rem <- function(object, newdata = NULL, se.fit=FALSE,...){
  if(is.null(newdata)) stats <- object$statistics #using the original model matrix
  if(!is.null(newdata)){#making a new model matrix
    newdata1 <- NULL #presetting to extract the correct data matrix
    if(inherits(newdata,"data.frame")) newdata1 <- newdata
    if(inherits(newdata,"dream_sequence")) newdata1 <- as.data.frame(newdata)
    if(is.null(newdata1)) base::stop("The `newdata` argument must be a `data.frame` or `dream_sequence` object. Please update this argument and reuse the function!")
    #getting only the stats for the appropriate model and checking if all variables are here
    stats <- model.matrix(object$specification,data=newdata1) #building the model matrix with the new data
    vars <- names(object$parameters)
    stats<-stats[,vars]
  }
  #the predicted hazard rate for the specific event
  lambda <- exp(stats%*%object$parameters)
  if(se.fit){ #if the standard error should be returned
    #the variance of the linear predictor
    var.linpred <- diag(stats%*%vcov(object)%*%t(stats))
    #per the delta method, V(exp[y_pred]) = (exp(y_pred))^2 * V(y_pred) [the derivative of exp(y_pred) = exp(y_pred)]
    var.predlambda <- (lambda*lambda)*var.linpred
    list(fit = lambda, se.fit = sqrt(var.predlambda))
  }else{ #if not
    lambda
  }
}



#' Model residuals for `dream_rem` relational event model fits
#'
#' This function returns a set of model residuals based upon the realized
#' events in a relational event sequence (Butts 2008). The residuals included are: unit deviance,
#' residual deviance, standardized residual deviance, and for each covariate,
#' the Schoenfeld residual (please see Schoenfeld 1982).
#'
#' @param object An object of class "dream_rem".
#' @param ... Additional arguments for other methods.
#'
#'
#' @return A `data.frame` object that contains the following results:
#' \itemize{
#'   \item \code{time} - The timing of each realized event.
#'   \item \code{unit.deviance} - The unit deviance for each realized event (the contribution of each realized event to the log-likelihood).
#'   \item \code{residual.deviance} - The residual deviance for each realized event (the square of the contribution of each realized event to the log-likelihood).
#'   \item \code{standardized.deviance} - The standardized residual deviance for each realized event (the square of the contribution of each realized event to the log-likelihood).
#'   \item \code{schof.resid} - The Schoenfeld residual for each included covariate.
#'}
#'
#' @details
#'
#' The residuals, based upon an estimated relational event model, included are: unit deviance,
#' residual deviance, standardized residual deviance, and for each covariate,
#' the Schoenfeld residual (please see Schoenfeld 1982).
#'
#' @references
#' Butts, Carter T. 2008. "A Relational Event Framework for Social Action." *Sociological Methodology* 38(1): 155-200.
#'
#' Schoenfeld, David. 1982 "Partial Residuals for The Proportional Hazards Regression Model." *Biometrika* 69(1): 239-241.
#'
#' @examples
#'#Creating a psuedo one-mode relational event sequence with ordinal timing
#'relational.seq <- simulate_rem_seq(n_actors = 8,
#'                                   n_events = 50,
#'                                   inertia = TRUE,
#'                                   inertia_p = 0.10,
#'                                   sender_outdegree = TRUE,
#'                                   sender_outdegree_p = 0.05)
#'
#'#Creating a post-processing event sequence for the above relational sequence
#'post.processing <-  create_res(type = "one-mode",
#'                               ordinal = TRUE,
#'                               riskset = "fixed",
#'                               time = relational.seq$eventID,
#'                               sender = as.character(relational.seq$sender),
#'                               receiver = as.character(relational.seq$target),
#'                               case_control=TRUE,
#'                               n_controls = 5)
#'
#'#Computing the sender-outdegree statistic for the above post-processing
#'#one-mode relational event sequence
#'post.processing <- dreamstats_degree(formation = "sender-outdegree",
#'                                   data = post.processing,
#'                                   halflife = 2)
#'
#'#Computing the inertia/repetition statistic for the above post-processing
#'#one-mode relational event sequence
#'post.processing <- dreamstats_repetition(data = post.processing,
#'                                         halflife = 2)
#'
#'#Fitting an ordinal timing relational event model to the above one-mode relational
#'#event sequence
#'rem <- estimate_rem(~ sender.outdegree + repetition,
#'                          data=post.processing)
#'summary(rem) #summary of the relational event model
#'
#'#the model residuals
#'residuals(rem)
#'@export
residuals.dream_rem <- function(object,...){
  if(!inherits(object, "dream_rem")) base::stop("Error: The `object` argument must be a `dream_rem` object.")
  m1 <- as.data.frame(object$rem.data) #the modeling dataframe
  m1 <- m1[m1$sampled == 1,] #extracting only those sampled events
  model.mat <- object$statistics
  coefnames <- names(coef(object)) #the coefficent names
  colnames(model.mat) <- coefnames
  check.if.need.to.drop.no.null <- aggregate(x=1-m1$observed,by=list(m1$time),FUN=sum)
  if(any(check.if.need.to.drop.no.null$x == 0)){
    drop.these <- check.if.need.to.drop.no.null[,1][which(check.if.need.to.drop.no.null$x == 0)]
    m1 <- subset(m1,!(m1$time %in% drop.these) )
  }
  m1 <- as.data.frame(cbind(model.mat, observed =m1$observed, rates = predict(object)[,1],
              time = m1$time))
  #computing the pij (relevent for the ordinal timing log likelihodo)
  m1a <- split(m1, f=m1$time)
  m1a <- lapply(m1a,function(j){j$pij <- j[,"rates"] / sum(j[,"rates"]); j})
  #adding the intervevent times
  ordinal<-object$rem.data$ordinal
  if(!ordinal){
    for(i in 1:length(m1a)){m1a[[i]]$interevent <- object$rem.data$interevent_times[i]}
  }

  if(ordinal) m1a <- lapply(m1a,function(res){ dev <- -2*log(res$pij)
                                                      dev[which(res$observed==1)]})
  if(!ordinal){
      if(object$rem.data$t>0){
        m1a <- m1a[-length(m1a)] #removing the set of null events
      }

      m1a <- lapply(m1a,function(res){
        idx <-which(res$observed == 1)
        dev <- log(res$pij) - sum(res$rates*res$interevent)
        -2*dev[idx]
      })
  }
  unit.deviance <- unlist(m1a) #the deviance residuals for the ordinal timing log likelihood
  deviance.resid <- sqrt(unit.deviance)
  standardized.deviance <- (deviance.resid - mean(deviance.resid))/sd(deviance.resid)

  m1a <- split(m1, f=m1$time)
  m1a <- lapply(m1a,function(j){j$pij <- j[,"rates"] / sum(j[,"rates"]); j})
  if(!ordinal){
    if(object$rem.data$t>0){
      m1a <- m1a[-length(m1a)] #removing the set of null events
    }
  }
  #estimating the schoenfeld residual for the kth covariate
  m1a <- lapply(m1a,function(res){
    resid.sch <- rep(0,length(coefnames))
    for(j in 1:length(coefnames)){
      effect <- coefnames[j]
      xik <- res[res$observed==1, effect ]
      all <- res[,effect ]
      expected <- all%*%res[,"rates"]/sum(res[,"rates"])
      resid.sch[j] <- (xik - expected)
    }
    resid.sch
  })

  schof.resid <- do.call(rbind, m1a)
  colnames(schof.resid) <- paste0(coefnames, ".schoenfeld")

  resid.data <- data.frame(time= as.numeric(names(m1a)),
                           unit.deviance = unit.deviance,
                           residual.deviance=deviance.resid,
                           standardized.deviance=standardized.deviance,
                           schof.resid)
  rownames(resid.data) <- NULL #making the rownames null
  resid.data

}


#' Estimate the proportion of dyads predicted by `dream_rem` relational event model fits
#'
#' This function returns the proportion of dyads, senders, and targets predicted
#' by a `dream_rem` relational event model object. To compute the proportion for each
#' category, the function finds, for each realized event time, the dyad (based upon
#' the user's provided risk set) that is most likely to occur (i.e., the dyad
#' with the largest hazard rate) (Butts 2008). The predicted dyad is then compared to the
#' realized event dyad.
#'
#' @param object An object of class "dream_rem".
#' @param realized How should cases be handled when the largest hazard for a specific event
#' is the same for the realized event and one or more control events? TRUE indicates that
#' the method will default to the realized event. FALSE indicates that the function will
#' take a sample from the set of realized and control events that have the largest event
#' hazard.
#' @param rseed The random seed for reproducibility of results. When more than 1 dyad has
#' the maximum hazard rate, the predicted dyad is selected at random with equal probabilty
#' for each dyad. For example, if three dyads have the same hazard value, then the probability
#' of selecting each dyad will be 1/3.
#' @param ... Additional arguments for other methods.
#'
#' @return A `dream_gof` S3 object that contains the following results:
#' \itemize{
#'   \item \code{props} - A `data.frame` that contains the proportion of correctly predicted dyads, senders, and targets.
#'   \item \code{edgelist} - A `data.frame` object that contains the predicted and realized relational event sequence.
#'   \item \code{estimated_rem} - The `dream_rem` S3 object provided to the function.
#'   \item \code{match_realized} - The  `realized` argument provided to the function.
#'}
#'
#' @details
#' To compute the proportion for each category, the function finds, for each realized event time, the dyad (based upon
#' the user's provided risk set) that is most likely to occur (i.e., the dyad
#' with the largest hazard rate) (Butts 2008). The predicted dyad is then compared to the
#' realized event dyad.
#'
#' @references
#' Butts, Carter T. 2008. "A Relational Event Framework for Social Action." *Sociological Methodology* 38(1): 155-200.
#'#'
#' @examples
#'#Creating a psuedo one-mode relational event sequence with ordinal timing
#'relational.seq <- simulate_rem_seq(n_actors = 8,
#'                                   n_events = 50,
#'                                   inertia = TRUE,
#'                                   inertia_p = 0.10,
#'                                   sender_outdegree = TRUE,
#'                                   sender_outdegree_p = 0.05)
#'
#'#Creating a post-processing event sequence for the above relational sequence
#'post.processing <-  create_res(type = "one-mode",
#'                               ordinal = TRUE,
#'                                riskset = "fixed",
#'                               time = relational.seq$eventID,
#'                               sender = as.character(relational.seq$sender),
#'                               receiver = as.character(relational.seq$target),
#'                               case_control=TRUE,
#'                               n_controls = 5)
#'
#'#Computing the sender-outdegree statistic for the above post-processing
#'#one-mode relational event sequence
#'post.processing <- dreamstats_degree(formation = "sender-outdegree",
#'                                   data = post.processing,
#'                                   halflife = 2)
#'
#'#Computing the inertia/repetition statistic for the above post-processing
#'#one-mode relational event sequence
#'post.processing <- dreamstats_repetition(data = post.processing,
#'                                         halflife = 2)
#'
#'#Fitting an ordinal timing relational event model to the above one-mode relational
#'#event sequence
#'rem <- estimate_rem(~ sender.outdegree + repetition,
#'                          data=post.processing)
#'summary(rem) #summary of the relational event model
#'
#'#the model goodness of fit
#'gof <- gof_rem(rem,rseed=3718)
#'summary(gof) #printing the summary of the goodness of fit assessment
#'plot(gof) #plotting the goodness of fit assessment
#'gof$edgelist #check the predicted edgelist for the dyadic events
#'@export
gof_rem <- function(object,realized=TRUE,rseed=NULL,...){
  #the proportion of correctly specified dyads
  #a properly fit model should be able to predict what dyadic
  #events will occur (fail), that is, replicate the original
  #event sequence edgelist
  if(!inherits(object, "dream_rem")) base::stop("Error: The `object` argument must be a `dream_rem` object.")
  set.seed(rseed)
  edges <- as.data.frame(object$rem.data) #the modeling dataframe
  check.if.need.to.drop.no.null <- aggregate(x=1-edges$observed,by=list(edges$time),FUN=sum)
  if(any(check.if.need.to.drop.no.null$x == 0)){
    drop.these <- check.if.need.to.drop.no.null[,1][which(check.if.need.to.drop.no.null$x == 0)]
    edges <- subset(edges,!(edges$time %in% drop.these) )
  }
  rates <- data.frame(ratei = predict(object),
                      time = edges$time,
                      realized = edges$observed,
                      sender=edges$sender,
                      receiver=edges$receiver) #extracting the predicted hazard rates
  if(object$rem.data$t > 0)  rates <- subset(rates, rates$time != object$rem.data$t )
  correct.dyads <- edges[edges$observed==1,c("time","sender","receiver")]
  ratesj <- split(rates,f=rates$time) #splitting as a list to start our guessing

  guess <- lapply(ratesj, function(data,check=realized){
    max.value <- max(data$ratei)
    ids <- which(data$ratei == max.value) #checking the guess
    real <- which(data$realized!=0) #the real event
    if(check & (real %in% ids) ) ids <- real
    #if we have dyads with the same probabilty, sample 1 at random
    if(length(ids) > 1) ids <- sample(ids, 1) #with equal probability
    guess <- data.frame(obs.time=data$time[ids],rate=data$ratei[ids],
                        predicted.sender=data$sender[ids],predicted.receiver=data$receiver[ids])
  })

  guess <- do.call(rbind,guess) #combining the results
  rownames(guess) <- NULL #making them null
  #now to check how well the model did in terms of picking the "best" dyad
  guess$predicted.dyad <- paste0(guess$predicted.sender,"_+_",guess$predicted.receiver)
  guess$observed.dyad  <- paste0(edges$sender[edges$observed==1],"_+_",edges$receiver[edges$observed==1])
  guess$correct.dyad <- ifelse(guess$predicted.dyad==guess$observed.dyad,1,0)
  guess$correct.sender <- ifelse(guess$predicted.sender==edges$sender[edges$observed==1],1,0)
  guess$correct.receiver <- ifelse(guess$predicted.receiver==edges$receiver[edges$observed==1],1,0)

  N <- nrow(guess) #the number of realized events
  Y <- sum(guess$correct.dyad) #the number that we predicted correctly
  Y_sender <- sum(guess$predicted.sender == edges$sender[edges$observed==1])
  Y_receiver<- sum(guess$predicted.receiver == edges$receiver[edges$observed==1])
  prop.correct <- Y/N
  prop.sender <- Y_sender/N
  prop.receiver <- Y_receiver/N
  gof <- list(props = data.frame(Formation = c("Dyadic", "Sender", "Receiver"),
                                 Events = N,
                                 Correct = c(Y,Y_sender,Y_receiver),
                                 Prop=c(prop.correct,prop.sender,prop.receiver)),
              edgelist = guess,
              match_realized=realized,
              estimated_rem=object)
  class(gof) <- "dream_gof"
  gof

}


#' Print Method for `dream_gof` Model
#'
#' @param x An object of class `dream_gof`.
#' @param digits The number of digits to print after the decimal point.
#' @param ... Additional arguments (currently unused).
#' @return No return value. Prints out the main results of a 'dream_gof' object.
#' @export
print.dream_gof <- function(x,digits=5,...){
  cat("Relational Event Model Goodness-of-Fit\n")
  cat("\nFitted REM Call:\n")
  print(x$estimated_rem$call)
  cat("\n")
  cat("Goodness-of-Fit Information: \n")
  rownames(x$props) <- NULL
  x$props$Prop <- round(x$props$Prop,digits)
  print(x$props)
  y<-ifelse(x$match_realized,"Default to realized event", "Sample from the set")
  cat("\nMethod for handling multiple hazards that match the realized event hazard:\n")
  cat("  -",y)
  invisible(x)
}


#' Summary Method for `dream_gof` Objects
#'
#' Summarizes the results of a goodness-of-fit assessment for
#' estimated relational event models.
#'
#' @param object An object of class `dream_gof`.
#' @param digits The number of digits to print after the decimal point.
#' @param ... Additional arguments (currently unused).
#' @return A list of summary statistics for the goodness-of-fit assessment.
#' @export
#'
summary.dream_gof <- function(object,digits=5,...){
  new.table <- object$props
  rownames(new.table) <- NULL
  new.table$Prop <- round(new.table$Prop,digits)
  res <- list(props = new.table,
              call=object$estimated_rem$call,
              nevents=object$estimated_rem$n.events,
              nullevents=object$estimated_rem$null.events,
              matching = object$match_realized)
  class(res) <- "summary.dream_gof"
  return(res)
}

#' Print Method for `summary.dream_gof` Model
#'
#' @param x An object of class `summary.dream_gof`
#' @param digits The number of digits to print after the decimal point.
#' @param ... Additional arguments (currently unused).
#' @return No return value. Prints out the main results of a `summary.dream_gof` summary object.
#' @export
print.summary.dream_gof <- function(x,digits=5,...){
  cat("Relational Event Model Goodness-of-Fit\n")
  cat("\nFitted REM Call:\n")
  print(x$call)
  cat("\n")
  cat("Relational Event Sequence Information:\n")
  cat("The number of realized events =",x$nevents,"\nThe number of control events =", x$nullevents,"\n")
  cat("\nGoodness-of-Fit Information Based on Predicted Event Hazards: \n")
  print(x$props)
  y<-ifelse(x$matching,"Default to realized event\n", "Sample from the set\n")
  cat("\nMethod for handling multiple hazards if they match the realized event hazard:\n")
  cat("  -",y)
}




#' Plot method for `dream_rem` Relational Event Model Fits
#'
#' Plot the residuals of a `dream_rem` relational event model fit. Currently,
#' the `plot` function returns plots of either: (1) standardized deviance
#' residual values (y) by realized event times (x) or, for each
#' covariate, (2) the Schoenfeld residual values (y) by realized event times (x). By default,
#' the function plots the standardized deviance residual values.
#'
#' @param x An object of class "dream_rem".
#' @param type If "std.deviance", the returned plot is the absolute value standardized deviance
#' residual values (y) by realized event times (x). If "schoenfeld", then the
#' returned plots are for each covariates, (2) the Schoenfeld residual values (y) by realized event times (x).
#' @param ... Additional arguments for other methods.
#' @importFrom graphics abline
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom graphics mtext
#' @details
#' Generate plots for `dream_rem` relational event model fits to plot
#' model diagnostics.
#'
#' @examples
#'#Creating a psuedo one-mode relational event sequence with ordinal timing
#'relational.seq <- simulate_rem_seq(n_actors = 8,
#'                                   n_events = 50,
#'                                   inertia = TRUE,
#'                                   inertia_p = 0.10,
#'                                   sender_outdegree = TRUE,
#'                                   sender_outdegree_p = 0.05)
#'
#'#Creating a post-processing event sequence for the above relational sequence
#'post.processing <-  create_res(type = "one-mode",
#'                               ordinal = TRUE,
#'                               riskset = "fixed",
#'                               time = relational.seq$eventID,
#'                               sender = as.character(relational.seq$sender),
#'                               receiver = as.character(relational.seq$target),
#'                               case_control=TRUE,
#'                               n_controls = 5)
#'
#'#Computing the sender-outdegree statistic for the above post-processing
#'#one-mode relational event sequence
#'post.processing <- dreamstats_degree(formation = "sender-outdegree",
#'                                   data = post.processing,
#'                                   halflife = 2)
#'
#'#Computing the inertia/repetition statistic for the above post-processing
#'#one-mode relational event sequence
#'post.processing <- dreamstats_repetition(data = post.processing,
#'                                         halflife = 2)
#'
#'#Fitting an ordinal timing relational event model to the above one-mode relational
#'#event sequence
#'rem <- estimate_rem(~ sender.outdegree + repetition,
#'                          data=post.processing)
#'summary(rem) #summary of the relational event model
#'
#'#plotting the absolute value standardized deviance residuals for the estimated model
#'plot(rem, type="std.deviance")

#'@export
plot.dream_rem <- function(x,type=c("std.deviance", "schoenfeld"),...) {
  if(!inherits(x, "dream_rem")) base::stop("Error: The `x` argument must be a `dream_rem` object.")
  type <- match.arg(type, c("std.deviance", "schoenfeld")) #this provides a baseline option of std.deviance
  if(!(type %in% c("std.deviance", "schoenfeld"))) base::stop("Error: The `type` argument must be of either `std.deviance` or `schoenfeld`.")

  res <- residuals(x) #the model residuals based upon the x object
  #the standarized deviance residual plot
  if(type=="std.deviance"){
  plot(y=abs(res$standardized.deviance),
         x= res$time,
         ylab = "Abs. Standardized Deviance Residuals",
         xlab = "Realized Event Times",
         main = "Absolute Standardized Deviance Residuals",
         pch = 20,
         cex.lab = 1.25,
         cex.main = 1.5,
         ...) #the plot
  abline(h = 2, lty = "solid",col="red")
  n <- length(res$standardized.deviance)
  k <- sum(abs(res$standardized.deviance) > 2) #the number of sd resid greater than 2
  mtext(paste0("number (and percent) of residuals with |st. dev.| > 2: ", k, " (", round(k/n,4)*100 ,"%)"),
    side = 3, line = 0.25, adj = 1, font=3,cex = 1.1,... )
  }
  if(type == "schoenfeld"){
  schoenfeld.res <- grep("schoenfeld", names(res),value = TRUE) #the schoenfeld residuals
  k <- length(schoenfeld.res) #the number of schoenfeld residuals
  #the .schoenfeld residual plots
  op <- par(no.readonly = TRUE)
  #making a square for each covariates
  par(mfrow = c(ceiling(k / 2), 2), mar = c(4, 4, 2, 1))
  for(j in 1:k) {
    y <- res[,schoenfeld.res[j]]
    name <- paste0("Covariate: ", sub("\\.schoenfeld$", "", schoenfeld.res[j]))
    plot(y,pch = 16,main = name, xlab = "Realized Event Times",ylab = "Schoenfeld Residual")
    lines(lowess(y), col = "blue", lwd = 2)
    abline(h = 0, lty = 2)
  }
  #resetting the plot parameters to the user
  par(op)
  }
}



#' Plot method for `dream_gof` Relational Event Model Goodness-of-Fits
#'
#' Creates barplots of a `dream_gof` to plot the proportion of
#' event dyads, senders, and receivers that are correctly specified.
#'
#' @param x An object of class `dream_gof`.
#' @param ... Additional arguments for the barplot method.
#' @importFrom graphics text
#' @importFrom graphics barplot
#' @details
#' Generate plots for `dream_gof` relational event model goodness-of-fits.
#' @examples
#'#Creating a psuedo one-mode relational event sequence with ordinal timing
#'relational.seq <- simulate_rem_seq(n_actors = 8,
#'                                   n_events = 50,
#'                                   inertia = TRUE,
#'                                   inertia_p = 0.10,
#'                                   sender_outdegree = TRUE,
#'                                   sender_outdegree_p = 0.05)
#'
#'#Creating a post-processing event sequence for the above relational sequence
#'post.processing <-  create_res(type = "one-mode",
#'                               ordinal = TRUE,
#'                                riskset = "fixed",
#'                               time = relational.seq$eventID,
#'                               sender = as.character(relational.seq$sender),
#'                               receiver = as.character(relational.seq$target),
#'                               case_control=TRUE,
#'                               n_controls = 5)
#'
#'#Computing the sender-outdegree statistic for the above post-processing
#'#one-mode relational event sequence
#'post.processing <- dreamstats_degree(formation = "sender-outdegree",
#'                                   data = post.processing,
#'                                   halflife = 2)
#'
#'#Computing the inertia/repetition statistic for the above post-processing
#'#one-mode relational event sequence
#'post.processing <- dreamstats_repetition(data = post.processing,
#'                                         halflife = 2)
#'
#'#Fitting an ordinal timing relational event model to the above one-mode relational
#'#event sequence
#'rem <- estimate_rem(~ sender.outdegree + repetition,
#'                          data=post.processing)
#'summary(rem) #summary of the relational event model
#'
#'#the model goodness of fit
#'gof <- gof_rem(rem,rseed=3718)
#'plot(gof) #plotting the goodness of fit assessment
#'@export
plot.dream_gof <- function(x,...){
  if(!inherits(x, "dream_gof")) base::stop("Error: The `x` argument must be a `dream_gof` object.")
  x$props <- x$props[order(x$props$Formation),] #barplot orders the formation types in the plot
  gofplot <- barplot(Correct ~ Formation,
                     data=x$props,
                     ylim=c(0,max(x$props$Events)),
                     col=c("#B22222", "#377EB8", "#4DAF4A"),
                     border="black",
                     xlab="Event Formation Type",
                     ylab="Number of Predicted Realized Events",
                     main="Goodness-of-Fit for Fitted Relational Event Model",
                     family = "serif",
                     cex.lab = 1,
                     cex.main = 1.3,
                     cex.axis = 1.0,
                     ...)
  text(x = gofplot,
       y = x$props$Correct + 2,
       labels = round(x$props$Prop, 4),
       family = "serif",
       font=2)

}





