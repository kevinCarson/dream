#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double remloglike(NumericVector beta, List event_stats) {
// A function to compute the log likelihood REM function

  double loglike = 0.0; // initializing the log likelihood to 0
  int K = beta.size(); // the number of parameters to estimate
  int N_events = event_stats.size(); // the numeric of variables

  for (int i = 0; i < N_events; i++) { // for all events clusters in the sequence

    NumericMatrix statsi = event_stats[i]; // the stats for event cluster i

    int nEvents =  statsi.nrow(); // the number of events in the cluster

    NumericVector labmdaj(nEvents,0.0); // the vector for lambdas

    for (int j = 0; j < K; j++) { // for all events clusters in the sequence
      labmdaj = labmdaj + beta[j]*statsi(_,j); // making all of the statistics
    }

    NumericVector dummy = statsi( _ , K+1); // the event dummy variable
    int Dummy = which_max(dummy); // assuming only one real event per cluster
    double lambdai = labmdaj[Dummy];//the true lambda
    NumericVector lambdaRS = labmdaj;
    double maxLambda = max(labmdaj);
    double lambdaRSum = maxLambda + log(sum(exp(labmdaj-maxLambda)));
    loglike = loglike + (lambdai-lambdaRSum);
  }
  return loglike; // returning the loglikelihood
}

// [[Rcpp::export]]
double remloglikeOPTIM(NumericVector beta, List event_stats) {
// A function to compute the log likelihood REM function for the optim function

  double loglike = 0.0; // initializing the log likelihood to 0
  int K = beta.size(); // the number of parameters to estimate
  int N_events = event_stats.size(); // the numeric of variables

  for (int i = 0; i < N_events; i++) { // for all events clusters in the sequence

    NumericMatrix statsi = event_stats[i]; // the stats for event cluster i

    int nEvents =  statsi.nrow(); // the number of events in the cluster

    NumericVector labmdaj(nEvents,0.0); // the vector for lambdas

    for (int j = 0; j < K; j++) { // for all events clusters in the sequence
      labmdaj = labmdaj + beta[j]*statsi(_,j); // making all of the statistics
    }

    NumericVector dummy = statsi( _ , K+1); // the event dummy variable
    int Dummy = which_max(dummy); // assuming only one real event per cluster
    double lambdai = labmdaj[Dummy];//the true lambda
    NumericVector lambdaRS = labmdaj;
    double maxLambda = max(labmdaj);
    double lambdaRSum = maxLambda + log(sum(exp(labmdaj-maxLambda)));
    loglike = loglike + (lambdai-lambdaRSum);
  }
  loglike = -loglike; // making the log liklelihood negative for the optim function!
  return loglike; // returning the loglikelihood
}

// [[Rcpp::export]]
NumericVector remGradient(NumericVector beta, List event_stats) {
// A function to compute the gradient vector for a REM

  int K = beta.size(); // the number of parameters to estimate
  NumericVector gradient(K); // an empty k x 1 vector
  int N_events = event_stats.size(); // the numeric of variables

  for (int i = 0; i < N_events; i++) { // for all events clusters in the sequence

    NumericMatrix statsi = event_stats[i]; // the stats for event cluster i

    int nEvents =  statsi.nrow(); // the number of events in the cluster

    NumericVector labmdaj(nEvents,0.0); // the vector for lambdas

    for (int j = 0; j < K; j++) { // for all events clusters in the sequence
      labmdaj = labmdaj + beta[j]*statsi(_,j); // making all of the statistics
    }
    NumericVector pij = exp(labmdaj)/sum(exp(labmdaj));
    NumericVector dummy = statsi( _ , K+1); // the event dummy variable
    int Dummy = which_max(dummy); // assuming only one real event per cluster
    for (int j = 0; j < K; j++) { // for all requested effects
      gradient[j] = gradient[j] + (statsi(Dummy,j)-sum(pij*statsi(_,j)));
    }
  }

  return gradient; // returning the loglikelihood


}

// [[Rcpp::export]]
NumericMatrix remHessian(NumericVector beta, List event_stats) {
 // A function to compute the hessian matrix for a REM

  int K = beta.size(); // the number of parameters to estimate
  NumericMatrix hessian(K,K); // an empty k x k matrix
  int N_events = event_stats.size(); // the numeric of variables

  for (int i = 0; i < N_events; i++) { // for all events clusters in the sequence

    NumericMatrix statsi = event_stats[i]; // the stats for event cluster i

    int nEvents =  statsi.nrow(); // the number of events in the cluster

    NumericVector labmdaj(nEvents,0.0); // the vector for lambdas

    for (int j = 0; j < K; j++) { // for all events clusters in the sequence
      labmdaj = labmdaj + beta[j]*statsi(_,j); // making all of the statistics
    }
    NumericVector pij = exp(labmdaj)/sum(exp(labmdaj));
    NumericMatrix crossprodmat(K,K); // an empty k x k matrix
    NumericMatrix startingmat(K,K); // an empty k x k matrix

    for (int j = 0; j < K; j++) { // for all requested effects (jth)
      // a double loop for all the cross products
      for (int k = 0; k < K; k++) { // for all requested effects (kth)
        crossprodmat(j,k)=sum(pij*statsi(_,j))*sum(pij*statsi(_,k));
        startingmat(j,k)=sum(pij*statsi(_,j)*statsi(_,k));
        //hessian = hessian +(startingmat -crossprodmat); //updating the hessian matrix
        hessian(j,k)= hessian(j,k)+startingmat(j,k)-crossprodmat(j,k);
    }
  }
  }

  for (int j = 0; j < K; j++) { // for all requested effects (jth)
    // a double loop for all the cross products
    for (int k = 0; k < K; k++) { // for all requested effects (kth)
      hessian(j,k)= -hessian(j,k);
    }
  }
  return hessian;

}

