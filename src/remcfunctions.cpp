// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <numbers>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat remhess_interval(arma::vec beta,
                           Rcpp::List event_stats,
                           arma::vec interevent,
                           bool tknown) {

  int N_events; //initializing the variable
  if(tknown){
    N_events = event_stats.size() - 1; // the number of event clusters
  }else{
    N_events = event_stats.size(); // the number of event clusters
  }
  int K = beta.size();
  arma::vec labmdaj;
  arma::mat statsi;
  arma::mat hessian(K,K); // an empty vector to store the gradient to be updated at each event
  for (int i = 0; i < N_events; i++) { // for all events clusters in the sequence
    arma::mat statsi = event_stats[i]; // the stats for event cluster i
    labmdaj = statsi * beta;
    labmdaj = exp(labmdaj)*interevent[i];
    hessian += -statsi.t() * diagmat(labmdaj) * statsi;
  }
  if(tknown){
    // adding the right censoring into the dataframe (if the timing is unknown, the contribution will be 0, log(exp(0)) = 0)
    arma::mat statsend = event_stats[event_stats.size() - 1]; // the stats for event cluster i
    labmdaj = statsend * beta;
    labmdaj = exp(labmdaj)*interevent[interevent.size()-1];
    hessian  += -statsend.t() * diagmat(labmdaj) * statsend;
  }
  return hessian;
}



// [[Rcpp::export]]
arma::vec remgrad_interval(arma::vec beta,
                                  Rcpp::List event_stats,
                                  arma::vec interevent,
                                  arma::vec obsindex,
                                  bool tknown) {
  int N_events; //initializing the variable
  if(tknown){
    N_events = event_stats.size() - 1; // the number of event clusters
  }else{
    N_events = event_stats.size(); // the number of event clusters
  }
  int K = beta.size(); // the number of parameters to estimate
  arma::rowvec gradient(K); // an empty vector to store the gradient to be updated at each event
  arma::rowvec realstats; // an empty vector to store the gradient to be updated at each event
  arma::vec labmdaj;
  arma::mat statsi;
  //arma::mat updatemat;
  for (int i = 0; i < N_events; i++) { // for all events clusters in the sequence
    arma::mat statsi = event_stats[i]; // the stats for event cluster i
    labmdaj = statsi * beta;
    realstats = statsi.row(obsindex[i]);
    labmdaj = exp(labmdaj)*interevent[i];
    statsi.each_col() %= labmdaj;
    gradient += realstats - sum(statsi, 0);  // adding the column sums
  }
  if(tknown){
  // adding the right censoring into the dataframe (if the timing is unknown, the contribution will be 0, log(exp(0)) = 0)
  arma::mat statsend = event_stats[event_stats.size() - 1]; // the stats for event cluster i
  labmdaj = statsend * beta;
  labmdaj = exp(labmdaj)*interevent[interevent.size()-1];
  statsend.each_col() %= labmdaj;
  gradient += -sum(statsend, 0);
  }
  return gradient.t(); // returning the loglikelihood
}


// [[Rcpp::export]]
double remloglik_interval(arma::vec beta,
                          Rcpp::List event_stats,
                          arma::vec interevent,
                          arma::vec obsindex,
                          bool tknown) {

  double loglike = 0.0; // presetting the value
  int N_events; //initializing the variable
  if(tknown){
    N_events = event_stats.size() - 1; // the number of event clusters
  }else{
    N_events = event_stats.size(); // the number of event clusters
  }
  arma::vec labmdaj;
  arma::vec labmdak;
  arma::mat statsi;
  double outersum;
  double innersum;
  //arma::mat updatemat;
  for (int i = 0; i < N_events; i++) { // for all events clusters in the sequence
    arma::mat statsi = event_stats[i]; // the stats for event cluster i
    labmdak = statsi * beta;
    labmdaj = exp(labmdak)*interevent[i];
    outersum = -sum(labmdaj);  // the summation over all the event hazards
    innersum = labmdak[obsindex[i]]; // the original value
    loglike += innersum + outersum; //updating the value here!
  }
  if(tknown){
  // adding the right censoring into the dataframe (if the timing is unknown, the contribution will be 0, log(exp(0)) = 0)
  arma::mat statsend = event_stats[event_stats.size() - 1]; // the stats for event cluster i
  labmdaj = statsend * beta;
  labmdaj = exp(labmdaj)*interevent[interevent.size()-1];
  loglike += -sum(labmdaj);
  }
  return -loglike; // returning the loglikelihood
}



// [[Rcpp::export]]
double remloglik_ordinal(arma::vec beta,
                         Rcpp::List event_stats,
                         arma::vec obsindex) {

  double loglike = 0.0; // presetting the value
  int N_events = event_stats.size(); // the number of event clusters
  double labmdai;
  arma::vec labmdaj;
  arma::mat statsi;
  double maxLambda;
  double lambdaRSum;
  //arma::mat updatemat;
  for (int i = 0; i < N_events; i++) { // for all events clusters in the sequence
    arma::mat statsi = event_stats[i]; // the stats for event cluster i
    labmdaj = statsi * beta;
    labmdai = labmdaj[obsindex[i]];
    //labmdaj = exp(labmdaj);
    maxLambda = max(labmdaj);
    lambdaRSum = maxLambda + log(sum(exp(labmdaj-maxLambda)));
    loglike += labmdai - lambdaRSum; //updating the value here!
  }
  return -loglike; // returning the loglikelihood
}



// [[Rcpp::export]]
arma::vec remgrad_ordinal(arma::vec beta,
                          Rcpp::List event_stats,
                          arma::vec obsindex) {

  int K = beta.size(); // the number of parameters to estimate
  int N_events = event_stats.size(); // the number of event clusters
  arma::rowvec gradient(K); // an empty vector to store the gradient to be updated at each event
  arma::rowvec realstats; // an empty vector to store the gradient to be updated at each event
  arma::vec labmdaj;
  arma::mat statsi;
  arma::vec pij;

  for (int i = 0; i < N_events; i++) { // for all events clusters in the sequence

    arma::mat statsi = event_stats[i]; // the stats for event cluster i
    labmdaj = statsi * beta;
    pij = exp(labmdaj)/sum(exp(labmdaj)); // the probabilty of each event
    realstats = statsi.row(obsindex[i]);
    statsi.each_col() %= pij;
    gradient += realstats - sum(statsi, 0);  // adding the column sums

  }
  return gradient.t(); // returning the computed gradient
}


// [[Rcpp::export]]
arma::mat remhess_ordinal(arma::vec beta,
                          Rcpp::List event_stats) {
  int N_events = event_stats.size(); // the number of event clusters
  int K = beta.size();
  arma::vec labmdaj;
  arma::mat statsi;
  arma::vec pij;
  arma::mat hessian(K,K); // an empty vector to store the gradient to be updated at each event
  arma::mat crossprodmat(K,K); // an empty k x k matrix
  arma::mat startingmat(K,K); // an empty k x k matrix
  for (int i = 0; i < N_events; i++) { // for all events clusters in the sequence
    arma::mat statsi = event_stats[i]; // the stats for event cluster i
    labmdaj = statsi * beta;
    pij = exp(labmdaj)/sum(exp(labmdaj)); // the probabilty of each event
    for (int j = 0; j < K; j++) { // for all requested effects (jth)
      // a double loop for all the cross products
      for (int k = 0; k < K; k++) { // for all requested effects (kth)
        crossprodmat(j,k)=sum(pij % statsi.col(j))*sum(pij %statsi.col(k));
        //startingmat(j,k)=sum(pij*statsi.col(j)*statsi.col(k));
        startingmat(j,k)=sum(pij % statsi.col(j) % statsi.col(k));
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

