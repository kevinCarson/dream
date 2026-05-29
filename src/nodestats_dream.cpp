#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector dream_nodestats(NumericVector time,
                              NumericVector time_start,
                              NumericVector time_end,
                              CharacterVector observed_actor,
                              CharacterVector actor_id,
                              NumericVector value){
  double n = time.length(); // the number of elements in the sampled vector
  double nt = actor_id.length(); // the number of observations with information
  NumericVector weights(n); // an empty vector to store the information
  String curactor;
  double starti;
  double stopi;
  double valuei;
  for(int i = 0; i < nt; i++){
    curactor = actor_id[i];
    starti = time_start[i];
    stopi = time_end[i];
    valuei = value[i];
    for(int j = 0; j < n; j++){
      if(curactor == observed_actor[j] &&
         time[j] >= starti &&
         stopi >= time[j]){
        weights[j] = valuei; // inputting the value for the node-level stat
      }
    }
  }
  return weights; // return the weights to the user
}










