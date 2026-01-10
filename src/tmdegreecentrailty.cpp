#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector tmdegcentraility(NumericMatrix tmnet) {
  double nactors = tmnet.nrow(); // the number of actors in the network
  double ntargets = tmnet.ncol(); // the number of actors in the network
  NumericVector degreestat(nactors); // a numeric vector to store the computed degree
  double ijk;
  // a loop to compute the degree statistic
  for(int i = 0; i < nactors; i++){ // for all actors
    for(int j = 0; j < nactors; j++){ // for all actors
      if(i != j){
            for(int k = 0; k < ntargets; k++){  // for all targets
              ijk =  tmnet(i,k) * tmnet(j,k);
              if(ijk > 0){
                degreestat[i] += 1; // updating the value based on the network nomination score
                break; // ending the search
              } // ending the if statement
            } // ending the k loop
      } // ending the if statement [i != j]
    } //ending the j loop
  }//ending the i loop
  return degreestat;
}


