#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector tmegodist(NumericMatrix net,
                           NumericVector mem) {
  
  int nactors = net.nrow(); // the number of actors in the network
  int nlev2s = net.ncol(); // the number of level 2 nominations in the network
  NumericVector dist(nactors); // an empty vector to store the network
  for(int i = 0; i < nactors; i++){ // for all k level 1 actors in the network
    NumericVector newmem(nactors); // a new membership vector when membership ()
    for(int k = 0; k < nactors; k++){ // for all k level 1 actors in the network
      if(mem[k] == mem[i]){ // checking if the values are the same (i.e., same membership)
        newmem[k] = 1; }
    }
    NumericVector lev2ties = net(i,_); // the level 2 nominations
    int n2ties = sum(lev2ties); // the number of level 2 nominations
    if(n2ties > 0){ // if there are more than 2 ties
      for(int j = 0; j < nlev2s; j++){ // for all j level 2 actors in the network
        if(net(i,j)>0){ // if the network nomination is greater than 0
          NumericVector level2aties = net(_,j); // the level 2 nominations
          if(sum(level2aties) > 1){ // if the sum is greater than 1 (more than just the single actor)
            //(sum((net[, a] == 1)) - 1) from R (denominator)
            //(sum(new[net[, a] == 1]) - 1) from R  (numerator)
            double pia = (sum(level2aties*newmem)-1)/(sum(level2aties) - 1); 
            dist[i] += (1 - std::abs(1 - pia)); 
          } // ending the if sum statement for > 
        } // ending the if statement 
      } // ending the j loop
    }else{
      dist[i] = 0; // they are pendants, so move the value to 0
    }
  }
  return dist;
}
