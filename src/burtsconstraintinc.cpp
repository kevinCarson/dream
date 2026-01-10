#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double burtconstraint(NumericMatrix net,
                        double nactors) {
  
  NumericMatrix pij(nactors, nactors); // Creates a matrix filled with zeros
  NumericMatrix cij(nactors, nactors); // Creates a matrix filled with zeros
  NumericMatrix pij1(nactors, nactors); // Creates a matrix filled with zeros
  double zijsum; // presetting the value to zero
  double maxzjk = 0; // presetting the value to zero
  
  for(int i = 0; i < nactors; i++){ // for all actors in the ego network
    zijsum = 0; // resetting the values 
    for(int j = 0; j < nactors; j++){ // for all alters in the network
      pij(i,j) = net(i,j) + net(j,i); // from burts formula
      zijsum += pij(i,j); //adding the value to the sum for all actors
      if(maxzjk < pij(i,j)){maxzjk = pij(i,j);} // if the value is greater than the current maximum update it!
    } // ending the j loop
    for(int j = 0; j < nactors; j++){// for all alters in the network
      pij1(i,j) = pij(i,j)/zijsum; // from burts formula
    } // ending the j loop
  } // ending the i loop
  
  //for all j alters'
  double constraint = 0; 
  double innerqsum = 0;
  for(int j = 1; j < nactors; j++ ){// for all alters in the network
    innerqsum = 0; // resetting the value
    for(int q = 1; q < nactors; q++ ){// for all alters in the network
      if(j != q){
        innerqsum += pij1(0,q)*pij1(q,j); // from burts formula
      }
    }
    constraint += (pij1(0,j) + innerqsum)*(pij1(0,j) + innerqsum); // from burts formula
  }
  
  return constraint;
}


