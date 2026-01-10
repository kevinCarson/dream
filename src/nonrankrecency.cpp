#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector computerecencynorank(NumericVector time, 
                              NumericVector sampledevent, 
                              NumericVector controlevent,
                              std::vector<std::string> dyad_id,
                              std::vector<std::string> sender,
                              std::vector<std::string> target,
                              bool i_neighborhood,
                              bool raw_diff,
                              double nopastEvents,
                              std::string appender) {
  
  
  std::unordered_map<std::string, double> dyad_times; // an unordered map to store the most recent event times
  int nevents = time.size(); // the number of events in the relational event sequence
  int sampleevent = sum(sampledevent); //the number of sampled events for which we need to compute the statistic for
  int cursampled = 0; // an empty starting point
  NumericVector statvector(sampleevent); //an empty starting vector to store the inertia scores
  std::string curdyad; // an empty container
  double sumweights = 0;
  //looping through all events to find the past event times
  for (int i = 0; i < nevents; i++) {
    
    if(i_neighborhood){ // if the computation should be based on the sender neighborhood
      curdyad = sender[i] + appender + target[i]; //appending the dyad id
    }else{// if the computation should be based on the target neighborhood
      curdyad = target[i] + appender + sender[i]; //appending the dyad id
    }
    
    // we dont need to compute the score for all events just those that are sampled!
    if(sampledevent[i] == 1){ // if the event is sampled, we need to compute the statistics
      sumweights = 0; // resetting the container to zero
      // based upon how risk sets are created, the same dyad SHOULD NOT APPEAR TWICE FOR A SINGLE Ti
      
      double mostrecent = dyad_times[curdyad];  // this extracts all of the past event times
      if(mostrecent == 0){ // if there are no past events, then it is undefined
        sumweights = nopastEvents;  // set the value to the nopastevents object
      }else{
        if(raw_diff){ // if the raw difference should be returned
          sumweights = (time[i] - mostrecent); 
        }else{ // if the inverse plus 1 shoudl be returned (check the manual)
          sumweights = 1/(time[i] - mostrecent + 1); 
        }
      }
      statvector[cursampled] = sumweights; // storing the values for the ith statistic
      cursampled += 1; // update this position indicator for storing the computed values
    }
    
    if(controlevent[i] != 1){ // if it is not a fake event, append the current time to the vector
      std::string dyadid = sender[i] + appender + target[i]; //appending the dyad id
      double& mostrecent = dyad_times[curdyad]; // the current dyad
      mostrecent = time[i]; // reinserting the most recent time value
    }
    
  }
  return(statvector);  //return the entire vector of values (we will use R to extract the appropriate values from the vector)
  
}

