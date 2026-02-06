#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector computeremprefattach(NumericVector sampledevent,
                                      NumericVector controlevent,
                                      NumericVector time,
                                      std::vector<std::string> sender,
                                      std::vector<std::string> target){
  //
  // std::vector<double> sampledevent
  // std::vector<double> controlevent
  // std::vector<std::string> sender
  // std::vector<std::string> target
  // bool dependency
  // double reltimespan
  //

  double ncount = 0;// for the denominator
  double pasteventtime = 0; // an empty object to store the last time update for the count
  std::unordered_map<std::string, double> sender_counts; //an empty list (container) to store the sender id and the past event counts
  std::unordered_map<std::string, double> target_counts; //an empty list (container) to store the target ids and the past event counts
  int nevents = time.size(); // the number of events in the relational event sequence
  int sampleevent = sum(sampledevent); //the number of sampled events for which we need to compute the statistic for
  double cursampled = 0; // an empty starting point
  NumericVector statvector(sampleevent); //an empty starting vector to store the inertia scores
  double sumweights=0;
  double denom=0;
  std::string pasttarget; // to discount in the internal computation
  std::string pastsender; // to discount in the internal computation

  //looping through all events to find the past event times
  for (int i = 0; i < nevents; i++) {
    std::string curreceiver = (target[i]); //the current event target

      // we dont need to compute the score for all events just those that are sampled!
      if(sampledevent[i] == 1){ // if the event is sampled, we need to compute the statistics
        double past_times_sender = sender_counts[curreceiver];  // this extracts all of the past event times (outdegree for target i)
        double past_times_target = target_counts[curreceiver];  // this extracts all of the past event times (indegree for target i)
        if(time[i] == pasteventtime){
          denom = (ncount-1) ;
          if(pasttarget == curreceiver){
            past_times_target = past_times_target - 1;
          }
          if(pastsender == curreceiver){
            past_times_sender = past_times_sender - 1;
          }
        }else{
          denom = ncount;
        }
        double ntiestarget = past_times_sender + past_times_target;
        sumweights = ntiestarget/(2*denom);
        statvector[cursampled] = sumweights; // storing the values for the ith statistic
        cursampled += 1; // update this position indicator for storing the computed values
      }

      if(controlevent[i] != 1){ // if it is not a fake event, append the current time to the vector
        double& sender_countsi = sender_counts[sender[i]]; //extracting the count for the actor (event sender)
        double& target_countsi = target_counts[target[i]]; //extracting the count for the actor (event target)
        sender_countsi += 1; //updating the score by one
        target_countsi += 1; //updating the score by one
        pasteventtime = time[i]; //updating the current event time
        ncount += 1; //updating the value by 1 for the number of past events
        pasttarget = target[i];  // the most recent target
        pastsender = sender[i];  // the most recent sender
      }
  }

  return statvector;
}


