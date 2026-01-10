#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector prefattachrelspanrem(NumericVector sampledevent, 
                                             NumericVector controlevent,
                                             NumericVector time, 
                                             std::vector<std::string> sender,
                                             std::vector<std::string> target,
                                             double reltimespan){
  //
  // std::vector<double> sampledevent
  // std::vector<double> controlevent
  // std::vector<std::string> sender
  // std::vector<std::string> target
  // bool dependency
  // double reltimespan
  //
  
  std::unordered_map<std::string, std::vector<double>> sender_times; //an empty list (container) to store the sender id and the past event counts
  std::unordered_map<std::string, std::vector<double>> target_times; //an empty list (container) to store the target ids and the past event counts
  int nevents = time.size(); // the number of events in the relational event sequence
  int sampleevent = sum(sampledevent); //the number of sampled events for which we need to compute the statistic for
  int cursampled = 0; // an empty starting point
  NumericVector statvector(sampleevent); //an empty starting vector to store the inertia scores
  double sumweights=0; 
  double denom=0; 
  std::string pasttarget; // to discount in the internal computation
  std::string pastsender; // to discount in the internal computation
  NumericVector pasteventtimes; 
  int pastindex = 0; 
  double senderout; 
  double targetout; 
  double mostrecenttime; 
  
  //looping through all events to find the past event times
  for (int i = 0; i < nevents; i++) {
    
    std::string curreceiver = (target[i]); //the current event target
    // we dont need to compute the score for all events just those that are sampled!
    if(sampledevent[i] == 1){ // if the event is sampled, we need to compute the statistics
      senderout = 0;
      targetout = 0;
      
      std::vector<double> past_times_sender = sender_times[curreceiver];  // this extracts all of the past event times (outdegree for target i)
      std::vector<double> past_times_target = target_times[curreceiver];  // this extracts all of the past event times (indegree for target i)

      for(int z = 0; z < past_times_sender.size(); z++){
        if(((time[i] - past_times_sender[z]) <= reltimespan) &&
           ((time[i] - past_times_sender[z]) != 0)){
          senderout += 1; //updating the count by 1
          }
      }
      
      for(int z = 0; z < past_times_target.size(); z++){
        if(((time[i] - past_times_target[z]) <= reltimespan) &&
           ((time[i] - past_times_target[z]) != 0)){
          targetout += 1; //updating the count by 1
        }
      }
      
      double ntiestarget = targetout + senderout; 
      /// we now need to compute denom with the past relevancy that matters
      denom = 0; // resetting the value to 2
      for(int j = pastindex; j < pasteventtimes.length(); j++){ // it will be bound with 2
        // this should turn to the minimum time that is good for the event count!
        if((time[i] - pasteventtimes[j]) <= reltimespan){ 
          pastindex = j; 
          break; 
        }
      }
      if(time[i] == mostrecenttime){
        denom = (pasteventtimes.length() -1) - pastindex; 
      }else{
        denom = pasteventtimes.length() - pastindex; 
      }
      sumweights = ntiestarget/(2*denom); 
      statvector[cursampled] = sumweights; // storing the values for the ith statistic
      cursampled += 1; // update this position indicator for storing the computed values
    }
    
    if(controlevent[i] != 1){ // if it is not a fake event, append the current time to the vector
      std::vector<double>& sender_timesi = sender_times[sender[i]]; //extracting the count for the actor (event sender)
      std::vector<double>& target_timesi = target_times[target[i]]; //extracting the count for the actor (event target)
      sender_timesi.push_back(time[i]); //updating the score by one 
      target_timesi.push_back(time[i]); //updating the score by one 
      pasttarget = target[i];  // the most recent target
      pastsender = sender[i];  // the most recent sender
      pasteventtimes.push_back(time[i]); // the past event times vector
      mostrecenttime = time[i]; // the most recent event time
    }
  }
  
  return statvector;
}


