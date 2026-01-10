#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector persistencerem(NumericVector time, 
                                 NumericVector sampledevent, 
                                 NumericVector controlevent,
                                 std::vector<std::string> dyad_id, 
                                 std::vector<std::string> actor,
                                 bool timedependency,
                                 double cuttime,
                                 double nopastEvents) { 
  

  
  
  int nevents = time.size(); // the number of events in the relational event sequence
  int sampleevent = sum(sampledevent); //the number of sampled events for which we need to compute the statistic for
  double cursampled = 0; // an empty starting point
  NumericVector statvector(sampleevent); //an empty starting vector to store the inertia scores
  
  // we need to create a dyadic list that stores for each event that past event times
  // and once they go to zero, we drop them from the set
  std::unordered_map<std::string, std::vector<double>> dyad_events; //an empty list (container) to store the dyad ids and the past event times
  std::unordered_map<std::string, std::vector<double>> actor_events; //an empty list (container) to store the dyad ids and the past event times

  
  //looping through all events to find the past event times
  for (int i = 0; i < nevents; i++) {
    std::string curdyad = dyad_id[i]; //the current event dyad
    std::string curactor = actor[i]; //the current event dyad
    
    std::vector<double>& past_times = dyad_events[curdyad];  // this extracts all of the past event times
    std::vector<double>& past_times_actor = actor_events[curactor];  // this extracts all of the past event times

    // we dont need to compute the score for all events just those that are sampled!
    if(sampledevent[i] == 1){ // if the event is sampled, we need to compute the statistics
      
      double sumweights = 0; // the event times
      
      double ni = past_times_actor.size(); // predefined integers
      double ndyad = 0; // predefined integers
 
      if(ni > 0){ //if there are past weight times then update!
        
        
        // need to check the event times
        
        ni = 0; //rewriting it
        ndyad = 0; //rewriting it
        
        
        for(int k = 0; k < past_times_actor.size(); k++){
          
          if(past_times_actor[k] < time[i]){
            
            if(timedependency){
              
              if(past_times_actor[k] > (time[i] - cuttime)){
                ni += 1; //updating the value by one if the value is sooner than the cutoff 
              }
              
            }else{
              ni += 1; //updating the value by 1 (i.e., the event time is good)
            }
            
          }
          
        }
        
        if(ni == 0){ // if after the update the times are no longer good
          
          sumweights = nopastEvents;
          
        }else{
          
          for(int k = 0; k < past_times.size(); k++){
            if(past_times[k] < time[i]){
              if(timedependency){
                if(past_times[k] > (time[i] - cuttime)){
                  ndyad += 1; //updating the value by one if the value is sooner than the cutoff 
                }
              }else{
                ndyad += 1; //updating the value by 1 (i.e., the event time is good)
              }
            }
          }
          
          sumweights = ndyad/ni; // from the persistence formula (the number of ties to a specific person) / total number of past ties
          
          
        }

      }else{
        sumweights = nopastEvents;
      }
      
      statvector[cursampled] = sumweights; // storing the values for the ith statistic
      cursampled += 1; // update this position indicator for storing the computed values
    }
    
    if(controlevent[i] == 1){ // if it is not a fake event, append the current time to the vector
      
      past_times.push_back(time[i]); // updating the network of past events at the dyad level
      past_times_actor.push_back(time[i]); // updating the network of past events at the actor level
      
    }
    
  }
  return(statvector);  //return the entire vector of values (we will use R to extract the appropriate values from the vector)
}
