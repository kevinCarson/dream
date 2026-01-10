#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector computerecencyrank(  NumericVector time, 
                                   NumericVector sampledevent, 
                                   NumericVector controlevent,
                                   std::vector<std::string> sender,
                                   std::vector<std::string> target,
                                   bool i_neighborhood,
                                   double nopastEvents,
                                   std::string appender) {
  
  

  // for rank, we need to do this:
  // 1. make a unordered_map<std::string, std::vector<double>> past targets
  // 2. update this after each event
  // 3. then check if that user is even in the list std::find(), if find is each to .length() then no past events
  // 4. if they are in the group (find the elemtns in the list [see how it sorts])
  // 5. take the inverse of the rank!  
  
  
  std::unordered_map<std::string, std::vector<double>> dyad_times; // an unordered map to store the most recent event times
  std::unordered_map<std::string, std::vector<double>> sender_times; // an unordered map to store the most recent event times
  int nevents = time.size(); // the number of events in the relational event sequence
  int sampleevent = sum(sampledevent); //the number of sampled events for which we need to compute the statistic for
  int cursampled = 0; // an empty starting point
  NumericVector statvector(sampleevent); //an empty starting vector to store the inertia scores
  double sumweights = 0;
  std::string cursender;  // the current event sender 
  std::string dyadi; // the dyadi value
  std::string curi; // the current i sender based on the neighborhood
  std::string curj; // the current j sender based on the neighborhood
  double mostrecenttime; // the most recent time of the last event (updating the network to do a check for contorl events at time i)
  double recentdyadtime; // the most recent time for the dyad i 
  double dyaddist; 

  //looping through all events to find the past event times
  for (int i = 0; i < nevents; i++) {
    
    // we dont need to compute the score for all events just those that are sampled!
    if(sampledevent[i] == 1){ // if the event is sampled, we need to compute the statistics
      sumweights = 0; // resetting the container to zero
      if(i_neighborhood){ // if the computation should be based on the sender neighborhood
        dyadi = sender[i] + appender + target[i]; // the dyad from the current 
        curi = sender[i]; // the current sending actor
      }else{// if the computation should be based on the target neighborhood
        dyadi = target[i] + appender + sender[i]; // the dyad from the current 
        curi = target[i];// the current sending actor
      }
      
      std::vector<double> curdyadtimes = dyad_times[dyadi]; // the current dyadic times
      std::vector<double> curactitimes = sender_times[curi]; // the current dyadic times
      
      if(curdyadtimes.size() == 0){ // if there are no past events in the dyad
        sumweights = nopastEvents; // since there are no past events within the dyad -> undefined receny rank
      }else{
        recentdyadtime = curdyadtimes[curdyadtimes.size() - 1]; // the most recent dyadic time for the user
        
        if((recentdyadtime == mostrecenttime) && curdyadtimes.size() == 1){
          sumweights = nopastEvents; // since there are no past events within the dyad -> undefined receny rank
        }else{
        
            if((recentdyadtime == mostrecenttime) && curdyadtimes.size() > 1){
              recentdyadtime = curdyadtimes[curdyadtimes.size() - 2]; // the most recent dyadic time for the user
            }else{
              recentdyadtime = curdyadtimes[curdyadtimes.size() - 1]; // the most recent dyadic time for the user
            }
            
            // searching for the recency rank 
            for(int z = 0; z < curactitimes.size(); z++){
              if(curactitimes[curactitimes.size() - z - 1] == recentdyadtime){
                dyaddist = z + 1; 
                break; 
              }
            }
            
            // time to do some back searching for the most recent event
            double checkend = curactitimes[curactitimes.size() - 1]; 
            double dist = (dyaddist); 
            if(checkend == mostrecenttime){
              dist = dist - 1; // the total number of events (minus 1)
            }
            sumweights = 1/(dist); // the reciprocal of the number
        } 
       }
      
      statvector[cursampled] = sumweights; // storing the values for the ith statistic
      cursampled += 1; // update this position indicator for storing the computed values
      
    }
      
    
    if(controlevent[i] != 1){ // if it is not a fake event, append the current time to the vector
      cursender = sender[i]; // the id for the sender
      dyadi = sender[i] + appender + target[i]; // creating the dyad id
      std::vector<double>& pastdyadytimes = dyad_times[dyadi]; // the past targets for the current sender
      std::vector<double>& pasttime   = sender_times[cursender]; // the past targets for the current sender
      pastdyadytimes.push_back(time[i]);  // adding the target to the current list
      pasttime.push_back(time[i]);  // adding the time to the current list
      mostrecenttime = time[i]; // the most recent time
    }
    
  }
  return(statvector);  //return the entire vector of values (we will use R to extract the appropriate values from the vector)
  
}

