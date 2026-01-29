#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector computeremweightsv2(NumericVector time,
                                 NumericVector sampledevent,
                                 NumericVector controlevent,
                                 CharacterVector dyad_id,
                                 CharacterVector dyad_idOpposite,
                                 double weightScheme,
                                 double counts,
                                 double cutweight,
                                 double halflife) {


  //The input variables for the function:
  //  time: the vector of event timings (i.e., the time the event occured)
  //  sampledevent: a binary vector (1 = sampled event [compute stats]; 0 = no)
  //  controlevent: a binary vector (1 = control event [don't update past timing]; 0 = no)
  //  dyad_id: the unique character vector for each dyad (network of past events)
  //  dyad_idC: the unique character vector for which we need to compute values for
  //  weightScheme: 1 = use lerner and lomi (2020) weight formula; 0 = lerner et al. (2013)
  //  counts: 1 = return the counts of events; 0 = return the exponential weights
  //  cutweight: the dyadic weight cutoff
  //  halflife: the halflife parameter for the exponential weighting formula


  // This is a C++ function that computes the abitrary statistic for
  // relational event sequences. Please see the R help page for the
  // respective formula for the statistic.
  //
  // This function assumes that the event sequence was "pre-prepped" in R
  // and sent to c++ to compute the statistic. As Lerner and Lomi (2020) argue
  // that in large relational event sequences, it is computationally faster to
  // update the network of past events (i.e., the past event times) than it is
  // to compute the relevant statistic for each event time point.
  //

  int nevents = time.size(); // the number of events in the relational event sequence
  int sampleevent = sum(sampledevent); //the number of sampled events for which we need to compute the statistic for
  double cursampled = 0; // an empty starting point
  NumericVector statvector(sampleevent); //an empty starting vector to store the inertia scores

  // we need to create a dyadic list that stores for each event that past event times
  // and once they go to zero, we drop them from the set
  std::unordered_map<std::string, std::vector<double>> dyad_events; //an empty list (container) to store the dyad ids and the past event times

  //looping through all events to find the past event times
  for (int i = 0; i < nevents; i++) {
    std::string curdyad = Rcpp::as<std::string>(dyad_id[i]); //the current event dyad
    std::vector<double>& past_times = dyad_events[curdyad];  // this extracts all of the past event times

    // we dont need to compute the score for all events just those that are sampled!
    if(sampledevent[i] == 1){ // if the event is sampled, we need to compute the statistics

      double sumweights = 0; // the event times
      int ni = past_times.size(); // the number of past event times

      if(ni > 0){ //if there are past weight times then update!
        NumericVector weights(past_times.size()); // the numeric vector of exponential weights

        for (int z = 0; z < ni; z++) {

          if(time[i] != past_times[z]){

            if(weightScheme == 1){ //if the weighting scheme should be of Lerner and Lomi (2020)
              weights[z] = exp((-(time[i] - past_times[z]) * log(2)/(halflife)));
            }else{//if the weighting scheme should be of Lerner et al. (2013)
              weights[z] = exp((-(time[i] - past_times[z]) * log(2)/(halflife))) * log(2)/(halflife);
            }
            //check for the dyadic cutoff weight
            if(weights[z] < cutweight){
              weights[z] = 0; // if the value is too small drop the scores
            }

          }else{
            weights[z] = 0; // the evnet time is the current time (this should not really happen empirically)
          }
        }

        if(counts == 1){
          //the below function will remove all elements in whi
          weights.erase(std::remove(weights.begin(), weights.end(), 0), weights.end());
          sumweights += weights.size(); // the number of past non-zero exponential weights
        }else{
          sumweights += sum(weights); // the sum of the exponential weights
        }

      }

      statvector[cursampled] = sumweights; // storing the values for the ith statistic
      cursampled += 1; // update this position indicator for storing the computed values
    }

    if(controlevent[i] != 1){ // if it is not a fake event, append the current time to the vector
      std::string curdyadC = Rcpp::as<std::string>(dyad_idOpposite[i]); //the current event dyad
      std::vector<double>& past_times = dyad_events[curdyadC];  // this extracts all of the past event times
      past_times.push_back(time[i]);
      if(cutweight < 0){ // updating the past edges!

        int npast = past_times.size(); // the number of past weights
        double pastweight; // the past cut weight
        int drop = -1;
        for (int z = 0; z < npast; z++) { // for all past events
          if(weightScheme == 1){ //if the weighting scheme should be of Lerner and Lomi (2020)
            pastweight= exp((-(time[i] - past_times[z]) * log(2)/(halflife)));
          }else{//if the weighting scheme should be of Lerner et al. (2013)
            pastweight= exp((-(time[i] - past_times[z]) * log(2)/(halflife))) * log(2)/(halflife);
          }
          if(pastweight >= cutweight){// if the weight is good
            drop = z; //
            break; // breaking out of the past weight check for loop
          }
        }
        if(drop == 0){ // no need to drop past event times
          past_times.push_back(time[i]); // appending the past times vector
        }else{
          if(drop == -1){ // need to remove all past weights as they are all zero
            past_times.clear();
            past_times.push_back(time[i]); // keep only the current event
          }else{ // drop those irrelevant temporal times
            past_times.erase(past_times.begin(), past_times.begin() + drop);
            past_times.push_back(time[i]); // appending the past times vector
          }
        }
      }else{
        past_times.push_back(time[i]); // appending the past times vector
      }
    }

  }
  return(statvector);  //return the entire vector of values (we will use R to extract the appropriate values from the vector)
}
