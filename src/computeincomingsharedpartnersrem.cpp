#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector computeincomingsharedparts(NumericVector time,
                                            NumericVector sampledevent,
                                            NumericVector controlevent,
                                            CharacterVector sender,
                                            CharacterVector target,
                                            CharacterVector dyad_id,
                                            double weightScheme,
                                            double counts,
                                            double cutweight,
                                            double halflife,
                                            std::string appender) {


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
  NumericVector statvector(sampleevent); //an empty starting vector to store the statistics scores

  // we need to create a dyadic list that stores for each event the past event times
  // and once they go to zero, we drop them from the set
  std::unordered_map<std::string, std::vector<std::string>> sender_neighbors; //an empty list (container) to store the past targets of each sender
  std::unordered_map<std::string, std::vector<std::string>> target_neighbors; //an empty list (container) to store the past senders of each target
  std::unordered_map<std::string, std::vector<double>> dyad_events; //an empty list (container) to store the dyad ids and the past event times
  double countsh = 0; // the count for the sender to the new target (h)
  double countrh = 0; // the count for the receiver to the new target (h)


  //looping through all events to find the past event times
  for (int i = 0; i < nevents; i++) {

    // we dont need to compute the score for all events just those that are sampled!
    if(sampledevent[i] == 1){ // if the event is sampled, we need to compute the statistics

      std::string curSender = Rcpp::as<std::string>(sender[i]); // the current event sender
      std::string curTarget = Rcpp::as<std::string>(target[i]); // the current event target
      // this extracts all of the past actors that the current sender sent a tie to
      std::vector<std::string>& sender_alters = target_neighbors[curSender];
      std::vector<std::string>& target_alters = target_neighbors[curTarget];

      std::vector<std::string> shared_actors; //an empty vector to stores the shared nodes
      //shared_actors.reserve(sender_neighbors.size()); // preallocating memory to avoid copying (for large vectors)
      int ntargets = sender_alters.size(); //the number of past actors for whom the actor is tied to
      std::unordered_set<std::string> target_set(target_alters.begin(), target_alters.end());

      //  for all actors in the targets set (that is the actors that the sender has sent a tie to)
      for (int k = 0; k < ntargets; k++) {
        if(target_set.count(sender_alters[k]) > 0){ //checking if the actors are in the dataset
          shared_actors.push_back(sender_alters[k]); //if they are present, add them to the list
        } //ending the if statement
      } //ending the loop
      // sorting the vector to remove duplicate h actors
      std::sort(shared_actors.begin(), shared_actors.end());
      // Remove duplicates h actors
      auto last = std::unique(shared_actors.begin(), shared_actors.end());
      shared_actors.erase(last, shared_actors.end()); //erasing the empty space

      double sumweights = 0; // the event exponential weights
      int nshared = shared_actors.size(); // the number of shared past actors (h)

      if(nshared > 0){ //if there are shared actors then we need to update the weights

        NumericVector weights(nshared); // an empty vector to store the weights for each actor

        for (int z = 0; z < nshared; z++) { //for each shared actor

          std::string curactor = shared_actors[z]; //the current shared actor
          std::string dyadA = curactor + appender +curTarget; //combining the strings (s,h)
          std::string dyadB = curactor + appender +curSender; //combining the strings (h,r)

          std::vector<double>& past_timesA = dyad_events[dyadA];  // this extracts all of the past event times
          std::vector<double>& past_timesB = dyad_events[dyadB];  // this extracts all of the past event times

          NumericVector weightsA(past_timesA.size()); // an empty vector to store the weights for dyad A (cur sender + shared)
          NumericVector weightsB(past_timesB.size()); // an empty vector to store the weights for dyad B (shared + cur sender)
          countsh = 0; // resetting the value to 0
          countrh = 0; // resetting the value to 0

          for (int g = 0; g < past_timesA.size(); g++) { // for all past event times

            if(past_timesA[g] < time[i]){

              if(weightScheme == 1){ //if the weighting scheme should be of Lerner and Lomi (2020)
                weightsA[g] = exp((-(time[i] - past_timesA[g]) * log(2)/(halflife)));
              }else{//if the weighting scheme should be of Lerner et al. (2013)
                weightsA[g] = exp((-(time[i] - past_timesA[g]) * log(2)/(halflife))) * log(2)/(halflife);
              }
              //check for the dyadic cutoff weight
              if(weightsA[g] < cutweight){
                weightsA[g] = 0; // if the value is too small drop the scores
              }
              if(counts == 1){ // if the user requested that the count be computed
                if(weightsA[g] > 0){// if the weight is non-zero
                  countsh += 1;  // update the count by 1
                }
              }
            }else{
              weightsA[g] = 0; // the value goes to 0
            }
          }

          for (int g = 0; g < past_timesB.size(); g++) { // for all past event times

            if(past_timesB[g] < time[i]){
              if(weightScheme == 1){ //if the weighting scheme should be of Lerner and Lomi (2020)
                weightsB[g] = exp((-(time[i] - past_timesB[g]) * log(2)/(halflife)));
              }else{//if the weighting scheme should be of Lerner et al. (2013)
                weightsB[g] = exp((-(time[i] - past_timesB[g]) * log(2)/(halflife))) * log(2)/(halflife);
              }
              //check for the dyadic cutoff weight
              if(weightsB[g] < cutweight){
                weightsB[g] = 0; // if the value is too small drop the scores
              }
              if(counts == 1){ // if the user requested that the count be computed
                if(weightsB[g] > 0){// if the weight is non-zero
                  countrh += 1;  // update the count by 1
                }
              }
            }else{
              weightsB[g] = 0; // if the value is too small drop the scores
            }
          }

          if(counts == 0){// if the exponential weights should be returned
            for(int g = 0; g < weightsA.size(); g++){ //for all shared actor times for the sender
              for(int k = 0;  k < weightsB.size(); k++){//for all shared actor times for the target
                weights[z]+= weightsA[g]*weightsB[k]; //taking the product of the two values
              }
            }

          }else{ // if the counts of events should be returned
            // per the formula, it should be the minimum number of events, thus which ever is less
            if(countrh == countsh){ // if they are the same
              weights[z] = countrh;  // use either or (the value will be the same)
            }
            if(countrh > countsh){ // if the receiver h count is greater
              weights[z] = countsh;  // return the sender h count
            }
            if(countsh > countrh){// if the sender h count is greater
              weights[z] = countrh; // return the receiver h count
            }

          }

        } // ending the if shared > 0 statement

        if(counts == 0){ // if the exponential weights should be returned
          sumweights = sqrt(sum(weights)); // from the formula take the square root of the sum
        }else{ // / if the counts of events should be returned
          sumweights = sum(weights); // from the formula take the sum of the minimum number
        }
      }

      statvector[cursampled] = sumweights; // storing the values for the ith statistic
      cursampled += 1; // update this position indicator for storing the computed values
    }

    if(controlevent[i] == 0){ // if it is not a fake event, append the current time to the vector

      //
      //    Updating the network of past events
      //

      std::string curSender = Rcpp::as<std::string>(sender[i]); // the current event sender
      std::string curTarget = Rcpp::as<std::string>(target[i]); // the current event target
      // extracting the vector of past targets the current sender has interacted with in the past
      std::vector<std::string>& past_sender_neighbors = sender_neighbors[curSender];
      // extracting the vector of past sender the current target has interacted with in the past
      std::vector<std::string>& past_target_neighbors = target_neighbors[curTarget];
      // extracting the vector of past event times for the sender where they have sent a tie
      past_sender_neighbors.push_back(Rcpp::as<std::string>(target[i])); //appending the current event target
      past_target_neighbors.push_back(Rcpp::as<std::string>(sender[i])); //appending the current event sender
      std::string curdyad = Rcpp::as<std::string>(dyad_id[i]); //the current event dyad
      std::vector<double>& past_times = dyad_events[curdyad];  // this extracts all of the past event times
      past_times.push_back(time[i]);


    }


  }
  return(statvector);  //return the entire vector of values (we will use R to extract the appropriate values from the vector)
}






