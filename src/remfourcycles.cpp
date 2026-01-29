#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector computefourcyclesrem(NumericVector time,
                                 NumericVector sampledevent,
                                 NumericVector controlevent,
                                 CharacterVector sender,
                                 CharacterVector target,
                                 CharacterVector dyad_id,
                                 double weightScheme,
                                 double counts,
                                 double cutweight,
                                 double halflife,
                                 std::string delim) {


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


  // This is a C++ function that computes the four-cycle statistic for
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
  double countAB = 0; // the count for the sender to the new target (h)
  double countAR = 0; // the count for the receiver to the new target (h)
  double countSB = 0; // the count for the receiver to the new target (h)


  //looping through all events to find the past event times
  for (int i = 0; i < nevents; i++) {

    // we dont need to compute the score for all events just those that are sampled!
    if(sampledevent[i] == 1){ // if the event is sampled, we need to compute the statistics

      std::string curSender = Rcpp::as<std::string>(sender[i]); // the current event sender
      std::string curTarget = Rcpp::as<std::string>(target[i]); // the current event target
      // this extracts all of the past actors that the current sender sent a tie to
      std::vector<std::string>& sender_alters = sender_neighbors[curSender];
      std::vector<std::string>& target_alters = target_neighbors[curTarget];

      //std::vector<std::string> shared_pairs_A; //an empty vector to stores the shared sender from the grouped dyad (s')
      //std::vector<std::string> shared_pairs_B; //an empty vector to stores the shared target from the grouped dyad (r')


      //shared_actors.reserve(sender_neighbors.size()); // preallocating memory to avoid copying (for large vectors)
      int ntargets = sender_alters.size(); //the number of past actors for whom the actor is tied to
      int nsenders = target_alters.size(); //the number of past actors for whom the actor is tied to

      std::vector<std::string> dyadic_pairs;//an empty vector to stores the shared sender from the grouped dyad (s')

      //  for all actors in the targets set (that is the actors that the sender has sent a tie to)
      for (int k = 0; k < ntargets; k++) {

        //we have to find the w(s',r') pairs that fit the criteria: w(s,r') and w(s',r) from the current pair w(s,r)
        std::string target_k = sender_alters[k]; // the current target for the current sender (r')

        for (int j = 0; j < nsenders; j++) {
          std::string sender_j = target_alters[j]; // the current target for the current sender (s')
          //checking if they have interacted with each other in the past

          if(target_k != curTarget && sender_j != curSender){ // making sure that the current dyad is not in the list! (a fail safe I hope)

          std::string dyad_kj = sender_j+delim+target_k; //combining the actors together
          std::vector<double> kj_past_times = dyad_events[dyad_kj]; // the past interaction times
          // auto itA = dyad_events.find(dy_sprime_r); do something here
            if(kj_past_times.size()>0){ //if they have interacted in the past
              dyadic_pairs.push_back(dyad_kj);
              //shared_pairs_A.push_back(sender_j);//combining the actors together
              //shared_pairs_B.push_back(target_k);//combining the actors together
            }
          }
        }
      }
      std::sort(dyadic_pairs.begin(), dyadic_pairs.end());
      // Remove duplicates h actors
      auto last = std::unique(dyadic_pairs.begin(), dyadic_pairs.end());
      dyadic_pairs.erase(last, dyadic_pairs.end()); //erasing the empty space

      double sumweights = 0; // the event exponential weights
      int nshared = dyadic_pairs.size(); // the number of shared past actors (h) _A and _B are the same length

      if(nshared > 0){ //if there are shared actors then we need to update the weights

        NumericVector weights(nshared); // an empty vector to store the weights for each dyadic pair

        for (int z = 0; z < nshared; z++) { //for each shared actor

          std::string curdyad = dyadic_pairs[z];
          std::size_t delimpos = curdyad.find(delim);
          std::string curh_sender = curdyad.substr(0,delimpos);
          std::string curh_target = curdyad.substr(delimpos + delim.length()); // 5 = length of "_DOG_"

          std::string dyadAB = curh_sender   + delim + curTarget; //curSender + "_DOG_" + curTarget; //combining the strings (s,h)
          std::string dyadAR = curSender + delim + curh_target; // //combining the strings (s,h)
          std::string dyadSB = curh_sender + delim + curh_target; //combining the strings (s,h)

          std::vector<double>& past_timesA = dyad_events[dyadAB];  // this extracts all of the past event times (for the real observed dyad)
          std::vector<double>& past_timesB = dyad_events[dyadAR];  // this extracts all of the past event times (for the observed sender and h target)
          std::vector<double>& past_timesC = dyad_events[dyadSB];  // this extracts all of the past event times (for the h sender and observed target)

          NumericVector weightsA(past_timesA.size()); // an empty vector to store the weights for dyad A (cur sender + shared)
          NumericVector weightsB(past_timesB.size()); // an empty vector to store the weights for dyad B (shared + cur sender)
          NumericVector weightsC(past_timesC.size()); // an empty vector to store the weights for dyad B (shared + cur sender)
          countAB = 0; // resetting the value to 0
          countAR = 0; // resetting the value to 0
          countSB = 0; // resetting the value to 0


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
                  countAB += 1;  // update the count by 1
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
                  countAR += 1;  // update the count by 1
                }
              }

            }else{
              weightsB[g] = 0; // if the value is too small drop the scores
            }
          }

          for (int g = 0; g < past_timesC.size(); g++) { // for all past event times

            if(past_timesC[g] < time[i]){

              if(weightScheme == 1){ //if the weighting scheme should be of Lerner and Lomi (2020)
                weightsC[g] = exp((-(time[i] - past_timesC[g]) * log(2)/(halflife)));
              }else{//if the weighting scheme should be of Lerner et al. (2013)
                weightsC[g] = exp((-(time[i] - past_timesC[g]) * log(2)/(halflife))) * log(2)/(halflife);
              }
              //check for the dyadic cutoff weight
              if(weightsC[g] < cutweight){
                weightsC[g] = 0; // if the value is too small drop the scores
              }
              if(counts == 1){ // if the user requested that the count be computed
                if(weightsC[g] > 0){// if the weight is non-zero
                  countSB += 1;  // update the count by 1
                }
              }
            }else{
              weightsC[g] = 0; // the value goes to 0
            }
          }


       if(counts == 0){// if the exponential weights should be returned

          for(int g = 0; g < weightsA.size(); g++){ //for all shared actor times for the sender

            for(int k = 0;  k < weightsB.size(); k++){//for all shared actor times for the target

              for(int j = 0;  j < weightsC.size(); j++){

              weights[z]+= weightsA[g]*weightsB[k]*weightsC[j]; //taking the product of the two values

              }
            }
          }


       }else{
         std::vector<double> pastcounts = {countAB,countAR,countSB};
         weights[z] = *std::min_element(pastcounts.begin(),pastcounts.end());
       }


       }

        if(counts == 0){ // if the exponential weights should be returned
          sumweights = std::cbrt(sum(weights)); // from the formula take the cube of the formula
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
