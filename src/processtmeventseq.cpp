#include <Rcpp.h>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
List processREMseqTM(std::vector<double> time,
                     std::vector<double> seqid,
                     std::vector<std::string> sender,
                     std::vector<std::string> target,
                     double pobserved = 1,
                     double ncontrols = 1,
                     double rseed = 9999) {

  double nevents = time.size(); // the number of observed evetns in the event seqeucne
  double howmany = std::round(nevents*pobserved); // this should return the number of events we need to sample
  //
  //
  //          Sampling from the observed event sequence
  //
  //
  std::vector<double> sampledevents;

  if(pobserved < 1){ // if the probability of sampling is less than 1 (i.e., we have to sample)
      //std::random_device rd; // for random uniform sampling from the observed sequence id
      //std::uniform_int_distribution<int> dist(0, nevents); //sampling from the sequence id
      //std::vector<double> sampledevents = di
      //sampling howmany observed events from the observed event sequence
      std::mt19937 gen(rseed);  // create a generator
      std::sample(seqid.begin(), seqid.end(),
                  std::back_inserter(sampledevents),
                  static_cast<size_t>(howmany),
                  gen);
  }else{
    sampledevents = seqid; // an empty vector to store the sampled events
  }

  List processedevents(howmany); //each list element will be a data frame! with 1 + ncontrols rows


  //
  //    We need to make the sender and target vectors unique here for RS
  //
  std::vector<std::string> obssender = sender;
  std::sort(sender.begin(), sender.end()); // sorting the sender pairs (we only need the unique actors for rs)
  auto last = std::unique(sender.begin(), sender.end()); // the last unique element
  sender.erase(last, sender.end()); //erasing the empty space (i.e., repeated actors)

  std::vector<std::string> obstarget = target;
  std::sort(target.begin(), target.end()); // sorting the sender pairs (we only need the unique actors for rs)
  auto last1 = std::unique(target.begin(), target.end()); // the last unique element
  target.erase(last1, target.end()); //erasing the empty space (i.e., repeated actors)

  std::random_device rd; //prepping the random device
  std::mt19937 gen(rseed); // setting the random seed
  std::uniform_int_distribution<int> sender_dist(0, sender.size()-1); //setting the distribution
  std::uniform_int_distribution<int> target_dist(0, target.size()-1); //setting the distribution

  for(int i = 0; i < howmany; i++){ // for all sampled events

    // now for each sampled event time point, we need to sample n null events (non-observed events at time t)
    double curevent = sampledevents[i] - 1; // the current observed events

    std::string cursender = obssender[curevent];
    std::string curtarget = obstarget[curevent];

    std::vector<double> curtime(ncontrols + 1); //a vector of length: n control + 1
    curtime[0] = time[curevent]; //the first place should be the observed event
    std::vector<double> curobserved(ncontrols + 1);
    curobserved[0] = 1; // this is the real event
    std::vector<double> curseqid(ncontrols + 1);
    curseqid[0] = seqid[curevent]; // the current event sequence
    std::vector<std::string> curfullsender(ncontrols + 1);
    curfullsender[0] = cursender; // the current event sender
    std::vector<std::string> curfulltarget(ncontrols + 1);
    curfulltarget[0] = curtarget; // the current event target

    // getting the list of null event senders
    std::vector<std::string> samsenders(ncontrols); // an empty vector to store the sampled senders
    std::vector<std::string> samtargets(ncontrols); // an empty vector to store the sampled targets
    std::vector<std::string> predyads; // an empty vector to store the previously sampled targets

    for(int j = 0; j < ncontrols; j++){

      double good = 0; // once good = 1 we accept the sampled dyad!

      while(good < 1){ // while good is still equal to zero

        int senderid = sender_dist(gen); // randomly sampling a sender
        int targetid = target_dist(gen); // randomly sampling a target
        std::string ssender = sender[senderid]; //the randomly sampled sender
        std::string starget = target[targetid]; //the randomly sampled target
        if(ssender == cursender && starget == curtarget){ // if it is the current dyad, sample a new dyad
          good = 0;
        }else{// if it is not the current dyad, check to make sure it is not an already sampled dyad
          if(j > 0){ //if we have already sampled an actor
            std::string checkID = (ssender + "__NIKOACAR__" + starget);
            auto check = std::find(predyads.begin(),predyads.end(),checkID);
            if(check == predyads.end()){ // if the dyad is not present then;
              samsenders[j] = ssender; //store the result
              samtargets[j] = starget; //store the result
              good = 1; // if this is not a repeated dyad, then we are good!
            }else{
              good = 0;
            }
          }else{
            good = 1; // if this is the first sample, we are good, store the results
            samsenders[j] = ssender; //store the result
            samtargets[j] = starget; //store the result
          }
        } // the end of the intial if statement
        } // the end of the while loop


      curtime[j + 1] = time[curevent]; //the first place should be the observed event
      curobserved[j + 1] = 0; // this is a null event
      curseqid[j + 1] = seqid[curevent]; // the current event sequence
      curfullsender[j + 1] = samsenders[j]; // the current event sender
      curfulltarget[j + 1] = samtargets[j]; // the current event target
      predyads.push_back(samsenders[j] + "__NIKOACAR__" + samtargets[j]);


    } // the end of the for loop for searching and sampling controls

    //storing the full resulting list that is:
    // the current dyad
    // the sampled dyad
    // the dummy vector
    // the time vector
    // creating a dataframe to store the results
    DataFrame processi =  DataFrame::create(Rcpp::Named("time") = curtime,
                                            Rcpp::Named("seqeuence_id") =curseqid,
                                            Rcpp::Named("sender") = curfullsender,
                                            Rcpp::Named("receiver") = curfulltarget,
                                            Rcpp::Named("observed") = curobserved);
    processedevents[i] = processi; // storing the dataframe

  }
  return(processedevents);
}


