#include <Rcpp.h>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
List processREMseqTM_varying(std::vector<double> time,
                             std::vector<double> seqid,
                             std::vector<std::string> sender,
                             std::vector<std::string> target,
                             double pobserved = 1,
                             double ncontrols = 1,
                             std::string appender = "__NIKOACAR3718__",
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
  std::unordered_set<std::string> targetset; // the past active targets at time t (for uniqueness)
  std::unordered_set<std::string> senderset; // the past active senders at time t (for uniqueness)
  std::vector<std::string> pasttargets; // the past active targets at time t (for sampling)
  std::vector<std::string> pastsenders; // the past active senders at time t (for sampling)
  std::size_t pastupdateindex = 0; // the most recent past time
  std::random_device rd; //prepping the random device
  std::mt19937 gen(rseed); // setting the random seed
  //
  //
  //         The First Sampled Event
  //
  //
  //double curevent = sampledevents[0] - 1; // the current observed events
  size_t curevent = static_cast<size_t>(sampledevents[0]) - 1; // the current observed events
  std::string cursender = sender[curevent]; // the current event sender
  std::string curtarget = target[curevent]; // the current event target
  std::string curdyad;

  if(sampledevents[0] == 1){
    std::vector<double> curtime(1); //a vector of length: n control + 1
    curtime[0] = time[curevent]; //the first place should be the observed event
    std::vector<double> curobserved(1);
    curobserved[0] = 1; // this is the real event
    std::vector<double> curseqid(1);
    curseqid[0] = seqid[curevent]; // the current event sequence
    std::vector<std::string> curfullsender(1);
    curfullsender[0] = cursender; // the current event sender
    std::vector<std::string> curfulltarget(1);
    curfulltarget[0] = curtarget; // the current event target
    DataFrame processi =  DataFrame::create(Rcpp::Named("time") = curtime,
                                            Rcpp::Named("seqeuence_id") =curseqid,
                                            Rcpp::Named("sender") = curfullsender,
                                            Rcpp::Named("receiver") = curfulltarget,
                                            Rcpp::Named("observed") = curobserved);
    processedevents[0] = processi; // storing the dataframe
  }else{

    double curtimei = time[curevent]; //a vector of length: n control + 1
    for(int j = 0; j < nevents; j++){ // for all events until the first one
      if(time[j] <= curtimei){ // if the time j is less or equal to the current event time
        if(targetset.insert(target[j]).second){ // if the current target can be inserted (i.e., not already in the set)
          pasttargets.push_back(target[j]); // adding the current target
        }
        if(senderset.insert(sender[j]).second){ // if the current sender can be inserted (i.e., not already in the set)
          pastsenders.push_back(sender[j]); // adding the current target
        }
      }else{
        pastupdateindex = j; // the most recent time that is beyond the current time point
        break; // escaping the loop as time is beyond ti
      }
    } // ending the searching loop to add nodes to their respective lists!

    if((pastsenders.size() * pasttargets.size()) > (ncontrols + 1)){ // if there are enough nodes for sampling!
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
      std::uniform_int_distribution<int> sender_dist(0, pastsenders.size()-1); //setting the distribution
      std::uniform_int_distribution<int> target_dist(0, pasttargets.size()-1); //setting the distribution

      for(int j = 0; j < ncontrols; j++){
        double good = 0; // once good = 1 we accept the sampled dyad!
        while(good < 1){ // while good is still equal to zero

          int senderid = sender_dist(gen); // randomly sampling a sender
          int targetid = target_dist(gen); // randomly sampling a target
          std::string ssender = pastsenders[senderid]; //the randomly sampled sender
          std::string starget = pasttargets[targetid]; //the randomly sampled target
          if(ssender == cursender && starget == curtarget){ // if it is the current dyad, sample a new dyad
            good = 0;
          }else{// if it is not the current dyad, check to make sure it is not an already sampled dyad
            if(j > 0){ //if we have already sampled an actor
              std::string checkID = (ssender + appender + starget);
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
        predyads.push_back(samsenders[j] + appender + samtargets[j]);


      } // the end of the for loop for searching and sampling controls

      DataFrame processi =  DataFrame::create(Rcpp::Named("time") = curtime,
                                              Rcpp::Named("seqeuence_id") =curseqid,
                                              Rcpp::Named("sender") = curfullsender,
                                              Rcpp::Named("receiver") = curfulltarget,
                                              Rcpp::Named("observed") = curobserved);
      processedevents[0] = processi; // storing the dataframe


    }else{ // if not, then we have to make the full event set
      std::vector<std::string> samsenders; // an empty vector to store the sampled senders
      std::vector<std::string> samtargets; // an empty vector to store the sampled targets
      std::vector<double> curtime; //a vector of length: n control + 1
      curtime.push_back(time[curevent]); //the first place should be the observed event
      std::vector<double> curobserved;
      curobserved.push_back(1); // this is the real event
      std::vector<double> curseqid;
      curseqid.push_back(seqid[curevent]);
      std::vector<std::string> curfullsender;
      curfullsender.push_back(cursender); // the current event sender
      std::vector<std::string> curfulltarget;
      curfulltarget.push_back(curtarget); // the current event sender
      for(int z = 0; z < pastsenders.size(); z++){
        for(int y = 0; y < pasttargets.size(); y++){
          if(!(pastsenders[z] == cursender && pasttargets[y] == curtarget)){
            curfullsender.push_back(pastsenders[z]); // adding the past z sender
            curfulltarget.push_back(pasttargets[y]); // adding the past y target
            curtime.push_back(time[curevent]); //the first place should be the observed event
            curobserved.push_back(0); // this is the real event
            curseqid.push_back(seqid[curevent]);
          }
        } // end of y loop
      } // end of z loop
      DataFrame processi =  DataFrame::create(Rcpp::Named("time") = curtime,
                                              Rcpp::Named("seqeuence_id") =curseqid,
                                              Rcpp::Named("sender") = curfullsender,
                                              Rcpp::Named("receiver") = curfulltarget,
                                              Rcpp::Named("observed") = curobserved);
      processedevents[0] = processi; // storing the dataframe
    }

  }


  for(int i = 1; i < howmany; i++){ // for all sampled events (we start at the second event here)

    // now for each sampled event time point, we need to sample n null events (non-observed events at time t)
    //double curevent = sampledevents[i] - 1; // the current observed events
    size_t curevent = static_cast<size_t>(sampledevents[i]) - 1; // the current observed events
    std::string cursender = sender[curevent]; // the current event sender
    std::string curtarget = target[curevent]; // the current event target
    // getting the list of null event senders
    std::vector<std::string> samsenders(ncontrols); // an empty vector to store the sampled senders
    std::vector<std::string> samtargets(ncontrols); // an empty vector to store the sampled targets
    std::vector<std::string> predyads; // an empty vector to store the previously sampled targets
    //
    //
    //          Updating the network of past events
    //
    //
    double curtimei = time[curevent]; //a vector of length: n control + 1
    for(int k = pastupdateindex; k < nevents; k++){ // for all events until the first one
      if(time[k] <= curtimei){ // if the time j is less or equal to the current event time
        if(targetset.insert(target[k]).second){ // if the current target can be inserted (i.e., not already in the set)
          pasttargets.push_back(target[k]); // adding the current target
        }
        if(senderset.insert(sender[k]).second){ // if the current sender can be inserted (i.e., not already in the set)
          pastsenders.push_back(sender[k]); // adding the current target
        }
      }else{
        pastupdateindex = k; // the most recent time that is beyond the current time point
        break; // escaping the loop as time is beyond ti
      }
    } // ending the searching loop to add nodes to their respective lists!


   if((pastsenders.size() * pasttargets.size()) > (ncontrols + 1)){ // if there are enough nodes for sampling!
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
      std::uniform_int_distribution<int> sender_dist(0, pastsenders.size()-1); //setting the distribution
      std::uniform_int_distribution<int> target_dist(0, pasttargets.size()-1); //setting the distribution

      for(int j = 0; j < ncontrols; j++){
        double good = 0; // once good = 1 we accept the sampled dyad!
        while(good < 1){ // while good is still equal to zero

          int senderid = sender_dist(gen); // randomly sampling a sender
          int targetid = target_dist(gen); // randomly sampling a target
          std::string ssender = pastsenders[senderid]; //the randomly sampled sender
          std::string starget = pasttargets[targetid]; //the randomly sampled target
          if(ssender == cursender && starget == curtarget){ // if it is the current dyad, sample a new dyad
            good = 0;
          }else{// if it is not the current dyad, check to make sure it is not an already sampled dyad
            if(j > 0){ //if we have already sampled an actor
              std::string checkID = (ssender + appender + starget);
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
        predyads.push_back(samsenders[j] + appender + samtargets[j]);


      } // the end of the for loop for searching and sampling controls

      DataFrame processi =  DataFrame::create(Rcpp::Named("time") = curtime,
                                              Rcpp::Named("seqeuence_id") =curseqid,
                                              Rcpp::Named("sender") = curfullsender,
                                              Rcpp::Named("receiver") = curfulltarget,
                                              Rcpp::Named("observed") = curobserved);
      processedevents[i] = processi; // storing the dataframe

    }else{ // if not, then we have to make the full event set
      std::vector<std::string> samsenders; // an empty vector to store the sampled senders
      std::vector<std::string> samtargets; // an empty vector to store the sampled targets
      std::vector<double> curtime; //a vector of length: n control + 1
      curtime.push_back(time[curevent]); //the first place should be the observed event
      std::vector<double> curobserved;
      curobserved.push_back(1); // this is the real event
      std::vector<double> curseqid;
      curseqid.push_back(seqid[curevent]);
      std::vector<std::string> curfullsender;
      curfullsender.push_back(cursender); // the current event sender
      std::vector<std::string> curfulltarget;
      curfulltarget.push_back(curtarget); // the current event sender
      for(int z = 0; z < pastsenders.size(); z++){
        for(int y = 0; y < pasttargets.size(); y++){
          if(!(pastsenders[z] == cursender && pasttargets[y] == curtarget)){
            curfullsender.push_back(pastsenders[z]); // adding the past z sender
            curfulltarget.push_back(pasttargets[y]); // adding the past y target
            curtime.push_back(time[curevent]); //the first place should be the observed event
            curobserved.push_back(0); // this is the real event
            curseqid.push_back(seqid[curevent]);
          }
        } // end of y loop
      } // end of z loop
      DataFrame processi =  DataFrame::create(Rcpp::Named("time") = curtime,
                                              Rcpp::Named("seqeuence_id") =curseqid,
                                              Rcpp::Named("sender") = curfullsender,
                                              Rcpp::Named("receiver") = curfulltarget,
                                              Rcpp::Named("observed") = curobserved);
      processedevents[i] = processi; // storing the dataframe
    }
  }
  return(processedevents);
}
