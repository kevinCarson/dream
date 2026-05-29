#include <Rcpp.h>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
List processREMseqOM(std::vector<double> time,
                     NumericVector seqid,
                     std::vector<std::string> sender,
                     std::vector<std::string> target,
                     double pobserved = 1,
                     double ncontrols = 1,
                     bool interval = false,
                     double t = 0) {

  double nevents = time.size(); // the number of observed evetns in the event seqeucne
  int howmany = std::round(nevents*pobserved); // this should return the number of events we need to sample
  RNGScope scope;
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
    sampledevents = as<std::vector<double>>(sample(seqid,howmany,false));
  }else{
    sampledevents = as<std::vector<double>>(seqid); // an empty vector to store the sampled events
  }


  List processedevents(howmany);//each list element will be a data frame! with 1 + ncontrols rows
  //
  //    We need to make the sender and target vectors unique here for RS
  //
  std::vector<std::string> actors = sender;
  actors.insert(actors.end(), target.begin(), target.end()); //since it is one mode, append the two vectors together (this may become too big)
  std::sort(actors.begin(), actors.end()); // sorting the sender pairs (we only need the unique actors for rs)
  auto last = std::unique(actors.begin(), actors.end()); // the last unique element
  actors.erase(last, actors.end()); //erasing the empty space (i.e., repeated actors)
  CharacterVector actors1 = Rcpp::wrap(actors);
  for(int i = 0; i < howmany; i++){ // for all sampled events

    // now for each sampled event time point, we need to sample n null events (non-observed events at time t)
    double curevent = sampledevents[i] - 1; // the current observed events

    std::string cursender = sender[curevent];
    std::string curtarget = target[curevent];

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

        CharacterVector sendsample = sample(actors1,1); //the randomly sampled sender
        CharacterVector recsample = sample(actors1,1); //the randomly sampled sender
        std::string ssender = as<std::string>(sendsample[0]); //the randomly sampled sender
        std::string starget =  as<std::string>(recsample[0]); //the randomly sampled target
        if(ssender != starget){
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
        }

      } // the end of the while loop


      curtime[j + 1] = time[curevent]; //the first place should be the observed event
      curobserved[j + 1] = 0; // this is a null event
      curseqid[j + 1] = seqid[curevent]; // the current event sequence
      curfullsender[j + 1] = samsenders[j]; // the current event sender
      curfulltarget[j + 1] = samtargets[j]; // the current event target
      predyads.push_back(samsenders[j] + "__NIKOACAR__" + samtargets[j]);//adding the new dyad to the vector!




    } // the end of the for loop for searching and sampling controls

    //storing the full resulting list that is:
    // the current dyad
    // the sampled dyad
    // the dummy vector
    // the time vector
    // creating a dataframe to store the results
    DataFrame processi =  DataFrame::create(Rcpp::Named("time") = curtime,
                                            Rcpp::Named("seqeuence_id") = curseqid,
                                            Rcpp::Named("sender") = curfullsender,
                                            Rcpp::Named("target") = curfulltarget,
                                            Rcpp::Named("observed") = curobserved);
    processedevents[i] = processi; // storing the dataframe

  }

  if(interval && t > 0){

    double curevent = seqid.size() + 1; // updating the count by 1

    std::vector<double> curtime(ncontrols); //a vector of length: n control + 1
    std::vector<double> curobserved(ncontrols);
    std::vector<double> curseqid(ncontrols);
    std::vector<std::string> curfullsender(ncontrols);
    std::vector<std::string> curfulltarget(ncontrols);

    // getting the list of null event senders
    std::vector<std::string> samsenders(ncontrols); // an empty vector to store the sampled senders
    std::vector<std::string> samtargets(ncontrols); // an empty vector to store the sampled targets
    std::vector<std::string> predyads; // an empty vector to store the previously sampled targets

    for(int j = 0; j < ncontrols; j++){

      double good = 0; // once good = 1 we accept the sampled dyad!

      while(good < 1){ // while good is still equal to zero

        CharacterVector sendsample = sample(actors1,1); //the randomly sampled sender
        CharacterVector recsample = sample(actors1,1); //the randomly sampled sender
        std::string ssender = as<std::string>(sendsample[0]); //the randomly sampled sender
        std::string starget =  as<std::string>(recsample[0]); //the randomly sampled target

        if(ssender != starget){
          // if it is not the current dyad, check to make sure it is not an already sampled dyad
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
        } // the end of the while loop
      }
      curtime[j] = t; //the first place should be the observed event
      curobserved[j] = 0; // this is a null event
      curseqid[j ] = curevent; // the current event sequence
      curfullsender[j] = samsenders[j]; // the current event sender
      curfulltarget[j] = samtargets[j]; // the current event target
      predyads.push_back(samsenders[j] + "__NIKOACAR__" + samtargets[j]);//adding the new dyad to the vector!
    } // the end of the for loop for searching and sampling controls

    //storing the full resulting list that is:
    // the current dyad
    // the sampled dyad
    // the dummy vector
    // the time vector
    // creating a dataframe to store the results
    DataFrame processi =  DataFrame::create(Rcpp::Named("time") = curtime,
                                            Rcpp::Named("seqeuence_id") = curseqid,
                                            Rcpp::Named("sender") = curfullsender,
                                            Rcpp::Named("target") = curfulltarget,
                                            Rcpp::Named("observed") = curobserved);

    processedevents.push_back(processi);
  }


  return(processedevents);
}

