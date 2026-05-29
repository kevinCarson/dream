#include <Rcpp.h>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
List timevaryingomrisksetwithsample(std::vector<double> time,
                                   NumericVector seqid,
                                   std::vector<std::string> sender,
                                   std::vector<std::string> target,
                                   std::vector<double> timestart,
                                   std::vector<double> timeend,
                                   std::vector<std::string> actors,
                                   double pobserved = 1,
                                   double ncontrols = 1,
                                   std::string appender = "__NIKOACAR3718__",
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
  std::vector<std::string> senderobs =sender;
  std::vector<std::string> targetobs =target;
  List processedevents(howmany);//each list element will be a data frame! with 1 + ncontrols rows
  int n = actors.size();

  for(int i = 0; i < howmany; i++){ // for all sampled events

    // now for each sampled event time point, we need to sample n null events (non-observed events at time t)
    double curevent = sampledevents[i] - 1; // the current observed events
    double curtime = time[curevent]; // getting the current event time
    std::string cursender = sender[curevent];
    std::string curtarget = target[curevent];
    std::vector<std::string> possibleactors;
    possibleactors.reserve(n); // preallocating the memory
    //searching through the possible range of actors
    for(int j = 0; j < n; j++){
      if(timeend[j] >= curtime && curtime >= timestart[j]){
        possibleactors.push_back(actors[j]);
      }
    }
    int outer = possibleactors.size();
    std::vector<std::string> samsenders(ncontrols); // an empty vector to store the sampled senders
    std::vector<std::string> samtargets(ncontrols); // an empty vector to store the sampled targets
    std::vector<std::string> predyads; // an empty vector to store the previously sampled targets

    if((outer * (outer-1)) > (ncontrols + 1)){

      std::vector<double> curtimevec(ncontrols + 1); //a vector of length: n control + 1
      curtimevec[0] = curtime; //the first place should be the observed event
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
      CharacterVector actors1 = Rcpp::wrap(possibleactors);

      for(int j = 0; j < ncontrols; j++){
        double good = 0; // once good = 1 we accept the sampled dyad!
        while(good < 1){ // while good is still equal to zero

          CharacterVector sendsample = sample(actors1,1); //the randomly sampled sender
          CharacterVector recsample = sample(actors1,1); //the randomly sampled sender
          std::string ssender = as<std::string>(sendsample[0]); //the randomly sampled sender
          std::string starget =  as<std::string>(recsample[0]); //the randomly sampled target

          if((ssender == cursender && starget == curtarget) &&
             (ssender != curtarget)){ // if it is the current dyad, sample a new dyad, or if this is a self-loop
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


        curtimevec[j + 1] = time[curevent]; //the first place should be the observed event
        curobserved[j + 1] = 0; // this is a null event
        curseqid[j + 1] = seqid[curevent]; // the current event sequence
        curfullsender[j + 1] = samsenders[j]; // the current event sender
        curfulltarget[j + 1] = samtargets[j]; // the current event target
        predyads.push_back(samsenders[j] + appender + samtargets[j]);


      } // the end of the for loop for searching and sampling controls

      DataFrame processi =  DataFrame::create(Rcpp::Named("time") = curtimevec,
                                              Rcpp::Named("seqeuence_id") =curseqid,
                                              Rcpp::Named("sender") = curfullsender,
                                              Rcpp::Named("receiver") = curfulltarget,
                                              Rcpp::Named("observed") = curobserved);
      processedevents[i] = processi; // storing the dataframe

    }else{
      std::vector<std::string> samsenders; // an empty vector to store the sampled senders
      std::vector<std::string> samtargets; // an empty vector to store the sampled targets
      std::vector<double> curtimevec; //a vector of length: n control + 1
      curtimevec.push_back(curtime); //the first place should be the observed event
      std::vector<double> curobserved;
      curobserved.push_back(1); // this is the real event
      std::vector<double> curseqid;
      curseqid.push_back(seqid[curevent]);
      std::vector<std::string> curfullsender;
      curfullsender.push_back(cursender); // the current event sender
      std::vector<std::string> curfulltarget;
      curfulltarget.push_back(curtarget); // the current event sender
      for(int z = 0; z < possibleactors.size(); z++){
        for(int y = 0; y < possibleactors.size(); y++){
          if(!(possibleactors[z] == cursender && possibleactors[y] == curtarget) &&
             (z != y)){
            curfullsender.push_back(possibleactors[z]); // adding the past z sender
            curfulltarget.push_back(possibleactors[y]); // adding the past y target
            curtimevec.push_back(curtime); //the first place should be the observed event
            curobserved.push_back(0); // this is the real event
            curseqid.push_back(seqid[curevent]);
          }
        } // end of y loop
      } // end of z loop
      DataFrame processi =  DataFrame::create(Rcpp::Named("time") = curtimevec,
                                              Rcpp::Named("seqeuence_id") =curseqid,
                                              Rcpp::Named("sender") = curfullsender,
                                              Rcpp::Named("receiver") = curfulltarget,
                                              Rcpp::Named("observed") = curobserved);

      processedevents[i] = processi; // storing the dataframe

    }


  }

  if(interval && t > 0){

    std::vector<std::string> samsenders(ncontrols); // an empty vector to store the sampled senders
    std::vector<std::string> samtargets(ncontrols); // an empty vector to store the sampled targets
    std::vector<std::string> predyads; // an empty vector to store the previously sampled targets
    double curtime = t; // getting the current event time
    std::vector<std::string> possibleactors;
    possibleactors.reserve(n); // preallocating the memory
    //searching through the possible range of actors
    for(int j = 0; j < n; j++){
      if(timeend[j] >= curtime && curtime >= timestart[j]){
        possibleactors.push_back(actors[j]);
      }
    }
    int outer = possibleactors.size();
    if((outer * (outer - 1)) > ncontrols){

      CharacterVector actors1 = Rcpp::wrap(possibleactors);
      std::vector<double> curtimevec; //a vector of length: n control + 1
      std::vector<double> curobserved;
      std::vector<double> curseqid;
      std::vector<std::string> curfullsender;
      std::vector<std::string> curfulltarget;

      for(int j = 0; j < ncontrols; j++){
        double good = 0; // once good = 1 we accept the sampled dyad!
        while(good < 1){ // while good is still equal to zero

          CharacterVector sendsample = sample(actors1,1); //the randomly sampled sender
          CharacterVector recsample = sample(actors1,1); //the randomly sampled sender
          std::string ssender = as<std::string>(sendsample[0]); //the randomly sampled sender
          std::string starget =  as<std::string>(recsample[0]); //the randomly sampled target
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

        curtimevec[j + 1] = curtime; //the first place should be the observed event
        curobserved[j + 1] = 0; // this is a null event
        curseqid[j + 1] = seqid.size() + 1; // the current event sequence
        curfullsender[j + 1] = samsenders[j]; // the current event sender
        curfulltarget[j + 1] = samtargets[j]; // the current event target
        predyads.push_back(samsenders[j] + appender + samtargets[j]);

      }
      }

    DataFrame processi =  DataFrame::create(Rcpp::Named("time") = curtimevec,
                                            Rcpp::Named("seqeuence_id") =curseqid,
                                            Rcpp::Named("sender") = curfullsender,
                                            Rcpp::Named("receiver") = curfulltarget,
                                            Rcpp::Named("observed") = curobserved);

      processedevents.push_back(processi);

    }else{
      int nt = outer * (outer - 1);
      std::vector<double> curtimevec(nt); //a vector of length: n control + 1
      std::fill(curtimevec.begin(), curtimevec.end(), curtime);
      std::vector<double> curobserved(nt);
      std::fill(curobserved.begin(), curobserved.end(), 0);
      std::vector<double> curseqid(nt);
      std::fill(curseqid.begin(), curseqid.end(), seqid.size() + 1);

      std::vector<std::string> sendersRS(nt);
      std::vector<std::string> receiversRS(nt);
      int index = 0;
      //creating the empty vectors to store the sender charcter vectors
      for(int k = 0; k < outer; k++){
        for(int j = 0; j < k; j++){
          sendersRS[index] = possibleactors[k];
          receiversRS[index] = possibleactors[j];
          index += 1;

          sendersRS[index] = possibleactors[j];
          receiversRS[index] = possibleactors[k];
          index += 1;
        }
      }
      DataFrame processi =  DataFrame::create(Rcpp::Named("time") = curtimevec,
                                              Rcpp::Named("seqeuence_id") =curseqid,
                                              Rcpp::Named("sender") = sendersRS,
                                              Rcpp::Named("receiver") = receiversRS,
                                              Rcpp::Named("observed") = curobserved);
      processedevents.push_back(processi);

    }
  }
  return processedevents;
}





// [[Rcpp::export]]
List timevaryingtmrisksetwithsample(std::vector<double> time,
                                    NumericVector seqid,
                                    std::vector<std::string> sender,
                                    std::vector<std::string> target,
                                    std::vector<double> timestartsender,
                                    std::vector<double> timeendsenders,
                                    std::vector<std::string> actorsenders,
                                    std::vector<double> timestarttarget,
                                    std::vector<double> timeendtargets,
                                    std::vector<std::string> actortargets,
                                    double pobserved = 1,
                                    double ncontrols = 1,
                                    std::string appender = "__NIKOACAR3718__",
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
  int nt = actortargets.size();
  int ns = actorsenders.size();


  for(int i = 0; i < howmany; i++){ // for all sampled events

    // now for each sampled event time point, we need to sample n null events (non-observed events at time t)
    double curevent = sampledevents[i] - 1; // the current observed events
    double curtime = time[curevent]; // getting the current event time
    std::string cursender = sender[curevent];
    std::string curtarget = target[curevent];
    std::vector<std::string> possiblesenders;
    possiblesenders.reserve(ns); // preallocating the memory
    //searching through the possible range of actors
    for(int j = 0; j < ns; j++){
      if(timeendsenders[j] >= curtime && curtime >= timestartsender[j]){
        possiblesenders.push_back(actorsenders[j]);
      }
    }

    std::vector<std::string> possibletargets;
    possibletargets.reserve(nt); // preallocating the memory
    //searching through the possible range of actors
    for(int j = 0; j < nt; j++){
      if(timeendtargets[j] >= curtime && curtime >= timestarttarget[j]){
        possibletargets.push_back(actortargets[j]);
      }
    }

    int outer = possiblesenders.size();
    int innner = possibletargets.size();

    std::vector<std::string> samsenders(ncontrols); // an empty vector to store the sampled senders
    std::vector<std::string> samtargets(ncontrols); // an empty vector to store the sampled targets
    std::vector<std::string> predyads; // an empty vector to store the previously sampled targets


    if(outer*innner > (ncontrols + 1)){ // if there are enough controls to sample!

      std::vector<double> curtimevec(ncontrols + 1); //a vector of length: n control + 1
      curtimevec[0] = curtime; //the first place should be the observed event
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
      CharacterVector actors1 = Rcpp::wrap(possiblesenders);
      CharacterVector actors2 = Rcpp::wrap(possibletargets);

      for(int j = 0; j < ncontrols; j++){
        double good = 0; // once good = 1 we accept the sampled dyad!
        while(good < 1){ // while good is still equal to zero

          CharacterVector sendsample = sample(actors1,1); //the randomly sampled sender
          CharacterVector recsample = sample(actors1,1); //the randomly sampled sender
          std::string ssender = as<std::string>(sendsample[0]); //the randomly sampled sender
          std::string starget =  as<std::string>(recsample[0]); //the randomly sampled target

          if((ssender == cursender && starget == curtarget) &&
             (ssender != curtarget)){ // if it is the current dyad, sample a new dyad, or if this is a self-loop
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


        curtimevec[j + 1] = curtime; //the first place should be the observed event
        curobserved[j + 1] = 0; // this is a null event
        curseqid[j + 1] = seqid[curevent]; // the current event sequence
        curfullsender[j + 1] = samsenders[j]; // the current event sender
        curfulltarget[j + 1] = samtargets[j]; // the current event target
        predyads.push_back(samsenders[j] + appender + samtargets[j]);


      } // the end of the for loop for searching and sampling controls

      DataFrame processi =  DataFrame::create(Rcpp::Named("time") = curtimevec,
                                              Rcpp::Named("seqeuence_id") =curseqid,
                                              Rcpp::Named("sender") = curfullsender,
                                              Rcpp::Named("receiver") = curfulltarget,
                                              Rcpp::Named("observed") = curobserved);
      processedevents[i] = processi; // storing the dataframe

    }else{

      std::vector<std::string> samsenders; // an empty vector to store the sampled senders
      std::vector<std::string> samtargets; // an empty vector to store the sampled targets
      std::vector<double> curtimevec; //a vector of length: n control + 1
      curtimevec.push_back(curtime); //the first place should be the observed event
      std::vector<double> curobserved;
      curobserved.push_back(1); // this is the real event
      std::vector<double> curseqid;
      curseqid.push_back(seqid[curevent]);
      std::vector<std::string> curfullsender;
      curfullsender.push_back(cursender); // the current event sender
      std::vector<std::string> curfulltarget;
      curfulltarget.push_back(curtarget); // the current event sender
      for(int z = 0; z < possiblesenders.size(); z++){
        for(int y = 0; y < possibletargets.size(); y++){
          if(!(possiblesenders[z] == cursender && possibletargets[y] == curtarget)){
            curfullsender.push_back(possiblesenders[z]); // adding the past z sender
            curfulltarget.push_back(possibletargets[y]); // adding the past y target
            curtimevec.push_back(curtime); //the first place should be the observed event
            curobserved.push_back(0); // this is the real event
            curseqid.push_back(seqid[curevent]);
          }
        } // end of y loop
      } // end of z loop
      DataFrame processi =  DataFrame::create(Rcpp::Named("time") = curtimevec,
                                              Rcpp::Named("seqeuence_id") =curseqid,
                                              Rcpp::Named("sender") = curfullsender,
                                              Rcpp::Named("receiver") = curfulltarget,
                                              Rcpp::Named("observed") = curobserved);

      processedevents[i] = processi; // storing the dataframe

    }
  }




  if(interval && t > 0){

    double curtime = t; // getting the current event time
    std::vector<std::string> possiblesenders;
    possiblesenders.reserve(ns); // preallocating the memory
    //searching through the possible range of actors
    for(int j = 0; j < ns; j++){
      if(timeendsenders[j] >= curtime && curtime >= timestartsender[j]){
        possiblesenders.push_back(actorsenders[j]);
      }
    }
    std::vector<std::string> possibletargets;
    possibletargets.reserve(nt); // preallocating the memory
    //searching through the possible range of actors
    for(int j = 0; j < nt; j++){
      if(timeendtargets[j] >= curtime && curtime >= timestarttarget[j]){
        possibletargets.push_back(actortargets[j]);
      }
    }

    int outer = possiblesenders.size();
    int innner = possibletargets.size();

    std::vector<std::string> samsenders(ncontrols); // an empty vector to store the sampled senders
    std::vector<std::string> samtargets(ncontrols); // an empty vector to store the sampled targets
    std::vector<std::string> predyads; // an empty vector to store the previously sampled targets



    if(outer*innner > (ncontrols + 1)){ // if there are enough controls to sample!

      CharacterVector actors1 = Rcpp::wrap(possiblesenders);
      CharacterVector actors2 = Rcpp::wrap(possibletargets);
      std::vector<double> curtimevec; //a vector of length: n control + 1
      std::vector<double> curobserved;
      std::vector<double> curseqid;
      std::vector<std::string> curfullsender;
      std::vector<std::string> curfulltarget;

      for(int j = 0; j < ncontrols; j++){
        double good = 0; // once good = 1 we accept the sampled dyad!
        while(good < 1){ // while good is still equal to zero

          CharacterVector sendsample = sample(actors1,1); //the randomly sampled sender
          CharacterVector recsample = sample(actors1,1); //the randomly sampled sender
          std::string ssender = as<std::string>(sendsample[0]); //the randomly sampled sender
          std::string starget =  as<std::string>(recsample[0]); //the randomly sampled target
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

          curtimevec[j + 1] = t; //the first place should be the observed event
          curobserved[j + 1] = 0; // this is a null event
          curseqid[j + 1] = seqid.size() + 1; // the current event sequence
          curfullsender[j + 1] = samsenders[j]; // the current event sender
          curfulltarget[j + 1] = samtargets[j]; // the current event target
          predyads.push_back(samsenders[j] + appender + samtargets[j]);

        }
      }

      DataFrame processi =  DataFrame::create(Rcpp::Named("time") = curtimevec,
                                              Rcpp::Named("seqeuence_id") =curseqid,
                                              Rcpp::Named("sender") = curfullsender,
                                              Rcpp::Named("receiver") = curfulltarget,
                                              Rcpp::Named("observed") = curobserved);

      processedevents.push_back(processi);


    }else{

      int nt = outer * innner;
      std::vector<double> curtimevec(nt); //a vector of length: n control + 1
      std::fill(curtimevec.begin(), curtimevec.end(), curtime);
      std::vector<double> curobserved(nt);
      std::fill(curobserved.begin(), curobserved.end(), 0);
      std::vector<double> curseqid(nt);
      std::fill(curseqid.begin(), curseqid.end(), seqid.size() + 1);
      std::vector<std::string> sendersRS(nt);
      std::vector<std::string> receiversRS(nt);
      int index = 0;
      //creating the empty vectors to store the sender charcter vectors
      for(int k = 0; k < outer; k++){
        for(int j = 0; j < innner; j++){
          sendersRS[index] = possiblesenders[k];
          receiversRS[index] = possibletargets[j];
          index += 1;
        }
      }
      DataFrame processi =  DataFrame::create(Rcpp::Named("time") = curtimevec,
                                              Rcpp::Named("seqeuence_id") =curseqid,
                                              Rcpp::Named("sender") = sendersRS,
                                              Rcpp::Named("receiver") = receiversRS,
                                              Rcpp::Named("observed") = curobserved);
      processedevents.push_back(processi);


    }

  }
  return processedevents;
}






