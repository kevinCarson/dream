#include <Rcpp.h>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
List processREMseqOM_varying(std::vector<double> time,
                             NumericVector seqid,
                             std::vector<std::string> sender,
                             std::vector<std::string> target,
                             double pobserved = 1,
                             double ncontrols = 1,
                             std::string appender = "__NIKOACAR3718__",
                             bool interval = false,
                             double t = 0) {

  double nevents = time.size(); // the number of observed evetns in the event seqeucne
  double howmany = std::round(nevents*pobserved); // this should return the number of events we need to sample
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

  List processedevents(howmany); //each list element will be a data frame! with 1 + ncontrols rows
  std::unordered_set<std::string> actorset; // the past active targets at time t (for uniqueness)
  std::vector<std::string> pastactors; // the past active targets at time t (for sampling)
  std::size_t pastupdateindex = 0; // the most recent past time
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
    std::vector<double> curtime(2); //a vector of length: n control + 1
    curtime[0] = time[curevent]; //the first place should be the observed event
    curtime[1] = time[curevent]; //the first place should be the observed event

    std::vector<double> curobserved(2);
    curobserved[0] = 1; // this is the real event
    curobserved[1] = 0; // this is the real event

    std::vector<double> curseqid(2);
    curseqid[0] = seqid[curevent]; // the current event sequence
    curseqid[1] = seqid[curevent]; // the current event sequence

    std::vector<std::string> curfullsender(2);
    curfullsender[0] = cursender; // the current event sender
    curfullsender[1] = curtarget; // the current event sender

    std::vector<std::string> curfulltarget(2);
    curfulltarget[0] = curtarget; // the current event target
    curfulltarget[1] = cursender; // the current event target


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
        if(actorset.insert(target[j]).second){ // if the current target can be inserted (i.e., not already in the set)
          pastactors.push_back(target[j]); // adding the current target
        }
        if(actorset.insert(sender[j]).second){ // if the current sender can be inserted (i.e., not already in the set)
          pastactors.push_back(sender[j]); // adding the current target
        }
      }else{
        pastupdateindex = j; // the most recent time that is beyond the current time point
        break; // escaping the loop as time is beyond ti
      }
    } // ending the searching loop to add nodes to their respective lists!

    if((pastactors.size() * (pastactors.size()-1)) > (ncontrols + 1)){ // if there are enough nodes for sampling!
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
      CharacterVector actors1 = Rcpp::wrap(pastactors);
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
      for(int z = 0; z < pastactors.size(); z++){
        for(int y = 0; y < pastactors.size(); y++){
          if(!(pastactors[z] == cursender && pastactors[y] == curtarget) &&
             (z != y)){
            curfullsender.push_back(pastactors[z]); // adding the past z sender
            curfulltarget.push_back(pastactors[y]); // adding the past y target
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
    for(int j = pastupdateindex; j < nevents; j++){ // for all events until the first one
      if(time[j] <= curtimei){ // if the time j is less or equal to the current event time
        if(actorset.insert(target[j]).second){ // if the current target can be inserted (i.e., not already in the set)
          pastactors.push_back(target[j]); // adding the current target
        }
        if(actorset.insert(sender[j]).second){ // if the current sender can be inserted (i.e., not already in the set)
          pastactors.push_back(sender[j]); // adding the current target
        }
      }else{
        pastupdateindex = j; // the most recent time that is beyond the current time point
        break; // escaping the loop as time is beyond ti
      }
    } // ending the searching loop to add nodes to their respective lists!


    if((pastactors.size() * (pastactors.size()-1)) > (ncontrols + 1)){ // if there are enough nodes for sampling!
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
      CharacterVector actors1 = Rcpp::wrap(pastactors);


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
      for(int z = 0; z < pastactors.size(); z++){
        for(int y = 0; y < pastactors.size(); y++){
          if(!(pastactors[z] == cursender && pastactors[y] == curtarget) &&
             (z != y)){
            curfullsender.push_back(pastactors[z]); // adding the past z sender
            curfulltarget.push_back(pastactors[y]); // adding the past y target
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



  if(interval && t > 0){


    // now for each sampled event time point, we need to sample n null events (non-observed events at time t)
    //double curevent = sampledevents[i] - 1; // the current observed events
    // getting the list of null event senders
    std::vector<std::string> samsenders(ncontrols); // an empty vector to store the sampled senders
    std::vector<std::string> samtargets(ncontrols); // an empty vector to store the sampled targets
    std::vector<std::string> predyads; // an empty vector to store the previously sampled targets
    //
    //
    //          Updating the network of past events
    //
    //
    double curtimei = t; //a vector of length: n control + 1
    for(int j = pastupdateindex; j < nevents; j++){ // for all events until the first one
      if(time[j] <= curtimei){ // if the time j is less or equal to the current event time
        if(actorset.insert(target[j]).second){ // if the current target can be inserted (i.e., not already in the set)
          pastactors.push_back(target[j]); // adding the current target
        }
        if(actorset.insert(sender[j]).second){ // if the current sender can be inserted (i.e., not already in the set)
          pastactors.push_back(sender[j]); // adding the current target
        }
      }else{
        pastupdateindex = j; // the most recent time that is beyond the current time point
        break; // escaping the loop as time is beyond ti
      }
    } // ending the searching loop to add nodes to their respective lists!

    double curevent = seqid.size() + 1; // updating the count by 1
    std::vector<double> curtime(ncontrols); //a vector of length: n control + 1
    std::vector<double> curobserved(ncontrols);
    std::vector<double> curseqid(ncontrols);
    std::vector<std::string> curfullsender(ncontrols);
    std::vector<std::string> curfulltarget(ncontrols);
    CharacterVector actors1 = Rcpp::wrap(pastactors);
    for(int j = 0; j < ncontrols; j++){
      double good = 0; // once good = 1 we accept the sampled dyad!
      while(good < 1){ // while good is still equal to zero

        CharacterVector sendsample = sample(actors1,1); //the randomly sampled sender
        CharacterVector recsample = sample(actors1,1); //the randomly sampled sender
        std::string ssender = as<std::string>(sendsample[0]); //the randomly sampled sender
        std::string starget =  as<std::string>(recsample[0]); //the randomly sampled target
        if(ssender != starget){ // if it is the current dyad, sample a new dyad, or if this is a self-loop
         // if it is not the current dyad, check to make sure it is not an already sampled dyad
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


      curtime[j] = t; //the first place should be the observed event
      curobserved[j] = 0; // this is a null event
      curseqid[j] = curevent; // the current event sequence
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
                                            Rcpp::Named("receiver") = curfulltarget,
                                            Rcpp::Named("observed") = curobserved);

    processedevents.push_back(processi);
  }


  return(processedevents);
}
