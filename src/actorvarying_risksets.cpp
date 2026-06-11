#include <Rcpp.h>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
List timevaryingomriskset(std::vector<double> time,
                          NumericVector seqid,
                          std::vector<std::string> sender,
                          std::vector<std::string> target,
                          std::vector<double> timestart,
                          std::vector<double> timeend,
                          std::vector<std::string> actors,
                          double pobserved = 1,
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
  std::string dyadupdate;
  std::string appender = "_OllieBug_";

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
    std::vector<std::string> sendersRS(outer*(outer-1));
    std::vector<std::string> receiversRS(outer*(outer-1));
    int index = 0;
    std::unordered_map<std::string, int> dyad_places;

    //creating the empty vectors to store the sender charcter vectors
    for(int k = 0; k < outer; k++){
      for(int j = 0; j < k; j++){
        sendersRS[index] = possibleactors[k];
        receiversRS[index] = possibleactors[j];
        dyadupdate = possibleactors[k] + appender + possibleactors[j];
        dyad_places[dyadupdate] = index;
        index +=1;

        sendersRS[index] = possibleactors[j];
        receiversRS[index] = possibleactors[k];
        dyadupdate = possibleactors[j] + appender + possibleactors[k];
        dyad_places[dyadupdate] = index;
        index +=1;
      }
    }
    int ncontrols = sendersRS.size();
    std::vector<double> curtimevec(ncontrols); //a vector of length: n control +
    std::fill(curtimevec.begin(), curtimevec.end(), curtime); // fill the vector with the current event time point
    std::vector<double> curobserved(ncontrols);
    std::fill(curobserved.begin(), curobserved.end(), 0); // fill the vector with the current event time point
    std::vector<double> curseqid(ncontrols);
    std::fill(curseqid.begin(), curseqid.end(),  seqid[curevent]); // fill the vector with the current event time point
    dyadupdate = cursender + appender + curtarget;
    auto here = dyad_places.find(dyadupdate);
    curobserved[here->second] = 1; //adding the value of 1 to the observed vector

    processedevents[i] =  DataFrame::create(Rcpp::Named("time") = curtimevec,
                                            Rcpp::Named("seqeuence_id") = curseqid,
                                            Rcpp::Named("sender") = sendersRS,
                                            Rcpp::Named("target") = receiversRS,
                                            Rcpp::Named("observed") = curobserved);


  }

  if(interval && t > 0){


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
    std::vector<std::string> sendersRS(outer*(outer-1));
    std::vector<std::string> receiversRS(outer*(outer-1));
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

    int ncontrols = sendersRS.size();
    std::vector<double> curtimevec(ncontrols); //a vector of length: n control +
    std::fill(curtimevec.begin(), curtimevec.end(), curtime); // fill the vector with the current event time point
    std::vector<double> curobserved(ncontrols);
    std::fill(curobserved.begin(), curobserved.end(), 0); // fill the vector with the current event time point
    std::vector<double> curseqid(ncontrols);
    std::fill(curseqid.begin(), curseqid.end(),  seqid.size() + 1); // fill the vector with the current event time point

    DataFrame processi =  DataFrame::create(Rcpp::Named("time") = curtimevec,
                                            Rcpp::Named("seqeuence_id") = curseqid,
                                            Rcpp::Named("sender") = sendersRS,
                                            Rcpp::Named("target") = receiversRS,
                                            Rcpp::Named("observed") = curobserved);


    processedevents.push_back(processi);
  }
  return processedevents;
}













// [[Rcpp::export]]
List timevaryingtmriskset(std::vector<double> time,
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
  std::string dyadupdate;
  std::string appender = "_OllieBug_";

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
    std::unordered_map<std::string, int> dyad_places;
    std::vector<std::string> sendersRS(outer*(innner));
    std::vector<std::string> receiversRS(outer*(innner));
    int index = 0;
    //creating the empty vectors to store the sender charcter vectors
    for(int k = 0; k < outer; k++){
      for(int j = 0; j < innner; j++){
        sendersRS[index] = possiblesenders[k];
        receiversRS[index] = possibletargets[j];
        dyadupdate = possiblesenders[k] + appender + possibletargets[j];
        dyad_places[dyadupdate] = index;
        index +=1;
      }
    }

    int ncontrols = sendersRS.size();
    std::vector<double> curtimevec(ncontrols); //a vector of length: n control +
    std::fill(curtimevec.begin(), curtimevec.end(), curtime); // fill the vector with the current event time point
    std::vector<double> curobserved(ncontrols);
    std::fill(curobserved.begin(), curobserved.end(), 0); // fill the vector with the current event time point
    std::vector<double> curseqid(ncontrols);
    std::fill(curseqid.begin(), curseqid.end(),  seqid[curevent]); // fill the vector with the current event time point
    dyadupdate = cursender + appender + curtarget;
    auto here = dyad_places.find(dyadupdate);
    curobserved[here->second] = 1; //adding the value of 1 to the observed vector

    processedevents[i] =  DataFrame::create(Rcpp::Named("time") = curtimevec,
                                            Rcpp::Named("seqeuence_id") = curseqid,
                                            Rcpp::Named("sender") = sendersRS,
                                            Rcpp::Named("target") = receiversRS,
                                            Rcpp::Named("observed") = curobserved);


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

    std::vector<std::string> sendersRS(outer*(innner));
    std::vector<std::string> receiversRS(outer*(innner));
    int index = 0;
    //creating the empty vectors to store the sender charcter vectors
    for(int k = 0; k < outer; k++){
      for(int j = 0; j < innner; j++){
        sendersRS[index] = possiblesenders[k];
        receiversRS[index] = possibletargets[j];
        index += 1;
      }
    }

    int ncontrols = sendersRS.size();
    std::vector<double> curtimevec(ncontrols); //a vector of length: n control +
    std::fill(curtimevec.begin(), curtimevec.end(), curtime); // fill the vector with the current event time point
    std::vector<double> curobserved(ncontrols);
    std::fill(curobserved.begin(), curobserved.end(), 0); // fill the vector with the current event time point
    std::vector<double> curseqid(ncontrols);
    std::fill(curseqid.begin(), curseqid.end(),  seqid.size() + 1); // fill the vector with the current event time point

    DataFrame processi =  DataFrame::create(Rcpp::Named("time") = curtimevec,
                                            Rcpp::Named("seqeuence_id") = curseqid,
                                            Rcpp::Named("sender") = sendersRS,
                                            Rcpp::Named("target") = receiversRS,
                                            Rcpp::Named("observed") = curobserved);


    processedevents.push_back(processi);
  }
  return processedevents;
}




