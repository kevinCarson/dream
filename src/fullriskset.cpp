#include <Rcpp.h>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
List fullrisksetom(std::vector<double> time,
                   NumericVector seqid,
                   std::vector<std::string> sender,
                   std::vector<std::string> target,
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

 //first, we need to get the unique actors in each set
  sender.insert(sender.end(), target.begin(), target.end());
  std::unordered_set<std::string> actors1(sender.begin(), sender.end());
  std::vector<std::string> actors(actors1.begin(), actors1.end());
  int outer = actors.size();
  std::vector<std::string> sendersRS(outer*(outer-1));
  std::vector<std::string> receiversRS(outer*(outer-1));
  int index = 0;
  std::string dyadupdate;
  std::string appender = "_OllieBug_";
  std::unordered_map<std::string, int> dyad_places;

  //creating the empty vectors to store the sender charcter vectors
  for(int i = 0; i < outer; i++){
   for(int j = 0; j < i; j++){
     sendersRS[index] = actors[i];
     receiversRS[index] = actors[j];
     dyadupdate = actors[i] + appender + actors[j];
     dyad_places[dyadupdate] = index;
     index += 1;
     sendersRS[index] = actors[j];
     receiversRS[index] = actors[i];
     dyadupdate = actors[j] + appender + actors[i];
     dyad_places[dyadupdate] = index;
     index += 1;
    }
  }

  double ncontrols = sendersRS.size();//the number of possible controls

  for(int i = 0; i < howmany; i++){ // for all sampled events

    // now for each sampled event time point, we need to sample n null events (non-observed events at time t)
    double curevent = sampledevents[i] - 1; // the current observed events

    std::string cursender = sender[curevent];
    std::string curtarget = target[curevent];
    std::vector<double> curtime(ncontrols); //a vector of length: n control + 1
    std::fill(curtime.begin(), curtime.end(), time[curevent]); // fill the vector with the current event time point
    std::vector<double> curobserved(ncontrols);
    std::fill(curobserved.begin(), curobserved.end(), 0); // fill the vector with the current event time point
    std::vector<double> curseqid(ncontrols);
    std::fill(curseqid.begin(), curseqid.end(),  seqid[curevent]); // fill the vector with the current event time point
    dyadupdate = cursender + appender + curtarget;
    auto here = dyad_places.find(dyadupdate);
    curobserved[here->second] = 1; //adding the value of 1 to the observed vector

    //storing the full resulting list that is:
    // the current dyad
    // the sampled dyad
    // the dummy vector
    // the time vector
    // creating a dataframe to store the results
    DataFrame processi =  DataFrame::create(Rcpp::Named("time") = curtime,
                                            Rcpp::Named("seqeuence_id") = curseqid,
                                            Rcpp::Named("sender") = sendersRS,
                                            Rcpp::Named("target") = receiversRS,
                                            Rcpp::Named("observed") = curobserved);
    processedevents[i] = processi; // storing the dataframe

  }

  if(interval && t > 0){
    std::vector<double> curtime(ncontrols); //a vector of length: n control + 1
    std::fill(curtime.begin(), curtime.end(), t); // fill the vector with the current event time point
    std::vector<double> curobserved(ncontrols);
    std::fill(curobserved.begin(), curobserved.end(), 0); // fill the vector with the current event time point
    std::vector<double> curseqid(ncontrols);
    std::fill(curseqid.begin(), curseqid.end(),  seqid.size() + 1); // fill the vector with the current event time point
    DataFrame processi =  DataFrame::create(Rcpp::Named("time") = curtime,
                                            Rcpp::Named("seqeuence_id") = curseqid,
                                            Rcpp::Named("sender") = sendersRS,
                                            Rcpp::Named("target") = receiversRS,
                                            Rcpp::Named("observed") = curobserved);
    processedevents.push_back(processi);
  }


  return processedevents;
}







// [[Rcpp::export]]
List fullrisksettm(std::vector<double> time,
                        NumericVector seqid,
                        std::vector<std::string> sender,
                        std::vector<std::string> target,
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
  //first, we need to get the unique actors in each set
  std::unordered_set<std::string> senders1(sender.begin(), sender.end());
  std::unordered_set<std::string> receivers1(target.begin(), target.end());
  std::vector<std::string> senders(senders1.begin(), senders1.end());
  std::vector<std::string> receivers(receivers1.begin(), receivers1.end());
  int outer = senders.size();
  int inner = receivers.size();
  std::vector<std::string> sendersRS(outer*inner);
  std::vector<std::string> receiversRS(outer*inner);
  int index = 0;
  std::string dyadupdate;
  std::string appender = "_OllieBug_";
  std::unordered_map<std::string, int> dyad_places;
  //creating the empty vectors to store the sender charcter vectors
  for(int i = 0; i < outer; i++){
    for(int j = 0; j < inner; j++){
      sendersRS[index] = senders[i];
      receiversRS[index] = receivers[j];
      dyadupdate = senders[i] + appender + receivers[j];
      dyad_places[dyadupdate] = index;
      index += 1;
    }
  }
  double ncontrols = sendersRS.size();//the number of possible controls

  for(int i = 0; i < howmany; i++){ // for all sampled events

    // now for each sampled event time point, we need to sample n null events (non-observed events at time t)
    double curevent = sampledevents[i] - 1; // the current observed events

    std::string cursender = sender[curevent];
    std::string curtarget = target[curevent];
    std::vector<double> curtime(ncontrols); //a vector of length: n control + 1
    std::fill(curtime.begin(), curtime.end(), time[curevent]); // fill the vector with the current event time point
    std::vector<double> curobserved(ncontrols);
    std::fill(curobserved.begin(), curobserved.end(), 0); // fill the vector with the current event time point
    std::vector<double> curseqid(ncontrols);
    std::fill(curseqid.begin(), curseqid.end(),  seqid[curevent]); // fill the vector with the current event time point
    dyadupdate = cursender + appender + curtarget;
    auto here = dyad_places.find(dyadupdate);
    curobserved[here->second] = 1; //adding the value of 1 to the observed vector

    //storing the full resulting list that is:
    // the current dyad
    // the sampled dyad
    // the dummy vector
    // the time vector
    // creating a dataframe to store the results
    DataFrame processi =  DataFrame::create(Rcpp::Named("time") = curtime,
                                            Rcpp::Named("seqeuence_id") = curseqid,
                                            Rcpp::Named("sender") = sendersRS,
                                            Rcpp::Named("target") = receiversRS,
                                            Rcpp::Named("observed") = curobserved);
    processedevents[i] = processi; // storing the dataframe

  }

  if(interval && t > 0){
    std::vector<double> curtime(ncontrols); //a vector of length: n control + 1
    std::fill(curtime.begin(), curtime.end(), t); // fill the vector with the current event time point
    std::vector<double> curobserved(ncontrols);
    std::fill(curobserved.begin(), curobserved.end(), 0); // fill the vector with the current event time point
    std::vector<double> curseqid(ncontrols);
    std::fill(curseqid.begin(), curseqid.end(),  seqid.size() + 1); // fill the vector with the current event time point
    DataFrame processi =  DataFrame::create(Rcpp::Named("time") = curtime,
                                            Rcpp::Named("seqeuence_id") = curseqid,
                                            Rcpp::Named("sender") = sendersRS,
                                            Rcpp::Named("target") = receiversRS,
                                            Rcpp::Named("observed") = curobserved);
    processedevents.push_back(processi);
  }


  return processedevents;

}










