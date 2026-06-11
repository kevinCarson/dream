#include <Rcpp.h>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
List cumulativeomriskset(std::vector<double> time,
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
  std::vector<std::string> pastactors; // the past active targets at time t (for sampling)
  std::size_t pastupdateindex = 0; // the most recent past time
  std::vector<std::string> sendersRS;
  std::vector<std::string> receiversRS;
  std::unordered_set<std::string> actorset; // the past active targets at time t (for uniqueness)
  std::unordered_map<std::string, int> dyad_places;
  int npast = 0;
  std::string curnew;
  int dif;
  int search;
  std::string dyadupdate;
  std::string appender = "_OllieBug_";
  int index = 0;
  for(int i = 0; i < howmany; i++){ // for all sampled events (we start at the second event here)

    // now for each sampled event time point, we need to sample n null events (non-observed events at time t)
    //double curevent = sampledevents[i] - 1; // the current observed events
    size_t curevent = static_cast<size_t>(sampledevents[i]) - 1; // the current observed events
    std::string cursender = sender[curevent]; // the current event sender
    std::string curtarget = target[curevent]; // the current event target
    // getting the list of null event senders
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


    if(pastactors.size() > npast){
      //creating the empty vectors to store the sender charcter vectors
      dif = pastactors.size() - npast; // the number of new actors
      search = pastactors.size();
      for(int k = 0; k < dif; k++){
        curnew = pastactors[pastactors.size() - k - 1]; // to account for the zero indexing
         for(int j = 0; j < search; j++){
            if(curnew != pastactors[j]){
                sendersRS.push_back(curnew);
                receiversRS.push_back(pastactors[j]);
                dyadupdate = curnew + appender + pastactors[j];
                dyad_places[dyadupdate] = index;
                index +=1;
                sendersRS.push_back(pastactors[j]);
                receiversRS.push_back(curnew);
                dyadupdate = pastactors[j] + appender + curnew;
                dyad_places[dyadupdate] = index;
                index +=1;
            }
          }
         search -= 1; //decreasing the count by 1 to not add this actor to the other actors
        }
      npast = pastactors.size(); //updating the index value
    }

    int ncontrols = sendersRS.size();
    std::vector<double> curtimevec(ncontrols); //a vector of length: n control +
    std::fill(curtimevec.begin(), curtimevec.end(), curtimei); // fill the vector with the current event time point
    std::vector<double> curobserved(ncontrols);
    std::fill(curobserved.begin(), curobserved.end(), 0); // fill the vector with the current event time point
    dyadupdate = cursender + appender + curtarget;
    auto here = dyad_places.find(dyadupdate);
    curobserved[here->second] = 1; //adding the value of 1 to the observed vector
    std::vector<double> curseqid(ncontrols);
    std::fill(curseqid.begin(), curseqid.end(),  seqid[curevent]); // fill the vector with the current event time point
    processedevents[i] =  DataFrame::create(Rcpp::Named("time") = curtimevec,
                                            Rcpp::Named("seqeuence_id") = curseqid,
                                            Rcpp::Named("sender") = sendersRS,
                                            Rcpp::Named("target") = receiversRS,
                                            Rcpp::Named("observed") = curobserved);

  }

  if(interval && t > 0){

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
    }

    if(pastactors.size() > npast){
      //creating the empty vectors to store the sender charcter vectors
      dif = pastactors.size() - npast; // the number of new actors
      search = pastactors.size();
      for(int k = 0; k < dif; k++){
        curnew = pastactors[pastactors.size() - k - 1]; // to account for the zero indexing
        for(int j = 0; j < search; j++){
          if(curnew != pastactors[j]){
            sendersRS.push_back(curnew);
            receiversRS.push_back(pastactors[j]);
            sendersRS.push_back(pastactors[j]);
            receiversRS.push_back(curnew);
          }
        }
        search -= 1; //decreasing the count by 1 to not add this actor to the other actors
      }
      npast = pastactors.size(); //updating the index value
    }

    int ncontrols = sendersRS.size();
    std::vector<double> curtimevec(ncontrols); //a vector of length: n control +
    std::fill(curtimevec.begin(), curtimevec.end(), curtimei); // fill the vector with the current event time point
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
List cumulativetmriskset(std::vector<double> time,
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
  std::size_t pastupdateindex = 0; // the most recent past time
  std::vector<std::string> sendersRS;
  std::vector<std::string> receiversRS;
  std::unordered_set<std::string> targetset; // the past active targets at time t (for uniqueness)
  std::unordered_set<std::string> senderset; // the past active senders at time t (for uniqueness)
  std::vector<std::string> pasttargets; // the past active targets at time t (for sampling)
  std::vector<std::string> pastsenders; // the past active senders at time t (for sampling)
  int npastsenders = 0;
  int npastreceivers = 0;
  std::string curnew;
  int dif;
  int search;
  int newsender=0;
  std::string dyadupdate;
  std::string appender = "_OllieBug_";
  int index = 0;
  std::unordered_map<std::string, int> dyad_places;

  for(int i = 0; i < howmany; i++){ // for all sampled events (we start at the second event here)

    // now for each sampled event time point, we need to sample n null events (non-observed events at time t)
    //double curevent = sampledevents[i] - 1; // the current observed events
    size_t curevent = static_cast<size_t>(sampledevents[i]) - 1; // the current observed events
    std::string cursender = sender[curevent]; // the current event sender
    std::string curtarget = target[curevent]; // the current event target
    // getting the list of null event senders
    //
    //
    //          Updating the network of past events
    //
    //
    double curtimei = time[curevent]; //a vector of length: n control + 1

    if(curevent == 0){ // that is, we sampled the first event

    std::vector<double> curtimevec(1); //a vector of length: n control +
    std::fill(curtimevec.begin(), curtimevec.end(), curtimei); // fill the vector with the current event time point
    std::vector<double> curobserved(1);
    std::fill(curobserved.begin(), curobserved.end(), 1); // fill the vector with the current event time point
    std::vector<double> curseqid(1);
    std::fill(curseqid.begin(), curseqid.end(),  seqid[curevent]); // fill the vector with the current event time point

      processedevents[i] =  DataFrame::create(Rcpp::Named("time") = curtimevec,
                                              Rcpp::Named("seqeuence_id") = curseqid,
                                              Rcpp::Named("sender") = cursender,
                                              Rcpp::Named("target") = curtarget,
                                              Rcpp::Named("observed") = curobserved);

    }else{

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



    if(pastsenders.size() > npastsenders){
      //creating the empty vectors to store the sender charcter vectors
      dif = pastsenders.size() - npastsenders; // the number of new actors
      newsender = dif;
      search = pasttargets.size();
      for(int k = 0; k < dif; k++){
        curnew = pastsenders[pastsenders.size() - k - 1]; // to account for the zero indexing
        for(int j = 0; j < search; j++){
            sendersRS.push_back(curnew);
            receiversRS.push_back(pasttargets[j]);
            dyadupdate = curnew + appender + pasttargets[j];
            dyad_places[dyadupdate] = index;
            index +=1;
        }
      }
      npastsenders = pastsenders.size(); //updating the index value
    }else{
      newsender=0; // fixing back the count the next loop
    }



    if(pasttargets.size() > npastreceivers){
      //creating the empty vectors to store the sender charcter vectors
      dif = pasttargets.size() - npastreceivers; // the number of new actors
      search = pastsenders.size() - newsender;
      for(int k = 0; k < dif; k++){
        curnew = pasttargets[pasttargets.size() - k - 1]; // to account for the zero indexing
        for(int j = 0; j < search; j++){
          sendersRS.push_back(pastsenders[j]);
          receiversRS.push_back(curnew);
          dyadupdate = pastsenders[j] + appender + curnew;
          dyad_places[dyadupdate] = index;
          index +=1;
        }
      }
      npastreceivers = pasttargets.size(); //updating the index value
    }





    int ncontrols = sendersRS.size();
    std::vector<double> curtimevec(ncontrols); //a vector of length: n control +
    std::fill(curtimevec.begin(), curtimevec.end(), curtimei); // fill the vector with the current event time point
    std::vector<double> curobserved(ncontrols);
    std::fill(curobserved.begin(), curobserved.end(), 0); // fill the vector with the current event time point
    dyadupdate = cursender + appender + curtarget;
    auto here = dyad_places.find(dyadupdate);
    curobserved[here->second] = 1; //adding the value of 1 to the observed vector
    std::vector<double> curseqid(ncontrols);
    std::fill(curseqid.begin(), curseqid.end(),  seqid[curevent]); // fill the vector with the current event time point
    processedevents[i] =  DataFrame::create(Rcpp::Named("time") = curtimevec,
                                            Rcpp::Named("seqeuence_id") = curseqid,
                                            Rcpp::Named("sender") = sendersRS,
                                            Rcpp::Named("target") = receiversRS,
                                            Rcpp::Named("observed") = curobserved);

  }




  }

  if(interval && t > 0){

    double curtimei = t; //a vector of length: n control + 1
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
    }

    if(pastsenders.size() > npastsenders){
      //creating the empty vectors to store the sender charcter vectors
      dif = pastsenders.size() - npastsenders; // the number of new actors
      newsender = dif;
      search = pasttargets.size();
      for(int k = 0; k < dif; k++){
        curnew = pastsenders[pastsenders.size() - k - 1]; // to account for the zero indexing
        for(int j = 0; j < search; j++){
          sendersRS.push_back(curnew);
          receiversRS.push_back(pasttargets[j]);
        }
      }
      npastsenders = pastsenders.size(); //updating the index value
    }else{
      newsender=0; // fixing back the count the next loop
    }



    if(pasttargets.size() > npastreceivers){
      //creating the empty vectors to store the sender charcter vectors
      dif = pasttargets.size() - npastreceivers; // the number of new actors
      search = pastsenders.size() - newsender;
      for(int k = 0; k < dif; k++){
        curnew = pasttargets[pasttargets.size() - k - 1]; // to account for the zero indexing
        for(int j = 0; j < search; j++){
          sendersRS.push_back(pastsenders[j]);
          receiversRS.push_back(curnew);
        }
      }
      npastreceivers = pasttargets.size(); //updating the index value
    }



    int ncontrols = sendersRS.size();
    std::vector<double> curtimevec(ncontrols); //a vector of length: n control +
    std::fill(curtimevec.begin(), curtimevec.end(), curtimei); // fill the vector with the current event time point
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












