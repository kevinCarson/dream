#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector pibcpp(NumericMatrix net,
                      NumericVector gmem,
                      bool symmetric,
                      std::string traidtype,
                      bool count, 
                      double isolate) {
  
  double nactors = gmem.size(); // the number of actors in the network
  NumericVector broker_scores(nactors); // an empty vector to store the brokerage scores
  NumericVector broker_sts_scores(nactors); // an empty vector to store the brokerage scores
  double redactor = 0; // presetting the variable
  
  if(count){
    
      if(!symmetric){ // if the networks are not symmetric
        
        //
        //    The below code does the following:
        //      1. Loops through all possible i, j dyads in the network
        //      2. Checks if the actors are tied and if they are not, if they have 
        //         different categorical memberships (gmem)
        //      3. Then, if the 2 is true, loops through all j actors to see if they 
        //         bridge the structural hole
        //      4. Given that a j actor is a bridge, the function checks if there are 
        //          any competiting actors that also bridge the hole (z actors)
        //      5. Based on 4 and 5, the strucutral scores for the j actor are updated!
        
        //rebuild to but the triad type inside the j loop
        
        
        for(int i = 0; i < nactors; i++){ // for all i actors in the network
        
          for(int k = 0; k < nactors; k++){ // for all k actors in the network
            
            // check that a cultural-structural hole appears between i and k (see Leal 2025)
            if ((i != k && i > k) && 
                (net(i, k) == 0 && net(k, i) == 0) && 
                (gmem[k] != gmem[i])) {
            
                // searching for all potential cultural hole brokers j
                for(int j = 0; j < nactors; j++){ // for all potential brokered actors
                  
                  // disaggreagating the type of triad formation for the brokerage role
    
                  if(traidtype == "OTS"){ // if the structural hole type is outgoing shared partners
                  
                      //check if they bridge the cultural hole between i and k 
                      if ((j != k && j != i) && 
                          (net(j, k) == 1) && 
                          (net(j, i) == 1)) {
                        
                        redactor = 0; //resetting the value
                        // if they are a structural broker, then go and check for all competiting actors
                        for(int z = 0; z < nactors; z++){ // for all third party brokers
                          // for all potential third party (competiting) brokers
                          if ((z != k && z != i && z != j) && 
                              (net(z, k) == 1) && 
                              (net(z, i) == 1)) {
                            redactor += 1; // updating the value by one
                          }
                        }
                        
                        // per Leal (2025), the formula for PIB  
                        broker_scores[j] += (1/(redactor + 1)); 
                        
                      } //ending the check for a outgoing shared partner cultural broker
                      
                  } //ending the outgoing shared partner if statement
                  
                  
                  if(traidtype == "ANY"){ // if the structural hole type is any type
                    
                    //check if they bridge the cultural hole between i and k 
                    if ((j != k && j != i) && 
                        (net(j, k) == 1 || net(k, j) == 1) &&
                        (net(j, i) == 1 || net(i, j) == 1) ) {
                      
                      redactor = 0; //resetting the value
                      // if they are a structural broker, then go and check for all competiting actors
                      for(int z = 0; z < nactors; z++){ // for all third party brokers
                        // for all potential third party (competiting) brokers
                        if ((z != k && z != i && z != j) && 
                            (net(z, k) == 1 || net(k, z) == 1) &&
                            (net(z, i) == 1 || net(i, z) == 1)) {
                          redactor += 1; // updating the value by one
                        }
                      }
                      
                      // per Leal (2025), the formula for PIB  
                      broker_scores[j] += (1/(redactor + 1)); 
                      
                    } //ending the check for a any shared partner cultural broker
                    
                    
                  } //ending the any type if statement           
                  
                  
                  if(traidtype == "MTS"){ // if the structural hole type is mutual two stars
    
                    //check if they bridge the cultural hole between i and k 
                    if ((j != k && j != i) && 
                        ((  net(k, j) == 1 && net(j, i) == 1) ||
                         (net(j, k) == 1 && net(i, j) == 1) )) {
                      redactor = 0; //resetting the value
                      // if they are a structural broker, then go and check for all competiting actors
                      for(int z = 0; z < nactors; z++){ // for all third party brokers
                        // for all potential third party (competiting) brokers
                        if ((z != k && z != i && z != j) && 
                            ((  net(k, z) == 1 && net(z, i) == 1) ||
                            (net(z, k) == 1 && net(i, z) == 1) )) {
                          redactor += 1; // updating the value by one
                        }
                      }
                      
                      // per Leal (2025), the formula for PIB  
                      broker_scores[j] += (1/(redactor + 1)); 
                      
                    } //ending the check for a mutual two stars cultural broker
                  } //ending the mutual two stars if statement  
                  
                  
                  if(traidtype == "ITS"){ // if the structural hole type is incoming two stars
                    
                    //check if they bridge the cultural hole between i and k incoming two stars
                    if ((j != k && j != i) && 
                        (net(k, j) == 1 && net(i, j) == 1)) {
                      
                      redactor = 0; //resetting the value
                      // if they are a structural broker, then go and check for all competiting actors
                      for(int z = 0; z < nactors; z++){ // for all third party brokers
                        // for all potential third party (competiting) brokers
                        if ((z != k && z != i && z != j) && 
                            (net(k, z) == 1 && net(i, z) == 1)) {
                          redactor += 1; // updating the value by one
                        }
                      }
                      // per Leal (2025), the formula for PIB  
                      broker_scores[j] += (1/(redactor + 1)); 
    
                    } //ending the check for a incoming two stars cultural broker
                    
                    
                  } //ending the incoming two stars if statement  
                  
                  
                } // ending the j loop for finding potential cultural hole brokers
            
            
            } // ending the if statement for checking if the dyad is a cultural hole
            
          } // ending of the k loop for searching for dyads
        } // ending of the i loop for searching for dyads
            
            //ending if the symmetric statement     
      
      }else{  // if the network is not symmetric
        
        
        for(int i = 0; i < nactors; i++){ // for all i actors in the network
          
          for(int k = 0; k < nactors; k++){ // for all k actors in the network
            
            // check that a cultural-structural hole appears between i and k (see Leal 2025)
            if ((i != k && i > k) && 
                (net(i, k) == 0) && 
                (gmem[k] != gmem[i])) {
              
                // searching for all potential cultural hole brokers j
                for(int j = 0; j < nactors; j++){ // for all potential brokered actors
                
                if ((j != k && j != i) && 
                    (net(j, k) == 1 && net(j, i) == 1)) {
                    
                    redactor = 0; //resetting the value
                  
                  // searching for all potential cultural hole brokers j
                  for(int z = 0; z < nactors; z++){ // for all potential brokered actors
                    
                    if ((z != k && z != i && z != j) & 
                        (net(z, k) == 1 && net(z, i) == 1)){
                      redactor += 1; // updating the value by 1
                      
                    } // ending the z if statement
                    
                  } // ending the z loop for competiting actors
    
                   // per Leal (2025), the formula for PIB  
                   broker_scores[j] += (1/(redactor + 1));
                }
                
                
                }  // ending of the j loop for searching for all potential PIB bridges
              
    
            } // ending the if statement for checking if the dyad is a cultural hole
            
            
          } // ending of the k loop for searching for dyads
          
        } // ending of the i loop for searching for dyads
        
      }
  
  }
  
  
  
  
  if(!count){
    
    if(!symmetric){ // if the networks are not symmetric
      
      //
      //    The below code does the following:
      //      1. Loops through all possible i, j dyads in the network
      //      2. Checks if the actors are tied and if they are not, if they have 
      //         different categorical memberships (gmem)
      //      3. Then, if the 2 is true, loops through all j actors to see if they 
      //         bridge the structural hole
      //      4. Given that a j actor is a bridge, the function checks if there are 
      //          any competiting actors that also bridge the hole (z actors)
      //      5. Based on 4 and 5, the strucutral scores for the j actor are updated!
      
      //rebuild to but the triad type inside the j loop
      
      
      for(int i = 0; i < nactors; i++){ // for all i actors in the network
        
        for(int k = 0; k < nactors; k++){ // for all k actors in the network
          
          // check that a cultural-structural hole appears between i and k (see Leal 2025)
          if ((i != k && i > k) && 
              (net(i, k) == 0 && net(k, i) == 0)) {
            
            for(int j = 0; j < nactors; j++){ // for all j actors in the network to be a broker
              
              // disaggreagating the type of triad formation for the brokerage role
              
              if(traidtype == "OTS"){ // if the structural hole type is outgoing shared partners
                //checking if j is a structural hole broker
                if ((j != k && j != i) && 
                    (net(j, k) == 1) && 
                    (net(j, i) == 1)) {
                  broker_sts_scores[j] += 1; // updating the value by 1
                  if(gmem[k] != gmem[i]){ // checking if it is a cultural hole
                    broker_scores[j] += 1; // updating the value by 1
                  }
                } // ending the internal if statement
              } // ending the OTS if statement
              
              
              if(traidtype == "ANY"){ // if the structural hole type is ANY
                //checking if j is a structural hole broker
                if ((j != k && j != i) && 
                    (net(j, k) == 1 || net(k, j) == 1) &&
                    (net(j, i) == 1 || net(i, j) == 1) ) {
                  broker_sts_scores[j] += 1; // updating the value by 1
                  if(gmem[k] != gmem[i]){ // checking if it is a cultural hole
                    broker_scores[j] += 1; // updating the value by 1
                  }
                } // ending the internal if statement
              } // ending the ANY if statement
              
              
              if(traidtype == "ITS"){ // if the structural hole type is incoming shared partner
                //checking if j is a structural hole broker
                if ((j != k && j != i) && 
                    (net(k, j) == 1 && net(i, j) == 1)) {
                  broker_sts_scores[j] += 1; // updating the value by 1
                  if(gmem[k] != gmem[i]){ // checking if it is a cultural hole
                    broker_scores[j] += 1; // updating the value by 1
                  }
                } // ending the internal if statement
              } // ending the ITS if statement
              
              
              if(traidtype == "MTS"){ // if the structural hole type is incoming shared partner
                //checking if j is a structural hole broker
                if ((j != k && j != i) && 
                    ((  net(k, j) == 1 && net(j, i) == 1) ||
                    (net(j, k) == 1 && net(i, j) == 1) )) {
                  broker_sts_scores[j] += 1; // updating the value by 1
                  if(gmem[k] != gmem[i]){ // checking if it is a cultural hole
                    broker_scores[j] += 1; // updating the value by 1
                  }
                } // ending the internal if statement
              } // ending the ITS if statement
            } // ending the search for cultural broker js
          } // ending the if statemetn for the dyad i,k
        } // ending of the k loop for searching for dyads
      } // ending of the i loop for searching for dyads
      
      //ending if the symmetric statement     
      
    }else{  // if the network is not symmetric
      
      for(int i = 0; i < nactors; i++){ // for all i actors in the network
        
        for(int k = 0; k < nactors; k++){ // for all k actors in the network
          
          // check that a cultural-structural hole appears between i and k (see Leal 2025)
          if ((i != k && i > k) && 
              (net(i, k) == 0)) {
            
            // searching for all potential cultural hole brokers j
            for(int j = 0; j < nactors; j++){ // for all potential brokered actors
              
              if ((j != k && j != i) && 
                  (net(j, k) == 1 && net(j, i) == 1)) {
                  broker_sts_scores[j] += 1; // updating the value by 1
                  if(gmem[k] != gmem[i]){ // checking if it is a cultural hole
                    broker_scores[j] += 1; // updating the value by 1
                  }
              } // ending the if statement
            }  // ending of the j loop for searching for all potential PIB bridges
          } // ending the if statement for checking if the dyad is a cultural hole
        } // ending of the k loop for searching for dyads
      } // ending of the i loop for searching for dyads
    }
    //updating the values to be the proportions
    for(int i = 0; i < nactors; i++){ // for all i actors in the network
      if(broker_sts_scores[i] == 0){
        broker_sts_scores[i] = isolate;
      }else{
        broker_scores[i] = broker_scores[i]/broker_sts_scores[i];
      }
    }
  }
  
  return broker_scores; // returning the brokerage scores to the user
}


