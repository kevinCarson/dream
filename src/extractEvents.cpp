#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<NumericMatrix> extractEventData(
    NumericMatrix stats,
    NumericVector outcome,
    IntegerVector event_cluster,
    CharacterVector names) {

  // a function to extract the network statistics per event cluster

  int N = stats.nrow(); // number of events in the sequence
  int K = stats.ncol();  // number of covariates requested by the user

  // 1. Stratify by event clusters
  std::map<int, std::vector<int>> strata;
  for (int i = 0; i < N; i++) { // for all events in the sequence
    strata[event_cluster[i]].push_back(i); // adding the respective place
  }

  // 2. Preallocate result vector of matrices
  std::vector<NumericMatrix> result; // a matrix to store the ith event cluster stats
  result.reserve(strata.size());

  // 3. Build each matrix: [net_stats | dummy | id]
  for (const auto& entry : strata) { // for each event sequence element
    const std::vector<int>& indices = entry.second; // the rows indices for the cluster
    int n_rows = indices.size(); // the number of events in this cluster
    NumericMatrix mat(n_rows, K + 2);  // matrix with n rows and k+2 columns

    for (int i = 0; i < n_rows; i++) { // for each row in the specific event cluster
      int idx = indices[i]; // getting the current index
      for (int j = 0; j < K; j++) { // for all parameters requested
        mat(i, j) = stats(idx, j); // extract the kth values for the ith person
      }
      mat(i, K) = outcome[idx];  // the outcome variable for the ith person
      mat(i, K + 1) = event_cluster[idx];     // the cluster id for the ith person
    }
    colnames(mat) = names;

    result.push_back(mat); // adding the matrix to the list!
  }

  return result;
}

// [[Rcpp::export]]
NumericVector checkVarianceData(List event_stats, // the event statistics matrix
                                int K){ // the number of estimated parameters
  int N_events = event_stats.size(); // the numeric of variables
  NumericVector varianceEst(K); // an empty vector of length k to store the variance estimates
  for (int i = 0; i < N_events; i++) { // for each event sequence cluster
    NumericMatrix statsi = event_stats[i]; //the current dataset
    for (int j = 0; j < K; j++) { // for each variable requested
      varianceEst[j] =  varianceEst[j] + var(statsi(_,j)); // adding the estimated variance to the vector
    }
  }
  return varianceEst; // returning the variance estimate
}





