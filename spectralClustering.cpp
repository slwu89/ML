#include <Rcpp.h>
using namespace Rcpp;


// Functions to generate similarity matrix S based on Gaussian RBF kernel function

// Evaluate Gaussian RBF kernel function over an entire matrix (rows are observations; columns are measured variables)
// [[Rcpp::export]]
NumericMatrix gaussian_kernCpp(NumericMatrix data, NumericVector gridA, NumericVector gridB, double sigma){
            int n_row = pow(gridA.size(),0.5);
            int n_iter = gridA.size();
            NumericMatrix output(n_row,n_row);
            for(int k=0; k<n_iter; ++k){
                int i = gridA[k];
                int j = gridB[k];
                NumericVector xi = data.row(i);
                NumericVector xj = data.row(j);
                double norm;
                norm = sum(pow((xi - xj),2));
                double ans;
                ans = norm / (2*pow(sigma,2));
                output(i,j) = exp(-(ans));
            }
            return(output);
}


// expand.grid functionality in C++
// [[Rcpp::export]]
NumericMatrix expand_gridCpp(IntegerVector vecA, IntegerVector vecB){
  int lengthA = vecA.size();
  int lengthB = vecB.size();
  int gridSize = lengthA*lengthB;
  NumericVector outA;
  NumericVector outB;
  NumericMatrix out(Dimension(gridSize,2));
  for(int i=0; i<lengthA; i++){
    for(int j=0; j<lengthB; j++){
      outA.push_back(vecA(i));
      outB.push_back(vecB(j));
    }
  }
  out(_,0) = outA;
  out(_,1) = outB;
  return(out);
}


// Generate similarity matrix with Gaussian RBF kernel in C++
// [[Rcpp::export]]
NumericMatrix sim_matRBF(NumericMatrix data, double sigma){
  //IntegerVector data_indicies = seq_len(data.nrow()) - 1;
  IntegerVector data_indicies;
  for(int i=0; i<data.nrow(); ++i){
    data_indicies.push_back(i);
  }
  NumericMatrix grid = expand_gridCpp(data_indicies,data_indicies);
  NumericVector gridA = grid(_,0);
  NumericVector gridB = grid(_,1);
  NumericMatrix out = gaussian_kernCpp(data, gridA, gridB, sigma);
  return(out);
}


// Functions to generate affinity matrix A
// A will be a matrix formed by passing a KNN filter over S

// helper function to return n highest values of a sorted vector
// [[Rcpp::export]]
NumericVector sort_n(NumericVector vector, int n_samp){
  NumericVector vector_clone = clone(vector);
  NumericVector sort_vec = vector_clone.sort();
  std::reverse(sort_vec.begin(),sort_vec.end());
  NumericVector subset_vec;
  NumericVector index;
  for(int i=0; i<n_samp; i++){
    index.push_back(i);
  }
  subset_vec = sort_vec[index];
  return(subset_vec);
}


// which function in C++ (given two vectors, produce a vector providing indicies of matching values)
// [[Rcpp::export]]
IntegerVector whichCpp(NumericVector vecA, NumericVector vecB){
  IntegerVector index;
  for(int i=0; i<vecA.size(); i++){ //iterate over smaller vector 
    for(int j=0; j<vecB.size(); j++){ //iterate over larger vector (search space) 
      if(vecB(j) == vecA(i)){
        index.push_back(j);
      }
    }
  }
  return(index);
}


// Generate affinity matrix A with KNN filter
// [[Rcpp::export]]
NumericMatrix aff_matKNN(NumericMatrix data, int neighbor){
  int n_points = data.nrow();
  NumericMatrix out(Dimension(n_points,n_points));
  for(int i=0; i<n_points; i++){
    NumericVector i_sim = data(i,_); // pull out ith row
    NumericVector n_sim = sort_n(i_sim,neighbor); // pass KNN filter over ith row

    // if 0 values, alert the user so they can adjust sigma
    if(is_true(any(n_sim == 0))){
      Rcout << "KNN filter for row " << i << " has encountered 0 value(s)" << std::endl;
      n_sim = n_sim[n_sim != 0];
    }
    
    // extract indicies of the K nearest neighbors after filtering
    IntegerVector n_index = whichCpp(n_sim,i_sim);
    
    // fill the affinity matrix A with the filtered vector(s)
    for(int j=0; j<n_index.size(); j++){
      out(i,n_index(j)) = n_sim(j);
      out(n_index(j),i) = n_sim(j);
    }
    
  }
  return(out);
}

