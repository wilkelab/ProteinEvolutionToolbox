#include <iostream>
#include <Rcpp.h>

#include <vector>
#include <math.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix matrix_scalar(NumericMatrix x, double window) {
  int nrow = x.nrow();

  NumericMatrix y(nrow, nrow);
  
  for(int i = 0; i < nrow; i++) {
    for(int j = 0; j < nrow; j++) {
      y(i, j) = x(i, j) * window;
    } 
  }

  return y;
}

// [[Rcpp::export]]
double Factorial(int x) {
  return (x == 1 ? x : x * Factorial(x - 1));
}

// [[Rcpp::export]]
NumericMatrix multiply_matrix(NumericMatrix x, NumericMatrix y) {
  int nrow = x.nrow();
  
  NumericMatrix z(nrow, nrow);
  fill(z.begin(), z.end(), 0);
  
  for(int i = 0; i < nrow; i++) {
    for(int j = 0; j < nrow; j++) {
      for(int k = 0; k < nrow; k++) {
        z(i, j) += x(i, k) * y(k, j);
      }
    } 
  }
  
  return z;
}

// [[Rcpp::export]]
NumericMatrix mat_exp(NumericMatrix mat) {
  int nrow = mat.nrow();
  
  NumericMatrix z(nrow, nrow);
  
  for(int k = 1; k < 100; k++) {
    NumericMatrix tmp_X(nrow, nrow);
    for(int i = 0; i < nrow; i++)
      for(int j = 0; j < nrow; j++)
        tmp_X(i, j) = mat(i, j); 
    
    NumericMatrix z_tmp(nrow, nrow);
    
    for(int i = 0; i < k; i++)
      tmp_X = multiply_matrix(tmp_X, mat);

    z_tmp = matrix_scalar(tmp_X, 1/Factorial(k));

    for(int i = 0; i < nrow*nrow; i++)
      z[i] = z[i] + z_tmp[i];
  }

  return z;
}