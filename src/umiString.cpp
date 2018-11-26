#include <iostream>
#include <string>
#include <vector>

#include <Rcpp.h>

// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp14)]]
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::IntegerMatrix hdist(std::vector<std::string> umi) {
  
  Rcpp::IntegerMatrix distmat(umi.size(), umi.size());
  distmat.attr("dimnames") = Rcpp::List::create(umi, umi);
  
  for(int i = 0; i < umi.size(); i++)
  {
    std::string umiI = umi[i]; // learn how to use pointers to reduce memory usage!!!
    for(int j = 0; j < umi.size(); j++)
    {
      int dist = 0; // initialise new distance counter for each i vs j comparison
      std::string umiJ = umi[j];
      // length error check
      if(umiI.length() != umiJ.length())
      {
        Rcpp::Rcerr << "ERROR: UMIs must be of equal length!" << std::endl;
        distmat(i, j) = NA_INTEGER;
      }else{
        
        for(int n = 0; n < umiI.length(); n++)
        {
          if(umiI[n] != umiJ[n])
          {
            dist++;
          }
        }
      }
      distmat(i, j) = dist;
    }
  }
  return distmat;
}
