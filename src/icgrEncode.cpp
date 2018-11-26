#include <iostream>
#include <math.h>
#include <string>
#include <vector>

#include <Rcpp.h>

// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp14)]]
using namespace Rcpp;

// [[Rcpp::export]]
int icgrEncode() {
  
  std::string seq = "CGTAACTAGT"; // same example sequence from paper and in worked example in the dev notebook
  
  std::vector<int> pix;
  pix.reserve(seq.length());
  std::vector<int> piy;
  piy.reserve(seq.length());
  
  for(int i = 0; i < seq.length(); i++)
  {
    int aix;
    int aiy;
    
    if(seq[i] == 'A'){
      aix = 1;
      aiy = 1;
    }
    else if (seq[i] == 'T'){
      aix = -1;
      aiy = 1;
    }
    else if (seq[i] == 'C'){
      aix = -1;
      aiy = -1;
    }
    else if (seq[i] == 'G'){
      aix = 1;
      aiy = -1;
    }
    else if (seq[i] == 'N'){
      aix = 0;
      aiy = 0;
      Rcpp::Rcerr << "Encountered N Base" << std::endl;
    }
    else{
      Rcpp::Rcerr << "Not a nucleobase!" << std::endl;
      break;
    }
    
    if(i == 0){
      pix[0] = aix;
      piy[0] = aiy;
    }else
    {
      pix[i] = (pix[i - 1]) + (std::pow(2, i)*aix); // need to change 2^i-1 to 2^i (c++ is 0 indexed, original formula is not.)
      piy[i] = (piy[i - 1]) + (std::pow(2, i)*aiy);
    }
  }
  
  for(int j = 0; j < seq.length(); j++)
  {
    Rcpp::Rcout << pix[j] << "\t" << piy[j] << std::endl;
  }
  
  return 0;
}
