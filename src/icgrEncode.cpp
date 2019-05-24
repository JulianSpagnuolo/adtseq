#include <iostream>
#include <math.h>
#include <string>
#include <vector>

#include <Rcpp.h>

// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp14)]]
using namespace Rcpp;

// Original python implementation:
// [Yin, C. (2018). Encoding DNA sequences by integer Chaos Game Representation, Journal of Computational Biology](https://arxiv.org/abs/1712.04546)
// https://github.com/cyinbox/iCGR
//' @title Integer Chaos Game Encoder
//' @name icgrEncode
//' @author Julian Spagnuolo
//' @description This function will encode a DNA sequence as an integer chaos game representation see: 
//' Yin, C., 2018. Encoding DNA sequences by integer Chaos Game Representation, Journal of Computational Biology; https://arxiv.org/abs/1712.04546; https://github.com/cyinbox/iCGR
//' for more detail.
//'
//' @param seq, a character vector of length = 1
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix icgrEncode(std::string seq) {
  
  // implement a return matrix.
  Rcpp::IntegerMatrix results(seq.length(), 2);
  /*
  std::vector<char> rnames;
  rnames.reserve(seq.length());
  std::vector<std::string> cnames;
  cnames.push_back("x");
  cnames.push_back("y");
  */
  
  /*
  std::vector<int> pix;
  pix.reserve(seq.length());
  std::vector<int> piy;
  piy.reserve(seq.length());
  */
  
  for(int i = 0; i < seq.length(); i++)
  {
    int aix;
    int aiy;
    
    //rnames[i] = seq[i]; // get rowname for pi.
    
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
      results(0,0) = aix;
      results(0,1) = aiy;
    }else
    {
      results(i,0) = (results(i - 1, 0)) + ((2^i)*aix); // need to change 2^i-1 to 2^i (c++ is 0 indexed, original formula is not.)
      results(i,1) = (results(i - 1, 1)) + ((2^i)*aiy);
    }
  }
  
  //results.attr("dimnames") = Rcpp::List::create(rnames, cnames);
  
  return results;
}
