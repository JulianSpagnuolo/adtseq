
#include <string>
#include <vector>
#include <Rcpp.h>

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RSeqAn)]]
// [[Rcpp::plugins(cpp14)]]

using namespace Rcpp;


// [[Rcpp::export]]
int hdist(std::vector<std::string> seqs) {
  
  for(int i = 0; i < seqs.size(); i++)
  {
    std::string seqi = seqs[i];
    
    for(int j = 0; j < seqs.size(); j++)
    {
      int dist = 0;
      std::string seqj = seqs[j];
      
      for(int u = 0; u < seqi.length(); u++)
      {
        if(seqi[u] != seqj[u])
        {
          dist++;
        }
      }
      Rcpp::Rcout << dist << std::endl;
    }
  }
  
  return 1;
}


