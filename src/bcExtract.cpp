#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/range/algorithm/count.hpp>
#include <seqan/basic.h>
#include <seqan/bam_io.h>
#include <seqan/stream.h>
#include <seqan/sequence.h>
#include <seqan/store.h>

#include <Rcpp.h>

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RSeqAn)]]
// [[Rcpp::plugins(cpp14)]]
using namespace Rcpp;
using namespace seqan;

//' Extraction of read names from specific cell barcodes.
//' 
//' @param bamIn The input BAM/SAM filepath, this is typically output from the dropseq-tools pipeline (string).
//' @param bamOut The output BAM/SAM filepath (string).
//' @param sumoutput Filepath for tab delimited text file output of cell-barcodes, umi, fastq read name (string).
//' @param cellBC Character string specifying the 12 nt cell barcode for which you want to extract reads for.
//' @export
// [[Rcpp::export]]
int bcExtract(std::string bamIn, std::string bamOut, std::string sumoutput, std::string cellBC) {
  
  
  Rcpp::Rcout << "Opening summary file stream" << std::endl;
  std::ofstream summary;
  summary.open(sumoutput);
  
  if(!summary.is_open())
  {
    Rcpp::Rcerr << "ERROR: Could not open " << sumoutput << std::endl;
    return 1;
  }
  //summary << "cellbc" << "\t" << "umi" << "\t" << "read2name";
  summary << "cellbc" << "\t" << "umi" << "\t" << "a.rd" << "\t" << "t.rd" << "\t" << "g.rd" << "\t" << "c.rd" << "\t" << "n.rd" << "\t" << "rd.length" << "\t" << "a.bc" << "\t" << "t.bc" << "\t" << "g.bc" << "\t" << "c.bc" << "\t" << "n.bc" << "\t" << "a.umi" << "\t" << "t.umi" << "\t" << "g.umi" << "\t" << "c.umi" << "\t" << "n.umi" << "\n";
  
  // Open input file, BamFileIn can read SAM and BAM files.
  seqan::String<char> bamFI = bamIn;
  seqan::BamFileIn bamFileIn;
  
  if (!open(bamFileIn, toCString(bamFI)))
  {
    Rcpp::Rcerr << "ERROR: Could not open " << bamIn << std::endl;
    return 1;
  }
  
  /*
   seqan::String<char> bamFO = bamOut;
   seqan::BamFileOut bamFileOut;
   if(!open(bamFileOut, toCString(bamFO)))
   {
   Rcpp::Rcerr << "ERROR: Could not open " << bamOut << std::endl;
   return 1;
   }
   */
  try
  {
    // Copy header.
    BamHeader header;
    readHeader(header, bamFileIn);
    //writeHeader(bamFileOut, header);
    
    // Copy records.
    BamAlignmentRecord record;
    
    while (!atEnd(bamFileIn))
    {
      readRecord(record, bamFileIn);
      
      /// qual-weighted hamming distance pattern matching. add refname to bam
      // initialise a BAM tag dictionary (seqan bs)
      BamTagsDict tags(record.tags);
      buildIndex(tags);
      
      seqan::CharString readName = toCString(record.qName);
      seqan::CharString cellbc;
      seqan::CharString umi;
      
      seqan::CharString seq = record.seq;
      std::string seqStr = toCString(seq);
      
      // counters for bases in the read
      int a_cnt_rd = boost::count(seqStr, 'A');
      int t_cnt_rd = boost::count(seqStr, 'T');
      int g_cnt_rd = boost::count(seqStr, 'G');
      int c_cnt_rd = boost::count(seqStr, 'C');
      int n_cnt_rd = boost::count(seqStr, 'N');
      int rd_length = length(seq);
      
      // counters for bases in barcode to determine cellBC diversity
      int a_cnt_bc;
      int t_cnt_bc;
      int g_cnt_bc;
      int c_cnt_bc;
      int n_cnt_bc;
      
      // counters for bases in umi to determine diversity
      int a_cnt_umi;
      int t_cnt_umi;
      int g_cnt_umi;
      int c_cnt_umi;
      int n_cnt_umi;
      
      unsigned uid; // get the uid output val.
      if(findTagKey(uid, tags, "XC"))
      {
        extractTagValue(cellbc, tags, uid);
        std::string cellbc_str = toCString(cellbc);
        a_cnt_bc = boost::count(cellbc_str, 'A');
        t_cnt_bc = boost::count(cellbc_str, 'T');
        g_cnt_bc = boost::count(cellbc_str, 'G');
        c_cnt_bc = boost::count(cellbc_str, 'C');
        n_cnt_bc = boost::count(cellbc_str, 'N');
      }else { 
        cellbc = "NO_CBC_FOUND";
        a_cnt_bc = 0;
        t_cnt_bc = 0;
        c_cnt_bc = 0;
        g_cnt_bc = 0;
        n_cnt_bc = 0;
      }
      if(findTagKey(uid, tags, "XM"))
      {
        extractTagValue(umi, tags, uid);
        std::string umi_str = toCString(umi);
        a_cnt_umi = boost::count(umi_str, 'A');
        t_cnt_umi = boost::count(umi_str, 'T');
        g_cnt_umi = boost::count(umi_str, 'G');
        c_cnt_umi = boost::count(umi_str, 'C');
        n_cnt_umi = boost::count(umi_str, 'N');
      }else { 
        umi = "NO_UMI_FOUND";
        a_cnt_umi = 0;
        t_cnt_umi = 0;
        g_cnt_umi = 0;
        c_cnt_umi = 0;
        n_cnt_umi = 0;
      }
      
      // some verbosity as sanity check
      Rcpp::Rcout << toCString(cellbc) << "\t" << toCString(umi) << std::endl;
      // write the cell_bc, umi, etc to file for analysis.
      summary << toCString(cellbc) << "\t" << toCString(umi) << "\t" << a_cnt_rd << "\t" << t_cnt_rd << "\t" << g_cnt_rd << "\t" << c_cnt_rd << "\t" << n_cnt_rd << "\t" << rd_length << "\t" << a_cnt_bc << "\t" << t_cnt_bc << "\t" << g_cnt_bc << "\t" << c_cnt_bc << "\t" << n_cnt_bc << "\t" << a_cnt_umi << "\t" << t_cnt_umi << "\t" << g_cnt_umi << "\t" << c_cnt_umi << "\t" << n_cnt_umi << "\t" << "\n";
      //summary << toCString(cellbc) << "\t" << toCString(umi) << "\t" << toCString(readName) << "\n";
      /// write bam out
      //writeRecord(bamFileOut, record);
      
      /*
       if(cellBC.compare(toCString(cellbc)) == 0)
       {
       Rcpp::Rcout << toCString(cellbc) << "\t" << toCString(umi) << "\t" << toCString(readName) << std::endl;
       // write the cell_bc, umi, etc to file for analysis.
       summary << toCString(cellbc) << "\t" << toCString(umi) << "\t" << toCString(readName) << "\n";
       /// write bam out
       writeRecord(bamFileOut, record);
       }
       */
    }
  }
  catch (Exception const & e)
  {
    Rcpp::Rcerr << "ERROR: " << e.what() << std::endl;
    return 1;
  }
  summary.close();
  return 0;
}
