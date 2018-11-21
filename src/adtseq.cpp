#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <regex>
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


std::vector<std::string> adtbc(std::string sequence, std::string phredq){
  
  // init container for results
  std::vector<std::string> res;
  res.reserve(2);
  
  // set the regex patterns for finding the ADT barcode and trim off the excess polyA
  std::regex adt_bc ("(^[ATGCN]{15}[TGCN]{1}[AN]{6,}).*");
  std::regex polyA ("([TGCN]{1}[AN]{6,})(.*)");
  
  // find the barcode
  std::smatch bc_polyA;
  std::regex_match(sequence, bc_polyA, adt_bc, std::regex_constants::match_default);
  
  // trim the end off phredq to match the barcode read.
  sequence = bc_polyA[1];
  phredq = phredq.substr(0, bc_polyA.length(1));
  
  // find and remove the polyA.
  std::smatch bc_seq;
  std::regex_search(sequence, bc_seq, polyA);
  
  res.push_back(sequence.substr(0, bc_seq.position(1)));
  res.push_back(phredq.substr(0, bc_seq.position(1)));
  
  return res;
}

class adtDictionary {
public:
  
  std::vector<std::string> id;
  std::vector<std::string> seq;
  std::vector<int> counts; // add count var for future proofing - could be useful to provide prelim stats.
  int hdist; // storage for best match hamming dist
  int mapq; // storage for best match map score
  std::string match_id; //storage for best match id
  
  // dictionary constructor using fasta file of ADT barcodes
  void createDictionary(std::string fastaFile)
  {
    
    std::string line;
    std::ifstream fastaIn (fastaFile.c_str());
    
    if(fastaIn.is_open())
    {
      std::regex idRegex ("(^>).*");
      std::smatch idmatch;
      
      std::regex seqRegex ("(^[ATCGN]).*");
      std::smatch seqmatch;
      
      while ( std::getline(fastaIn, line) )
      {
        if(std::regex_match(line, idmatch, idRegex))
        {
          id.push_back (line.substr(1,line.length()));
        }
        else
        {
          if(std::regex_match(line, seqmatch, seqRegex))
          {
            seq.push_back (line);
          }
        }
      }
      fastaIn.close();
    }
  };
  
  // hamming dist function for finding matches
  void match(std::string read, std::string phredq_str, int max_dist)
  {
    // reset the result holders in the adtDictionary
    match_id = "NO_MATCH"; // no match yet.
    hdist = 15; // highest hamming distance possible
    mapq = -255; // worst possible score
    
    char read_char[read.length()];
    std::strcpy(read_char, read.c_str());
    
    for(int i = 0; i < id.size(); i++)
    {
      char adt_char[seq[i].length()];
      std::strcpy(adt_char, seq[i].c_str());
      
      int map_score = 0;
      int count = 0;
      
      // get hamming distance & map score between read and each adt barcode
      for(int u = 0; u < seq[i].length(); u++)
      {
        int score;
        int phredq;
        phredq = (int(phredq_str[u]) - 33);
        score = (2 + (4 * (phredq / 40)));
        
        if(read_char[u] != adt_char[u])
        {
          count++;
          map_score -= score;
          
          if(read_char[u] == 'N')
          {
            map_score -= 1;
          }
        }
        if(read_char[u] == adt_char[u])
        {
          map_score += score;
        }
      }
      
      if(count <= hdist & count <= max_dist & map_score >= mapq)
      {
        mapq = map_score;
        hdist = count;
        match_id = id[i];
      }
    }
  }
  
};

// [[Rcpp::export]]
int adtseq(std::string bamFileName, std::string bamOut, std::string adtFasta, int max_dist, std::string sumoutput) {
  
  
  Rcpp::Rcout << "Opening summary file stream" << std::endl;
  std::ofstream summary;
  summary.open(sumoutput);
  
  if(!summary.is_open())
  {
    Rcpp::Rcerr << "ERROR: Could not open " << sumoutput << std::endl;
    return 1;
  }
  summary << "cellbc" << "\t" << "umi" << "\t" << "antibody" << "\t" << "hdist" << "\t" << "mapq" << "\n";
  // this method loads and holds the dictionary in memory (should not be an issue since it is small compared to BAM and BAM is streamed record by record)
  adtDictionary adtDict; // initialise the barcode dictionary
  adtDict.createDictionary(adtFasta); // fill the dictionary from adtFasta
  
  
  // Open input file, BamFileIn can read SAM and BAM files.
  seqan::String<char> bamFI = bamFileName;
  seqan::BamFileIn bamFileIn;
  
  if (!open(bamFileIn, toCString(bamFI)))
  {
    Rcpp::Rcerr << "ERROR: Could not open " << bamFileName << std::endl;
    return 1;
  }
  
  seqan::String<char> bamFO = bamOut;
  seqan::BamFileOut bamFileOut;
  if(!open(bamFileOut, toCString(bamFO)))
  {
    Rcpp::Rcerr << "ERROR: Could not open " << bamOut << std::endl;
    return 1;
  }
  
  try
  {
    // Copy header.
    BamHeader header;
    readHeader(header, bamFileIn);
    writeHeader(bamFileOut, header);
    
    // Copy records.
    BamAlignmentRecord record;
    
    while (!atEnd(bamFileIn))
    {
      readRecord(record, bamFileIn);
      
      // access the current records sequence and phred quality members
      seqan::CharString current_seq;
      seqan::CharString current_q;
      current_seq = record.seq;
      current_q = record.qual;
      
      // convert them to std c++ strings for regex
      std::string pre_seq;
      std::string pre_q;
      pre_seq = toCString(current_seq);
      pre_q = toCString(current_q);
      
      /// regex matching
      std::vector<std::string> barcodes;
      barcodes = adtbc(pre_seq, pre_q);
      
      /// replace old seq with regexed and trimmed seq, trim the qual line
      record.seq = barcodes[0];
      record.qual = barcodes[1];
      current_seq = record.seq;
      
      /// qual-weighted hamming distance pattern matching. add refname to bam
      // initialise a BAM tag dictionary (seqan bs)
      BamTagsDict tags(record.tags);
      buildIndex(tags);
      
      
      seqan::CharString cellbc;
      seqan::CharString umi;
      unsigned uid; // get the uid output val.
      if(findTagKey(uid, tags, "XC"))
      {
        extractTagValue(cellbc, tags, uid);
      }else { cellbc = "NO_CBC_FOUND"; }
      if(findTagKey(uid, tags, "XM"))
      {
        extractTagValue(umi, tags, uid);
      }else { umi = "NO_UMI_FOUND"; }
      
      // find the matching ADT barcode and write the id to a new tag "XA".
      adtDict.match(toCString(barcodes[0]), toCString(barcodes[1]), max_dist);
      appendTagValue(tags, "XA", toCString(adtDict.match_id));
      appendTagValue(tags, "XN", adtDict.hdist);
      appendTagValue(tags, "AS", adtDict.mapq);
      
      record.tags = host(tags);
      
      // some verbosity as sanity check
      Rcpp::Rcout << toCString(cellbc) << "\t" << toCString(umi) << "\t" << adtDict.match_id << "\t" << adtDict.hdist << "\t" << adtDict.mapq << std::endl;
      // write the cell_bc, umi, etc to file for analysis.
      summary << toCString(cellbc) << "\t" << toCString(umi) << "\t" << adtDict.match_id << "\t" << adtDict.hdist << "\t" << adtDict.mapq << "\n";
      
      
      /// write bam out
      writeRecord(bamFileOut, record);
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
