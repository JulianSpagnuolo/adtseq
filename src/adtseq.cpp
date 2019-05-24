#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <regex>
#include <seqan/basic.h>
#include <seqan/seq_io.h>
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


std::vector<std::string> adtbc(std::string sequence, std::string phredq, std::string adt_panel, int bc_length){
  
  // init container for results
  std::vector<std::string> res;
  res.reserve(2);
  
  // set the regex patterns for finding the ADT barcode and trim off the excess polyA
  // TODO: Add new regex rules for TotalSeq A and B panels.
  // TODO: Add regex rules for custom panels.
  std::regex adt_bc;
  std::regex polyA;
  std::stringstream custom_rule;
  std::string rule_str;
  if(adt_panel == "A")
  {
    adt_bc = ("(^[ATGCN]{15}[TGCN]{1}[AN]{6,}).*");
    polyA = ("([TGCN]{1}[AN]{6,})(.*)");
  }
  else if(adt_panel == "B")
  {
    adt_bc = ("(^[ATGCN]{15}[ATGCN]{9}(GCTTTAAGGCCGGTCCTAGCAA){1}[AN]{6,}).*");
    polyA = ("([ATGCN]{9}(GCTTTAAGGCCGGTCCTAGCAA){1}[AN]{6,})(.*)");
  }
  else if(adt_panel == "C")
  {
    adt_bc = ("(^[ATGCN]{15}[ATGCN]{9}(CCCATATAAGAAA){1}[AN]{6,}).*");
    polyA = ("([ATGCN]{9}(CCCATATAAGAAA){1}[AN]{6,})(.*)");
  }
  else if(adt_panel == "custom")
  {
    custom_rule << "(^[ATGCN]{" << bc_length << "}).";
    rule_str = custom_rule.str();
    adt_bc = (rule_str);
    polyA = ("([TGCN]{1}[AN]{6,})(.*)");
  }
  else{
    Rcpp::Rcout << "Please Define adt_panel correctly, ie A, B, C or custom" << std::endl;
  }

  
  
  // find the barcode
  if(adt_panel == "custom")
  {
    std::smatch bc_polyA;
    std::regex_match(sequence, bc_polyA, adt_bc, std::regex_constants::match_default);
    
    // trim the end off phredq to match the barcode read.
    sequence = bc_polyA[1];
    phredq = phredq.substr(0, bc_polyA.length(1));
    
    // find and remove the polyA.
    //std::smatch bc_seq;
    //std::regex_search(sequence, bc_seq, polyA);
    res.push_back(bc_polyA[1]);
    res.push_back(phredq);
  }
  else{
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
  }
  
  
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
    seqan::String<char> seqFileName = fastaFile;
    
    seqan::CharString fid;
    seqan::CharString dnaseq;
    
    seqan::SeqFileIn FaFileIn(toCString(seqFileName));
    if (!open(FaFileIn, toCString(seqFileName)))
    {
      Rcpp::Rcerr << "ERROR: Could not open " << fastaFile << std::endl;
    }
    
    try
    {
      
      while(!atEnd(FaFileIn))
      {
        readRecord(fid, dnaseq, FaFileIn);
        
        adtDictionary::id.push_back(toCString(fid));
        adtDictionary::seq.push_back(toCString(dnaseq));
      }
    }
    catch (Exception const & e)
    {
      Rcpp::Rcerr << "ERROR: " << e.what() << std::endl;
    }
     /*
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
    */
  };
  
  // hamming dist function for finding matches
  void match(std::string read, std::string phredq_str, int max_dist, int bc_length)
  {
    // reset the result holders in the adtDictionary
    match_id = "NO_MATCH"; // no match yet.
    hdist = bc_length; // highest hamming distance possible
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
//' @title Antibody Derived Tag Extraction
//' @name adtseq
//' @author Julian Spagnuolo
//' @description This function will read a processed/tagged DropSeq (not tested on 10x) BAM/SAM file, tagging each record with antibody barcode name (from the provided ADT-tag fasta file)
//' using a combination of hamming distance (less than max_dist) and best Phred quality weighted mapping score calculated as the sum of scores for each base, 
//' subtracting the score at mismatches or adding at matches. This is performed similarly to the method used by Bowtie2 to determine alignment scores.
//' At present Bowtie2’s defaults for MN (minimum score) and MX (maximum score), 2 and 6, respectively, are used.
//' The function will place the name of the identified antibody in a new tag “XA”, the mapping score under “AS”, and the hamming distance under “XN”.
//' Additionally, the function will output a second tab-delimited text file containing the Cell Barcodes, UMI and identified antibody along with the hamming distance and mapping scores for each bam/sam record.
//' The resulting BAM/SAM file can be directly used in latter stages of the dropseq-tools pipeline which attempt to correct for synthesis errors in the cell barcodes and UMIs.
//' However, should users wish, they may retrieve count data directly from the text file and omit the latter steps (not recommended).
//'
//' @param bamFileName, An processed BAM or SAM file that has been tagged with the cell barcode and UMI using the Dropseq-tools TagBamWithReadSequenceExtended.
//' @param bamOut, output BAM/SAM file, will autodetect format by the suffix i.e. *.bam for BAM, *.sam for SAM
//' @param adtFasta, a fasta file containing the sequences for the antibody derived tags
//' @param max_dist, Integer, the maximum allowed hamming distance a read can be from the expected antibody derived tag sequence in the fasta file
//' @param sumoutput, output summary txt file
//' @param adt_panel, character indicating which TotalSeq panel (i.e. one of "A","B" or "C") is being used or "custom" if custom home-made barcodes are used
//' @param bc_length, character indicating how long the antibody barcodes are (only required if adt_panel == "custom").
//' @export
// [[Rcpp::export]]
int adtseq(std::string bamFileName, std::string bamOut, std::string adtFasta, int max_dist, std::string sumoutput, std::string adt_panel = "A", int bc_length = 15) {
  
  
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
      if(adt_panel != "custom")
      {
        barcodes = adtbc(pre_seq, pre_q, adt_panel, 15);
      }else{
        barcodes = adtbc(pre_seq, pre_q, adt_panel, bc_length);
      }
      
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
      adtDict.match(toCString(barcodes[0]), toCString(barcodes[1]), max_dist, bc_length);
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
