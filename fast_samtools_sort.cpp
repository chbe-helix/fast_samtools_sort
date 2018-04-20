/*
 * Copyright 2018, Christopher Bennett <Christopher.Bennett@UTSouthwestern.edu> and Daehwan Kim <infphilo@gmail.com>
 *
 * This file is part of fast-samtools-sort.
 *
 * fast-samtools-sort is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * fast-samtools-sort is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with fast-samtools-sort.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <string.h>
#include <stdexcept>
#include <fstream>
#include <string>
#include <cassert>
#include <vector>
#include <map>
#include <algorithm>
#include <memory>

struct SamRecord {
  size_t read_id;
  size_t pos;
  std::string line;
};

struct SamRecord_cmp {
  bool operator() (const SamRecord& a, const SamRecord& b) const {
    if(a.pos != b.pos) return a.pos < b.pos;
    return a.read_id < b.read_id; // Preserve reads
  }
};

int simple_samtools_sort(const char* bam_fname) {
  std::vector<SamRecord> samRecords;
  std::map<std::string, size_t> contig2pos;
  
  // Read BAM file
  std::string cmd = "samtools view -h ";
  cmd += bam_fname;
  std::shared_ptr<FILE> pipe(popen(cmd.c_str(), "r"), pclose);
  if(!pipe) throw std::runtime_error("popen() failed!");
  size_t size_sofar = 0;
  char buffer[2048];
  std::string line;
  std::vector<std::string> headers;
  while(!feof(pipe.get())) {
    buffer[0] = 0;
    if(fgets(buffer, 2048, pipe.get()) == nullptr) break;
    if(strlen(buffer) == 0) continue;

    // Is the current line header?
    bool header = (buffer[0] == '@');
    // Is the current line sequence info?
    bool sequence = false;
    if(header) {
      headers.push_back(buffer);
    } else {
      line = buffer;
    }

    // Split and parse fields
    int field_num = 0;
    char* pch = strtok(buffer, "\t");
    std::string contig_name = "";
    while(pch != nullptr) {
      if(header) {
	if(field_num == 0) {
	  sequence = (strcmp(pch, "@SQ") == 0);
	}
	if(sequence) {
	  if(field_num == 1) { // e.g. SN:14
	    if(strlen(pch) <= 3) {
	      throw std::runtime_error("");
	    }
	    std::string contig_name = pch + 3;
	    contig2pos[contig_name] = size_sofar;
	  } else if(field_num == 2) { // e.g. LN:107043718
	    if(strlen(pch) <= 3) {
	      throw std::runtime_error("");
	    }
	    size_t contig_len = atoi(pch + 3);
	    size_sofar += contig_len;
	    break;
	  }
	}
      } else {
	if(field_num == 2) { // chromosome or contig
	  contig_name = pch;
	} else if(field_num == 3) { // position
	  SamRecord samRecord;
	  samRecord.read_id = samRecords.size();
	  if(contig_name[0] == '*') {
	    samRecord.pos = std::numeric_limits<size_t>::max();
	  } else {
	    samRecord.pos = contig2pos[contig_name] + atoi(pch);
	  }
	  samRecord.line = line;
	  samRecords.push_back(samRecord);
	  break;
	}      
      }
      pch = strtok(NULL, "\t");
      field_num++;
    }
  }


  // DK - debugging purposes
#if 1
  std::cout << "Number of sam records: " << samRecords.size() << std::endl;
  // Show the first 20 entries
  for(size_t i = 0; i < std::min<size_t>(20, samRecords.size()); i++) {
    const SamRecord& samRecord = samRecords[i];
    std::cout << "ReadID: " << samRecord.read_id << " Pos: " << samRecord.pos << "\t" << samRecord.line;
  }
#endif

  // Sort
  std::sort(samRecords.begin(), samRecords.end(), SamRecord_cmp());

  // DK - debugging purposes
#if 1
  // Show the first 20 entries
  std::cout << std::endl << std::endl;
  std::cout << "After sorting:" << std::endl;
  for(size_t i = 0; i < std::min<size_t>(20, samRecords.size()); i++) {
    const SamRecord& samRecord = samRecords[i];
    std::cout << "ReadID: " << samRecord.read_id << " Pos: " << samRecord.pos << "\t" << samRecord.line;
  }
#endif

  // Write BAM file
  cmd = "samtools view -bS - > ";
  cmd += bam_fname;
  cmd += ".sorted";

  std::shared_ptr<FILE> pipe2(popen(cmd.c_str(), "w"), pclose);
  if(!pipe2) throw std::runtime_error("popen() failed!");
  for(size_t i = 0; i < headers.size(); i++) {
    fputs(headers[i].c_str(), pipe2.get());
  }
  for(size_t i = 0; i < samRecords.size(); i++) {
    const SamRecord& samRecord = samRecords[i];
    fputs(samRecord.line.c_str(), pipe2.get());
  }

  return 0;
}


int main(int argc,char **argv) {
  simple_samtools_sort("RNA_5M.bam");
  
  return 0;
}

