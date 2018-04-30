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
#include <set>
#include <algorithm>
#include <memory>
#include <chrono>

#include "tinythread.h"

// Program options
static std::string opt_infname = "";
static size_t opt_threads = 1;
static size_t opt_memory = size_t(1) << 31; // 2 GB
static size_t opt_compression = 5; // DK - we need to figure out what the default value is
static std::string opt_outfname = "-";
static bool opt_verbose = false;

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

void thread_worker(void *vp) {
}

void thread_test() {
  std::vector<tthread::thread*> threads(opt_threads);
  std::vector<size_t> thread_ids(opt_threads);
  for(size_t i = 0; i < opt_threads; i++) {
    threads[i] = new tthread::thread(thread_worker, (void*)&thread_ids[i]);
  }
  
  for(size_t i = 0; i < opt_threads; i++) {
    threads[i]->join();
  }

  for(size_t i = 0; i < opt_threads; i++) {
    delete threads[i];
  }
}

void print_usage(std::ostream& out) {
  out << "fast-samtools-sort version " << FAST_SAMTOOLS_SORT_VERSION << " by Chris Bennett (Christopher.Bennett@UTSouthwestern.edu) and Daehwan Kim (infphilo@gmail.com)" << std::endl;
  std::string tool_name = "fast-samtools-sort";
  out << "Usage: " << std::endl
      << "  " << tool_name << " [options] [in.bam]" << std::endl
      << "Options:" << std::endl
      << "  -l INT          Compression level from 0 (no compression, fastest) to 9 (highest compression, slowest) (Default: ?)" << std::endl
      << "  -m INT[G/M/K]   Maximum memory in total, shared by threads (Default: ?G)" << std::endl
      << "  -o STR          Output filename" << std::endl
      << "  -@/--threads    Number of threads (Default: 1)" << std::endl
      << "  -v/--verbose    Verbose" << std::endl;
}

int main(int argc, char** argv) {
  if(argc == 1) {
    print_usage(std::cerr);
    return 0;
  }

  // Parse options
  opt_infname = argv[1];
  {
    std::ifstream f(opt_infname);
    if(!f.good()) {
      std::cerr << "Error: " << opt_infname << " does not exist." << std::endl;
      return 0;
    }
  }

  opt_outfname = opt_infname + ".sorted";
  
  std::set<std::string> uint_options {"-l", "-@", "--threads"};
  std::set<std::string> str_options {"-m", "-o"};
  std::set<std::string> arg_needed_options = uint_options; arg_needed_options.insert(str_options.begin(), str_options.end());
  int curr_argc = 2;
  while(curr_argc < argc) {
    std::string option = argv[curr_argc];
    std::string str_value = "";
    size_t uint_value = 0;
    if(arg_needed_options.find(option) != arg_needed_options.end()) {
      curr_argc++;
      if(curr_argc >= argc) {
	std::cerr << "Error: option, " << option << ", needs an argument." << std::endl << std::endl;
	return 0;
      }
      str_value = argv[curr_argc];
      if(uint_options.find(option) != uint_options.end()) {
	if(!std::all_of(str_value.begin(), str_value.end(), ::isdigit)) {
	  std::cerr << "Error: option, " << option << ", needs an integer argument" << option << std::endl << std::endl;
	  return 0;
	}
	uint_value = strtol(str_value.c_str(), nullptr, 10);
      }
    }
    if(option == "-l") {
      opt_compression = std::min<size_t>(9, uint_value);
    } else if(option == "-m") {
      size_t multi = 1;
      char last = str_value[str_value.length() - 1];
      if(last == 'K' || last == 'k') {
	multi = 1 << 10;
      } else if(last == 'M' || last == 'm') {
	multi = 1 << 20;
      } else if(last == 'G' || last == 'g') {
	multi = 1 << 30;
      }
      uint_value = strtol(str_value.c_str(), nullptr, 10);
      opt_memory = uint_value * multi;
    } else if(option == "-o") {
      opt_outfname = str_value;
    } else if(option == "-@" || option == "--threads") {
      opt_threads = std::max<size_t>(1, uint_value);
    } else if(option == "-v" || option == "--verbose") {
      opt_verbose = true;
    } else {
      std::cerr << "Error: unrecognized option, " << option << std::endl << std::endl;
      return 0;
    }
    curr_argc++;
  }

  if(opt_verbose) {
    size_t out_memory = opt_memory;
    char out_memory_suffix = ' ';
    if(out_memory > (1 << 14)) {
      out_memory >>= 10;
      out_memory_suffix = 'K';
      if(out_memory > (1 << 14)) {
	out_memory >>= 10;
	out_memory_suffix = 'M';
	if(out_memory > (1 << 14)) {
	  out_memory >>= 10;
	  out_memory_suffix = 'G';
	}
      }
    }
    std::cerr << "fast-samtools-sort is executed with the following options." << std::endl
	      << " " << out_memory << out_memory_suffix << " memory" << std::endl
	      << " " << opt_threads << (opt_threads == 1 ? " thread" : " threads") << std::endl;
  }

  auto program_begin = std::chrono::system_clock::now();

  simple_samtools_sort(opt_infname.c_str());

  auto program_end = std::chrono::system_clock::now();

  if(opt_verbose) {
    auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(program_end - program_begin).count();
    std::cerr << "Elapsed: " << milliseconds / 1000.0 << " seconds." << std::endl;
  }
  
  return 0;
}

