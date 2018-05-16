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
#include <iomanip>
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
static size_t opt_memory_per_thread = opt_memory / opt_threads;
static size_t opt_compression = 6; // DK - we need to figure out what the default value is CB - if it's gzip compression the default is 6
static std::string opt_outfname = "";
static bool opt_verbose = false;
static bool opt_sambamba = false;
static bool opt_sam = false; // CB Edit SAM

/**
 * Use std::chrono to keep track of elapsed time between creation and
 * destruction. If verbose is true, Timer will print a message showing
 * elapsed time to the given output stream upon destruction.
 */
class Timer {
public:
  Timer(std::ostream& out = std::cerr, const std::string& msg = "", bool verbose = true) :
    _t(std::chrono::system_clock::now()), _out(out), _msg(msg), _verbose(verbose) { }
  
  /// Optionally print message
  ~Timer() {
    if(_verbose) write(_out);
  }
  
  /// Return elapsed time since Timer object was created
  double elapsed() const {
    double milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - _t).count();
    return milliseconds;
  }
  
  void write(std::ostream& out) {
    double passed = elapsed() / 1000.0;
    // Print the message supplied at construction time followed
    // by time elapsed formatted HH:MM:SS 
    int hours   = ((int)passed / 60) / 60;
    int minutes = ((int)passed / 60) % 60;
    int seconds = ((int)passed % 60);
    int milliseconds = int(passed * 1000) % 1000;
    out << _msg << " " << std::setfill ('0') << std::setw (2) << hours << ":"
	<< std::setfill ('0') << std::setw (2) << minutes << ":"
	<< std::setfill ('0') << std::setw (2) << seconds << "."
      	<< std::setfill ('0') << std::setw (3) << milliseconds << std::endl;
  }
	
private:
  std::chrono::time_point<std::chrono::system_clock> _t;
  std::ostream&  _out;
  std::string    _msg;
  bool           _verbose;
};

struct SamRecord {
  size_t read_id;
  size_t pos;
  char* line;
};

struct fileLines {
	size_t numLines = 0;
	bool bypass = 0;
};

struct table_records {
	size_t num_lines;
	size_t num_char;
};

struct SamRecord_cmp {
  bool operator() (const SamRecord& a, const SamRecord& b) const {
    if(a.pos != b.pos) return a.pos < b.pos;
    return a.read_id < b.read_id; // Preserve reads
  }
};

// Contig to position table that does use char* instead of std::string as key,
//    thus avoiding numerous memory allocations/deallocations
class Contig2Pos {
public:
  ~Contig2Pos() {
    for(auto itr = _map.begin(); itr != _map.end(); itr++) {
      assert(itr->first != nullptr);
      delete []itr->first;
    }
  }
  
  void add(const char* str, size_t pos) {
    assert(_map.find(str) == _map.end());
    char* str_copy = new char[strlen(str) + 1];
    strcpy(str_copy, str);
    _map[str_copy] = pos;
  }

  size_t operator[](const char* str) {
    assert(_map.find(str) != _map.end());
    return _map[str];
  }

  struct Contig2PosCmp {
    bool operator()(const char* str1, const char* str2) const {
      return strcmp(str1, str2) < 0;
    }
  };
	   
private:
  std::map<const char*, size_t, Contig2PosCmp> _map;
};

static tthread::mutex thread_mutex;

struct ThreadParam {
  std::string fname_base;
  size_t* next_block;
  size_t num_block;
  Contig2Pos* contig2pos;
  std::vector<std::string>* headers;
  std::vector<fileLines>* file_lines;
  std::vector<table_records>* table;

  size_t thread_id;
  size_t num_threads;
};

// CB todo clean code

void fieldSplitter(int &pass,
	std::vector<table_records> *table,
	size_t *table_size,
	std::vector<std::string> &headers,
	Contig2Pos &contig2pos,
	std::string &cmd,
	size_t *aligned_file_num,
	std::ofstream* vec_pipes = nullptr,
	std::vector<SamRecord>* samRecords = nullptr,
	char* sam_cur = nullptr){

	char buffer[2048];
	char line[2048];
	size_t size_sofar = 0;
	size_t unalign_itr = 0;
	const size_t interval = 1 << 10;

	std::shared_ptr<FILE> pipe(popen(cmd.c_str(), "r"), pclose);
	if(!pipe) throw std::runtime_error("popen() failed!");
	while(!feof(pipe.get())) {
		buffer[0] = 0;
		if(fgets(buffer, sizeof(buffer), pipe.get()) == nullptr) break;
		if(strlen(buffer) == 0) continue;

		if (pass == 2){
			if(buffer[0] == '@') continue;
			strcpy(line, buffer);
		} else if (pass == 3){
			strcpy(sam_cur, buffer);
		}

		// Is the current line header?
		bool header = (buffer[0] == '@');
		// Is the current line sequence info?
		bool sequence = false;
		if(header && pass == 1) {
			headers.push_back(buffer);
		} else if (pass == 1){
			strcpy(line, buffer);
			if(table->size() == 0) {
				*table_size = (size_sofar + interval - 1) / interval + 1;
				table->resize(*table_size);
				table->reserve(*table_size + 10); // Added 10 here to keep from reallocating memory if the table grows due to unaligned reads
				for(size_t i = 0; i < table->size(); i++) {
					(*table)[i].num_char = 0;
	  				(*table)[i].num_lines = 0;
					}
				};
	 		}

		// Split and parse fields
		int field_num = 0;
		char* pch_next = nullptr;
		char* pch = strtok_r(buffer, "\t", &pch_next);
		char* contig_name = nullptr;
		while(pch != nullptr) {
			if(header && pass == 1) {
				if(field_num == 0) {
					sequence = (strcmp(pch, "@SQ") == 0);
				}
				if(sequence) {
					if(field_num == 1) { // e.g. SN:14
						if(strlen(pch) <= 3) {
							throw std::runtime_error("");
						}
					contig_name = pch + 3;
					contig2pos.add(contig_name, size_sofar);
					} else if(field_num == 2) { // e.g. LN:107043718
						if(strlen(pch) <= 3) {
							throw std::runtime_error("");
						}
						char* end;
						size_t contig_len = strtol(pch + 3, &end, 10);
						size_sofar += contig_len;
						break;
					}
				}
			} else {
				if(field_num == 2) { // chromosome or contig
					contig_name = pch;
				} else if(field_num == 3 && pass == 1) { // position
					if(contig_name[0] == '*') {
						if(((*table)[table->size() - 1].num_char + (strlen(line) + 1)) > opt_memory_per_thread){
							table_records tbl;
							tbl.num_char = 0;
							tbl.num_lines = 0;
							table->push_back(tbl);
						}
						(*table)[table->size() - 1].num_char += (strlen(line) + 1);
						(*table)[table->size() - 1].num_lines++;
					} else {
						size_t pos = contig2pos[contig_name] + strtol(pch, nullptr, 10);
						(*table)[(pos / interval)].num_char += (strlen(line) + 1);
						(*table)[(pos / interval)].num_lines++;
					}
					break;
				} else if(field_num == 3 && pass == 2) { // position
					if(contig_name[0] == '*') {
						if((*table)[*table_size - 1 + unalign_itr].num_lines == 0){
							unalign_itr++;
						}
						vec_pipes[*aligned_file_num + unalign_itr] << line;
						(*table)[*table_size - 1 + unalign_itr].num_lines--;
					} else {
						size_t pos = contig2pos[contig_name] + strtol(pch, nullptr, 10);
						vec_pipes[(*table)[(pos / interval)].num_char] << line;
					}
					break;
				} else if(field_num == 3 && pass == 3) { // position
					SamRecord samRecord;
					samRecord.read_id = samRecords->size();
					if(contig_name[0] == '*') {
						samRecord.pos = std::numeric_limits<size_t>::max();
					} else {
						samRecord.pos = contig2pos[contig_name] + strtol(pch, nullptr, 10);
					}
					samRecord.line = sam_cur;
					sam_cur += (strlen(sam_cur) + 1);
					samRecords[(*table)[samRecord.pos / interval].num_lines].push_back(samRecord);
					break;
				}
			}
			pch = strtok_r(pch_next, "\t", &pch_next);
			field_num++;
		}
	}
};

void thread_worker(void* vp) {
  const ThreadParam& threadParam = *(ThreadParam*)vp;
  Contig2Pos& contig2pos = *threadParam.contig2pos;
  std::vector<std::string>& headers = *threadParam.headers;
  size_t thread_id = threadParam.thread_id;
  std::vector<fileLines>& arrFileLines = *threadParam.file_lines;
  std::vector<table_records>& table = *threadParam.table;
    
  char* sam = new char[opt_memory_per_thread];
  
  while(*threadParam.next_block < threadParam.num_block) {
    //
    size_t cur_block = threadParam.num_block;
    thread_mutex.lock();
    cur_block = *threadParam.next_block;
    *threadParam.next_block += 1;

    if(opt_verbose) {
      std::cerr << "Thread #" << thread_id << " is processing block #" << cur_block << "." << std::endl;
    }
    
    thread_mutex.unlock();
    if(cur_block >= threadParam.num_block) break;

    std::string in_fname = threadParam.fname_base + ".tmp." + std::to_string(cur_block);
    std::string cmd = "cat " + in_fname;

    // CB todo get bucket sort going here
    if(!arrFileLines[cur_block].bypass){
    	std::vector<SamRecord> samRecords[(arrFileLines[cur_block].numLines / 10000) + 1];
    	for(size_t itr = 0; itr < (arrFileLines[cur_block].numLines / 10000) + 1 ; itr++){
    		samRecords[itr].reserve(10000);
    	}

    	{
    		Timer t(std::cerr, "\tThread #0 reading SAM", opt_verbose && thread_id == 0);

    		// Read SAM file
    	    int pass = 3;
    	    fieldSplitter(pass,
    	    		&table,
    	    		nullptr,
    	    		headers,
    	    		contig2pos,
    	    		cmd,
    	    		nullptr,
    	    		nullptr,
    	    		samRecords,
    	    		sam);
    	}

    	remove(in_fname.c_str()); // Remove the input file

    	if(opt_verbose && thread_id == 0) {
    		#if 0
    		std::cout << "Number of sam records: " << samRecords.size() << std::endl;
    		// Show the first 10 entries
    		for(size_t i = 0; i < std::min<size_t>(10, samRecords.size()); i++) {
    			const SamRecord& samRecord = samRecords[i];
    			std::cout << "ReadID: " << samRecord.read_id << " Pos: " << samRecord.pos << "\t" << samRecord.line;
    		}
    		#endif
    		}
    	// Sort
    	{
    		Timer t(std::cerr, "\tThread #0 sorting", opt_verbose && thread_id == 0);
    		for(size_t i = 0; i < ((arrFileLines[cur_block].numLines / 10000) + 1); i++){
    			std::sort(samRecords[i].begin(), samRecords[i].end(), SamRecord_cmp());
    		}
    	}
    	if(opt_verbose && thread_id == 0) {
    		#if 0
    		// Show the first 10 entries
    		std::cout << std::endl << std::endl;
    		std::cout << "After sorting:" << std::endl;
    		for(size_t i = 0; i < std::min<size_t>(10, samRecords.size()); i++) {
    			const SamRecord& samRecord = samRecords[i];
    			std::cout << "ReadID: " << samRecord.read_id << " Pos: " << samRecord.pos << "\t" << samRecord.line;
    		}
    		#endif
    		}
    	// Write BAM file aligned
    	{
    		Timer t(std::cerr, "\tThread #0 writing into BAM", opt_verbose && thread_id == 0);
    		cmd = (opt_sambamba ? "sambamba" : "samtools");
    		cmd += " view ";
    		if(opt_sambamba) {
    			cmd += " -f bam -S /dev/stdin -o ";
    		} else {
    			cmd += " -bS - > ";
    		}
    		cmd = cmd + threadParam.fname_base + ".tmp.sorted." + std::to_string(cur_block);
    		std::shared_ptr<FILE> pipe2(popen(cmd.c_str(), "w"), pclose);
    		if(!pipe2) throw std::runtime_error("popen() failed!");
    		for(size_t i = 0; i < headers.size(); i++) {
    			fputs(headers[i].c_str(), pipe2.get());
    		}
    		for(size_t i = 0; i < ((arrFileLines[cur_block].numLines / 10000) + 1); i++){
    			for(size_t j = 0; j < samRecords[i].size(); j++) {
    				const SamRecord& samRecord = samRecords[i][j];
    				fputs(samRecord.line, pipe2.get());
    		}
    		}
    	}

    } else if(arrFileLines[cur_block].bypass)
    // Write BAM file unaligned
    {
    	{
    		Timer t(std::cerr, "\tWriting unaligned reads: ", opt_verbose && thread_id == 0);
    	    std::shared_ptr<FILE> pipe(popen(cmd.c_str(), "r"), pclose);
    	    if(!pipe) throw std::runtime_error("popen() failed!");
    	    char buffer[2048];

    	    cmd = (opt_sambamba ? "sambamba" : "samtools");
    		cmd += " view ";
    		if(opt_sambamba) {
    			cmd += " -f bam -S /dev/stdin -o ";
    		} else {
    			cmd += " -bS - > ";
    		}
    		cmd = cmd + threadParam.fname_base + ".tmp.sorted." + std::to_string(cur_block);
    		std::shared_ptr<FILE> pipe2(popen(cmd.c_str(), "w"), pclose);
    		if(!pipe2) throw std::runtime_error("popen() failed!");
    		for(size_t i = 0; i < headers.size(); i++) {
    			fputs(headers[i].c_str(), pipe2.get());
    		}
    		while(!feof(pipe.get())) {
    		   	buffer[0] = 0;
    		   	if(fgets(buffer, sizeof(buffer), pipe.get()) == nullptr) break;
    		   	if(strlen(buffer) == 0) continue;
    			fputs(buffer, pipe2.get());
    		}
    	}
    	remove(in_fname.c_str());
    }
  }
  delete []sam;
}

int fast_samtools_sort(const std::string& in_fname, const std::string& out_fname) {
  std::vector<std::string> headers;
  Contig2Pos contig2pos;
  std::vector<table_records> table;
  size_t table_size = 0;
  size_t aligned_file_num = 0;
  std::string cmd = (opt_sambamba ? "sambamba" : "samtools");
  cmd += (opt_sam && opt_sambamba ? " view -h -S " : " view -h ");
  cmd += (opt_sambamba ? "--nthreads " : "--threads ");
  cmd += std::to_string(opt_threads) + " ";
  cmd += in_fname;

  int pass = 1;
  // Read BAM file
  // First pass
  {
    Timer t(std::cerr, "\t1st pass) Reading BAM/SAM file: " + cmd, opt_verbose);

    fieldSplitter(pass,
    		&table,
    		&table_size,
    		headers,
    		contig2pos,
    		cmd,
    		&aligned_file_num
    		);
  }  

  // Determine number of files and lines per file
  size_t sam_size = 0, file_num = 1;
  fileLines lines_per_file;
  std::vector<fileLines> arrFileLines;
  assert(table_size != 0);

  for(size_t itr = 0; itr < (table_size - 1); itr++) {
    assert(sam_size <= opt_memory_per_thread);
    if(sam_size + table[itr].num_char > opt_memory_per_thread) {
      arrFileLines.push_back(lines_per_file);
      file_num++;
      sam_size = table[itr].num_char;
      lines_per_file.numLines = table[itr].num_lines;
    } else {
      sam_size += table[itr].num_char;
      lines_per_file.numLines += table[itr].num_lines;
    }
    table[itr].num_lines = (lines_per_file.numLines / 10000);
    table[itr].num_char = file_num - 1;
  }
  arrFileLines.push_back(lines_per_file);

  aligned_file_num = file_num;
  file_num += (table.size() - (table_size - 1));

  for(size_t itr = (table_size - 1); itr < table.size(); itr++){
	  lines_per_file.numLines = table[itr].num_lines;
	  lines_per_file.bypass = 1;
	  arrFileLines.push_back(lines_per_file);
  }

  // Second pass
  {
    Timer t(std::cerr, "\t2nd pass) Reading BAM/SAM file: " + cmd, opt_verbose);
    std::ofstream vec_pipes[file_num];
    for(size_t i = 0; i < file_num; i++) {
      std::string fname = in_fname + ".tmp." + std::to_string(i);
      vec_pipes[i].open(fname);
    }

    pass = 2;
    fieldSplitter(pass,
    		&table,
    		&table_size,
    		headers,
    		contig2pos,
    		cmd,
    		&aligned_file_num,
    		vec_pipes);

    for(size_t i = 0; i < file_num; i++) {
      vec_pipes[i].close();
    }
  }

  // Sort blocks using multiple threads
  size_t next_block = 0; // will need to change this number to reflect the new unaligned files
  {
    Timer t(std::cerr, "\tSorting SAM blocks: ", opt_verbose);
    std::vector<tthread::thread*> threads(opt_threads);
    std::vector<ThreadParam> threadParams(opt_threads);
    for(size_t i = 0; i < opt_threads; i++) {
      threadParams[i].fname_base  = in_fname;
      threadParams[i].next_block  = &next_block;
      threadParams[i].num_block   = file_num; // Make this aligned_file_num
      threadParams[i].contig2pos  = &contig2pos;
      threadParams[i].headers     = &headers;
      threadParams[i].thread_id   = i;
      threadParams[i].num_threads = opt_threads;
      threadParams[i].file_lines = &arrFileLines;
      threadParams[i].table = &table;
      threads[i] = new tthread::thread(thread_worker, (void*)&threadParams[i]);
    }
    
    for(size_t i = 0; i < opt_threads; i++) {
      threads[i]->join();
    }
    
    for(size_t i = 0; i < opt_threads; i++) {
      delete threads[i];
    }
  }

  // DK -> CB
  // CB todo test implementation on large BAM
  // the number of unalined reads is usually very large, perhaps over hundreds of millions
  //   and there is no need to sort unaligned reads.
  // thus, we may want to write unalinged reads directly into a BAM file using multiple threads
  // CB: I think I got it. See Above in threadworker


  // Use samtools's cat to concatenate BAM blocks
  //  Note: sambamba hasn't implemented "cat" function
  {
    std::string cmd = "samtools cat -o " + out_fname;
    std::vector<std::string> block_fnames;
    for(size_t i = 0; i < file_num; i++) {
      std::string block_fname = in_fname + ".tmp.sorted." + std::to_string(i);
      block_fnames.push_back(block_fname);
      cmd += (" " + block_fname);
    }
    Timer t(std::cerr, "\tConcatenating BAM blocks: " + cmd, opt_verbose);
    int return_value = system(cmd.c_str());
    if(return_value != 0) {
      std::cerr << "BAM concatenation failed." << "\n\t" << cmd << std::endl;
    }
    // Remove block BAM files
    for(size_t i = 0; i < block_fnames.size(); i++) {
      remove(block_fnames[i].c_str());
    }
  }

  return 0;
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
      << "  -S/--SAM        Input File format is SAM (only needed if using Sambamba)" <<std::endl // CB Edit
      << "  -@/--threads    Number of threads (Default: 1)" << std::endl
      << "  -v/--verbose    Verbose" << std::endl;
}

int main(int argc, char** argv) {
  if(argc == 1) {
    print_usage(std::cerr);
    return 0;
  }

  // Parse options
  std::set<std::string> uint_options {"-l", "-@", "--threads"};
  std::set<std::string> str_options  {"-m", "-o"};
  std::set<std::string> arg_needed_options = uint_options;
  arg_needed_options.insert(str_options.begin(), str_options.end());
  int curr_argc = 1;
  while(curr_argc < argc) {
    if(curr_argc + 1 == argc) {
      opt_infname = argv[curr_argc];
      break;
    }
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
    } else if(option == "--sambamba") {
      opt_sambamba = true;
    } else if(option == "-S" || option == "--SAM"){
    	opt_sam = true;
    } else {
      std::cerr << "Error: unrecognized option, " << option << std::endl << std::endl;
      return 0;
    }
    curr_argc++;
  }

  opt_memory_per_thread = opt_memory / opt_threads;

  // Check if the input BAM file exists.
  {
    std::ifstream f(opt_infname);
    if(!f.good()) {
      std::cerr << "Error: " << opt_infname << " does not exist." << std::endl;
      return 0;
    }
  }
  // Update the output BAM file name if it is empty.
  if(opt_outfname == "") {
    opt_outfname = opt_infname + ".sorted";
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

    std::cerr << "\tEquivalent samtools command: time samtools sort --threads " << opt_threads
	      << " -m " << opt_memory_per_thread << " " << opt_infname
	      << " -o " << opt_outfname << std::endl;
    std::cerr << "\t           sambamba command: time sambamba sort --nthreads " << opt_threads
	      << " -m " << opt_memory << " " << opt_infname
	      << " -o " << opt_outfname << std::endl;
  }

  {
    Timer t(std::cerr, "Overall:", opt_verbose);
    fast_samtools_sort(opt_infname,
		       opt_outfname);
  }

  return 0;
}

