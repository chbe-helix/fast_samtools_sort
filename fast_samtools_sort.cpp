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

struct SamRecord_cmp {
  bool operator() (const SamRecord& a, const SamRecord& b) const {
    if(a.pos != b.pos) return a.pos < b.pos;
    return a.read_id < b.read_id; // Preserve reads
  }
};

//CB additions ############################################
// heap sort
void makeHeap(std::vector<SamRecord>& array, int i, int n){
	int largest = i;
	int l = 2*i + 1;
	int r = 2*i + 2;

	if(l < n && array[l].pos > array[largest].pos){
		largest = l;
	}

	if(r < n && array[r].pos > array[largest].pos){
		largest = r;
	}

	if(largest != i){
		std::swap(array[i], array[largest]);
		makeHeap(array, largest, n);
	}
};

void sortHeap(std::vector<SamRecord>& array){
	int n = array.size();

	for (int i = n/2; i >= 1; i--){
		makeHeap(array, i, n);
	}

	for (int i = n-1; i >= 0; i--){
		std::swap(array[0], array[i]);
		makeHeap(array, 0, i);
	}
};

// merge sort
void merge(std::vector<SamRecord>& array, std::vector<SamRecord>& array1, std::vector<SamRecord>& array2){
	array.clear();

	size_t i, j;
	for(i = 0, j = 0; i < array1.size() && j < array2.size(); ){
		if(array1[i].pos <= array2[j].pos){
			array.push_back(array1[i]);
			i++;
		} else if(array1[i].pos > array2[j].pos){
			array.push_back(array2[j]);
			j++;
		}
	}
	while(i < array1.size()){
		array.push_back(array1[i]);
		i++;
	}
	while(j < array2.size()){
		array.push_back(array2[j]);
		j++;
	}
};

void mergeSort(std::vector<SamRecord>& array){
	if(1 < array.size()){
		std::vector<SamRecord> array1(array.begin(), array.begin() + array.size() / 2);
		mergeSort(array1);
		std::vector<SamRecord> array2(array.begin() + array.size() / 2, array.end());
		mergeSort(array2);
		merge(array, array1, array2);
	}
};
//end CB additions ###########################################

static tthread::mutex thread_mutex;

struct ThreadParam {
  std::string fname_base;
  size_t* next_block;
  size_t num_block;
  std::map<std::string, size_t>* contig2pos;
  std::vector<std::string>* headers;
};

void thread_worker(void* vp) {
  const ThreadParam& threadParam = *(ThreadParam*)vp;
  char* sam = new char[opt_memory_per_thread];
  
  while(*threadParam.next_block < threadParam.num_block) {
    size_t next_block = threadParam.num_block;
    thread_mutex.lock();
    next_block = *threadParam.next_block;
    *threadParam.next_block += 1;
    thread_mutex.unlock();
    if(next_block >= threadParam.num_block) break;
  
    std::vector<SamRecord> samRecords;
    std::map<std::string, size_t>& contig2pos = *threadParam.contig2pos;
    std::vector<std::string>& headers = *threadParam.headers;
    
    // Read SAM file
    std::string cmd = "cat " + threadParam.fname_base + ".tmp" + std::to_string(next_block);
    std::shared_ptr<FILE> pipe(popen(cmd.c_str(), "r"), pclose);
    if(!pipe) throw std::runtime_error("popen() failed!");
    char buffer[2048];
    
    char* sam_cur = sam;
    while(!feof(pipe.get())) {
      buffer[0] = 0;
      if(fgets(buffer, sizeof(buffer), pipe.get()) == nullptr) break;
      if(strlen(buffer) == 0) continue;
      
      // Split and parse fields
      int field_num = 0;
      char* pch = strtok(buffer, "\t");
      std::string contig_name = "";
      while(pch != nullptr) {
	if(field_num == 2) { // chromosome or contig
	  contig_name = pch;
	} else if(field_num == 3) { // position
	  SamRecord samRecord;
	  samRecord.read_id = samRecords.size();
	  if(contig_name[0] == '*') {
	    samRecord.pos = std::numeric_limits<size_t>::max();
	  // table[table_size - 1] += strlen(sam_cur);
	  } else {
	    samRecord.pos = contig2pos[contig_name] + strtol(pch, nullptr, 10);
	    // table[(samRecord.pos / interval) + 1] += strlen(sam_cur);
	  }
	  samRecord.line = sam_cur;
	  sam_cur += (strlen(sam_cur) + 1);
	  samRecords.push_back(samRecord);
	  break;
	}
      }
      pch = strtok(NULL, "\t");
      field_num++;
    }
    
#if 0
    std::cout << "Number of sam records: " << samRecords.size() << std::endl;
    // Show the first 20 entries
    for(size_t i = 0; i < std::min<size_t>(10, samRecords.size()); i++) {
      const SamRecord& samRecord = samRecords[i];
      std::cout << "ReadID: " << samRecord.read_id << " Pos: " << samRecord.pos << "\t" << samRecord.line;
    }
#endif
    
    // this is the field we need to change #########################
    // Sort
    {
      Timer t(std::cerr, "\tSorting", opt_verbose);
      std::sort(samRecords.begin(), samRecords.end(), SamRecord_cmp());
      //std::make_heap(samRecords.begin(), samRecords.end(), SamRecord_cmp());
      //std::sort_heap(samRecords.begin(), samRecords.end(), SamRecord_cmp());
      //sortHeap(samRecords);
      //mergeSort(samRecords);
    }
    
    if(opt_verbose) {
#if 0
      // Show the first 20 entries
      std::cout << std::endl << std::endl;
      std::cout << "After sorting:" << std::endl;
      for(size_t i = 0; i < std::min<size_t>(10, samRecords.size()); i++) {
	const SamRecord& samRecord = samRecords[i];
	std::cout << "ReadID: " << samRecord.read_id << " Pos: " << samRecord.pos << "\t" << samRecord.line;
      }
#endif
    }
    
    // Write BAM file
    cmd = (opt_sambamba ? "sambamba" : "samtools");
    cmd += " view ";
    cmd += (opt_sambamba ? "--nthreads " : "--threads ");
    cmd += std::to_string(opt_threads);
    if(opt_sambamba) {
      cmd += " -f bam -S /dev/stdin -o ";
    } else {
      cmd += " -bS - > ";
    }
    cmd = cmd + threadParam.fname_base + ".tmp.sorted." + std::to_string(opt_threads);
    
    std::shared_ptr<FILE> pipe2(popen(cmd.c_str(), "w"), pclose);
    if(!pipe2) throw std::runtime_error("popen() failed!");
    for(size_t i = 0; i < headers.size(); i++) {
    fputs(headers[i].c_str(), pipe2.get());
    }
    for(size_t i = 0; i < samRecords.size(); i++) {
      const SamRecord& samRecord = samRecords[i];
      fputs(samRecord.line, pipe2.get());
    }
  }
  
  delete []sam;
}

int simple_samtools_sort(const std::string& in_fname,
			 const std::string& out_fname) {
  std::vector<SamRecord> samRecords;
  std::map<std::string, size_t> contig2pos;

  // Read BAM file
  std::string cmd = (opt_sambamba ? "sambamba" : "samtools");
  cmd += (opt_sam && opt_sambamba ? " view -h -S " : " view -h ");
  cmd += (opt_sambamba ? "--nthreads " : "--threads ");
  cmd += std::to_string(opt_threads) + " ";
  cmd += in_fname;
  std::shared_ptr<FILE> pipe(popen(cmd.c_str(), "r"), pclose);
  if(!pipe) throw std::runtime_error("popen() failed!");
  size_t size_sofar = 0;
  char buffer[2048], line[2048];

  size_t* table = nullptr;
  const size_t interval = 1000000;
  size_t table_size = 0;

  // First pass
  std::vector<std::string> headers;
  {
    Timer t(std::cerr, "\tReading BAM/SAM file: " + cmd, opt_verbose);
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
	strcpy(line, buffer);
	if(table == nullptr) {
	  table_size = (size_sofar + interval - 1) / interval + 1;
	  table = new size_t[table_size];
	  for(size_t i = 0; i < table_size; i++) {
	    table[i] = 0;
	  }
	};
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
	      char* end;
	      size_t contig_len = strtol(pch + 3, &end, 10); //was atoi(pch + 3)
	      size_sofar += contig_len;
	      break;
	    }
	  }
	} else {
	  if(field_num == 2) { // chromosome or contig
	    contig_name = pch;
	  } else if(field_num == 3) { // position
		  size_t pos = 0;
	    if(contig_name[0] == '*') {
	      pos = std::numeric_limits<size_t>::max();
	      table[table_size - 1] += strlen(line);
	    } else {
	      pos = contig2pos[contig_name] + strtol(pch, nullptr, 10);
	      table[(pos / interval)] += strlen(line);
	    }
	    break;
	  }
	}
	pch = strtok(NULL, "\t");
	field_num++;
      }
    }
  }  

  // DK - debugging purposes
  for(size_t i = 0; i < table_size; i++) {
	  std::cout << i << "M: " << table[i] << std::endl;
  }

  // Write a list of blocks in SAM
#if 0
  cmd = (opt_sambamba ? "sambamba" : "samtools");
  cmd += " view ";
  cmd += (opt_sambamba ? "--nthreads " : "--threads ");
  cmd += std::to_string(opt_threads);
  if(opt_sambamba) {
    cmd += " -f bam -S /dev/stdin -o ";
  } else {
    cmd += " -bS - > ";    
  }
  cmd += out_fname;


  {
    Timer t(std::cerr, "\tWriting into a BAM file: " + cmd, opt_verbose);
    std::shared_ptr<FILE> pipe2(popen(cmd.c_str(), "w"), pclose);
    if(!pipe2) throw std::runtime_error("popen() failed!");
    for(size_t i = 0; i < headers.size(); i++) {
      fputs(headers[i].c_str(), pipe2.get());
    }
    size_t sam2_pos = 0;
    while(true) {
      fputs(sam2 + sam2_pos, pipe2.get());
      sam2_pos = sam2_pos + strlen(sam2 + sam2_pos) + 1;
    }
  }

#endif

  // Sort blocks using multiple threads
  {
    Timer t(std::cerr, "\tSorting SAM blocks: ", opt_verbose);
    std::vector<tthread::thread*> threads(opt_threads);
    std::vector<ThreadParam> threadParams(opt_threads);
    for(size_t i = 0; i < opt_threads; i++) {
      threads[i] = new tthread::thread(thread_worker, (void*)&threadParams[i]);
    }
    
    for(size_t i = 0; i < opt_threads; i++) {
      threads[i]->join();
    }
    
    for(size_t i = 0; i < opt_threads; i++) {
      delete threads[i];
    }
  }

  // Concatenate BAM blocks
  {
    Timer t(std::cerr, "\tConcatenating BAM blocks: ", opt_verbose);
    std::string cmd = (opt_sambamba ? "sambamba" : "samtools");
    cmd += " cat -h ";
    // cmd += header_fname;
    cmd += " -o ";
    // cmd += out_fname;
    // cmd += in_fnames;
    std::shared_ptr<FILE> pipe(popen(cmd.c_str(), "r"), pclose);
    if(!pipe) throw std::runtime_error("popen() failed!");
    while(!feof(pipe.get())) {
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
  }

  {
    Timer t(std::cerr, "Overall:", opt_verbose);
    simple_samtools_sort(opt_infname,
			 opt_outfname);
  }

  return 0;
}

