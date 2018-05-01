# fast-samtools-sort
fast-samtools-sort is a fast BAM read sorter implemented on top of Samtools or Sambamba. This tool allows for faster sorting of BAM files than the default sorter supplied in Samtools or Sambamba.

Authors: Christopher Bennett <Christopher.Bennett@UTSouthwestern.edu> and Daehwan Kim <infphilo@gmail.com>

## Installation
### Dependencies
* C++ compiler (gcc >= 5.4.0)
* samtools >= 1.3 or sambamba >= 0.6.6 

### Install
For this example you clone fast_samtools_sort into a directory within your user directory called fast-samtools, make the executable, and add the directory to your PATH.

```sh
git clone https://github.com/chbe-helix/fast_samtools_sort.git ~/fast-samtools/
cd ~/fast-samtools/

make
echo "export=PATH=~/fast-samtools/:$PATH" >> ~/.bashrc
```

## Usage
```sh
fast-samtools-sort [-l level] [-m maxMem] [-o out.bam] [-O format] [-n] [-t tag] [-T tmpprefix] [-@ threads] [in.sam|in.bam|in.cram]
```

Sort alignments by alignment position in genome.

The sorted output is written to .bam.sorted file, or to the specified file (out.bam) when -o is used.

Options | Description
--------- | --------------------------
-l INT | Compression level from 0 (no compression, fastest) to 9 (highest compression, slowest) (Default: 6?)
-m INT[G/M/K] | Maximum memory in total, shared by threads (Default: 2G?)
-o STR | Output filename (Default: $file-name.bam.sorted)
-@/--threads INT | Number of threads to use (Default: 1)
-v/--verbose | Verbose
