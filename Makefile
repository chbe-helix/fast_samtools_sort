#
# Copyright 2018, Christopher Bennett <Christopher.Bennett@UTSouthwestern.edu> and Daehwan Kim <infphilo@gmail.com>
#
# This file is part of fast-samtools-sort.
#
# CMA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CMA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT.  If not, see <http://www.gnu.org/licenses/>.
#
#
# Makefile for fast-samtools-sort
#

INC =
GCC_PREFIX = $(shell dirname `which gcc`)
GCC_SUFFIX =
CC = $(GCC_PREFIX)/gcc$(GCC_SUFFIX)
CPP = $(GCC_PREFIX)/g++$(GCC_SUFFIX)
CXX = $(CPP)
HEADERS = $(wildcard *.h)
FAST_SAMTOOLS_SORT_MM = 1

LINUX = 1

# Detect Cygwin or MinGW
WINDOWS = 0
CYGWIN = 0
MINGW = 0
ifneq (,$(findstring CYGWIN,$(shell uname)))
	WINDOWS = 1 
	CYGWIN = 1
	LINUX = 0
else
	ifneq (,$(findstring MINGW,$(shell uname)))
		WINDOWS = 1
		MINGW = 1
		LINUX = 0
	endif
endif

MACOS = 0
ifneq (,$(findstring Darwin,$(shell uname)))
	MACOS = 1
	LINUX = 0
endif

PTHREAD_PKG =
PTHREAD_LIB = 

ifeq (1,$(MINGW))
	PTHREAD_LIB = 
else
	PTHREAD_LIB = -lpthread
endif

LIBS = $(PTHREAD_LIB)

SHARED_CPPS = tinythread.cpp

VERSION = $(shell cat VERSION)

# Convert BITS=?? to a -m flag
BITS=32
ifeq (x86_64,$(shell uname -m))
BITS=64
endif
# msys will always be 32 bit so look at the cpu arch instead.
ifneq (,$(findstring AMD64,$(PROCESSOR_ARCHITEW6432)))
	ifeq (1,$(MINGW))
		BITS=64
	endif
endif
BITS_FLAG =

ifeq (32,$(BITS))
	BITS_FLAG = -m32
endif

ifeq (64,$(BITS))
	BITS_FLAG = -m64
endif
SSE_FLAG=-msse2

DEBUG_FLAGS    = -O0 -g3 $(BIToS_FLAG) $(SSE_FLAG) -std=c++11
DEBUG_DEFS     = -DCOMPILER_OPTIONS="\"$(DEBUG_FLAGS) $(EXTRA_FLAGS)\""
RELEASE_FLAGS  = -O3 $(BITS_FLAG) $(SSE_FLAG) -funroll-loops -g3 -std=c++11
RELEASE_DEFS   = -DCOMPILER_OPTIONS="\"$(RELEASE_FLAGS) $(EXTRA_FLAGS)\""
NOASSERT_FLAGS = -DNDEBUG
FILE_FLAGS     = -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE

CMA_BIN_LIST = fast-samtools-sort
CMA_BIN_LIST_AUX = fast-samtools-sort-debug

GENERAL_LIST = $(wildcard scripts/*.sh) \
	$(wildcard scripts/*.pl) \
	$(wildcard *.py) \
	$(wildcard hisatgenotype_modules/*.py) \
	$(wildcard hisatgenotype_scripts/*.py) \
	doc/manual.inc.html \
	doc/README \
	doc/style.css \
	$(wildcard example/index/*.ht2) \
	$(wildcard example/reads/*.fa) \
	example/reference/22_20-21M.fa \
	example/reference/22_20-21M.snp \
	$(PTHREAD_PKG) \
	cma \
	cma-bin \
	AUTHORS \
	LICENSE \
	NEWS \
	MANUAL \
	MANUAL.markdown \
	TUTORIAL \
	VERSION

ifeq (1,$(WINDOWS))
	CMA_BIN_LIST := $(CMA_BIN_LIST) hisat2.bat hisat2-build.bat hisat2-inspect.bat 
endif

# This is helpful on Windows under MinGW/MSYS, where Make might go for
# the Windows FIND tool instead.
FIND=$(shell which find)

SRC_PKG_LIST = $(wildcard *.h) \
	$(wildcard *.hh) \
	$(wildcard *.c) \
	$(wildcard *.cpp) \
	doc/strip_markdown.pl \
	Makefile \
	$(GENERAL_LIST)

BIN_PKG_LIST = $(GENERAL_LIST)

.PHONY: all allall both both-debug

all: $(CMA_BIN_LIST)

allall: $(CMA_BIN_LIST) $(CMA_BIN_LIST_AUX)

both: fast-samtools-sort

both-debug: fast-samtools-sort-debug

DEFS=-fno-strict-aliasing \
     -DFAST_SAMTOOLS_SORT_VERSION="\"`cat VERSION`\"" \
     -DBUILD_HOST="\"`hostname`\"" \
     -DBUILD_TIME="\"`date`\"" \
     -DCOMPILER_VERSION="\"`$(CXX) -v 2>&1 | tail -1`\"" \
     $(FILE_FLAGS) \
     $(PREF_DEF) \
     $(MM_DEF) \
     $(SHMEM_DEF)

fast-samtools-sort: fast_samtools_sort.cpp $(SHARED_CPPS) $(HEADERS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) $(NOASSERT_FLAGS) -Wall \
	$(INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(HISAT_CPPS_MAIN) \
	$(LIBS) $(SEARCH_LIBS)

fast-samtools-sort-debug: fast_samtools_sort.cpp $(SHARED_CPPS) $(HEADERS)
	$(CXX) $(DEBUG_FLAGS) \
	$(DEBUG_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -Wall \
	$(INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(HISAT_CPPS_MAIN) \
	$(LIBS) $(SEARCH_LIBS)


cma: ;

cma.bat:
	echo "@echo off" > hisat2.bat
	echo "perl %~dp0/cma %*" >> hisat2.bat

hisat2-build.bat:
	echo "@echo off" > hisat2-build.bat
	echo "python %~dp0/hisat2-build %*" >> hisat2-build.bat

hisat2-inspect.bat:
	echo "@echo off" > hisat2-inspect.bat
	echo "python %~dp0/hisat2-inspect %*" >> hisat2-inspect.bat


.PHONY: cma-src
cma-src-pkg: $(SRC_PKG_LIST)
	chmod a+x scripts/*.sh scripts/*.pl
	mkdir .src.tmp
	mkdir .src.tmp/hisat2-$(VERSION)
	zip tmp.zip $(SRC_PKG_LIST)
	mv tmp.zip .src.tmp/hisat2-$(VERSION)
	cd .src.tmp/hisat2-$(VERSION) ; unzip tmp.zip ; rm -f tmp.zip
	cd .src.tmp ; zip -r hisat2-$(VERSION)-source.zip hisat2-$(VERSION)
	cp .src.tmp/hisat2-$(VERSION)-source.zip .
	rm -rf .src.tmp

.PHONY: cma-bin
cma-bin-pkg: $(BIN_PKG_LIST) $(HISAT2_BIN_LIST) $(HISAT2_BIN_LIST_AUX)
	chmod a+x scripts/*.sh scripts/*.pl
	rm -rf .bin.tmp
	mkdir .bin.tmp
	mkdir .bin.tmp/hisat2-$(VERSION)
	if [ -f hisat2.exe ] ; then \
		zip tmp.zip $(BIN_PKG_LIST) $(addsuffix .exe,$(HISAT2_BIN_LIST) $(HISAT2_BIN_LIST_AUX)) ; \
	else \
		zip tmp.zip $(BIN_PKG_LIST) $(HISAT2_BIN_LIST) $(HISAT2_BIN_LIST_AUX) ; \
	fi
	mv tmp.zip .bin.tmp/hisat2-$(VERSION)
	cd .bin.tmp/hisat2-$(VERSION) ; unzip tmp.zip ; rm -f tmp.zip
	cd .bin.tmp ; zip -r hisat2-$(VERSION)-$(BITS).zip hisat2-$(VERSION)
	cp .bin.tmp/hisat2-$(VERSION)-$(BITS).zip .
	rm -rf .bin.tmp

.PHONY: doc
doc: doc/manual.inc.html MANUAL

doc/manual.inc.html: MANUAL.markdown
	pandoc -T "HISAT2 Manual" -o $@ \
	 --from markdown --to HTML --toc $^
	perl -i -ne \
	 '$$w=0 if m|^</body>|;print if $$w;$$w=1 if m|^<body>|;' $@

MANUAL: MANUAL.markdown
	perl doc/strip_markdown.pl < $^ > $@

.PHONY: clean
clean:
	rm -f $(CMA_BIN_LIST) $(CMA_BIN_LIST_AUX) \
	$(addsuffix .exe,$(CMA_BIN_LIST) $(CMA_BIN_LIST_AUX)) \
	cma-src.zip cma-bin.zip
	rm -f core.* .tmp.head
	rm -rf *.dSYM

.PHONY: push-doc
push-doc: doc/manual.inc.html
	scp doc/*.*html doc/indexes.txt salz-dmz:/ccb/salz7-data/www/ccb.jhu.edu/html/software/hisat2/
