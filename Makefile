#~ LHAPDF_INC_PATH := -I/cvmfs/sft.cern.ch/lcg/views/LCG_95/x86_64-centos7-gcc7-opt/include
#~ LHAPDF_LIB_PATH := -L/cvmfs/sft.cern.ch/lcg/views/LCG_95/x86_64-centos7-gcc7-opt/lib
#~ LHAPDF_INC_PATH := -I/$(LHAPATH)/include
#~ LHAPDF_LIB_PATH := -L/$(LHAPATH)/lib
LHAPDF_INC_PATH := -I/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-omkpbe2/include
LHAPDF_LIB_PATH := -L/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-omkpbe2/lib
LHAPDF_LIB := -lLHAPDF

EXTRA_CFLAGS := $(LHAPDF_INC_PATH)
EXTRA_LDFLAGS := $(LHAPDF_LIB_PATH) $(LHAPDF_LIB)

CFLAGS := -g -O3 $(EXTRA_CFLAGS)
CFLAGS += $(shell root-config --cflags --libs)
LDFLAGS := $(shell root-config --cflags --libs) -L. $(EXTRA_LDFLAGS)

CC := g++
LD := g++

DIRS := include

PROGRAM := guent.her
TARGETS := $(PROGRAM)

SOURCES := $(wildcard *.cc)
SOURCES += $(foreach dir, $(DIRS), $(wildcard $(dir)/*.cc))

OBJECTS := $(SOURCES:.cc=.o)

.PHONY: all clean
all: $(TARGETS)

$(PROGRAM): $(OBJECTS)
		$(LD) $^ $(CFLAGS) $(LDFLAGS) -o $@
		
%.o : %.cc
		$(CC) $(CFLAGS) -c -o $@ $<
		
clean:
		@rm -f $(PROGRAM) $(OBJECTS)
