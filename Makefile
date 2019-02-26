CFLAGS := $(shell root-config --cflags --libs)
CFLAGS += -g -O3
LDFLAGS := $(shell root-config --cflags --libs)

CC := g++
LD := g++

PROGRAM := a.out
TARGETS := $(PROGRAM)
SOURCES := $(wildcard *.cc)
OBJECTS := $(SOURCES:.cc=.o)

.PHONY: all clean
all: $(TARGETS)

$(PROGRAM): $(OBJECTS)
		$(LD) $^ $(CFLAGS) $(LDFLAGS) -o $@
		
%.o : %.cc
		$(CC) $(CFLAGS) -c -o $@ $<
		
clean:
		@rm -f $(PROGRAM) $(OBJECTS)
