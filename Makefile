CFLAGS := $(shell root-config --cflags --libs)
CFLAGS += -g -O3
LDFLAGS := $(shell root-config --cflags --libs)

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
