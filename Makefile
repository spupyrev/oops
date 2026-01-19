##  Makefile for Standard, Debug, Release (with LTO), and Profile versions
##
##    "make"      standard (optimized, debug info, asserts on)
##    "make r"    release (asserts off, LTO on)
##    "make p"    profiled build
##    "make d"    debug (no optimizations)

EXEC      = oops
DEPDIR    = breakid glucose satsuma

PWD       = $(CURDIR)
MROOT     = $(PWD)/src

# MAKEFLAGS += -j

CSRCS      = $(wildcard $(MROOT)/*.cpp)
CHDRS      = $(wildcard $(MROOT)/*.h)
DSRCS      = $(foreach dir,$(DEPDIR), \
              $(filter-out $(MROOT)/$(dir)/main*.cpp, \
              $(wildcard $(MROOT)/$(dir)/*.cpp)))

COBJS      = $(CSRCS:.cpp=.o) $(DSRCS:.cpp=.o)
PCOBJS     = $(addsuffix p,$(COBJS))
DCOBJS     = $(addsuffix d,$(COBJS))
RCOBJS     = $(addsuffix r,$(COBJS))

# Compiler / linker
CXX    ?= g++
CFLAGS ?= -Wall -Wextra -Wno-pedantic -Wno-shadow -Wformat=2 -Wundef -Wno-unused-parameter \
          -Wno-deprecated-declarations -Wsign-compare -std=c++17 \
LFLAGS ?= -Wall
LDLIBS += -lz

CFLAGS += -I$(MROOT)

.PHONY : s p d r clean
.DEFAULT_GOAL := s

s:   $(EXEC)
r:   $(EXEC)_release
p:   $(EXEC)_profile
d:   $(EXEC)_debug

## Compile options
%.o:  CFLAGS += -O3 -g
%.or: CFLAGS += -O3 -g -DNDEBUG -flto
%.op: CFLAGS += -O2 -pg -g -DNDEBUG -fno-omit-frame-pointer
%.od: CFLAGS += -O0 -g -DDEBUG -fno-omit-frame-pointer

## Link options
$(EXEC):            LFLAGS += -g
$(EXEC)_release:    LFLAGS += -flto
$(EXEC)_profile:    LFLAGS += -g -pg
$(EXEC)_debug:      LFLAGS += -g

## Dependencies
$(EXEC):            $(COBJS)
$(EXEC)_release:    $(RCOBJS)
$(EXEC)_profile:    $(PCOBJS)
$(EXEC)_debug:      $(DCOBJS)

## Build rule (auto dependency generation)
%.o %.op %.od %.or: %.cpp
	@echo Compiling: "$(subst $(MROOT)/,,$@)"
	@$(CXX) $(CFLAGS) -MMD -MP -c -o $@ $<

## Link rule
$(EXEC) $(EXEC)_profile $(EXEC)_debug $(EXEC)_release:
	@echo Linking: "$@ [$(LFLAGS) $(LDLIBS)]"
	@$(CXX) $^ $(LFLAGS) $(LDLIBS) -o $@

## Clean
clean:
	@rm -f $(EXEC) $(EXEC)_* \
	  $(COBJS) $(PCOBJS) $(DCOBJS) $(RCOBJS) \
	  $(COBJS:.o=.d) $(PCOBJS:.op=.d) \
	  $(DCOBJS:.od=.d) $(RCOBJS:.or=.d) \
	  log_*

## Auto-include generated dependencies
-include $(COBJS:.o=.d) \
         $(PCOBJS:.op=.d) \
         $(DCOBJS:.od=.d) \
         $(RCOBJS:.or=.d)
