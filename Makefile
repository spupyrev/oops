##  Makefile for Standard, Debug, Release, and Release-static versions
##
##    "make"    for the standard version (optimized, but with debug information and assertions active)
##    "make r"  for a release version
##    "make rs" for a statically linked release version
##    "make d"  for a debug version (no optimizations)

EXEC      = oops
DEPDIR    = breakid glucose satsuma
MROOT = $(PWD)/src

PWD       = $(shell pwd)
## run in parallel
# MAKEFLAGS += -j

CSRCS      = $(wildcard $(MROOT)/*.cpp)
CHDRS      = $(wildcard $(MROOT)/*.h)
DSRCS      = $(foreach dir, $(DEPDIR), $(filter-out $(MROOT)/$(dir)/main*.cpp, $(wildcard $(MROOT)/$(dir)/*.cpp)))


COBJS      = $(CSRCS:.cpp=.o) $(DSRCS:.cpp=.o)
LCOBJS     = $(addsuffix l,  $(COBJS))
PCOBJS     = $(addsuffix p,  $(COBJS))
DCOBJS     = $(addsuffix d,  $(COBJS))
RCOBJS     = $(addsuffix r,  $(COBJS))

# Use clang, if installed
# CXX       = ${LLVM_HOME}/clang++
# CFLAGS    ?= -Wall -Wno-deprecated-declarations -Wsign-compare -std=c++17 -fno-omit-frame-pointer
# LFLAGS    ?= -Wall -fuse-ld=lld --ld-path=${LLVM_HOME}/ld.lld
# Otherwise, fallback to gcc
CXX       = g++
CFLAGS    ?= -Wall -Wno-deprecated-declarations -Wsign-compare -std=c++17
LFLAGS    ?= -Wall

## standard flags
CFLAGS    += -I$(MROOT)
LFLAGS    += -lz

.PHONY : s lto p d r rs clean

s:	 $(EXEC)
lto: $(EXEC)_lto
p:	 $(EXEC)_profile
d:	 $(EXEC)_debug
r:	 $(EXEC)_release
rs:	 $(EXEC)_static

## Compile options
%.o:			CFLAGS += -O3 -g -D DEBUG
%.ol:			CFLAGS += -O3 -g -D NDEBUG -flto
%.op:			CFLAGS += -O3 -pg -g -D NDEBUG
%.od:			CFLAGS += -O0 -g -D DEBUG
%.or:			CFLAGS += -O3 -g -D NDEBUG

## Link options
$(EXEC):		    LFLAGS += -g
$(EXEC)_lto:		LFLAGS += -g -flto
$(EXEC)_profile:	LFLAGS += -g -pg
$(EXEC)_debug:		LFLAGS += -g
$(EXEC)_release:	LFLAGS +=
$(EXEC)_static:		LFLAGS += --static

## Dependencies
$(EXEC):		    $(COBJS)
$(EXEC)_lto:	    $(LCOBJS)
$(EXEC)_profile:	$(PCOBJS)
$(EXEC)_debug:		$(DCOBJS)
$(EXEC)_release:	$(RCOBJS)
$(EXEC)_static:		$(RCOBJS)


## Build rule
%.o %.ol %.op %.od %.or:	%.cpp
	@echo Compiling: "$(subst $(MROOT)/,,$@)"
	@$(CXX) $(CFLAGS) -c -o $@ $<

## Link rule
$(EXEC) $(EXEC)_lto $(EXEC)_profile $(EXEC)_debug $(EXEC)_release $(EXEC)_static:
	@echo Linking: "$@ [$(LFLAGS)]"
	@$(CXX) $^ $(LFLAGS) -o ${EXEC}

## Clean rule
clean:
	@rm -f $(EXEC) $(EXEC)_*  $(COBJS) $(LCOBJS) $(PCOBJS) $(DCOBJS) $(RCOBJS) \
	  log_* depend.mk

## Make dependencies
depend.mk: $(CSRCS) $(CHDRS)
	@echo Making dependencies
	@$(CXX) $(CFLAGS) -I$(MROOT) \
	   $(CSRCS) -MM | sed 's|\(.*\):|$(PWD)/\1 $(PWD)/\1l $(PWD)/\1r $(PWD)/\1d $(PWD)/\1p:|' > depend.mk
	@for dir in $(DEPDIR); do \
	      if [ -r $(MROOT)/$${dir}/depend.mk ]; then \
		  echo Depends on: $${dir}; \
		  cat $(MROOT)/$${dir}/depend.mk >> depend.mk; \
	      fi; \
	  done
