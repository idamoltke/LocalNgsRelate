#modied from htslib makefile
FLAGS=-O3

CXXFLAGS += $(FLAGS)

CXXSRC = $(wildcard *.cpp)
OBJ = $(CXXSRC:.cpp=.o)

all: LocalNgsRelate

PACKAGE_VERSION  = 0.123

# Adjust $(HTSSRC) to point to your top-level htslib directory
ifdef HTSSRC
$(info HTSSRC defined)
HTS_INCDIR=$(realpath $(HTSSRC))
HTS_LIBDIR=$(realpath $(HTSSRC))/libhts.a
else
$(info HTSSRC not defined, assuming systemwide installation -lhts)
endif


ifneq "$(wildcard ../../.git)" ""
PACKAGE_VERSION := $(shell git describe --always --dirty)
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))
endif

version.h:
	echo '#define LocalNgsRelate_version "$(PACKAGE_VERSION)"' > $@

.PHONY: misc clean test

-include $(OBJ:.o=.d)

ifdef HTSSRC
%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS) $*.cpp -I$(HTS_INCDIR) -D__WITH_BCF__ ## -DDB_EMIS -DDB_MP
	$(CXX) -MM $(CXXFLAGS) $*.cpp -I$(HTS_INCDIR) -D__WITH_BCF__ >$*.d

LocalNgsRelate: version.h $(OBJ)
	$(CXX) $(FLAGS)  -o LocalNgsRelate *.o $(HTS_LIBDIR) -lz -lm -lbz2 -llzma -lpthread -lcurl -D__WITH_BCF__
else
%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  -D__WITH_BCF__ $*.cpp
	$(CXX) -MM $(CXXFLAGS)  -D__WITH_BCF__ $*.cpp >$*.d

LocalNgsRelate: version.h $(OBJ)
	$(CXX) $(FLAGS)  -o LocalNgsRelate *.o -lz -lpthread -lhts  -lbz2 -llzma 
endif

clean:
	rm  -f *.o *.d LocalNgsRelate *~ version.h

.PHONY: all clean install install-all install-misc misc test

test:
	echo "Only subset of analyses is being tested"
	cd test;./testAll.sh ../ngsRelate
force:
