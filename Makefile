RELVERSION  = $(shell cat .release)

ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

MakefileFullPath = $(abspath $(lastword $(MAKEFILE_LIST)))
MakefileDirFullPath = $(shell dirname $(MakefileFullPath))
INSTALLDIR = $(MakefileDirFullPath)/install.$(RELVERSION)/

CXX  = g++
CXX += -I./

CXXFLAGS  = -g -Wall -fPIC -Wno-deprecated
CXXFLAGS += $(ROOTCFLAGS)
CXXFLAGS += $(ROOTLIBS)
CXXFLAGS += $(ROOTGLIBS)
CXXFLAGS += -std=c++14
CXXFLAGS += -fconcepts

OUTLIB = ./obj/

#----------------------------------------------------#

all: makedir convert2root

makedir:
	mkdir -p $(OUTLIB);

.PHONY: printmakehelp_and_reminder
printmakehelp_and_reminder: convert2root.cpp Makefile
	$(info  /*****************************************************************/)
	$(info  * task --> printmakehelp_and_reminder: convert2root.cpp Makefile *)
	$(info  * $$@ ----> $@                            *)
	$(info  * $$< --------------------------------> $<          *)
	$(info  * $$^ --------------------------------> $^ *)
	$(info  /*****************************************************************/)

convert2root: convert2root.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS)

install: makedir obj/libconvert2root.so
	mkdir -p $(INSTALLDIR);
	cp $(OUTLIB)libconvert2root.so $(INSTALLDIR)libconvert2root.so
	cp src/*.hh $(INSTALLDIR).

cleaninstall:
	rm -rf $(INSTALLDIR)

clean:
	rm -f convert2root
	rm -f *~
	rm -f .*~
	rm -f $(OUTLIB)*.o
	rm -f $(OUTLIB)*.so
