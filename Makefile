ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs) -lNew -lMinuit libutils.a

CXXFLAGS += -I$(ROOTSYS)/include $(shell root-config --cflags)
EXTRALIBS = $(shell root-config --libs) $(shell root-config --glibs) -lGeomPainter -lGeom

CC       = gcc
CXX      = g++ 
CXXFLAGS = -W -Wall -ansi -pedantic
CXXFLAGS += $(ROOTCFLAGS)
OPTIM    = -O2 -fexpensive-optimizations -funroll-loops
LIBS     = $(ROOTLIBS) 
GLIBS    = $(ROOTGLIBS)

nuSpectrum: main.o classNuSpectrum.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) main.o classNuSpectrum.o -o nuSpectrum $(LIBS) $(GLIBS)

main.o: main.cc classNuSpectrum.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c main.cc -o main.o $(LIBS) $(GLIBS)

classNuSpectrum.o: classNuSpectrum.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c classNuSpectrum.cc -o classNuSpectrum.o $(LIBS) $(GLIBS)

.PHONY: clean mrproper

clean:
	rm -f nuSpectrum
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)
