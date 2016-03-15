CXX=g++
CXXINCLUDE=$(MKLROOT)/include
MKLLIBDIR=$(MKLROOT)/lib/intel64
CXXOPTIONS=-fopenmp -m64 -I$(MKLROOT)/include
LIBS=-Wl,--start-group $(MKLLIBDIR)/libmkl_intel_lp64.a $(MKLLIBDIR)/libmkl_core.a $(MKLLIBDIR)/libmkl_gnu_thread.a $(HOME)/ARPACK/libarpack_linux.a -Wl,--end-group -ldl -lpthread -lgfortran -lm

all: eigs_sps eigs_ges dsort

eigs_sps: eigs_sps.o eigs_sps.cpp
	$(CXX) -g -c -o eigs_sps.o eigs_sps.cpp $(CXXOPTIONS) $(LIBS)
	ar -rv libarpack.a eigs_sps.o
	ranlib libarpack.a

eigs_ges: eigs_ges.o eigs_ges.cpp
	$(CXX) -g -c -o eigs_ges.o eigs_ges.cpp $(CXXOPTIONS) $(LIBS)
	ar -rv libarpack.a eigs_ges.o
	ranlib libarpack.a

dsort: dsort.o dsort.cpp
	$(CXX) -g -c -o dsort.o dsort.cpp $(CXXOPTIONS) $(LIBS)
	ar -rv libarpack.a dsort.o
	ranlib libarpack.a

clean:
	rm eigs_sps.o eigs_ges.o dsort.o
