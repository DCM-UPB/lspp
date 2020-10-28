-include $(TOPDIR)/make.in
-include ../make.in
.PHONY: all clean 
all: libtoolbox.a histogram ndhistogram
libdimred.o: libdimred.cpp dimreduce.hpp tbdefs.hpp
	$(CXX) -c -o libdimred.o $(CCINC) libdimred.cpp $(CCOPTS) -DNOARPACK
libioatoms.o: libioatoms.cpp ioatoms.hpp tbdefs.hpp
	$(CXX) -c -o libioatoms.o $(CCINC) libioatoms.cpp $(CCOPTS) 
libtb.o: libtb.cpp tbdefs.hpp
	$(CXX) -c -o libtb.o $(CCINC) libtb.cpp $(CCOPTS) 
libminsearch.o: libminsearch.cpp tbdefs.hpp minsearch.hpp
	$(CXX) -c -o libminsearch.o $(CCINC) libminsearch.cpp $(CCOPTS) 
libioparser.o: libioparser.cpp ioparser.hpp tbdefs.hpp
	$(CXX) -c -o libioparser.o $(CCINC) libioparser.cpp $(CCOPTS)
libclparser.o: libclparser.cpp clparser.hpp tbdefs.hpp
	$(CXX) -c -o libclparser.o $(CCINC) libclparser.cpp $(CCOPTS)
librndgen.o: librndgen.cpp rndgen.hpp
	$(CXX) -c -o librndgen.o $(CCINC) librndgen.cpp $(CCOPTS) 
libfmblas.o: matrix-full-blas.hpp libfmblas.cpp
	$(CXX) -c -o libfmblas.o $(CCINC) libfmblas.cpp $(CCOPTS) 
liblinalg.o: liblinalg.cpp linalg.hpp
	$(CXX) -c -o liblinalg.o $(CCINC) liblinalg.cpp $(CCOPTS) 
libunits.o: libunits.cpp conv-units.hpp
	$(CXX) -c -o libunits.o $(CCINC) libunits.cpp $(CCOPTS) 
libtoolbox.a: libtb.o libminsearch.o libioparser.o libclparser.o librndgen.o libfmblas.o liblinalg.o libioatoms.o libunits.o libdimred.o
	ar rcs libtoolbox.a libtb.o libminsearch.o libioparser.o libclparser.o librndgen.o libfmblas.o liblinalg.o libioatoms.o libunits.o libdimred.o
autocorr: autocorr.cpp tools-autocorr.hpp tbdefs.hpp libtoolbox.a
	$(CXX) -o autocorr $(CCINC) $(CCOPTS) $(LDOPTS) autocorr.cpp libtoolbox.a -lfftw3 -lm 
	cp autocorr ../bin
crosscorr: crosscorr.cpp tools-autocorr.hpp tbdefs.hpp libtoolbox.a
	$(CXX) -o crosscorr $(CCINC) $(CCOPTS) $(LDOPTS) crosscorr.cpp libtoolbox.a -lfftw3 -lm 
	cp crosscorr ../bin
histogram: histogram.cpp tools-histogram.hpp tbdefs.hpp libtoolbox.a
	$(CXX) -o histogram $(CCINC) $(CCOPTS) $(LDOPTS) histogram.cpp libtoolbox.a -lm
ndhistogram: ndhistogram.cpp tools-histogram.hpp tbdefs.hpp libtoolbox.a
	$(CXX) -o ndhistogram $(CCINC) $(CCOPTS) $(LDOPTS) ndhistogram.cpp libtoolbox.a -lm
trajworks.o: trajworks.cpp ioatoms.hpp libtoolbox.a
	$(CXX) -c $(CCINC) $(CCOPTS) -I/usr/include trajworks.cpp -D__USEREGEX
trajworks: trajworks.o libtoolbox.a
	$(CXX) -o trajworks $(CCINC) $(CCOPTS) $(LDOPTS) -I/usr/include trajworks.o libtoolbox.a -lm -lfftw3 -llapack -lboost_regex
dimred: dimred.cpp libtoolbox.a
	$(CXX) -o dimred  $(CCINC) $(CCOPTS) $(LDOPTS) dimred.cpp libtoolbox.a -llapack 
clean:
	rm -f *~ *.o
