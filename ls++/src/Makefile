-include $(TOPDIR)/make.in
-include ../make.in
.PHONY: all clean
all: lsdm lsdm-cheb libls++.a
	mv lsdm lsdm-cheb ../bin
clean:
	rm -f *~  *.o lsdm
libls++.a: libls.o libls-cheb.o
	ar rcs $@ $^
libls.o: libls.cpp libls.hpp libls-io.hpp
	$(CXX) -c $(CCINC) libls.cpp $(CCOPTS)
libls-cheb.o: libls-cheb.cpp libls-cheb.hpp libls.hpp libls-io.hpp
	$(CXX) -c $(CCINC) libls-cheb.cpp $(CCOPTS)
lsdm: lsdm.cpp libls.hpp libls.o
	$(CXX) -o lsdm $(CCINC) libls.o lsdm.cpp  $(CCOPTS) $(LDOPTS)
# lsdm-cheb: libls.o libls-cheb.o lsdm-cheb.cpp
# 	$(CXX) -o lsdm-cheb $(CCINC) lsdm-cheb.cpp libls.o libls-cheb.o $(CCOPTS) $(LDOPTS)
lsdm-cheb: lsdm-cheb.o libls++.a
	$(CXX) -o $@ $(CCINC) $^ $(CCOPTS) $(LDOPTS)
lsdm-cheb.o: lsdm-cheb.cpp
	$(CXX) -c $(CCINC) $^ $(CCOPTS)
#$(CXX) -o lsdm-cheb $(CCINC) libls.o libls-cheb.o lsdm-cheb.cpp $(CCOPTS)
#$(LIBTB) /usr/lib/sse2/libcblas.a /usr/lib/sse2/libatlas.a
#lsdm-cheb-dim: lsdm-cheb-dim.o libls++.a
#	$(CXX) -o $@ $(CCINC) $^ $(CCOPTS) $(LDOPTS)
#lsdm-cheb-dim.o: lsdm-cheb-dim.cpp
#	$(CXX) -c $(CCINC) $^ $(CCOPTS)
