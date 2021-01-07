#TODO
#Add target that builds the http server

#CC=clang
#CXX=clang++
AR=ar
LD=clang++
DYN_SUFFIX=.dylib
DYN_OPT=-dynamiclib
PREFIX=/usr/local/
CXXFLAGS+= -g -fPIC -O3 -mavx -std=c++17
LDFLAGS+= -lcfitsio -lfftw3 -lcapnp-rpc -lcapnp -lcurl -lkj-async -lkj -pthread


.PHONY: all clean
all : modelEval driver

clean : 
	rm -f driver modelEval
	rm -rf build/*.o
	rm -rf build/distribute.capnp.h
	rm -rf build/distribute.capnp.c++
	rmdir build

modelEval : build/modelEval.o build/geometry.o build/diskPhysics.o build/HTTPRequests.o build/ParameterSet.o build/distribute.capnp.o
	$(CXX) -o modelEval $^ $(LDFLAGS)

driver: build/driver.o build/ParameterSet.o  build/distribute.capnp.o
	$(CXX) -o driver $^ $(LDFLAGS)

test : build/modelTest.o build/geometry.o build/diskPhysics.o
	$(CXX) -o modelTest $^ $(LDFLAGS)

build/modelEval.o : modelEval.cpp geometry.h grid.h HTTPRequests.h image.h worker.h build/distribute.capnp.h | build
	$(CXX) $(CXXFLAGS) modelEval.cpp -c -o build/modelEval.o
build/modelTest.o : modelEval.cpp geometry.h grid.h image.h | build
	$(CXX) $(CXXFLAGS) modelTest.cpp -c -o build/modelTest.o
build/sampler.o : sampler.cpp geometry.h grid.h image.h mcmc.h | build
	$(CXX) $(CXXFLAGS) sampler.cpp -c -o build/sampler.o
build/geometry.o : geometry.cpp diskPhysics.h geometry.h | build
	$(CXX) $(CXXFLAGS) $(INCFLAGS) geometry.cpp -c -o build/geometry.o
build/image.o : image.cpp image.h diskPhysics.h | build
	$(CXX) $(CXXFLAGS) $(INCFLAGS) image.cpp -c -o build/image.o
build/diskPhysics.o : diskPhysics.cpp diskPhysics.h | build
	$(CXX) $(CXXFLAGS) $(INCFLAGS) diskPhysics.cpp -c -o build/diskPhysics.o
build/ParameterSet.o : ParameterSet.cpp ParameterSet.h | build
	$(CXX) $(CXXFLAGS) $(INCFLAGS) ParameterSet.cpp -c -o build/ParameterSet.o
build/distribute.capnp.h : distribute.capnp Makefile | build
	capnp compile -oc++ distribute.capnp
	mv distribute.capnp.h build/
build/HTTPRequests.o : HTTPRequests.cpp HTTPRequests.h | build
	$(CXX) $(CXXFLAGS) $(INCFLAGS) HTTPRequests.cpp -c -o build/HTTPRequests.o

build/distribute.capnp.c++ : build/distribute.capnp.h Makefile | build
	mv distribute.capnp.c++ build/

build/distribute.capnp.o : build/distribute.capnp.c++ Makefile | build
	$(CXX) $(CXXFLAGS) -c build/distribute.capnp.c++ -o build/distribute.capnp.o

build/driver.o : driver.cpp driver.h build/distribute.capnp.h ParameterSet.h radmcInterface.h | build
	$(CXX) $(CXXFLAGS) -c driver.cpp -o build/driver.o

build :
	mkdir build

install :
	cp disk.so $(PREFIX)/include/
