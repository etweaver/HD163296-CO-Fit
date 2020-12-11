all : http_server

clean : 
	rm -f http_server *.o

BOOST_PREFIX:=/home/vagrant/install_prefix
CXXFLAGS+=-O3 -I$(BOOST_PREFIX)/include
LDFLAGS+=-L$(BOOST_PREFIX)/lib -lboost_thread -lboost_system -pthread -Wl,-rpath,$(BOOST_PREFIX)/lib

http_server : Makefile connection.o main.o mime_types.o reply.o request_handler.o request_parser.o server.o
	$(CXX) $(LDFLAGS) -pthread connection.o main.o mime_types.o reply.o request_handler.o request_parser.o server.o -o http_server

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@
