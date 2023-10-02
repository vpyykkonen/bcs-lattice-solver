CC = g++
#CC = clang++

CFLAGS = -std=c++17 -Wall -g -O3 -fopenmp -march=native #-Wextra -MD -MP

INCLUDES = -I./libs/eigen-3.3.9/ -I/usr/include/hdf5/serial


LIBS =  -lhdf5_cpp -lhdf5

LFLAGS = -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial -L/usr/lib/x86_64-linux-gnu/hdf5/serial

SRCS = $(wildcard *.cpp)

SRCS = BCSModel.cpp Lattice.cpp fd_dist.cpp ScfSolver.cpp ScfMethod.cpp config_parser.cpp file_io.cpp

SRC1 = main.cpp

OBJS = $(SRCS:.cpp=.o)
OBJ1 = $(SRC1:.cpp=.o)


ALL_OBJS = $(OBJS) $(OBJ1) 
depends = $(all_obj:%.o=%.d)


MAIN1 = test

.PHONY: depend clean


all:	$(MAIN1) 
		@echo Everything went better than expected

serial: $(MAIN1)


$(MAIN1): $(OBJ1) $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -o $(MAIN1) $(OBJ1) $(OBJS) $(LIBS)


%.o: %.cpp Makefile
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS)   -c $< -o $@ $(LIBS)

clean:
	$(RM) *.o *~ *.d $(MAIN1) $(MAIN2) $(MAIN3) $(MAIN4)

-include $(DEPENDS)

#depend: $(SRCS)
# 		makedepend $(INCLUDES) $^
