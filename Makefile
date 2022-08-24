CC=g++
DIC=$(PWD)
#CFLAGS=-std=c++11
CFLAGS=-std=c++11 -fopenmp
#CFLAGS=-c -Wall -g  -I $(DIC)
#LDFLAGS= -I $(DIC)/armadillo/include/ -DARMA_DONT_USE_WRAPPER -llapack -lblas -lgslcblas  -lgsl
LDFLAGS= -I $(DIC)/armadillo/include/ -DARMA_DONT_USE_WRAPPER -lgslcblas  -lgsl
SOURCES1=MsCaviar.cpp MsPostCal.cpp MsUtil.cpp
EXECUTABLE1=MsCAVIAR


all: $(SOURCES1) $(EXECUTABLE1)

$(EXECUTABLE1): $(SOURCES1)
	$(CC) $(CFLAGS) $(SOURCES1)   $(LDFLAGS) -o $@ \
	/home/tmirmira/miniconda3/envs/env_popcorn2/lib/libblas.so \
	/home/tmirmira/miniconda3/envs/env_popcorn2/lib/liblapack.so
clean:
	rm MsCAVIAR
