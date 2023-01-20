# Compiler
CC=gcc

# Options
SOURCES=src/lib/twister.c src/lib/lbfgs.c src/plm.c src/inference.c src/weights.c src/main.c
GCCFLAGS=-std=c99 -lm -O3 -msse4.2
CLANGFLAGS=-lm -Wall -Ofast -msse4.2

all:
	$(CC) $(SOURCES) -o bin/plmc $(GCCFLAGS)

all32:
	$(CC) $(SOURCES) -o bin/plmc $(GCCFLAGS) -D USE_FLOAT

all-dev:
	$(CC) $(SOURCES) -o bin/plmc $(GCCFLAGS) -Wall

all-openmp:
	$(CC) $(SOURCES) -o bin/plmc -fopenmp $(GCCFLAGS)

all-openmp32:
	$(CC) $(SOURCES) -o bin/plmc -fopenmp $(GCCFLAGS) -D USE_FLOAT

all-mac:
	clang $(SOURCES) -o bin/plmc $(CLANGFLAGS)

all-mac32:
	clang $(SOURCES) -o bin/plmc $(CLANGFLAGS) -D USE_FLOAT

# If using homebrew for openMP (libomp)
all-mac-openmp: LIBOMP_PREFIX=$(shell brew --prefix libomp)
all-mac-openmp:
	clang $(SOURCES) -o bin/plmc -Xpreprocessor -fopenmp $(CLANGFLAGS) -lomp -L$(LIBOMP_PREFIX)/lib/ -I$(LIBOMP_PREFIX)/include/

clean:
	rm -rf bin/*
