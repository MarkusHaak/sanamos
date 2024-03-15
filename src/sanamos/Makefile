CC=gcc
# Using --std=gnu11 and -lrt for clock_gettime
# -lm to make <math.h> work
# -fopenmp for openmp
# -funsigned-char because default is not defined for char!!
CFLAGS= --std=gnu11 -O3 -march=native -Wall -pedantic -lrt -fopenmp -lm -funsigned-char 
LDFLAGS= -fPIC -shared
DEPS = sa.h
OBJ = main.o sa.o sa.so

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: sa.so main

sa.so: sa.c $(DEPS)
	gcc $(LDFLAGS) $(CFLAGS) sa.c -o sa.so

main: $(OBJ)
	gcc -o $@ $^ $(CFLAGS)

.PHONY: clean
clean:
	-rm -f main $(OBJ)
