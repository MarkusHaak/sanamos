CC=gcc
# Using --std=gnu11 and -lrt for clock_gettime
# -lm to make <math.h> work
# -fopenmp for openmp
# -funsigned-char because default is not defined for char!!
CFLAGS= --std=gnu11 -Wall -pedantic -lrt -fopenmp -lm -funsigned-char 
LDFLAGS= -fPIC -shared
DEPS = suffix_array.h
OBJ = main.o suffix_array.o sa.so

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: suffix_array.so main

suffix_array.so: suffix_array.c $(DEPS)
	gcc $(LDFLAGS) $(CFLAGS) suffix_array.c -o sa.so

main: $(OBJ)
	gcc -o $@ $^ $(CFLAGS)

.PHONY: clean
clean:
	-rm -f main $(OBJ)
