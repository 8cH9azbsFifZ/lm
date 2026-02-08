CFLAGS=-Wall -pedantic -funroll-loops -std=c99
LINKFLAGS=-lgsl -lgslcblas -lm -O3

.PHONY: all clean

all:
	gcc $(CFLAGS) $(LINKFLAGS) -o static static.c

clean:
	rm -f static
