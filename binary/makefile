CFLAGS=-Wall -pedantic -funroll-loops -std=c99
LINKFLAGS=-lgsl -lgslcblas -lm -O3
all:
	rm -rf $(CLEANFILES)
	gcc $(CFLAGS) $(LINKFLAGS) -o static static.c 
