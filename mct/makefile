CFLAGS=-Wall -pedantic -funroll-loops -std=c99 -O3

.PHONY: all clean

all: nep nep_2006

# korrigiertes Schema (Standard)
nep: nep.c
	gcc $(CFLAGS) -DQUAD=1 -o nep nep.c -lm

# Originalschema der Diplomarbeit 2006
nep_2006: nep.c
	gcc $(CFLAGS) -DQUAD=0 -o nep_2006 nep.c -lm

clean:
	rm -f nep nep_2006
