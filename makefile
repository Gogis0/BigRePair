# compilation flags
CFLAGS=-g -Wall -std=c99
CC=gcc 

# executables in this directory
EXECS=procdic postproc iprocdic ipostproc

# targets not producing a file declared phony
.PHONY: all ctph myrepair large_repair clean

all: $(EXECS) ctph myrepair large_repair

# general rule for the targets in this directory 
%: %.c
	gcc $(CFLAGS) -o $@ $<

ctph:
	make -C ctph

myrepair:
	make -C myrepair

large_repair:
	make -C large_repair

clean:
	rm -f $(EXECS)
	make -C ctph clean
	make -C myrepair clean
	make -C large_repair clean

