# compilation flags
CFLAGS=-g -Wall -std=c99
CC=gcc 

# executables in this directory
EXECS=procdic postproc iprocdic ipostproc

# targets not producing a file declared phony
.PHONY: all ctph repair large_repair clean

all: $(EXECS) ctph repair large_repair

# general rule for the targets in this directory 
%: %.c
	gcc $(CFLAGS) -o $@ $< -D Unique=0x78000000

ctph:
	make -C ctph

repair:
	make -C repair

large_repair:
	make -C large_repair

clean:
	rm -f $(EXECS)
	make -C ctph clean
	make -C repair clean
	make -C large_repair clean

