# compilation flags
CFLAGS=-g -Wall -std=c99

# executables in this directory
EXECS=procdic procdic_new postproc 

# targets not producing a file declared phony
.PHONY: all ctph myrepair large_repair clean


all: $(EXECS) ctph myrepair large_repair

procdic: procdic.c
	gcc $(CFLAGS) -o $@ $<

procdic_new: procdic_new.c
	gcc $(CFLAGS) -o $@ $<

postproc: postproc.c
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

