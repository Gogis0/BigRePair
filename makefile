
# targets not producing a file declared phony
.PHONY: all ctph myrepair large_repair clean

all: procdic postproc ctph myrepair large_repair

procdic: procdic.c
	gcc -O9 -o procdic procdic.c

postproc: postproc.c
	gcc -O9 -o postproc postproc.c

ctph:
	make -C ctph

myrepair:
	make -C myrepair

large_repair:
	make -C large_repair

clean:
	rm -f procdic postproc
	make -C ctph clean
	make -C myrepair clean
	make -C large_repair clean

