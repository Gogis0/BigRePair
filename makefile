
# targets not producing a file declared phony
.PHONY: all ctph myrepair clean

all: procdic postproc ctph myrepair

procdic: procdic.c
	gcc -O9 -o procdic procdic.c

postproc: postproc.c
	gcc -O9 -o postproc postproc.c

ctph:
	make -C ctph

myrepair:
	make -C myrepair

clean:
	rm -f procdic postproc
	make -C ctph clean
	make -C myrepair

