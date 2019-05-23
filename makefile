
all: procdic postproc ctph myrepair

procdic: procdic.c
	gcc -O9 -o procdic procdic.c

postproc: postproc.c
	gcc -O9 -o postproc postproc.c

ctph:
	make -C ctph

myrepair:
	make -C myrepair


