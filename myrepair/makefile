
all: repair despair 

repair: repair.o array.o hash.o heap.o records.o basics.o makefile
	gcc -O9 -o repair repair.o array.o hash.o heap.o records.o basics.o

repair.o: repair.c array.h hash.h heap.h records.h basics.h makefile
	gcc -O9 -c repair.c

despair: despair.o basics.o makefile
	gcc -O9 -o despair despair.o basics.o

despair.o: despair.c basics.h makefile
	gcc -O9 -c despair.c

array.o: array.c array.h hash.h heap.h records.h basics.h makefile
	gcc -O9 -c array.c

hash.o: hash.c array.h hash.h heap.h records.h basics.h makefile
	gcc -O9 -c hash.c

heap.o: heap.c array.h hash.h heap.h records.h basics.h makefile
	gcc -O9 -c heap.c

records.o: records.c array.h hash.h heap.h records.h basics.h makefile
	gcc -O9 -c records.c

basics.o: basics.c basics.h makefile
	gcc -O9 -c basics.c

clean:
	rm -f repair despair *.o
