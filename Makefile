.PHONY: clean

fqgrep: fqgrep.o bitap.o bm.o
	gcc -g -lz -o fqgrep fqgrep.o bitap.o bm.o

fqgrep.o: fqgrep.c kseq.h
	gcc -g -I. -c fqgrep.c

bitap.o: bitap.c bitap.h
	gcc -g -I. -c bitap.c

bm.o: bm.c bm.h
	gcc -g -I. -c bm.c

clean:
	rm fqgrep *.o
