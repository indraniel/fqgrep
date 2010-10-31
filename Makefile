.PHONY: clean

fqgrep: fqgrep.o bm.o
	gcc -g -lz -ltre -o fqgrep fqgrep.o bm.o

fqgrep.o: fqgrep.c kseq.h
	gcc -g -I. -c fqgrep.c

bm.o: bm.c bm.h
	gcc -g -I. -c bm.c

clean:
	rm fqgrep *.o
