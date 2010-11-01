.PHONY: clean

fqgrep: libfqgrep.a
	gcc -static -g -L. -o fqgrep fqgrep.o -lfqgrep -lz -ltre

libfqgrep.a: fqgrep.o bm.o
	ar rc libfqgrep.a fqgrep.o bm.o
	ranlib libfqgrep.a

fqgrep.o: fqgrep.c kseq.h
	gcc -g -I. -c fqgrep.c

bm.o: bm.c bm.h
	gcc -g -I. -c bm.c

clean:
	rm fqgrep *.o *.a
