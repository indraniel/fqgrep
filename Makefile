.PHONY: clean

fqgrep: fqgrep.o libbitap.o bm.o
	gcc -g -lz -o fqgrep fqgrep.o libbitap.o bm.o

fqgrep.o: fqgrep.c kseq.h
	gcc -g -I. -c fqgrep.c

libbitap.o: libbitap.c libbitap.h
	gcc -g -I. -c libbitap.c

bm.o: bm.c bm.h
	gcc -g -I. -c bm.c

clean:
	rm fqgrep *.o
