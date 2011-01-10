.PHONY: clean macports genome clean-genome

fqgrep: fqgrep.o bm.o
	gcc -g -lz -ltre -o fqgrep fqgrep.o bm.o

macports: fqgrep.o bm.o
	gcc -g -L. -L/opt/local/lib -o fqgrep fqgrep.o bm.o -lz -ltre

genome: libfqgrep.a
	gcc -static -g -L. -o fqgrep fqgrep.o -lfqgrep -lz -ltre

libfqgrep.a: fqgrep.o bm.o
	ar rc libfqgrep.a fqgrep.o bm.o
	ranlib libfqgrep.a

fqgrep.o: fqgrep.c kseq.h
	gcc -g -I. -I /opt/local/include -c fqgrep.c

bm.o: bm.c bm.h
	gcc -g -I. -c bm.c

clean:
	rm fqgrep *.o

clean-genome:
	rm fqgrep *.o *.a
