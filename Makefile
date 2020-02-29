CC=gcc
CFLAGS=-Wall
LINKERFLAG=-lm

all: lang2.c
	$(CC) $(CFLAGS) lang2.c -O2 -o lang $(LINKERFLAG)

dbg: lang2.c
	$(CC) -g $(CFLAGS) lang2.c -o lang $(LINKERFLAG)

clean:
	rm -rf lang
