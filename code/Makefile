CFLAGS = -Ofast -Wall -g
CC = gcc
LIBS = -lm
LDFLAGS = $(LIBS)

all: cart cltcart

main.o: main.c
	$(CC) -c $(CFLAGS) -DOFFSET=0.005 -o $@ $^

cart: main.o cart.o
	$(CC) $(CFLAGS) $(LIBS) -lfftw3 -o $@ $^

cltcart: main.o cltcart.o

icltcart: main.o icltcart.o

clean:
	rm -fv *.o cart cltcart
