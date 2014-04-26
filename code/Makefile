CFLAGS = -Ofast -Wall
CC = gcc
LIBS = -lm
LDFLAGS = $(LIBS)

all: cart cltcart

cart-main.o: main.c
	$(CC) -c $(CFLAGS) -DOFFSET=0.005 -o $@ $^

cltcart-main.o: main.c
	$(CC) -c $(CFLAGS) -DOFFSET=0.005 -o $@ $^

cart: cart-main.o cart.o
	$(CC) $(CFLAGS) $(LIBS) -lfftw3 -o $@ $^

cartv: cart-main.o cartv.o
	$(CC) $(CFLAGS) $(LIBS) -lfftw3 -o $@ $^

cltcart: cltcart-main.o cltcart.o

clean:
	rm -fv *.o