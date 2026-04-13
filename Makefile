CC = gcc
CFLAGS = -Wall -Wextra -std=c99
LDFLAGS = -lm

.PHONY: all test clean

all: test_matrix

test_matrix: test_matrix.c matrix.h matrix.c
	$(CC) $(CFLAGS) -o $@ test_matrix.c $(LDFLAGS)

test: test_matrix
	./test_matrix

clean:
	rm -f test_matrix *.o
