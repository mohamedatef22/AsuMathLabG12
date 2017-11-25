# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
CFLAGS=-c -Wall

all: execute

execute: main.o  Matrix.o
	$(CC) main.o  Matrix.o -o execute

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp



hello.o: hello.cpp
	$(CC) $(CFLAGS) Matrix.cpp

clean:
	rm -rf *o execute
