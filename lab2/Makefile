CC = gcc

CFLAGS = -Wall -mmmx -msse -msse2 -mssse3 -msse4 -O2 

# Utility routines
#UTIL_OBJS += main.o 

#

all: main

#OBJ = $(UTIL_OBJS)
$(OBJ) : %.o : %.c
	$(CC) -c $(CFLAGS) -o $@ $<

main.o: main.c
	$(CC)  -c main.c -o main.o $(CFLAGS)  

main: $(OBJ) main.o
	$(CC)  $(OBJ) main.o -o main

clean:
	rm -f *.o 

