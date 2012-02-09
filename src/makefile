#makefile
CFLAGS = -g  -std=c99

all:
	mpicc -g  -std=c99 driver.c par_mult.c seq_mult.c -DNO_TIMING  -lm -Lhmetis -lhmetis -Wall -o run
	gcc -std=c99 generate.c -Wall -o generate
clean:
	rm -rf *.o 


