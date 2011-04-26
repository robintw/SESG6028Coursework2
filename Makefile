all:
	mpicc transpose_v2.c array_alloc.c timer.c -o transpose_v2 -I. -g -Wall -pedantic -std=c89
	mpicc transpose_v1.c array_alloc.c timer.c -o transpose_v1 -I. -g -Wall -pedantic -std=c89

