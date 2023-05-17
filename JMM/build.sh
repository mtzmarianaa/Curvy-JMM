#!/usr/bin/env bash

CFLAGS="-g -Wall -Werror -std=c99 -fsanitize=undefined"

gcc $CFLAGS -c eik_grid.c -o eik_grid.o
gcc $CFLAGS -c files_methods.c -o files_methods.o
gcc $CFLAGS -c neighbors.c -o neighbors.o
gcc $CFLAGS -c mesh2D.c -o mesh2D.o
gcc $CFLAGS -c linAlg.c -o linAlg.o
gcc $CFLAGS -c priority_queue.c -o priority_queue.o
gcc $CFLAGS -c test_eik_grid.c -o test_eik_grid.o
gcc $CFLAGS -o test_eik_grid test_eik_grid.o mesh2D.o eik_grid.o files_methods.o neighbors.o linAlg.o priority_queue.o -ljson-c -lm
