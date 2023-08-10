#!/usr/bin/env bash

CFLAGS="-g -Wall -Werror -std=c99 -fsanitize=undefined"

gcc $CFLAGS -c eik_grid.c -o eik_grid.o
gcc $CFLAGS -c files_methods.c -o files_methods.o
gcc $CFLAGS -c neighbors.c -o neighbors.o
gcc $CFLAGS -c mesh2D.c -o mesh2D.o
gcc $CFLAGS -c linAlg.c -o linAlg.o
gcc $CFLAGS -c marcher_T2.c -o marcher_T2.o
gcc $CFLAGS -c priority_queue.c -o priority_queue.o
gcc $CFLAGS -c test_eik_gridEasy.c -o test_eik_gridEasy.o
gcc $CFLAGS -c opti_method.c -o opti_method.o
gcc $CFLAGS -o test_eik_gridEasy test_eik_gridEasy.o mesh2D.o eik_grid.o marcher_T2.o  files_methods.o neighbors.o linAlg.o priority_queue.o opti_method.o -ljson-c -lm
