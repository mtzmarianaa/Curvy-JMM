#pragma once

static void swap_double(double *a, double *b);

static void swap_int(int *a, int *b);

static void heapify(double eik_queue[], int index_queue[], int size, int i);

static void insert(double eik_queue[], int index_queue[], double newNum, int newIndex);

static void insert_end(double eik_queue[], int index_queue[], double newNum, int newIndex);

static void delete_findValue(double eik_queue[], int index_queue[], double num);

static void delete_findIndex(double eik_queue[], int index_queue[], int ind);

static void deleteRoot(double eik_queue[], int index_queue[]);

static void printeik_queue(double eik_queue[], int index_queue[], int size);

static void update(double eik_queue[], int index_queue[], double new_valConsidered, int index);

static double get_valueAtIndex(double eik_queue[], int index_queue[], int index, int size);

