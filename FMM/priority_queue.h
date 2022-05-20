#pragma once

typedef struct Priority_queue p_queue;

void priority_queue_init( p_queue *p_queueImp  );

void priority_queue_deinit( p_queue *p_queueImp );

void grow_queue( p_queue *p_queueImp );

static void swap_double(double *a, double *b);

static void swap_int(int *a, int *b);

static void heapify(p_queue *p_queueImp, int i);

static void insert(p_queue *p_queueImp, double newNum, int newIndex);

static void insert_end(p_queue *p_queueImp, double newNum, int newIndex);

static void delete_findValue(p_queue *p_queueImp, double num);

static void delete_findIndex(p_queue *p_queueImp, int ind);

static void deleteRoot(p_queue *p_queueImp);

static void printeik_queue(p_queue *p_queueImp);

/*

static void update(double eik_queue[], int index_queue[], double new_valConsidered, int index);

static double get_valueAtIndex(double eik_queue[], int index_queue[], int index, int size);

*/