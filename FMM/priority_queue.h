#pragma once

typedef struct Priority_queue p_queue;

void priority_queue_alloc(p_queue **Priority_queue );

void priority_queue_dealloc(p_queue **Priority_queue );

void priority_queue_init( p_queue *p_queueImp  );

void grow_queue( p_queue *p_queueImp );

void swap_double(double *a, double *b);

void swap_int(int *a, int *b);

void heapify(p_queue *p_queueImp, int i);

void insert(p_queue *p_queueImp, double newNum, int newIndex);

void insert_end(p_queue *p_queueImp, double newNum, int newIndex);

void delete_findValue(p_queue *p_queueImp, double num);

void delete_findIndex(p_queue *p_queueImp, int ind);

void deleteRoot(p_queue *p_queueImp);

void printeik_queue(p_queue *p_queueImp);

void update(p_queue *p_queueImp, double new_valConsidered, int index);

double get_valueAtIndex(p_queue *p_queueImp, int index);

int getSize(p_queue *Priority_queue);

int getIndicesInQueue(p_queue *Priority_queue);
