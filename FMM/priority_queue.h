#pragma once

static void swap(int *a, int *b);

static void heapify(int eik_queue[], int index_queue[], int size, int i);

static void insert(int eik_queue[], int index_queue[], int newNum, int newIndex);

static void insert_end(int eik_queue[], int index_queue[], int newNum, int newIndex);

static void delete_findValue(int eik_queue[], int index_queue[], int num);

static void delete_findIndex(int eik_queue[], int index_queue[], int ind);

static void deleteRoot(int eik_queue[], int index_queue[]);

static void printeik_queue(int eik_queue[], int index_queue[], int size);

static void update(int eik_queue[], int index_queue[], int new_valConsidered, int index);

static int get_valueAtIndex(int eik_queue[], int index_queue[], int index, int size);