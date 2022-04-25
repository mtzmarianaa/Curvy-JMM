#pragma once


static void FMM_2D( double x_min, double y_min, int start[2], double *distance, int *Q, int M, int N, double h);

static double speed(double x, double y);

static void ActualSolution(double x_min, double y_min, int start[2], double *trueSolution, int M, int N, double h);

static void printQGridFromQueue(int *Q, int M, int N);

static void printGridFromDistance(double *distance, int M, int N);

static void generateDataForPlot(double *x, int M, int N, char data_name[10]);