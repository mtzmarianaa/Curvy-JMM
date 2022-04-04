/*
FMM (Fast Marching Method) naive implementation in 2D
Following (kind of) the idea of Dijkstra's algorithm 
which I also implemented.
 */

#include <stdio.h>
//#include <math.h> 

#define INFINITY 1e5
// Define the sizes of the grids
#define M 5
#define N 5

void FMM_2D( double x_min, double x_max, double y_min, double y_max, int start[2], double distance[M][N] );

double speed(double point[2]);


int main(){
    double x_min, x_max, y_min, y_max, distance[M][N];
    int start[2];
    x_min = 0;
    x_max = 10;
    y_min = 0;
    y_max = 10;
    start[0] = 0;
    start[1] = 1;

    FMM_2D( x_min, x_max, y_min, y_max, start, distance );

    return 0;
}

double speed(double point[2]){
    /*
    Function that defines the 1/f function used in the eikonal equation |u| = 1/f
    */
    return 1.0;
}


void FMM_2D( double x_min, double x_max, double y_min, double y_max, int start[2], double distance[M][N] ){
    /*
     Naive implementation of the fast marching method in a 2D grid. In the priority queue:
         - 0: far
         1: trial
         2: valid
     Hence the objective is to have the priority queue with just 2's. 
     IN:  x_min : from where in the x direction
             - x_max : to where in the x direction
             - y_min : from where in the y direction
             - y_max : to where in the y direction
             - start : two coordinates (from 0 to M-1 and from 0 to N-1) on where to start
        
    Now this is a 2D grid and we're going to "flatten" it in the row by row way, starting from the top and moving all the way 
    to the bottom (from southwst to northeast). Hence our vector is going to be a MxN vector (which represents the grid)
    */
     // Initialization: alive points (the starting point(s)), trial points (narrow band points), and far away points
     int Q[M*N], i;
     double minDistance;

     for (i = 0; i< M*N; i++)
     {
         Q[i] = 0;
     }
     // Add the starting point
     Q[start[0]*N + start[1] ] = 2;
     distance[start[0]][start[1]] = 0;
     // Initialize the points that are connected to that starting point (at most 4 points are connected to this starting point)
     // We need to know if we have all those 4 points or now (we don't want -1 as an index)
     if (start[0] == 0 ){ // most south of the grid
         if(start[1] == 0){ // if this happens we're in the most southwest corner
             Q[1] = 1;
             Q[N] = 1;
         }
         else if (start[1] == (N-1) ){ // if this happens we're in the most southeast corner
             Q[N-2] = 1;
             Q[ 2*N - 1 ] = 1;
         }
         else { // if this happens we're just in a non corner south edge of the grid
             Q[start[1] + 1] = 1;
             Q[start[1] - 1] = 1;
             Q[ N + start[1] ] = 1;
         }
     }
     if (start[0] = (M-1) ){ // most north of the grid
             if(start[1] == 0){ // most northwest point
                 Q[N*(M-1) + 1] = 1;
                 Q[(N-1)*(M-1)] = 1;
             }
             else if (start[1] = (N-1)){ // most northeast point
                 Q[N*M -2] = 1;
                 Q[ () ]
             }
     }

     minDistance = 0.0;
     
     }