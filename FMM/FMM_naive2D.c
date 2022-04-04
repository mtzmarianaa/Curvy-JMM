/*
FMM (Fast Marching Method) naive implementation in 2D
Following (kind of) the idea of Dijkstra's algorithm 
which I also implemented.
 */

#include <stdio.h>
#include <math.h> 
// Define the sizes of the grids
#define M 5
#define N 5

void FMM_2D( double x_min, double x_max, double y_min, double y_max, int start[2], double distance[M][N] );

double speed(double x, double y);

void printQGridFromQueue(int Q[M*N]);

void printGridFromDistance(double distance[M][N]);


int main(){
    double x_min, x_max, y_min, y_max, distance[M][N];
    int start[2];
    x_min = 0.0;
    x_max = 10.0;
    y_min = 0.0;
    y_max = 10.0;
    start[0] = 2;
    start[1] = 2;

    FMM_2D( x_min, x_max, y_min, y_max, start, distance );

    return 0;
}

double speed(double x, double y){
    /*
    Function that defines the 1/f function used in the eikonal equation |u| = 1/f
    */
    return 1.0;
}

void printQGridFromQueue(int Q[M*N]){
    for(int i = N-1; i>-1; i --){
        for(int j = 0; j<M; j++){
            printf("%d ", Q[i*N + j] );
        }
        printf("\n");
    }
}

void printGridFromDistance(double distance[M][N]){
    for(int j = 0; j<N; j ++){
        for(int i = 0; i<M; i++){
            printf("%lf ", distance[i][j]);
        }
        printf("\n");
    }
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
   // INITIALIZATION PART
     // alive points (the starting point(s)), trial points (narrow band points), and far away points
     int Q[M*N], i, j;
     double minDistance, x_linspace[N], y_linspace[M], h_x, h_y;

     // stepsize in x direction is h_x, stepsize in y direction is h_y
     h_x = (x_max - x_min)/N;
     h_y = (y_max - y_min)/N;

     // Start the linspace in the x and y directions
     for(i = 0; i<N; i++){
         x_linspace[i] = x_min + i*h_x;
     }
     for(i = 0; i<M; i++){
         y_linspace[i] = y_min + i*h_y;
     }
    // Initialize the queue
     for (i = 0; i< M*N; i++)
     {
         Q[i] = 0;
     }

     // Initialize the distances with infinity
     for(j = 0; j< M; j ++){
         for(i = 0; i<N; i ++){
             distance[j][i] = INFINITY;
         }
     }


     // Add the starting point
     Q[start[0]*N + start[1] ] = 2;
     distance[start[0]][start[1]] = 0;

     // Initialize the points that are connected to that starting point (at most 4 points are connected to this starting point)
     // We need to know if we have all those 4 points or now (we don't want -1 as an index)
     // We also know that these are the only connected points to our start point so we can initialize here
     // the distance of those points.
     if (start[0] !=0  ){ // we're not in the most southern part, we do have a neighbour south
         Q[N*(start[0]-1) + start[1] ] = 1;
         distance[start[0]][start[1] - 1 ] = speed( x_linspace[start[0]], y_linspace[start[1]]-1 )*h_y;
     }
     if ( start[1] != 0 ){ // we're not in a west edge, we can have a neighbour to the left
         Q[ N*start[0] + start[1] -1 ] = 1;
         distance[start[0]-1][start[1]] = speed( x_linspace[start[0] -1 ], y_linspace[ start[1] ] )*h_x;
     }
     if (start[1] != (N-1) ){ // we're not in a east edge, we can have a neighbour to the right
         Q[N*start[0] + start[1] + 1] = 1;
         distance[ start[0] + 1 ][start[1]] = speed( x_linspace[start[0] + 1], y_linspace[start[1]] )*h_x;
     }
     if (start[0] != (M-1)){ // we're not in a north edge, we can have a neighbour north
         Q[ N*(start[0]+1) + start[1] ] = 1;
         distance[start[0] ][start[1]+1] = speed( x_linspace[start[0]], y_linspace[start[1] + 1] )*h_y;
     }
     //printQGridFromQueue(Q); // in case we need this
     //printGridFromDistance(distance);



     // ITERATION PART



     minDistance = 0.0;
     
     }