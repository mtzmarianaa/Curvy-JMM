#include "fmm_2d.h"
#include "priority_queue.c"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

static double speed(double x, double y){
    /*
    Function that defines the 1/f function used in the eikonal equation |u| = 1/f
    */
    return 1;
}

static void ActualSolution(double x_min, double y_min, int start[2], double *trueSolution, int M, int N, double h){
    /*
    This is the true solution of constant speed of sound.
    */
    double x_linspace[N], y_linspace[M], x_0, y_0;
    int i;

    for(i = 0; i<N; i++){
         x_linspace[i] = x_min + i*h;
     }
     for(i = 0; i<M; i++){
         y_linspace[i] = y_min + i*h;
     }

     x_0 = x_linspace[start[0]];
     y_0 = y_linspace[start[1]];

     for (i = 0; i < N*M; i++){
         trueSolution[i] = sqrt(  pow( x_linspace[i%N] - x_0, 2  ) +  pow( y_linspace[i/N] - y_0, 2  )  );
     }

}

static double twoPointUpdate(double u1, double u2, double h, int coordinate, int N){
    /*
    Function that calculate a two point Eikonial update
    */
   int i, j;
   double d;
   i = coordinate%N;
   j = (coordinate - i)/N;
   d = 0.5*(u1+u2) + 0.5*sqrt( pow(u1 + u2, 2) - 2*(pow(u1, 2) + pow(u2, 2) -  pow(h,2)/speed(i, j) )  );
   return d;
}

static void printQGridFromQueue(int *Q, int M, int N){
    for(int i = N-1; i>-1; i --){
        for(int j = 0; j<M; j++){
            printf("%d ", Q[i*N + j] );
        }
        printf("\n");
    }
    printf("\n");
}

static void printGridFromDistance(double *distance, int M, int N){
    for(int i = N-1; i>-1; i --){
        for(int j = 0; j<M; j++){
            printf("%fl ", distance[i*N + j] );
        }
        printf("\n");
    }
    printf("\n");
}

static void generateDataForPlot(double *x, int M, int N, char data_name[10]){
    char ext[] = ".bin";
	FILE *fp = fopen(strcat(data_name, ext), "wb");
	fwrite(x, sizeof(double), M*N, fp);
	fclose(fp);

}


static void FMM_2D( double x_min, double y_min, int start[2], double *distance, double *eik_queue, int *index_queue, int *current_states, int M, int N, double h){
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
     int i, j, next_valid, coordinate;
     double minDistance, x_linspace[N], y_linspace[M], u1, two_point1, two_point2, one_point;


     // Start the linspace in the x and y directions
     for(i = 0; i<N; i++){
         x_linspace[i] = x_min + i*h;
     }
     for(i = 0; i<M; i++){
         y_linspace[i] = y_min + i*h;
     }
    // Initialize the queue and the distances with infinity
     for (i = 0; i< M*N; i++)
     {
         current_states[i] = 0; // all the nodes are far
         index_queue[i] = i; // all the nodes go into the priority queue (indeces)
         insert(eik_queue, index_queue, INFINITY, i);
         distance[i] = INFINITY;
     }


     // Add the starting point
     current_states[start[0]*N + start[1] ] = 2;
     update(eik_queue, index_queue, 0.0, (start[0]*N + start[1])); // add the point
     deleteRoot(eik_queue, index_queue); // and delete the root, it has been accepted
     distance[start[0]*N + start[1] ] = 0;

     // Initialize the points that are connected to that starting point (at most 4 points are connected to this starting point)
     // We need to know if we have all those 4 points or now (we don't want -1 as an index)
     // We also know that these are the only connected points to our start point so we can initialize here
     // the distance of those points.
     if (start[0] !=0  ){ // we're not in the most southern part, we do have a neighbour south
         current_states[N*(start[0]-1) + start[1] ] = 1; // update current states
         update(eik_queue, index_queue, speed( x_linspace[start[0]], y_linspace[start[1]]-1 )*h, (N*(start[0]-1) + start[1] )); // update the current value in the priority eikonal queue
         distance[N*(start[0]-1) + start[1] ] = speed( x_linspace[start[0]], y_linspace[start[1]]-1 )*h;
     }
     if ( start[1] != 0 ){ // we're not in a west edge, we can have a neighbour to the left
         current_states[ N*start[0] + start[1] -1 ] = 1;
         update(eik_queue, index_queue, speed( x_linspace[start[0] -1 ], y_linspace[ start[1] ] )*h, (N*start[0] + start[1] -1));
         distance[ N*start[0] + start[1] -1 ] = speed( x_linspace[start[0] -1 ], y_linspace[ start[1] ] )*h;
     }
     if (start[1] != (N-1) ){ // we're not in a east edge, we can have a neighbour to the right
         current_states[N*start[0] + start[1] + 1] = 1;
         update(eik_queue, index_queue, speed( x_linspace[start[0] + 1], y_linspace[start[1]] )*h, (N*start[0] + start[1] + 1));
         distance[N*start[0] + start[1] + 1] = speed( x_linspace[start[0] + 1], y_linspace[start[1]] )*h;
     }
     if (start[0] != (M-1)){ // we're not in a north edge, we can have a neighbour north
         current_states[ N*(start[0]+1) + start[1] ] = 1;
         update(eik_queue, index_queue, speed( x_linspace[start[0]], y_linspace[start[1] + 1] )*h, (N*(start[0]+1) + start[1]));
         distance[ N*(start[0]+1) + start[1] ] = speed( x_linspace[start[0]], y_linspace[start[1] + 1] )*h;
     }
     //printQGridFromQueue(Q, M, N); // in case we need this
     //printGridFromDistance(distance, M, N);

     // ITERATION PART

     // We iterate for each point in our grid
     for(i = 0; i < M*N; i ++)
     {
          minDistance = INFINITY;
          // Find the next optimal point (the one with minimum distance)
          minDistance = eik_queue[0]; // thanks to the binary tree we know its the first element
          next_valid = index_queue[0]; // the same with the index queue

          // Now that we've found the next point to put on valid, we update
          current_states[next_valid] = 2; // set this as valid
          deleteRoot(eik_queue, index_queue); // delete the root of the priority queue (binary tree)
          //printQGridFromQueue(Q); // Check
          // Now we need to add its neighbours, mark them as trial but we know where their neighbours are
          // but we just need to mark them as trial if they are currently marked as far
          // ADD NON VALID NEIGHBOURS TO PRIORITY QUEUE
          if( next_valid >= N && current_states[next_valid - N] == 0 ){ // if this happens then the next valid point is not in the most southern part of the grid
              current_states[next_valid - N] = 1;
          }
          if ( next_valid%N != 0 && current_states[next_valid -1] == 0  ){ // if this happens then the next valid point is not in the west edge
              current_states[next_valid -1] = 1;
          }
          if ( next_valid%N != N-1 && current_states[next_valid + 1] == 0 ){ // if this happens then the next valid point is not in the east edge
              current_states[next_valid + 1] = 1;
          }
          if(next_valid/N < M-1 && current_states[next_valid + N] == 0 ){ // if this happens then the next valid point is not in the northern edge of the grid
              current_states[next_valid + N] = 1;
          }
          //printQGridFromQueue(Q, M, N); // Check

          // UPDATE POINT IN THE 9 GRID STENCIL USING EITHER 2 POINT UPDATES OF 1 POINT UPDATES
          u1 = distance[next_valid];

          // next valid point is not a southern edge + this neighbour is not valid currently
          if( next_valid >= N && current_states[next_valid - N] != 2  ){
              coordinate = next_valid - N;
              if ( coordinate%N != 0 ){ // its not in the west edge as well
                  two_point1 = twoPointUpdate(minDistance, distance[next_valid - N - 1], h, coordinate, N);
              }
              else {
                  two_point1 = INFINITY;
              }
              if (coordinate%N != N-1 ){
                  two_point2 = twoPointUpdate(minDistance, distance[next_valid - N + 1], h, coordinate, N);
              }
              else{
                  two_point2 = INFINITY;
              }
              one_point = distance[next_valid] + h*speed( coordinate%N,  (coordinate - coordinate%N)/N );
              if(distance[next_valid-N] > two_point1){
                  distance[next_valid-N] = two_point1;
              }
              if(distance[next_valid-N] > two_point2){
                  distance[next_valid-N] = two_point2;
              }
              if(distance[next_valid-N] > one_point){
                  distance[next_valid-N] = one_point;
              }
              update(eik_queue, index_queue, distance[next_valid-N], next_valid-N);
          }
          
          // next valid point is not a western edge + this neighbour is not valid currently
          if( next_valid%N != 0 && current_states[next_valid -1] != 2  ){
              coordinate = next_valid - 1;
              if( coordinate>= N ){
                  two_point1 = twoPointUpdate(minDistance, distance[next_valid - N - 1], h, coordinate, N);
              }
              else{
                  two_point1 = INFINITY;
              }
              if ( coordinate/N < M-1 ){
                  two_point2 = twoPointUpdate(minDistance, distance[next_valid + N - 1], h, coordinate, N);
              }
              else{
                  two_point2 = INFINITY;
              }
              one_point = distance[next_valid] + h*speed( coordinate%N,  (coordinate - coordinate%N)/N );
              if(distance[next_valid-1] > two_point1){
                  distance[next_valid-1] = two_point1;
              }
              if(distance[next_valid-1] > two_point2){
                  distance[next_valid-1] = two_point2;
              }
              if(distance[next_valid-1] > one_point){
                  distance[next_valid-1] = one_point;
              }
              update(eik_queue, index_queue, distance[next_valid-1], next_valid-1);
          }

          // next valid point is not a eastern edge + this neighbour is not valid currently
          if( next_valid%N != N-1 && current_states[next_valid + 1] != 2  ){
              coordinate = next_valid + 1;
              if ( coordinate>= N  ){
                  two_point1 = twoPointUpdate(minDistance, distance[next_valid - N + 1], h, coordinate, N);
              }
              else{
                  two_point1 = INFINITY;
              }
              if( coordinate/N < M-1 ){
                  two_point2 = twoPointUpdate(minDistance, distance[next_valid + N + 1], h, coordinate, N);
              }
              else{
                  two_point2 = INFINITY;
              }
              one_point = distance[next_valid] + h*speed( coordinate%N,  (coordinate - coordinate%N)/N );
              if(distance[next_valid+1] > two_point1){
                  distance[next_valid+1] = two_point1;
              }
              if(distance[next_valid+1] > two_point2){
                  distance[next_valid+1] = two_point2;
              }
              if(distance[next_valid+1] > one_point){
                  distance[next_valid+1] = one_point;
              }
              update(eik_queue, index_queue, distance[next_valid+1], (next_valid+1) );
          }

          // next valid point is not a northern edge + this neighbour is not valid currently
          if( next_valid/N < M-1 && current_states[next_valid + N] != 2   ){
              coordinate = next_valid - N;
              if (coordinate%N != 0 ){
                  two_point1 = twoPointUpdate(minDistance, distance[next_valid + N - 1], h, coordinate, N);
              }
              else{
                  two_point1 = INFINITY;
              }
              if (coordinate%N != N-1){
                  two_point2 = twoPointUpdate(minDistance, distance[next_valid + N + 1], h, coordinate, N);
              }
              else{
                  two_point2 = INFINITY;
              }
              one_point = distance[next_valid] + h*speed( coordinate%N,  (coordinate - coordinate%N)/N );
              if(distance[next_valid+N] > two_point1){
                  distance[next_valid+N] = two_point1;
              }
              if(distance[next_valid+N] > two_point2){
                  distance[next_valid+N] = two_point2;
              }
              if(distance[next_valid+N] > one_point){
                  distance[next_valid+N] = one_point;
              }
              update(eik_queue, index_queue, distance[next_valid+N], (next_valid+N));
          }

          //printGridFromDistance(distance, M, N);


     }

}