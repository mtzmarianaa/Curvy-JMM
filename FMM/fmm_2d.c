#include "fmm_2d.h"
#include "eik_grid.c"

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

static double onePointUpdate(  eik_gridS *eikonal_g, int coordinate  ){
    int x_coor, y_coor;
    x_coor = coordinate%eikonal_g->N;
    y_coor = coordinate/eikonal_g->N;
    return speed( eikonal_g->x_linspace[x_coor], eikonal_g->y_linspace[y_coor] )*eikonal_g->h;
}

static void generateDataForPlot(double *x, int M, int N, char data_name[10]){
    char ext[] = ".bin";
	FILE *fp = fopen(strcat(data_name, ext), "wb");
	fwrite(x, sizeof(double), M*N, fp);
	fclose(fp);

}


void FMM_2D( eik_gridS *eikonal_g ){
    /*
     Naive implementation of the fast marching method in a 2D grid. In the priority queue:
         - 0: far
         - 1: trial
         - 2: valid
     Hence the objective is to have the priority queue with just 2's. 
     Remember that eik_gridS has:
         double x0;               x coordinate of the initial point
         double y0;               y coordinate of the initial point
         int start[2];            indices of the initial point
         int M;                   length of the grid
         int N;                   width of the grid
         double h;                step size
         double *x_linspace;      x coordinates of the nodes (length N)
         double *y_linspace;      y coordinates of the nodes (length M)
         double *eik_gridVals;    Eikonal values of the nodes on the grid (flatten grid)
         p_queue *p_queueG;       both Eikonal values of the nodes on the queue and their indices
         int *current_states;     array with 0,1 or 2 depending on the current state of all the nodes
        
    Now this is a 2D grid and we're going to "flatten" it in the row by row way, starting from the top and moving all the way 
    to the bottom (from southwst to northeast). Hence our vector is going to be a MxN vector (which represents the grid)
    */
   // INITIALIZATION PART IS GIVEN WHEN INITIALIZING THE eik_gridS struct
     // alive points (the starting point(s)), trial points (narrow band points), and far away points
     int i, next_valid, coordinate;
     double minDistance, u1, two_point1, two_point2, one_point;

     // Add the starting point
     eikonal_g->current_states[ eikonal_g->start[0]*eikonal_g->N + eikonal_g->start[1] ] = 2; // set the starting point as valid
     insert( eikonal_g->p_queueG, 0.0,  eikonal_g->start[0]*eikonal_g->N + eikonal_g->start[1]  ); // set the starting point's Eikonal value to 0
     deleteRoot(eikonal_g->p_queueG); // delete the root because we have set this point to valid
     eikonal_g->eik_gridVals[ eikonal_g->start[0]*eikonal_g->N + eikonal_g->start[1] ] = 0.0;

     // Initialize the points that are connected to that starting point (at most 4 points are connected to this starting point)
     // We need to know if we have all those 4 points or now (we don't want -1 as an index)
     // We also know that these are the only connected points to our start point so we can initialize here
     // the distance of those points.
     if (eikonal_g->start[0] !=0  ){ // we're not in the most southern part, we do have a neighbour south
         eikonal_g->current_states[eikonal_g->N *(eikonal_g->start[0]-1) + eikonal_g->start[1] ] = 1; // update current states
         insert( eikonal_g->p_queueG,  speed( eikonal_g->x_linspace[eikonal_g->start[0]], eikonal_g->y_linspace[eikonal_g->start[1]]-1 )*eikonal_g->h, (eikonal_g->N*(eikonal_g->start[0]-1) + eikonal_g->start[1] )  );
     }
     if ( eikonal_g->start[1] != 0 ){ // we're not in a west edge, we can have a neighbour to the left
         eikonal_g->current_states[eikonal_g->N *eikonal_g->start[0] + eikonal_g->start[1] -1 ] = 1; // update current states
         insert( eikonal_g->p_queueG,  speed( eikonal_g->x_linspace[eikonal_g->start[0] - 1], eikonal_g->y_linspace[eikonal_g->start[1]] )*eikonal_g->h, (eikonal_g->N*(eikonal_g->start[0]) + eikonal_g->start[1] -1 )  );
         //current_states[ N*start[0] + start[1] -1 ] = 1;
         //update(eik_queue, index_queue, speed( x_linspace[start[0] -1 ], y_linspace[ start[1] ] )*h, (N*start[0] + start[1] -1));
         //distance[ N*start[0] + start[1] -1 ] = speed( x_linspace[start[0] -1 ], y_linspace[ start[1] ] )*h;
     }
     if (eikonal_g->start[1] != (eikonal_g->N-1) ){ // we're not in a east edge, we can have a neighbour to the right
         eikonal_g->current_states[eikonal_g->N *eikonal_g->start[0] + eikonal_g->start[1] + 1 ] = 1; // update current states
         insert( eikonal_g->p_queueG,  speed( eikonal_g->x_linspace[eikonal_g->start[0] + 1], eikonal_g->y_linspace[eikonal_g->start[1]] )*eikonal_g->h, (eikonal_g->N*(eikonal_g->start[0]) + eikonal_g->start[1] +1 )  );
         //current_states[N*start[0] + start[1] + 1] = 1;
         //update(eik_queue, index_queue, speed( x_linspace[start[0] + 1], y_linspace[start[1]] )*h, (N*start[0] + start[1] + 1));
         //distance[N*start[0] + start[1] + 1] = speed( x_linspace[start[0] + 1], y_linspace[start[1]] )*h;
     }
     if (eikonal_g->start[0] != (eikonal_g->M-1)){ // we're not in a north edge, we can have a neighbour north
         
         eikonal_g->current_states[eikonal_g->N *(eikonal_g->start[0]+1) + eikonal_g->start[1] ] = 1; // update current states
         insert( eikonal_g->p_queueG,  speed( eikonal_g->x_linspace[eikonal_g->start[0]], eikonal_g->y_linspace[eikonal_g->start[1]+1] )*eikonal_g->h, (eikonal_g->N*(eikonal_g->start[0]+1) + eikonal_g->start[1] )  );
         //current_states[ N*(start[0]+1) + start[1] ] = 1;
         //update(eik_queue, index_queue, speed( x_linspace[start[0]], y_linspace[start[1] + 1] )*h, (N*(start[0]+1) + start[1]));
         //distance[ N*(start[0]+1) + start[1] ] = speed( x_linspace[start[0]], y_linspace[start[1] + 1] )*h;
     }

     printeik_queue(eikonal_g->p_queueG);
     print_eikonal_grid(eikonal_g);
     print_currentStates(eikonal_g);

     
     // ITERATION PART

     // We iterate for each point in our grid
     for(i = 0; i < eikonal_g->M*eikonal_g->N; i ++)
     {
          minDistance = INFINITY;
          // Find the next optimal point (the one with minimum distance)
          minDistance = eikonal_g->p_queueG->queue_vals[0]; // thanks to the binary tree we know its the first element
          next_valid = eikonal_g->p_queueG->queue_index[0]; // the same with the index queue

          // Now that we've found the next point to put on valid, we update
          eikonal_g->current_states[next_valid] = 2; // set this as valid
          deleteRoot(eikonal_g->p_queueG); // delete the root of the priority queue (binary tree)
          // Now we need to add its neighbours, mark them as trial but we know where their neighbours are
          // but we just need to mark them as trial if they are currently marked as far
          // ADD FAR NEIGHBOURS TO PRIORITY QUEUE
          if( next_valid >= eikonal_g->N && eikonal_g->current_states[next_valid - eikonal_g->N] == 0 ){ // if this happens then the next valid point is not in the most southern part of the grid
              eikonal_g->current_states[next_valid - eikonal_g->N] = 1;
              insert(eikonal_g->p_queueG, onePointUpdate( eikonal_g,  next_valid - eikonal_g->N ) , next_valid - eikonal_g->N);
          }
          if ( next_valid%eikonal_g->N != 0 && eikonal_g->current_states[next_valid -1] == 0  ){ // if this happens then the next valid point is not in the west edge
              eikonal_g->current_states[next_valid -1] = 1;
              insert(eikonal_g->p_queueG, onePointUpdate( eikonal_g,  next_valid -1 ) , next_valid -1);
          }
          if ( next_valid%eikonal_g->N != eikonal_g->N-1 && eikonal_g->current_states[next_valid + 1] == 0 ){ // if this happens then the next valid point is not in the east edge
              eikonal_g->current_states[next_valid + 1] = 1;
              insert(eikonal_g->p_queueG, onePointUpdate( eikonal_g,  next_valid +1 ) , next_valid +1);
          }
          if(next_valid/eikonal_g->N < eikonal_g->M-1 && eikonal_g->current_states[next_valid + eikonal_g->N] == 0 ){ // if this happens then the next valid point is not in the northern edge of the grid
              eikonal_g->current_states[next_valid + eikonal_g->N] = 1;
              insert(eikonal_g->p_queueG, onePointUpdate( eikonal_g,  next_valid + eikonal_g->N ) , next_valid + eikonal_g->N);
          }
          //printeik_queue(eikonal_g->p_queueG);
          //print_eikonal_grid(eikonal_g);
          //print_currentStates(eikonal_g);

          // UPDATE POINT IN THE 9 GRID STENCIL USING EITHER 2 POINT UPDATES OF 1 POINT UPDATES
         
          u1 = eikonal_g->eik_gridVals[next_valid];

          // next valid point is not a southern edge + this neighbour is not valid currently
          if( next_valid >= eikonal_g->N && eikonal_g->current_states[next_valid - eikonal_g->N] != 2  ){
              coordinate = next_valid - eikonal_g->N;
              if ( coordinate%eikonal_g->N != 0 ){ // its not in the west edge as well
                  two_point1 = twoPointUpdate(minDistance, eikonal_g->eik_gridVals[next_valid - eikonal_g->N - 1], eikonal_g->h, coordinate, eikonal_g->N);
              }
              else {
                  two_point1 = INFINITY;
              }
              if (coordinate%eikonal_g->N != eikonal_g->N-1 ){
                  two_point2 = twoPointUpdate(minDistance, eikonal_g->eik_gridVals[next_valid - eikonal_g->N + 1], eikonal_g->h, coordinate, eikonal_g->N);
              }
              else{
                  two_point2 = INFINITY;
              }
              one_point = eikonal_g->eik_gridVals[next_valid] + eikonal_g->h*speed( coordinate%eikonal_g->N,  (coordinate - coordinate%eikonal_g->N)/eikonal_g->N );
              if(eikonal_g->eik_gridVals[next_valid-eikonal_g->N] > two_point1){
                  eikonal_g->eik_gridVals[next_valid-eikonal_g->N] = two_point1;
              }
              if(eikonal_g->eik_gridVals[next_valid-eikonal_g->N] > two_point2){
                  eikonal_g->eik_gridVals[next_valid-eikonal_g->N] = two_point2;
              }
              if(eikonal_g->eik_gridVals[next_valid-eikonal_g->N] > one_point){
                  eikonal_g->eik_gridVals[next_valid-eikonal_g->N] = one_point;
              }
              update(eikonal_g->p_queueG, eikonal_g->eik_gridVals[next_valid-eikonal_g->N], next_valid-eikonal_g->N);
          }
          
          // next valid point is not a western edge + this neighbour is not valid currently
          if( next_valid%eikonal_g->N != 0 && eikonal_g->current_states[next_valid -1] != 2  ){
              coordinate = next_valid - 1;
              if( coordinate>= eikonal_g->N ){
                  two_point1 = twoPointUpdate(minDistance, eikonal_g->eik_gridVals[next_valid - eikonal_g->N - 1], eikonal_g->h, coordinate, eikonal_g->N);
              }
              else{
                  two_point1 = INFINITY;
              }
              if ( coordinate/eikonal_g->N < eikonal_g->M-1 ){
                  two_point2 = twoPointUpdate(minDistance, eikonal_g->eik_gridVals[next_valid + eikonal_g->N - 1], eikonal_g->h, coordinate, eikonal_g->N);
              }
              else{
                  two_point2 = INFINITY;
              }
              one_point = eikonal_g->eik_gridVals[next_valid] + eikonal_g->h*speed( coordinate%eikonal_g->N,  (coordinate - coordinate%eikonal_g->N)/eikonal_g->N );
              if(eikonal_g->eik_gridVals[next_valid-1] > two_point1){
                  eikonal_g->eik_gridVals[next_valid-1] = two_point1;
              }
              if(eikonal_g->eik_gridVals[next_valid-1] > two_point2){
                  eikonal_g->eik_gridVals[next_valid-1] = two_point2;
              }
              if(eikonal_g->eik_gridVals[next_valid-1] > one_point){
                  eikonal_g->eik_gridVals[next_valid-1] = one_point;
              }
              update(eikonal_g->p_queueG, eikonal_g->eik_gridVals[next_valid-1], next_valid-1);
          }

          // next valid point is not a eastern edge + this neighbour is not valid currently
          if( next_valid%eikonal_g->N != eikonal_g->N-1 && eikonal_g->current_states[next_valid + 1] != 2  ){
              coordinate = next_valid + 1;
              if ( coordinate>= eikonal_g->N  ){
                  two_point1 = twoPointUpdate(minDistance, eikonal_g->eik_gridVals[next_valid - eikonal_g->N + 1], eikonal_g->h, coordinate, eikonal_g->N);
              }
              else{
                  two_point1 = INFINITY;
              }
              if( coordinate/eikonal_g->N < eikonal_g->M-1 ){
                  two_point2 = twoPointUpdate(minDistance, eikonal_g->eik_gridVals[next_valid + eikonal_g->N + 1], eikonal_g->h, coordinate, eikonal_g->N);
              }
              else{
                  two_point2 = INFINITY;
              }
              one_point = eikonal_g->eik_gridVals[next_valid] + eikonal_g->h*speed( coordinate%eikonal_g->N,  (coordinate - coordinate%eikonal_g->N)/eikonal_g->N );
              if(eikonal_g->eik_gridVals[next_valid+1] > two_point1){
                  eikonal_g->eik_gridVals[next_valid+1] = two_point1;
              }
              if(eikonal_g->eik_gridVals[next_valid+1] > two_point2){
                  eikonal_g->eik_gridVals[next_valid+1] = two_point2;
              }
              if(eikonal_g->eik_gridVals[next_valid+1] > one_point){
                  eikonal_g->eik_gridVals[next_valid+1] = one_point;
              }
              update(eikonal_g->p_queueG, eikonal_g->eik_gridVals[next_valid+1], (next_valid+1) );
          }

          // next valid point is not a northern edge + this neighbour is not valid currently
          if( next_valid/eikonal_g->N < eikonal_g->M-1 && eikonal_g->current_states[next_valid + eikonal_g->N] != 2   ){
              coordinate = next_valid - eikonal_g->N;
              if (coordinate%eikonal_g->N != 0 ){
                  two_point1 = twoPointUpdate(minDistance, eikonal_g->eik_gridVals[next_valid + eikonal_g->N - 1], eikonal_g->h, coordinate, eikonal_g->N);
              }
              else{
                  two_point1 = INFINITY;
              }
              if (coordinate%eikonal_g->N != eikonal_g->N-1){
                  two_point2 = twoPointUpdate(minDistance, eikonal_g->eik_gridVals[next_valid + eikonal_g->N + 1], eikonal_g->h, coordinate, eikonal_g->N);
              }
              else{
                  two_point2 = INFINITY;
              }
              one_point = eikonal_g->eik_gridVals[next_valid] + eikonal_g->h*speed( coordinate%eikonal_g->N,  (coordinate - coordinate%eikonal_g->N)/eikonal_g->N );
              if(eikonal_g->eik_gridVals[next_valid+eikonal_g->N] > two_point1){
                  eikonal_g->eik_gridVals[next_valid+eikonal_g->N] = two_point1;
              }
              if(eikonal_g->eik_gridVals[next_valid+eikonal_g->N] > two_point2){
                  eikonal_g->eik_gridVals[next_valid+eikonal_g->N] = two_point2;
              }
              if(eikonal_g->eik_gridVals[next_valid+eikonal_g->N] > one_point){
                  eikonal_g->eik_gridVals[next_valid+eikonal_g->N] = one_point;
              }
              update(eikonal_g->p_queueG, eikonal_g->eik_gridVals[next_valid+eikonal_g->N], (next_valid+eikonal_g->N));
          }

          //printGridFromDistance(distance, M, N);


     }



}