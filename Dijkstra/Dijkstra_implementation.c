/*
DIJKSTRA'S ALGORITHM IMPLEMENTATION IN C
Chech the comments done in Dijkstra_implementation.py
I followed the same idea coded there but this is in C
 */

#include <stdio.h>
#define INFINITY 9999
#define MAX 5

int Dijkstra( int vertices[MAX][MAX], int edges[MAX][MAX], int start );

int main(){
    // Same example as in the python file

    int vertices[MAX][MAX], start, i;
    int edges[MAX][MAX], distance[MAX];

    start = 0;

    // We fill these matrices

    vertices[0][0] = 0;
    vertices[0][1] = 1;
    vertices[0][2] = 1;
    vertices[0][3] = 0;
    vertices[0][4] = 0;

    vertices[1][0] = 0;
    vertices[1][1] = 0;
    vertices[1][2] = 1;
    vertices[1][3] = 0;
    vertices[1][4] = 0;

    vertices[2][0] = 0;
    vertices[2][1] = 0;
    vertices[2][2] = 0;
    vertices[2][3] = 1;
    vertices[2][4] = 1;

    vertices[3][0] = 0;
    vertices[3][1] = 1;
    vertices[3][2] = 0;
    vertices[3][3] = 0;
    vertices[3][4] = 1;

    vertices[4][0] = 0;
    vertices[4][1] = 0;
    vertices[4][2] = 0;
    vertices[4][3] = 0;
    vertices[4][4] = 0;

    edges[0][0] = 0;
    edges[0][1] = 100;
    edges[0][2] = 30;
    edges[0][3] = 0;
    edges[0][4] = 0;

    edges[1][0] = 0;
    edges[1][1] = 0;
    edges[1][2] = 20;
    edges[1][3] = 0;
    edges[1][4] = 0;

    edges[2][0] = 0;
    edges[2][1] = 0;
    edges[2][2] = 0;
    edges[2][3] = 10;
    edges[2][4] = 60;

    edges[3][0] = 0;
    edges[3][1] = 15;
    edges[3][2] = 0;
    edges[3][3] = 0;
    edges[3][4] = 50;

    edges[4][0] = 0;
    edges[4][1] = 0;
    edges[4][2] = 0;
    edges[4][3] = 0;
    edges[4][4] = 0;

    distance = Dijkstra(vertices[MAX][MAX], edges[MAX][MAX], start )

    for (i = 0, i< MAX, i++)
    {
        printf("Distance from %d, to %d is %lf", start, i, distance[i]);
    }

    return 0;
}

int Dijkstra( int vertices[MAX][MAX], int edges[MAX][MAX], int start ){
    /*
     Dijkstra's algorithm to find the shortest path between all the nodes in a 
    grapha and the initial node (start).
    Arguments:
       vertices : matrix that represents if two nodes are connected or not SQUARE MATRIX
       edges : matrix that represents the weights of the edges between the nodes SQUARE MATRIX AND SAME SIZE AS VERTICES
       start: an integer, in which node to start the path
    Returns:
       distance : the distance of the shortest path up to that node
    */
   // Initialization
   int Q[MAX], distance[MAX], previous[MAX], minDistance, i, j, u, v;

   for(i = 0, i < MAX, i++)
   {
       Q[i] = 0;
       distance[i] = INFINITY;
       previous[i] = start;
   }

   Q[start] = 1;
   distance[start] = 0;

   // Initialize the distance of the nodes that are connected to the starting node
   for(j = 0, j < MAX, j++)
   {
       if(vertices[start][j] == 1)
       {
           distance[j] = edges[start][j];
       }
   }

   // If the nodes are not connected then their distance is infinity in the original edge matrix
   // We put those infinity entries in case the user put zeroes

   for(j = 0, j < MAX, j++)
   {
       for(i = 0, i < MAX, i++)
       {
           if(vertices[i][j] == 0)
           {
               edges[i][j] = INFINITY;
           }
       }
   }

   // We iterate for each vertex in our graph

   for(v = 0, v < MAX, v++)
   {
       minDistance = INFINITY;
       // Find the next optimal node in the path
       for(u = 0, u < MAX, u++)
       {
           if(distance[u] < minDistance && Q[u] == 0)
           {
               minDistance = distance[u];
               next_node_in_path = u;
           }
       }
       // We update, we found the next node to visit
       previous[next_node_in_path] = v;
       Q[next_node_in_path] = 1;

       // Update all the other current weights of the nodes in the priority queue
       // But we must update just for the nodes that are actually the neighbours of next_node_in_path

       for(i = 0, i< MAX, i++)
       {
           if(Q[i] == 0 && (minDistance + edges[next_node_in_path][i] < distance[i]) )
           {
               // If this happens then it is a good idea to visit next_in_path and then i
               distance[i] = minDistance + edges[next_node_in_path][i];
               previous[i] = next_node_in_path;
           } 
       }

       return(distance)

   }


}

        