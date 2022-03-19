#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Implementation of Dijkstra's algorithm

"""

# To look for an explanation read the pdf that is in the Notizen folder.
# Dijkstra's algorithm is used to find the shortes path in a weighted graph
# given an initial starting point.
# We're going to assume that the graph is given in two matrices (this is because
# we're tying to avoid using objects). The first matrix represents if two nodes
# are connected or not. (This in the FMM is going to be much easier hopefully). 
# So in this first graph 1 in the ij entry means that the two nodes are connected
# and 0 represents that the nodes i and j are not connected. For example,
# in the first row of such matrix, if there is a 0 in the second entry it means
# that the first node has no path to the second node. Notice that this matrix 
# doesn't have to be symmetric because such graphs are directed.
# The second graph represents the weights of the edges between the nodes.
# OBJECTIVE : USE THE LEAST AMOUNT OF LIBRARIES IN PYTHON/ALREADY PROGRAMMED FUNCTIONS

infinity = float('inf') # We need to define "infinity" and initialize the distance with this

def Dijstra(vertices, edges, start):
    '''
    Dijkstra's algorithm to find the shortest path between all the nodes in a 
    grapha and the initial node (start).
    Arguments:
       vertices : matrix that represents if two nodes are connected or not SQUARE MATRIX
       edges : matrix that represents the weights of the edges between the nodes SQUARE MATRIX AND SAME SIZE AS VERTICES
       start: an integer, in which node to start the path
    Returns:
       previous : list the node that goes before that node in the shortes path
       distance : the distance of the shortest path up to that node
    '''
    
    N = len( vertices[0] )  # This is the number of vertices we are dealing with
    
    # Initialization
    
    Q = [0] * N # if the node has or hasn't been visited, this is the priority queue, 0 if the node hasn't 
    # been visited (thus its still in the queue) and 1 if it has been visited (thus not in the queue)
    Q[start] = 1 # But the starting node is not in the priority queue
    distance = [infinity] * N # All the distances are initialized at infinity
    distance[start] = 0 # Except the distance from the start node to the start node
    previous = [start]* N # list of the node that goes before, they all start in the starting node
    
    # Initialize the distance of the nodes that are connected to the starting node
    for j in range(N):
        if vertices[start][j] == 1:
            distance[j] = edges[start][j]
    
    for vertex in range(N):
        # Here we iterate for each vertex in our graph
        minDistance = infinity # We start assuming that all the distances are infinite, we need to find the minimum
        # This is because we don't want to use the find min function in Python
        # We need to extract U, the node in Q that has the minimum distance (current distance) and hasnt been visited
        # yet (i.e.) is still in the priority queue/still has a temporary tag
        
        # FIND THE NEXT OPTIMAL NODE IN THE PATH
        for u in range(N):
            if distance[u] < minDistance and Q[u] == 0:
                minDistance = distance[u] # we found a better distance
                next_node_in_path = u
        # We update, we found the next node to visit
        previous[next_node_in_path] = vertex # thus the next node to visit after the current node is the optimal one
        Q[next_node_in_path] = 1 # we've visited this optimal node
        
        # WE UPDATE ALL THE OTHER CURRENT WEIGHTS OF THE NODES IN THE PRIORITY QUEUE
        # But we must update just for the nodes that are actually the neighbours of next_node_in_path
        
        for n in range(N):
            if Q[n] == 0 and ( minDistance + edges[next_node_in_path][n] < distance[n]):
                # If this happens then its is its worth it to visit next_in_path and then n
                distance[n] = minDistance + edges[next_node_in_path][n]
                previous[n] = next_node_in_path
    
    return( distance, previous, Q )
        
                
                
        
    

