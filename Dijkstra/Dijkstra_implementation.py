#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Implementation of Dijkstra's algorithm

"""

# To look for an explanation read the pdf that is in the Notizen folder.
# Dijkstra's algorithm is used to find the shortes path in a weighted graph
# given an initial starting point.
# We're going to assume that the graph is given in two matrices (this is because
# we're tying to avoid using objects). The first graph represents if two nodes
# are connected or not. (This in the FMM is going to be much easier hopefully). 
# So in this first graph 1 in the ij entry means that the two nodes are connected
# and 0 represents that the nodes i and j are not connected.
# The second graph represents the weights of the edges between the nodes.

