# Project with FMM (Fast Marching Method)

### First stage: getting familiarized with FMM
#### Spring Break 2022
Since I've never done something like this, the strategy to follow is to implement Dijkstra's algorithm and then try to implement the "vanilla" version of the FMM. So far:

   - Implemented Dijkstra's in Python (so that I have a template when I want to implement it in C), outputs are arrays as suggested
   - Implemented Dijkstra's in C with a lot of help and a lot of suggestions

### Second stage: using C + implementing a naive version of the FMM
   - Gave up using GDB for debugging in C + VSCode + Mac
   - "Standarized" the LDDB json files and tasks.json files, now I can pseudo debug and compile + execute in the same task
   - Started implementing FMM in C with help
   - Finished FMM initialization step
   - Finished FMM naive implementation (need to implement it with the correct data structures)
   - Corrected the FMM implementation
   - Added a file to plot the "solutions" I get using Python
   - Implemented the priority queue as a binary tree with a "parallel" array that has the information about the indices

### Third stage: building everything in order to have a working method in a triangular mesh structure + optimiation pov
   - Built the test geometry in python with meshpy and export it
   - Built structs coords, neighbors, facets, faces
   - Built the mesh structure
   - Built the secant method in C for this particular use

### Fourth stage: include higher order approximations
   - Include a higher order approximation to T(x_lambda) and T(x_mu) using cubic Hermite interpolation
   - Change the optimization problems (solve with secant and with projected gradient descent as well)
   - Compare this with the linear approximation we used to have
