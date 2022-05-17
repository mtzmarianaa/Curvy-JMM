#include "fmm_2d.c"

int main(){
	int M, N, n;
    double x_min, y_min, h;
    int start[2];

    n = pow(2, 8);
    M = 2*n+1;
    N = 2*n+1;
    h = 1.0/n;

    x_min = -1.0;
    y_min = -1.0;
    start[0] = n;
    start[1] = n;

	double *distance = malloc(M*N*sizeof(double));
    double *trueSolution = malloc(M*N*sizeof(double));
    double *eik_queue = malloc(M*N*sizeof(double));
    int *index_queue = malloc(M*N*sizeof(int));
    int *current_states = malloc(M*N*sizeof(int));
    int *Q = malloc(M*N*sizeof(int));
	if (distance == NULL) { // we just need to worry about either distance or trueSolution to be null because up until rn they're the same
		printf("oh no!\n");
		exit(EXIT_FAILURE);
	}
    if (Q == NULL){
        printf("oh no! \n");
        exit(EXIT_FAILURE);
    }

    // Use naive FMM to solve
    FMM_2D( x_min, y_min, start, distance, eik_queue, index_queue, current_states, M, N, h);

    // Compute actual solution
    ActualSolution(x_min, y_min, start, trueSolution, M, N, h);

    // Generate data from both sources so that we can plot them

    char data_name1[] = "from_FMM_h001k8_bTree";

    generateDataForPlot(distance, M, N, data_name1);

    char data_name2[] = "exact_h001k8_bTree";

    generateDataForPlot(trueSolution, M, N, data_name2);

    //printQGridFromQueue(Q, M, N);
    //printGridFromDistance(distance, M, N);

	free(distance);
    free(trueSolution);
    free(Q);

    return EXIT_SUCCESS;
}