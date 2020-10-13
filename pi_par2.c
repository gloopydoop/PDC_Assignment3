
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define SEED     921

int main(int argc, char* argv[])
{
    int local_count = 0; 
    int flip = 1 << 24;
    int rank, num_ranks, i, iter, provided;
    double x, y, z, pi;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    double start_time, stop_time, elapsed_time;
    start_time = MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

    srand(rank*100); // Important: Multiply SEED by "rank" when you introduce MPI!
    
    flip = flip/num_ranks;
    // Calculate PI following a Monte Carlo method
    for (int iter = 0; iter < flip; iter++)
    {
        // Generate random (X,Y) points
        x = (double)random() / (double)RAND_MAX;
        y = (double)random() / (double)RAND_MAX;
        z = sqrt((x*x) + (y*y));
        
        // Check if point is in unit circle
        if (z <= 1.0)
        {
            local_count++;
        }
    }

	// here is the tree bit
	// obviously the algorithm was stolen from some stack exchange
	//
	//
	for (int tree_depth =  pow(2,ceil(log(num_ranks)/log(2))-1); tree_depth > 0; tree_depth >>= 1)
     {
         // printf("tree depth %d\n", tree_depth);
             if ((rank >= tree_depth) && (rank < (tree_depth << 1)))
             {
                // printf("%d -> %d\n", k, k - tree_depth);
		 MPI_Send(&local_count, 1, MPI_INT, (rank-tree_depth), 0, MPI_COMM_WORLD); 
             }
            
             
             else if ((rank < tree_depth) && (rank + tree_depth) < num_ranks)
             {
		 int extra_count;
    		 MPI_Recv(&extra_count, 1, MPI_INT,rank+tree_depth, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
		 local_count += extra_count;
             }
     
     } 

	if(rank==0)
	{
	pi = ((double) local_count / (double) (flip * num_ranks)) * 4.0;
	}
    stop_time = MPI_Wtime();
    elapsed_time = stop_time - start_time;
    //printf("rank %d: %d / %d = %f\n",rank, local_count, flip, (double)local_count / (double)flip);
    if (rank == 0) {
	    printf("%d \t %f \t %f\n", num_ranks,elapsed_time,pi);
    }
    MPI_Finalize();
    return 0;
}

