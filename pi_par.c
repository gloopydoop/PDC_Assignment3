
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>


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


    	if (rank == 0)
	{
    		for (i=1; i < num_ranks; i++)
    		{
	        int extra_count;
		MPI_Recv(&extra_count,1,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		local_count += extra_count;
		}
	    
		pi = ((double) local_count / (double) (flip * num_ranks)) * 4.0;
	}
	else
	{
	    	MPI_Send(&local_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
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

