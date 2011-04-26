/* For printf() */
#include<stdio.h>

/* For constants for EXIT_SUCCESS and EXIT_FAILURE codes (set right for each OS) */
#include<stdlib.h>

/* Ian Bush's routines to allocate arrays */
#include<array_alloc.h>

/* Ian Bush's timing routines */
#include<timer.h>

/* To get sqrt() */
#include<math.h>

/* To get strcmp() */
#include<string.h>

/* For all of the MPI functions and constants */
#include<mpi.h>

/* Not a Number is defined in the IEEE standard, and one guarenteed way to get it is to divide
0 by 0 */
#define NAN (0.0/0.0)

int main(int argc, char ** argv)
{
	/* ----------- VARIABLE DECLARATIONS ----------- */

	/* Size of matrix - the matrix will be n x n */
	int n;

	/* Number of processes and rank of this process*/
	int nprocs, rank;
	
	/* Size and periodicity of cartesian grid */
	int dim_size[2];
	int periods[2];
	
	/* Number of processes in x and y directions */
	int npx, npy;
	
	/* Offset for this process in x and y directions */
	int x_offset, y_offset;
	
	/* The standard (normal) number of rows and cols
	per core - not valid for the last in a row or col */
	int std_rows_in_core, std_cols_in_core;
	
	/* Coordinates in cartesian grid, and num of rows and cols
	this process has */
	int coords[2];
	int rows_in_core;
	int cols_in_core;
	
	/* Array to store the chunk of the matrix that we have */
	double **array;
	double **new_array;
	
	/* Loop variables */
	int i, j;
	
	/* Cartesian communicator */
	MPI_Comm cart_comm;
	
	/* Hold the co-ordinates and rank of the process we're sending to/receiving from */
	int coords2[2];
	int rank2;
	
	/* Request objects returned from Isend and Irecv calls */
	MPI_Request send_request;
	MPI_Request recv_request;
	
	/* Stores the number of requests this core will need to keep track of
	and creates a pointer to eventually hold the array of these requests */
	int n_requests;
	MPI_Request *all_requests;
	int request_counter = 0;
	
	/* Hold the tag used for the current communication */
	int tag;

	/* Holds the derived type used to actually transpose the data */
	MPI_Datatype vector_type;
	
	/* Determines whether the final result should be printed or not */
	int print = 0;
	
	/* Holds start and end times for the transpose */
	double tstart, tend, walltime;
	
	/* ----------- PROCESS CMD LINE ARGUMENTS ----------- */
	
	/* If there isn't one command-line argument (argc also counts the program name)
	then print an error message and exit.
	Only do this if we're on process 0 - otherwise nprocs error messages will be printed!*/
	if (argc < 2 && rank == 0)
	{
		printf("You must specify the size of the matrix (n) as a command-line argument\n");
		printf("Exiting\n");
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		return EXIT_FAILURE;
	}
	
	/* As we've got here, there must be at least one command-line argument, so try and convert it to an
	integer */
	n = atoi(argv[1]);
	
	/* Process other command line arguments */
	/* Process command line arguments checking them and running the appropriate function */
	for (i = 2; i < argc; i++)  /* Skip argv[0] (program name) and argv[1] (size of matrix)*/
    {
    	/* if strcmp returns 0 then the strings are identical */
        if (strcmp(argv[i], "print") == 0)
        {
            print = 1;
        }

    }	
	
	/* If atoi returns zero it is either because it couldn't convert the string to an integer,
	or the integer given was 0. In either case, we must print an error message because it is a
	nonsensical size for the matrix */
	if (n == 0)
	{
		printf("You must an integer > 0 as the size of the matrix\n");
		printf("Exiting\n");
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		return EXIT_FAILURE;
	}
		
	/* ----------- MPI SETUP ----------- */

	/* Initialise MPI */
	MPI_Init(&argc, &argv);

	/* Get the rank for this process, and the number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	/* The method we'll use here will only work on a square grid of processors.
	Therefore, find the square number below the given number of processors, and
	only use that number of processors */
	
	npx = (int) sqrt(nprocs);
	npy = npx;

	/* Create the cartesian communicator */
	dim_size[0] = npy;
	dim_size[1] = npx;
	periods[0] = 1;
	periods[1] = 1;
	
	MPI_Cart_create(MPI_COMM_WORLD, 2, dim_size, periods, 1, &cart_comm);
	
	if (rank >= (npx * npx))
	{
		/* We're not using these processors because they don't factorise into an square
		so just tell these processors to stop doing anything and return */
		MPI_Finalize();
		return EXIT_SUCCESS;
	}
	else
	{	
		/* Get our co-ordinates within the cartesian communicator */
		MPI_Cart_coords(cart_comm, rank, 2, coords);
			
		/* Calculate how many rows and columns there will be in this core */
		rows_in_core = ceil(n / (float) npx);
		cols_in_core = ceil(n / (float) npy);
		
		std_cols_in_core = cols_in_core;
		std_rows_in_core = rows_in_core;
			
		if (coords[0] == (npy - 1))
		{
			/* We're at the far end of a row */
			cols_in_core = n - (cols_in_core * (npy - 1));
		}
		if (coords[1] == (npx - 1))
		{
			/* We're at the bottom of a col */
			rows_in_core = n - (rows_in_core * (npx - 1));
		}
			
		/* ----------- INITIALISE MATRIX CHUNKS ----------- */
		
		/* Allocate an array to store this chunk of the matrix.
		If the allocation fails, print an error */
		array = alloc_2d_double(rows_in_core, cols_in_core);
		new_array = alloc_2d_double(rows_in_core, cols_in_core);
		if (array == NULL || new_array == NULL)
		{
			printf("Problem with array allocation.\nExiting\n");
			MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
			return EXIT_FAILURE;
		}
		
		/* Initialise the output array to NaNs so that it's easy to tell if something goes wrong */
		for (i = 0; i < rows_in_core; i++)
		{
			for (j = 0; j < cols_in_core; j++)
			{
				new_array[i][j] = NAN;
			}
		}
		
		/* Calculate the offset of the top left-hand corner of our chunk of the matrix from the
		top left-hand corner of the whole matrix */
		x_offset = coords[0] * std_cols_in_core;
		y_offset = coords[1] * std_rows_in_core;
		
		/* Fill the matrix with values so that the matrix ends up going (for a 4x4 matrix):
		1  2  3  4
		5  6  7  8
		9  10 11 12
		13 14 15 16
		*/
		for (i = 0; i < rows_in_core; i++)
		{
			for (j = 0; j < cols_in_core; j++)
			{
				array[i][j] = (float) ( (i + y_offset) * n) + ( (j + x_offset) + 1);
			}
		}
		
		/* malloc some memory for the array to hold all of the requests from the sends and receives
		we've done */
		n_requests = rows_in_core + cols_in_core;
		all_requests = malloc(n_requests * sizeof(MPI_Request));
			
		/* ----------- DO TRANSPOSE ----------- */
		
		/* Record the start time, making sure all processes have got there at the same time using
		a barrier */
		MPI_Barrier(cart_comm);
		tstart = timer();
		
		/* Find the opposite co-ordinates (as we know it's a square) */
		coords2[0] = coords[1];
		coords2[1] = coords[0];
		
		/* Get the rank for the opposite co-ordinates */
		MPI_Cart_rank(cart_comm, coords2, &rank2);

		
		/* Create new derived type to receive as, which will actually do the transposing */
		MPI_Type_vector(rows_in_core, 1, cols_in_core, MPI_DOUBLE, &vector_type);
		MPI_Type_commit(&vector_type);

		
		/* Loop through the rows and send each row */
		for (i = 0; i < rows_in_core; i++)
		{
			tag = (coords[0] + 1) * (coords[1] + 1) + (i + 1);
			MPI_Isend(&array[i][0], cols_in_core, MPI_DOUBLE, rank2, tag, cart_comm, &send_request);
			all_requests[request_counter] = send_request;
			request_counter++;
		}
		
		/* Loop through the columns and receive into a column (using the derived datatype to convert
		from the row that has been sent) */
		for (i = 0; i < cols_in_core; i++)
		{
			tag = (coords[0] + 1) * (coords[1] + 1) + (i + 1);
			MPI_Irecv(&new_array[0][i], 1, vector_type, rank2, tag, cart_comm, &recv_request);
			all_requests[request_counter] = recv_request;
			request_counter++;
		}
		
		/* Wait for all send/receive requests to complete */
		MPI_Waitall(request_counter, all_requests, MPI_STATUSES_IGNORE);
		
		/* Store end time (using barrier before-hand to make the calculated times for all of the
		processes almost exactly equal, meaning we can just take the time from one process
		for output */
		MPI_Barrier(cart_comm);
		tend = timer();
			
		/* ----------- PRINT OUTPUT AND FINALISE ENVIRONMENT ----------- */
		
		walltime = tend - tstart;
		
		if (print == 1)
		{
			for (i = 0; i < rows_in_core; i++)
			{
				printf("Result at Process (%d, %d): ", coords[0], coords[1]);
				for (j = 0; j < cols_in_core; j++)
				{			
					printf("%f ", new_array[i][j]);
				}
				printf("\n");
			}
		}
		
		if (rank == 0)
		{
			printf("%d\t%d\t%f\n", nprocs, n, walltime);
		}
		
		free(all_requests);
		
		/* Free the memory we grabbed earlier for the data arrays */
		free_2d_double(array);
		free_2d_double(new_array);
		
		
		/* Close down the MPI environment */
		MPI_Finalize();
		
		return EXIT_SUCCESS;
	}
}