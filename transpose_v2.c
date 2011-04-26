#include<stdio.h>
#include<stdlib.h>
#include<array_alloc.h>
#include<math.h>
#include<string.h>
#include<timer.h>
#include<mpi.h>

#define NAN (0.0/0.0)

struct global_index
{
  int x;
  int y;
};

struct local_index
{
	int procX;
	int procY;
	int arrayX;
	int arrayY;
};

struct global_index local_to_global(int procX, int procY, int arrayX, int arrayY, int cols, int rows)
{
	struct global_index result;
	/* printf("Inside GI: pX: %d, pY: %d, aX: %d, aY: %d, cols: %d, rows: %d\n", procX, procY, arrayX, arrayY, cols, rows); */
	
	result.x = arrayX + (procX * cols);
	result.y = arrayY + (procY * rows);
	
	return result;
}

struct local_index global_to_local(int x, int y, int cols, int rows)
{
	struct local_index result;
	
	result.arrayX = x % cols;
	result.procX = (int) x / (int) cols;
	
	result.arrayY = y % rows;
	result.procY = (int) y / (int) rows;
	
	return result;
}

/* Factorisation function taken from Lab Sheet 9 - written (I assume) by Ian Bush */
/* (Reformatted to fit with my coding style) */
void factorise(int n, int *x, int*y)
{
    /* Factorise n into two numbers which are as close to equal as possible */
    *x = sqrt( n ) + 0.5;
    *y = n / *x;

    /* Loop will definitely terminate when x == 1 */
    while ( *x * *y != n )
    {
		(*x)--;
		*y = n / *x;
    }
}

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
	
	int coords2[2];
	int rank2;

	int global_x1, global_x2;
	int global_y1, global_y2;

	struct global_index gi;
	struct local_index li;
	
	MPI_Request send_request;
	MPI_Request recv_request;
	
	/* Stores the number of requests this core will need to keep track of
	and creates a pointer to eventually hold the array of these requests */
	int n_requests;
	MPI_Request *all_requests;
	
	int request_counter;
	
	int tag;
	
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
	
	/* ----------- MPI Setup ----------- */

	/* Initialise MPI */
	MPI_Init(&argc, &argv);

	/* Get the rank for this process, and the number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	/* Work out how big the grid should be, given the number of processes */
	factorise( nprocs, &npx, &npy );
	
	/* Create the cartesian communicator */
	dim_size[0] = npy;
	dim_size[1] = npx;
	periods[0] = 1;
	periods[1] = 1;
	
	MPI_Cart_create(MPI_COMM_WORLD, 2, dim_size, periods, 1, &cart_comm);
	
	/* Get our co-ordinates within that communicator */
	MPI_Cart_coords(cart_comm, rank, 2, coords);
		
	rows_in_core = ceil(n / (float) npx);
	cols_in_core = ceil(n / (float) npy);
	
	std_rows_in_core = rows_in_core;
	std_cols_in_core = cols_in_core;
		
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
	
	/* Calculate the number of individual communications that are needed to transpose the matrix.
	We know that this will be the number of requests we need to store, so malloc some memory for the
	array of the right size. Each process has its own requests array - so it's simply 2 * ncols * nrows
	*/
	n_requests = 2 * rows_in_core * cols_in_core;
	all_requests = malloc(n_requests * sizeof(MPI_Request));
	
	/* ----------- INITIALISE MATRIX CHUNKS ----------- */
	
	/* Allocate an array to store this chunk of the matrix.
	If the allocation fails, print an error */
	array = alloc_2d_double(rows_in_core, cols_in_core);
	new_array = alloc_2d_double(rows_in_core, cols_in_core);
	if (array == NULL || new_array == NULL)
	{
		printf("Problem with array allocation.\nExiting\n");
		return 1;
	}
	
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
	
	for (i = 0; i < rows_in_core; i++)
	{
		/*printf("Process (%d, %d): ", coords[0], coords[1]);*/
		for (j = 0; j < cols_in_core; j++)
		{
			array[i][j] = (float) ( (i + y_offset) * n) + ( (j + x_offset) + 1);
			
			/*printf("%f ", array[i][j]);*/
		}
		/*printf("\n");*/
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	/* ----------- DO TRANSPOSE ----------- */
	/* Record the start time, making sure all processes have got there at the same time using
	a barrier */
	MPI_Barrier(cart_comm);
	tstart = timer();
		
	request_counter = 0;
	
	for (i = 0; i < rows_in_core; i++)
	{
		for (j = 0; j < cols_in_core; j++)
		{
			/* For each item in this chunk of the array:
			
			Send it to the opposite index
			Receive from the opposite index */
						
			/* Calculate the opposite index */
			gi = local_to_global(coords[0], coords[1], j, i, std_cols_in_core, std_rows_in_core);
						
			/* The two global indices for swapping */
			global_x1 = gi.x;
			global_y1 = gi.y;
			
			global_x2 = global_y1;
			global_y2 = global_x1;
			

			/* We know which rank one of them is (as we're running as it!)
			but we need to calculate the rank of the other which involves:
			1. Getting the local index of the swapped global index
			2. Converting the X, Y process co-ordinates to a rank */
			li = global_to_local(global_x2, global_y2, std_cols_in_core, std_rows_in_core);
						
			coords2[0] = li.procX;
			coords2[1] = li.procY;
						
			MPI_Cart_rank(cart_comm, coords2, &rank2);
			
			tag = (global_x2 + 1) * (global_y2 + 1);
						
			/*printf("Swapping\n");
			printf("Value = %f\n", array[i][j]);
			printf("(%f) From global (%d, %d) - Processor (%d, %d) with local (%d, %d)\n", array[i][j], global_x1, global_y1, coords[0], coords[1], j, i);
			printf("(%f) To global (%d, %d)- Processor (%d, %d) with local (%d, %d)\n", array[i][j], global_x2, global_y2, coords2[0], coords2[1], li.arrayX, li.arrayY);
			printf("(%f) Going to rank %d\n", array[i][j], rank2); */
			
			MPI_Isend(&array[i][j], 1, MPI_DOUBLE, rank2, (global_y2 * n) + global_x2, MPI_COMM_WORLD, &send_request);

			MPI_Irecv(&new_array[i][j], 1, MPI_DOUBLE, rank2, (global_x2 * n) + global_y2, MPI_COMM_WORLD, &recv_request);
			all_requests[request_counter] = send_request;
			request_counter++;
			all_requests[request_counter] = recv_request;
			request_counter++;
		}
	}
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
	
	/* Free the memory we grabbed for the requests array */
	free(all_requests);
	
	/* Free the memory we grabbed earlier for the data arrays */
	free_2d_double(array);
	free_2d_double(new_array);
	
	/* Close down the MPI environment */
	MPI_Finalize();
	
	return EXIT_SUCCESS;
}