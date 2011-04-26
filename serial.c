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


/* Not a Number is defined in the IEEE standard, and one guarenteed way to get it is to divide
0 by 0 */
#define NAN (0.0/0.0)

int main(int argc, char ** argv)
{
	/* ----------- VARIABLE DECLARATIONS ----------- */

	/* Size of matrix - the matrix will be n x n */
	int n;

	/* Array to store the chunk of the matrix that we have */
	double **array;
	double **new_array;
	
	/* Loop variables */
	int i, j;
	
	/* Determines whether the final result should be printed or not */
	int print = 0;
	
	/* Holds start and end times for the transpose */
	double tstart, tend, walltime;
	
	/* ----------- PROCESS CMD LINE ARGUMENTS ----------- */
	
	/* If there isn't one command-line argument (argc also counts the program name)
	then print an error message and exit.
	Only do this if we're on process 0 - otherwise nprocs error messages will be printed!*/
	if (argc < 2)
	{
		printf("You must specify the size of the matrix (n) as a command-line argument\n");
		printf("Exiting\n");
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
		return EXIT_FAILURE;
	}
	
	array = alloc_2d_double(n, n);
	new_array = alloc_2d_double(n, n);

	
	for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				array[i][j] = (float) ( i * n) + ( j + 1);
			}
		}

	/* ----------- DO TRANSPOSE -------- */
	tstart = timer();
	
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			new_array[i][j] = array[j][i];
		}
	}
			
	tend = timer();
	
		/* ----------- PRINT OUTPUT AND FINALISE ENVIRONMENT ----------- */
		
		walltime = tend - tstart;
		
		if (print == 1)
		{
			for (i = 0; i < n; i++)
			{
				for (j = 0; j < n; j++)
				{			
					printf("%f ", new_array[i][j]);
				}
				printf("\n");
			}
		}

		printf("%d\t%d\t%f\n", 1, n, walltime);

		
		/* Free the memory we grabbed earlier for the data arrays */
		free_2d_double(array);
		free_2d_double(new_array);
		
		return EXIT_SUCCESS;
}
