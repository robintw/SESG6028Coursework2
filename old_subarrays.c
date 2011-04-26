		/* Code to use subarrays to do it
		
		tag = (coords[0] + 1) * (coords[1] + 1);
		
		sizes[0] = rows_in_core;
		sizes[1] = cols_in_core;
		
		subsizes[0] = rows_in_core;
		subsizes[1] = cols_in_core;
		
		starts[0] = 0;
		starts[1] = 0;
		
		MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE, &send_type);
		MPI_Type_commit(&send_type);
		
		MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &recv_type);
		MPI_Type_commit(&recv_type);
		
		MPI_Send(&array[0][0], 1, send_type, rank2, tag, cart_comm);
		
		MPI_Recv(&new_array[0][0], 1, recv_type, rank2, tag, cart_comm, &status);
		
		*/