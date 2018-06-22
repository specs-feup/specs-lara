#include <stdio.h>
#include <iostream>
#include "mpi.h"
#define N 1011

int numTasks;
int numWorkers;
int taskId;

void clavaMpiWorker() {
	
    MPI_Status status;
    int mpi_loop_num_elems;

    MPI_Recv(&mpi_loop_num_elems, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

    // Inputs
    double a[mpi_loop_num_elems];
    MPI_Recv(&a, mpi_loop_num_elems, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);

    double b[mpi_loop_num_elems];
    MPI_Recv(&b, mpi_loop_num_elems, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);

    // Declare outputs
    double c[mpi_loop_num_elems];

    for(int i=0; i<mpi_loop_num_elems; i++) {
        c[i] = a[i] + b[i];
    }

    // Send outputs
    MPI_Send(&c, mpi_loop_num_elems, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
}

void foo(double a[N], double b[N], double c[N]) {
	
	{
		// Loop to parallelize
		#pragma clava parallel
		//if(taskId == 0) {
		// Master
		// split iterations of the loop
		int clava_mpi_loop_limit = N;
		// A better distribution calculation could be used
		int clava_mpi_num_iter = N / numWorkers;
		int clava_mpi_num_iter_last = clava_mpi_num_iter + N % numWorkers;
		// int clava_mpi_num_iter_last = clava_mpi_num_iter + (clava_mpi_loop_limit - (clava_mpi_num_iter * numWorkers));
	
		// send number of iterations
		for(int i=0; i<numWorkers-1; i++) {
			MPI_Send(&clava_mpi_num_iter, 1, MPI_INT, i+1, 1, MPI_COMM_WORLD);
		}
		MPI_Send(&clava_mpi_num_iter_last, 1, MPI_INT, numWorkers, 1, MPI_COMM_WORLD);
	
		// send input a - elements: iteration_space
		for(int i=0; i<numWorkers-1; i++) {
			MPI_Send(&a[i*clava_mpi_num_iter], clava_mpi_num_iter, MPI_DOUBLE, i + 1,1, MPI_COMM_WORLD);
		}
		MPI_Send(&a[(numWorkers-1)*clava_mpi_num_iter], clava_mpi_num_iter_last, MPI_DOUBLE, numWorkers,1, MPI_COMM_WORLD);
	
		// send input b - elements: iteration_space
		for(int i=0; i<numWorkers-1; i++) {
			MPI_Send(&b[i*clava_mpi_num_iter], clava_mpi_num_iter, MPI_DOUBLE, i + 1,1, MPI_COMM_WORLD);
		}
		MPI_Send(&b[(numWorkers-1)*clava_mpi_num_iter], clava_mpi_num_iter_last, MPI_DOUBLE, numWorkers,1, MPI_COMM_WORLD);
	
	
		
		MPI_Status status;
		
		// receive output c - elements: iteration_space
		for(int i=0; i<numWorkers-1; i++) {
			MPI_Recv(&c[i*clava_mpi_num_iter], clava_mpi_num_iter, MPI_DOUBLE, i + 1, 1, MPI_COMM_WORLD, &status);
		}
		MPI_Recv(&c[(numWorkers-1)*clava_mpi_num_iter], clava_mpi_num_iter_last, MPI_DOUBLE, numWorkers, 1, MPI_COMM_WORLD, &status);		
	}
	
}

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskId);
    MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
    numWorkers = numTasks - 1;

    if(numWorkers == 0) {
        std::cerr << "This program does not support working with a single process." << std::endl;
        return 1;
    }

	if(taskId > 0) {
		clavaMpiWorker();
		MPI_Finalize();
		return 0;
	}

    double a[N], b[N], c[N];

    for(int i=0; i<N; i++) {
        a[i] = i;
        b[i] = i + 1;
    }

	foo(a, b, c);


    // test output
    double acc = 0;
    for(int i=0; i<N; i++) {
        acc += c[i];
    }

	printf("Result: %f\n", acc);

    MPI_Finalize();
	
    return 0;
    
}
