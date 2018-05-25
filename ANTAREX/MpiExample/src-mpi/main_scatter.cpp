#include <vector>
#include <iostream>
#include "mpi.h"

const int N = 1011;

int clava_mpi_num_tasks;
int clava_mpi_num_workers;
int clava_mpi_world_rank;

void clava_mpi_calc_scatter(int loop_limit, int num_workers, int* send_counts, int* displs, int* buffer_size)
{
    int sum = 0;
    int quo = loop_limit / num_workers;
    int rem = loop_limit % num_workers;
    
    *buffer_size = rem == 0 ? quo : quo + 1;

    // calculate send counts and displacements
    for (int i = 0; i < num_workers; i++) {

        send_counts[i] = quo;
        if (rem > 0) {
            send_counts[i]++;
            rem--;
        }

        displs[i] = sum;
        sum += send_counts[i];
    }
    
    // print calculated send counts and displacements for each process
    if (0 == clava_mpi_world_rank) {
        std::cout << "Buffer size: " << *buffer_size << std::endl;
        for (int i = 0; i < num_workers; i++) {
            printf("sendcounts[%d] = %d\tdispls[%d] = %d\n", i, send_counts[i], i, displs[i]);
        }
    }
}

int main(int argc, char **argv)
{

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &clava_mpi_world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &clava_mpi_num_tasks);
    clava_mpi_num_workers = clava_mpi_num_tasks; // in this version, the master also performs work

    if(clava_mpi_num_workers == 0) {
        std::cerr << "This program does not support working with a single process." << std::endl;
        return 1;
    }

    int a[N], b[N], c[N];

    for(int i=0; i<N; i++) {
        a[i] = i;
        b[i] = i + 1;
    }

    // Loop to parallelize
#pragma clava parallel

    // split iterations of the loop
    int clava_mpi_loop_limit = N;
    int* clava_mpi_send_counts = (int*) malloc(sizeof(*clava_mpi_send_counts) * clava_mpi_num_workers);
    int* clava_mpi_displs = (int*) malloc(sizeof(*clava_mpi_displs) * clava_mpi_num_workers);
    int clava_mpi_buffer_size;

    clava_mpi_calc_scatter(clava_mpi_loop_limit, clava_mpi_num_workers,
                           clava_mpi_send_counts, clava_mpi_displs, &clava_mpi_buffer_size);

    // Buffers for each subprocess
    int *clava_mpi_sub_a = (int *) malloc(sizeof(*clava_mpi_sub_a) * clava_mpi_buffer_size);
    int *clava_mpi_sub_b = (int *) malloc(sizeof(*clava_mpi_sub_b) * clava_mpi_buffer_size);
    int *clava_mpi_sub_c = (int *) malloc(sizeof(*clava_mpi_sub_c) * clava_mpi_buffer_size);

    // Scatter a
    MPI_Scatterv(a,
                 clava_mpi_send_counts,
                 clava_mpi_displs,
                 MPI_INT,
                 clava_mpi_sub_a,
                 clava_mpi_buffer_size,
                 MPI_INT,
                 0,
                 MPI_COMM_WORLD);

    // Scatter b
    MPI_Scatterv(b,
                 clava_mpi_send_counts,
                 clava_mpi_displs,
                 MPI_INT,
                 clava_mpi_sub_b,
                 clava_mpi_buffer_size,
                 MPI_INT,
                 0,
                 MPI_COMM_WORLD);

    // Compute c
    for(int i=0; i < clava_mpi_send_counts[clava_mpi_world_rank]; i++) {
        clava_mpi_sub_c[i] = clava_mpi_sub_a[i] + clava_mpi_sub_b[i];
    }

    // Gather c
    MPI_Gatherv(clava_mpi_sub_c,
                clava_mpi_send_counts[clava_mpi_world_rank],
                MPI_INT,
                c,
                clava_mpi_send_counts,
                clava_mpi_displs,
                MPI_INT,
                0,
                MPI_COMM_WORLD);

    // test output
    int acc = 0;
    for(int i=0; i<N; i++) {
        acc += c[i];
    }

    printf("[%d] Result: %d\n", clava_mpi_world_rank, acc);

    MPI_Finalize();
    
    free(clava_mpi_sub_c);
    free(clava_mpi_sub_b);
    free(clava_mpi_sub_a);
    free(clava_mpi_send_counts);
    free(clava_mpi_displs);
    
    return 0;

}
