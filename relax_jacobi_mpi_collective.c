/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil -*- */
/*
 * Copyright (c) 2016-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @Author: George Bosilca
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER$
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <assert.h>

#include "jacobi.h"
#include "header.h"

/**
 * This is a jacobi_mpi_collective 
 */


struct relaxation_jacobi_mpi_collective_hidden_params_s {
    struct relaxation_params_s super;
    double* data[2];
    uint32_t idx;
    // additional params
    uint32_t res;
    uint32_t P, Q;
    int err_msg, rank, size;
    MPI_Datatype north_south_type, east_west_type;
    // collective params
    MPI_Comm topocomm;
    int north, south, west, east;
};

void compute_inner_jacobi_collective(double* n, double* o, struct relaxation_jacobi_mpi_collective_hidden_params_s* rp){
    /* The size[x,y] account for the boundary */
    int i, j;
    for(i = 2; i < (rp->super.sizey-2); i++ ) {
        for(j = 2; j < (rp->super.sizex-2); j++ ) {
            n[i*rp->super.sizex+j] = 0.25 * (o[ i     * rp->super.sizex + (j-1) ]+  // left
                                       o[ i     * rp->super.sizex + (j+1) ]+  // right
                                       o[ (i-1) * rp->super.sizex + j     ]+  // top
                                       o[ (i+1) * rp->super.sizex + j     ]); // bottom
        }
    }
}

void compute_outer_jacobi_collective(double* n, double* o, struct relaxation_jacobi_mpi_collective_hidden_params_s* rp, int direction){
    int i, j;
    /* The size[x,y] account for the boundary */
    if (direction == 0){
        // printf("Iteration %d, Rank %d: top row (receive from north)\n", rp->idx, rp->rank);
        i = 1;
        for (j=1; j<rp->super.sizex-1;j++){
            n[i*rp->super.sizex+j] = 0.25 * (o[ i     * rp->super.sizex + (j-1) ]+  // left
                                           o[ i     * rp->super.sizex + (j+1) ]+  // right
                                           o[ (i-1) * rp->super.sizex + j     ]+  // top
                                           o[ (i+1) * rp->super.sizex + j     ]); // bottom
        }

    } else if (direction == 1){
        // printf("Iteration %d, Rank %d: bottom row (receive from south)\n", rp->idx, rp->rank);
        i = rp->super.sizey-2;
        for (j=1; j<rp->super.sizex-1;j++){
            n[i*rp->super.sizex+j] = 0.25 * (o[ i     * rp->super.sizex + (j-1) ]+  // left
                                           o[ i     * rp->super.sizex + (j+1) ]+  // right
                                           o[ (i-1) * rp->super.sizex + j     ]+  // top
                                           o[ (i+1) * rp->super.sizex + j     ]); // bottom
        }
    } else if (direction == 2){
        // printf("Iteration %d, Rank %d: right column (receive from east)\n", rp->idx, rp->rank);
        j = rp->super.sizex -2;
        for (i=1; i<rp->super.sizey-1;i++){
            n[i*rp->super.sizex+j] = 0.25 * (o[ i     * rp->super.sizex + (j-1) ]+  // left
                                           o[ i     * rp->super.sizex + (j+1) ]+  // right
                                           o[ (i-1) * rp->super.sizex + j     ]+  // top
                                           o[ (i+1) * rp->super.sizex + j     ]); // bottom
        }
    } else if (direction == 3){
        // printf("Iteration %d, Rank %d: left column (receive from west)\n", rp->idx, rp->rank);
        j = 1;
        for (i=1; i<rp->super.sizey-1;i++){
            n[i*rp->super.sizex+j] = 0.25 * (o[ i     * rp->super.sizex + (j-1) ]+  // left
                                           o[ i     * rp->super.sizex + (j+1) ]+  // right
                                           o[ (i-1) * rp->super.sizex + j     ]+  // top
                                           o[ (i+1) * rp->super.sizex + j     ]); // bottom
        }
    } else{
        printf("Incorrect direction specified\n");;
    }
    // printf("Iteration %d, Rank %d: jacobi completed \n", rp->idx, rp->rank);
}

const struct relaxation_function_class_s _relaxation_jacobi_mpi_collective;

static struct relaxation_params_s*
jacobi_mpi_collective_relaxation_init(struct hw_params_s* hw_params){

    struct relaxation_jacobi_mpi_collective_hidden_params_s* rp;
    uint32_t np = hw_params->resolution;

    rp = (struct relaxation_jacobi_mpi_collective_hidden_params_s*)
            malloc(sizeof(struct relaxation_jacobi_mpi_collective_hidden_params_s));

    if( NULL == rp ) {
        fprintf(stderr, "Cannot allocate memory for the relaxation structure\n");
        return NULL;
    }

    // initialize MPI 
    rp->err_msg = MPI_Init(NULL, NULL);                        // initialize MPI
    rp->err_msg = MPI_Comm_rank(MPI_COMM_WORLD, &rp->rank);         // Get the individual process ID (rank)
    rp->err_msg = MPI_Comm_size(MPI_COMM_WORLD, &rp->size);         // Get the number of processes (world)

    // Set up the processor grid P x Q
    rp->P = hw_params->num_proc_p;
    rp->Q = rp->size/rp->P; 
    assert(rp->size = (rp->P*rp->Q));

    // Set the size of the matrix
    rp->res = np;
    rp->super.sizex = np/rp->Q + 2;
    rp->super.sizey = np/rp->P + 2;
    rp->super.rel_class = &_relaxation_jacobi_mpi_collective;

    // Initialize iteration index and data buffers
    rp->idx = 0;
    rp->data[0] = calloc( (rp->super.sizey)*(rp->super.sizex), sizeof(double) );
    rp->data[1] = (double*) malloc((rp->super.sizey)*(rp->super.sizex)*sizeof(double));

    if( (NULL == rp->data[0]) && (NULL == rp->data[1]) ) {
        fprintf(stderr, "Cannot allocate the memory for the matrices\n");
        goto fail_and_return;
    }

    // Create a cartesian topology 
    int ndims = 2, reorder=0;
    int dims[2] = {rp->P,rp->Q}, periods[2] = {0,0}, coords[2];
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &rp->topocomm);
    
    // Obtain my local coordinate in process grid 
    MPI_Cart_coords(rp->topocomm, rp->rank, 2, coords);
    // printf("coord = (%d, %d)\n", coords[0], coords[1] );
    // int rx = coords[0], ry = coords[1];

    // Determine my neighbors
    MPI_Cart_shift(rp->topocomm, 0, 1, &rp->north, &rp->south);
    MPI_Cart_shift(rp->topocomm, 1, 1, &rp->west, &rp->east);
    // printf("Rank %d: north %d, south %d, east %d, west %d\n", rp->rank, rp->north, rp->south, rp->east, rp->west);

    // Determine the coordinates of current process in processor grid P x Q 
    int rx = rp->rank % rp->Q, ry = rp->rank / rp->Q; 
    int lm = rp->super.sizey, ln = rp->super.sizex;
    int colpos = rx*(np/rp->Q), rowpos = ry*(np/rp->P);

    // Split matrix into submatrices for all processes
    relaxation_matrix_set_ex(hw_params, rp->data[0], ln, lm, colpos, rowpos);
    // printf("Iteration %d, Rank %d:\n", rp->idx, rp->rank);
    // print_matrix(stdout, rp->data[0], rp->super.sizex, rp->super.sizey);
    
    /* Copy the boundary conditions on all matrices */
    memcpy(rp->data[1], rp->data[0], (rp->super.sizey)*(rp->super.sizex)*sizeof(double));

    // Create derived MPI_Datatype for north_south and east_west 
    MPI_Type_contiguous(rp->super.sizex-2, MPI_DOUBLE, &rp->north_south_type);
    MPI_Type_commit(&rp->north_south_type);
    MPI_Type_vector(rp->super.sizey-2,1,rp->super.sizex,MPI_DOUBLE, &rp->east_west_type);
    MPI_Type_commit(&rp->east_west_type);


    return (struct relaxation_params_s*)rp;
 fail_and_return:
    if( NULL != rp->data[0] ) free(rp->data[0]);
    if( NULL != rp->data[1] ) free(rp->data[1]);
    free(rp);
    return NULL;
}

static int jacobi_mpi_collective_relaxation_fini(relaxation_params_t** prp){

    struct relaxation_jacobi_mpi_collective_hidden_params_s* rp = (struct relaxation_jacobi_mpi_collective_hidden_params_s*)*prp;
    
    MPI_Type_free(&rp->north_south_type);
    MPI_Type_free(&rp->east_west_type);
    // collective params

    rp->err_msg = MPI_Finalize();

    if( NULL != rp->data[0] ) free(rp->data[0]);
    if( NULL != rp->data[1] ) free(rp->data[1]);
    free(rp);
    *prp = NULL;

    return 0;
}

static int jacobi_mpi_collective_relaxation_coarsen(relaxation_params_t* grp, double* dst, uint32_t dstx, uint32_t dsty){
    struct relaxation_jacobi_mpi_collective_hidden_params_s* rp = (struct relaxation_jacobi_mpi_collective_hidden_params_s*)grp;

    return coarsen(rp->data[(rp->idx + 1) % 2], rp->super.sizex, rp->super.sizey,
                   dst, dstx, dsty);
}

static int jacobi_mpi_collective_relaxation_print(relaxation_params_t* grp, FILE* fp){
    struct relaxation_jacobi_mpi_collective_hidden_params_s* rp = (struct relaxation_jacobi_mpi_collective_hidden_params_s*)grp;
    fprintf( fp, "\n# Iteration %d\n", rp->idx);
    print_matrix(fp, rp->data[(rp->idx + 1) % 2], rp->super.sizex, rp->super.sizey);
    fprintf( fp, "\n\n");
    return 0;
}

/**
 * One step of a simple jacobi_mpi_collective relaxation.
 */
static double jacobi_mpi_collective_relaxation_apply(relaxation_params_t* grp){
    struct relaxation_jacobi_mpi_collective_hidden_params_s* rp = (struct relaxation_jacobi_mpi_collective_hidden_params_s*)grp;
    double diff, sum = 0.0, *n, *o;
    int i, j;

    n = rp->data[(rp->idx + 0) % 2];
    o = rp->data[(rp->idx + 1) % 2];

    // printf("Iteration %d: Rank %d \n", rp->idx, rp->rank);
    // print_matrix(stdout, o, rp->super.sizex, rp->super.sizey);

    // Determine the size of local matrix for processor
    int bx = rp->super.sizex;
    int by = rp->super.sizey; 
        
    // allocate memory for send/recv buffer and copy data 
    // buffer ordering ( north, south, west, east)
    double* sendbuf = (double*) calloc(2*(bx-2)+2*(by-2), sizeof(double));
    double* recvbuf = (double*) calloc(2*(bx-2)+2*(by-2), sizeof(double));

    // printf("Iteration %d, Rank %d sendbuf:\n", rp->idx, rp->rank);
    // printf("send north\n");
    for (i=0; i<(bx-2); i++) {
        sendbuf[i] = o[1*bx+1+i];
    }               // north
    for (i=0; i<(bx-2); i++) {
        sendbuf[(bx-2)+i] = o[(by-2)*bx+1+i];
    }   //south
    for (i=0; i<(by-2); i++) {
        sendbuf[2*(bx-2)+i] = o[1*bx+1+(i*bx)];
    }                   // west
    for (i=0; i<(by-2); i++) {
        sendbuf[2*(bx-2)+(by-2)+i] = o[1*bx+(bx-2)+i*bx];
    }         // east
    
    int counts[4] = {bx-2, bx-2, by-2, by-2};
    int displs[4] = {0, bx-2, 2*(bx-2), 2*(bx-2)+(by-2)};

    MPI_Request req;
    MPI_Status status;
    MPI_Ineighbor_alltoallv(sendbuf, counts, displs, MPI_DOUBLE, 
        recvbuf, counts, displs, MPI_DOUBLE, rp->topocomm, &req);
    // can do a compute_inner_jacobi here for overlap computation and communication
    compute_inner_jacobi_collective(n,o,rp);
    MPI_Wait(&req, &status);
    // MPI_Neighbor_alltoallv(sendbuf, counts, displs, MPI_DOUBLE, 
    //     recvbuf, counts, displs, MPI_DOUBLE, rp->topocomm);

    if (rp->north >= 0){ // north
        for (i=0; i<(bx-2); i++) {o[0*bx+1+i] = recvbuf[i];} 
    }
         
    if (rp->south >= 0){ // south
        for (i=0; i<(bx-2); i++) {o[(by-1)*bx+1+i] = recvbuf[(bx-2)+i];}  
    }
    
    if (rp->west >=0){ //west
        for (i=0; i<(by-2); i++) {o[1*bx+0+i*bx] = recvbuf[2*(bx-2)+i];}  
    }    
     
    if (rp->east >= 0) { //east
        for (i=0; i<(by-2); i++) {o[1*bx+(bx-1)+i*bx] = recvbuf[2*(bx-2)+(by-2)+i];} 
    }                

    compute_outer_jacobi_collective(n, o, rp, 0);
    compute_outer_jacobi_collective(n, o, rp, 1);
    compute_outer_jacobi_collective(n, o, rp, 2);
    compute_outer_jacobi_collective(n, o, rp, 3);

    // Calculate the residual of local matrix after all communication and computation are done
    for( i = 1; i < (rp->super.sizey-1); i++ ) {
        for( j = 1; j < (rp->super.sizex-1); j++ ) {
            diff = n[i*rp->super.sizex+j] - o[i*rp->super.sizex+j];
            sum += diff * diff;
        }
    }
    

    rp->idx++;

    double global_sum = 0.0;
    MPI_Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return global_sum;

}

static double* jacobi_mpi_collective_relaxation_get_data(relaxation_params_t* grp){   
    struct relaxation_jacobi_mpi_collective_hidden_params_s* rp = (struct relaxation_jacobi_mpi_collective_hidden_params_s*)grp;
    
    return rp->data[(rp->idx + 1) % 2];
}


const struct relaxation_function_class_s _relaxation_jacobi_mpi_collective =
    {
        .type = RELAXATION_JACOBI_MPI_COLLECTIVE,
        .method_name = "MPI Collective for Jacobi relaxation",
		._init = jacobi_mpi_collective_relaxation_init,
        ._fini = jacobi_mpi_collective_relaxation_fini,
        ._coarsen = jacobi_mpi_collective_relaxation_coarsen,
        ._print = jacobi_mpi_collective_relaxation_print,
        ._relaxation = jacobi_mpi_collective_relaxation_apply,
        ._get_data = jacobi_mpi_collective_relaxation_get_data,
    };
