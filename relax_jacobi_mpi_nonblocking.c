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
 * This is a jacobi_mpi_nonblocking 
 */

struct relaxation_jacobi_mpi_nonblocking_hidden_params_s {
    struct relaxation_params_s super;
    double* data[2];
    uint32_t idx;
    // additional params
    uint32_t res;
    uint32_t P, Q;
    int err_msg, rank, size;
    double diff, sum, global_sum;
    MPI_Datatype north_south_type, east_west_type;
};

void compute_inner_jacobi(double* n, double* o, struct relaxation_jacobi_mpi_nonblocking_hidden_params_s* rp){
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

void compute_outer_jacobi(double* n, double* o, struct relaxation_jacobi_mpi_nonblocking_hidden_params_s* rp, int direction){
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


const struct relaxation_function_class_s _relaxation_jacobi_mpi_nonblocking;

static struct relaxation_params_s*
jacobi_mpi_nonblocking_relaxation_init(struct hw_params_s* hw_params){

    struct relaxation_jacobi_mpi_nonblocking_hidden_params_s* rp;
    uint32_t np = hw_params->resolution;

    rp = (struct relaxation_jacobi_mpi_nonblocking_hidden_params_s*)malloc(sizeof(struct relaxation_jacobi_mpi_nonblocking_hidden_params_s));
    if( NULL == rp ) {
        fprintf(stderr, "Cannot allocate memory for the relaxation structure\n");
        return NULL;
    }

    // initialize MPI 
    rp->err_msg = MPI_Init(NULL, NULL);                             // initialize MPI
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
    rp->super.rel_class = &_relaxation_jacobi_mpi_nonblocking;

    // Initialize iteration index and data buffers
    rp->idx = 0;
    rp->data[0] = calloc( (rp->super.sizey) * (rp->super.sizex), sizeof(double) );
    rp->data[1] = (double*)malloc( (rp->super.sizey)*(rp->super.sizex)*sizeof(double) );

    if( (NULL == rp->data[0]) && (NULL == rp->data[1]) ) {
        fprintf(stderr, "Cannot allocate the memory for the matrices\n");
        goto fail_and_return;
    }

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
    // printf("Iteration %d, Rank %d:\n", rp->idx, rp->rank);
    // print_matrix(stdout, rp->data[1], rp->super.sizex, rp->super.sizey);

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

static int jacobi_mpi_nonblocking_relaxation_fini(relaxation_params_t** prp){

    struct relaxation_jacobi_mpi_nonblocking_hidden_params_s* rp = (struct relaxation_jacobi_mpi_nonblocking_hidden_params_s*)*prp;
    
    MPI_Type_free(&rp->north_south_type);
    MPI_Type_free(&rp->east_west_type);
    rp->err_msg = MPI_Finalize();

    if( NULL != rp->data[0] ) free(rp->data[0]);
    if( NULL != rp->data[1] ) free(rp->data[1]);
    free(rp);
    *prp = NULL;

    return 0;
}

static int jacobi_mpi_nonblocking_relaxation_coarsen(relaxation_params_t* grp, double* dst, uint32_t dstx, uint32_t dsty){
    struct relaxation_jacobi_mpi_nonblocking_hidden_params_s* rp = (struct relaxation_jacobi_mpi_nonblocking_hidden_params_s*)grp;

    return coarsen(rp->data[(rp->idx + 1) % 2], rp->super.sizex, rp->super.sizey,
                   dst, dstx, dsty);
}

static int jacobi_mpi_nonblocking_relaxation_print(relaxation_params_t* grp, FILE* fp){
    struct relaxation_jacobi_mpi_nonblocking_hidden_params_s* rp = (struct relaxation_jacobi_mpi_nonblocking_hidden_params_s*)grp;
    fprintf( fp, "\n# Iteration %d\n", rp->idx);
    print_matrix(fp, rp->data[(rp->idx + 1) % 2], rp->super.sizex, rp->super.sizey);
    fprintf( fp, "\n\n");
    return 0;
}

/**
 * One step of a simple jacobi_mpi_nonblocking relaxation.
 */
static double jacobi_mpi_nonblocking_relaxation_apply(relaxation_params_t* grp){
    struct relaxation_jacobi_mpi_nonblocking_hidden_params_s* rp = (struct relaxation_jacobi_mpi_nonblocking_hidden_params_s*)grp;

    double diff, sum = 0.0, global_sum = 0.0;
    double *n, *o;
    int i, j;

    n = rp->data[(rp->idx + 0) % 2];
    o = rp->data[(rp->idx + 1) % 2];

    // printf("Iteration %d: Rank %d \n", rp->idx, rp->rank);
    // print_matrix(stdout, o, rp->super.sizex, rp->super.sizey);

    // Determine the coordinates of processor in processor grid P x Q 
    int px = rp->Q, py = rp->P;
    int rx = rp->rank % px, ry = rp->rank / px; 
    // printf("Iteration %d, Rank %d: [%d,%d]\n", rp->idx, rp->rank, ry, rx);

    // Determine the size of local matrix for processor
    int bx = rp->super.sizex;
    int by = rp->super.sizey; 

    // Determine my four neighbors in the processor grid 
    int north = (ry-1)*px+rx; if(ry-1 < 0)   north = MPI_PROC_NULL;
    int south = (ry+1)*px+rx; if(ry+1 >= py) south = MPI_PROC_NULL;
    int west= ry*px+rx-1;     if(rx-1 < 0)   west = MPI_PROC_NULL;
    int east = ry*px+rx+1;    if(rx+1 >= px) east = MPI_PROC_NULL;

    // Communication exchange between neighbors (blocking) within submatrix
    typedef enum{RECV_N, RECV_S, RECV_E, RECV_W, SEND_N, SEND_S, SEND_E, SEND_W, MAX_REQ} MyRequests;
    MPI_Request reqs[MAX_REQ];
    MPI_Status stats[MAX_REQ];
    reqs[SEND_N] = MPI_REQUEST_NULL, reqs[SEND_S] = MPI_REQUEST_NULL, reqs[SEND_E] = MPI_REQUEST_NULL, reqs[SEND_W] = MPI_REQUEST_NULL;

    int tag = 0;

    // north (send to north and receive from north)
    if (north != MPI_PROC_NULL){
        compute_inner_jacobi(n, o, rp);
        MPI_Wait(&reqs[SEND_N], &stats[SEND_N]);
        // MPI_Sendrecv(&o[1*bx+1], 1, rp->north_south_type, north, tag, &o[0*bx+1], 1, rp->north_south_type, north, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);   
        MPI_Isend(&o[1*bx+1], 1, rp->north_south_type, north, tag, MPI_COMM_WORLD, &reqs[SEND_N]);
        // compute_inner_jacobi(n, o, rp);
        // MPI_Wait(&reqs[SEND_N], &stats[SEND_N]);
        // MPI_Recv(&o[0*bx+1], 1, rp->north_south_type, north, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
        MPI_Irecv(&o[0*bx+1], 1, rp->north_south_type, north, tag, MPI_COMM_WORLD, &reqs[RECV_N]);   
        MPI_Wait(&reqs[RECV_N], &stats[RECV_N]);
        // compute_outer_jacobi(n, o, rp, 0);
        // printf("Iteration %d, Rank %d: Receive from north\n", rp->idx, rp->rank);
        // print_matrix(stdout, o, rp->super.sizex, rp->super.sizey);
    }
    // south
    if (south != MPI_PROC_NULL){
        compute_inner_jacobi(n, o, rp);
        // MPI_Sendrecv(&o[(by-2)*bx+1], 1, rp->north_south_type, south, tag, &o[(by-1)*bx+1], 1, rp->north_south_type, south, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Wait(&reqs[SEND_S], &stats[SEND_S]);
        MPI_Isend(&o[(by-2)*bx+1], 1, rp->north_south_type, south, tag, MPI_COMM_WORLD, &reqs[SEND_S]);
        // compute_inner_jacobi(n, o, rp);
        // MPI_Wait(&reqs[SEND_S], &stats[SEND_S]);
        // MPI_Recv(&o[(by-1)*bx+1], 1, rp->north_south_type, south, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Irecv(&o[(by-1)*bx+1], 1, rp->north_south_type, south, tag, MPI_COMM_WORLD, &reqs[RECV_S]);
        MPI_Wait(&reqs[RECV_S], &stats[RECV_S]);
        // compute_outer_jacobi(n, o, rp, 1);
        // printf("Iteration %d, Rank %d: Receive from south\n", rp->idx, rp->rank);
        // print_matrix(stdout, o, rp->super.sizex, rp->super.sizey);
    }
    // east
    if (east != MPI_PROC_NULL){
        compute_inner_jacobi(n, o, rp);
        // MPI_Sendrecv(&o[1*bx+(bx-2)], 1, rp->east_west_type, east, tag, &o[1*bx+(bx-1)], 1, rp->east_west_type, east, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Wait(&reqs[SEND_E], &stats[SEND_E]);
        MPI_Isend(&o[1*bx+(bx-2)], 1, rp->east_west_type, east, tag, MPI_COMM_WORLD, &reqs[SEND_E]);
        // compute_inner_jacobi(n, o, rp);
        // MPI_Wait(&reqs[SEND_E], &stats[SEND_E]);
        // MPI_Recv(&o[1*bx+(bx-1)], 1, rp->east_west_type, east, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Irecv(&o[1*bx+(bx-1)], 1, rp->east_west_type, east, tag, MPI_COMM_WORLD, &reqs[RECV_E]);
        MPI_Wait(&reqs[RECV_E], &stats[RECV_E]);
        // compute_outer_jacobi(n, o, rp, 2);
        // printf("Iteration %d, Rank %d: Receive from east\n", rp->idx, rp->rank);
        // print_matrix(stdout, o, rp->super.sizex, rp->super.sizey);
    }
    // west
    if (west != MPI_PROC_NULL){
        compute_inner_jacobi(n, o, rp);
        // MPI_Sendrecv(&o[1*bx+1], 1, rp->east_west_type, west, tag, &o[1*bx+0], 1, rp->east_west_type, west, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);        
        MPI_Wait(&reqs[SEND_W], &stats[SEND_W]);
        MPI_Isend(&o[1*bx+1], 1, rp->east_west_type, west, tag, MPI_COMM_WORLD, &reqs[SEND_W]);
        // compute_inner_jacobi(n, o, rp);
        // MPI_Wait(&reqs[SEND_W], &stats[SEND_W]);
        // MPI_Recv(&o[1*bx+0], 1, rp->east_west_type, west, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Irecv(&o[1*bx+0], 1, rp->east_west_type, west, tag, MPI_COMM_WORLD, &reqs[RECV_W]);
        MPI_Wait(&reqs[RECV_W], &stats[RECV_W]);
        // compute_outer_jacobi(n, o, rp, 3);
        // printf("Iteration %d, Rank %d: Receive from west\n", rp->idx, rp->rank);
        // print_matrix(stdout, o, rp->super.sizex, rp->super.sizey);
    }

    // Acount for when size == 1
    if (rp->size == 1){     
        compute_inner_jacobi(n,o,rp);
        // compute_outer_jacobi(n, o, rp, 0);
        // compute_outer_jacobi(n, o, rp, 1);
        // compute_outer_jacobi(n, o, rp, 2);
        // compute_outer_jacobi(n, o, rp, 3);
    }

    // compute_inner_jacobi(n,o,rp);
    compute_outer_jacobi(n, o, rp, 0);
    compute_outer_jacobi(n, o, rp, 1);
    compute_outer_jacobi(n, o, rp, 2);
    compute_outer_jacobi(n, o, rp, 3);

    // Calculate the residual of local matrix after all communication and computation are done
    for( i = 1; i < (rp->super.sizey-1); i++ ) {
        for( j = 1; j < (rp->super.sizex-1); j++ ) {
            diff = n[i*rp->super.sizex+j] - o[i*rp->super.sizex+j];
            sum += diff * diff;
        }
    }

    // printf("Iteration %d: Rank %d \n", rp->idx, rp->rank);
    // print_matrix(stdout, o, rp->super.sizex, rp->super.sizey);

    rp->idx++;

    // printf("Iteration %d, Rank %d: sum = %.6f\n", rp->idx, rp->rank, sum);

    MPI_Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // printf("Iteration %d, Rank %d: global sum = %.6f\n", rp->idx, rp->rank, global_sum);
    // printf("Iteration %d: Rank %d \n", rp->idx, rp->rank);
    // print_matrix(stdout, o, rp->super.sizex, rp->super.sizey);
    return global_sum;

}

static double* jacobi_mpi_nonblocking_relaxation_get_data(relaxation_params_t* grp){   
    struct relaxation_jacobi_mpi_nonblocking_hidden_params_s* rp = (struct relaxation_jacobi_mpi_nonblocking_hidden_params_s*)grp;
    
    return rp->data[(rp->idx + 1) % 2];

}


const struct relaxation_function_class_s _relaxation_jacobi_mpi_nonblocking =
    {
        .type = RELAXATION_JACOBI_MPI_NONBLOCKING,
		.method_name = "MPI Non-blocking for Jacobi relaxation",
		._init = jacobi_mpi_nonblocking_relaxation_init,
        ._fini = jacobi_mpi_nonblocking_relaxation_fini,
        ._coarsen = jacobi_mpi_nonblocking_relaxation_coarsen,
        ._print = jacobi_mpi_nonblocking_relaxation_print,
        ._relaxation = jacobi_mpi_nonblocking_relaxation_apply,
        ._get_data = jacobi_mpi_nonblocking_relaxation_get_data,
    };

