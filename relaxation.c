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
// name mangling 
#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

// // name mangling 
// #if defined(__cplusplus) || defined(c_plusplus)
// extern "C" {
// #endif
    #include "jacobi.h"
    #include "header.h"
// #if defined(__cplusplus) || defined(c_plusplus)
// }
// #endif

extern const struct relaxation_function_class_s _relaxation_jacobi;
extern const struct relaxation_function_class_s _relaxation_template;
extern const struct relaxation_function_class_s _relaxation_jacobi_mpi_blocking;
extern const struct relaxation_function_class_s _relaxation_jacobi_mpi_nonblocking;
extern const struct relaxation_function_class_s _relaxation_jacobi_mpi_collective;
extern const struct relaxation_function_class_s _relaxation_jacobi_mpi_rma;
extern const struct relaxation_function_class_s _relaxation_jacobi_openmp;

// // name mangling 
// #if defined(__cplusplus) || defined(c_plusplus)
// extern "C" {
// #endif

extern const struct relaxation_function_class_s _relaxation_jacobi_cuda;



const struct relaxation_function_class_s * const _relaxation_classes[] =
    {&_relaxation_jacobi, &_relaxation_template, 
     &_relaxation_jacobi_mpi_blocking, &_relaxation_jacobi_mpi_nonblocking, 
     &_relaxation_jacobi_mpi_collective, &_relaxation_jacobi_mpi_rma,
     &_relaxation_jacobi_openmp, &_relaxation_jacobi_cuda, NULL};

// #if defined(__cplusplus) || defined(c_plusplus)
// }
// #endif

relaxation_params_t* relaxation_init(struct hw_params_s* hw_params)
{
    relaxation_params_t* rp = NULL;

    for( int i = 0; NULL != _relaxation_classes[i]; i++ ) {
        
        if( _relaxation_classes[i]->type != hw_params->alg_type ){
           // printf("i = %d, type = %u\n", i, _relaxation_classes[i]->type);
            continue;
        }
            

        // printf("i = %d\n", i );
        if( NULL == (rp =  _relaxation_classes[i]->_init(hw_params)) ) {
            return NULL;
        }

        return rp;
    }
    fprintf(stderr, "Unable to find a suitable Jacobi for the requested algorithm %d\n",
            hw_params->alg_type);
    return NULL;
}

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif


#if 0
/**
 * One step of a red-black solver.
 *
 * The difference with a Jacobi is that the updates are happening
 * in place.
 */
double redblack_relaxation(relaxation_params_t* rp,
                           double* o,
                           uint32_t sizex,
                           uint32_t sizey)
{
    double nval, diff, sum = 0.0;
    int i, j;

    for( i = 1; i < sizex; i++ ) {
        for( j = 1; j < sizex; j++ ) {
            nval = 0.25 * (o[ i    *sizey + (j-1) ]+  // left
                           o[ i    *sizey + (j+1) ]+  // right
                           o[ (i-1)*sizey + j     ]+  // top
                           o[ (i+1)*sizey + j     ]); // bottom
            diff = nval - o[i*sizey+j];
            sum += diff * diff;
            o[i*sizey+j] = nval;
        }
    }
    return sum;
}
#endif

