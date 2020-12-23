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

#include "jacobi.h"
#include "header.h"

struct relaxation_jacobi_hidden_params_s {
    struct relaxation_params_s super;
    double* data[2];
    uint32_t idx;
};
const struct relaxation_function_class_s _relaxation_jacobi;

static struct relaxation_params_s*
jacobi_relaxation_init(struct hw_params_s* hw_params)
{
    struct relaxation_jacobi_hidden_params_s* rp;
    uint32_t np = hw_params->resolution + 2;

    rp = (struct relaxation_jacobi_hidden_params_s*)malloc(sizeof(struct relaxation_jacobi_hidden_params_s));
    if( NULL == rp ) {
        fprintf(stderr, "Cannot allocate memory for the relaxation structure\n");
        return NULL;
    }

    rp->super.sizex = np;
    rp->super.sizey = np;
    rp->super.rel_class = &_relaxation_jacobi;

    rp->idx = 0;
    rp->data[0] = calloc(np*np, sizeof(double));
    rp->data[1] = (double*)malloc(np*np*sizeof(double));

    if( (NULL == rp->data[0]) && (NULL == rp->data[1]) ) {
        fprintf(stderr, "Cannot allocate the memory for the matrices\n");
        goto fail_and_return;
    }
    /* we can use the square version of the initializer */
    relaxation_matrix_set(hw_params, rp->data[0], np, 0, 0);
    // print_matrix(stdout, rp->data[0], np, np);
    // printf("\n");

    /* Copy the boundary conditions on all matrices */
    memcpy(rp->data[1], rp->data[0], np*np*sizeof(double));

    return (struct relaxation_params_s*)rp;
 fail_and_return:
    if( NULL != rp->data[0] ) free(rp->data[0]);
    if( NULL != rp->data[1] ) free(rp->data[1]);
    free(rp);
    return NULL;
}

static int jacobi_relaxation_fini(relaxation_params_t** prp)
{
    struct relaxation_jacobi_hidden_params_s* rp = (struct relaxation_jacobi_hidden_params_s*)*prp;
    if( NULL != rp->data[0] ) free(rp->data[0]);
    if( NULL != rp->data[1] ) free(rp->data[1]);
    free(rp);
    *prp = NULL;
    return 0;
}

static int jacobi_relaxation_coarsen(relaxation_params_t* grp, double* dst, uint32_t dstx, uint32_t dsty)
{
    struct relaxation_jacobi_hidden_params_s* rp = (struct relaxation_jacobi_hidden_params_s*)grp;

    return coarsen(rp->data[(rp->idx + 1) % 2], rp->super.sizex, rp->super.sizey,
                   dst, dstx, dsty);
}

static int jacobi_relaxation_print(relaxation_params_t* grp, FILE* fp)
{
    struct relaxation_jacobi_hidden_params_s* rp = (struct relaxation_jacobi_hidden_params_s*)grp;
    fprintf( fp, "\n# Iteration %d\n", rp->idx);
    print_matrix(fp, rp->data[(rp->idx + 1) % 2], rp->super.sizex, rp->super.sizey);
    fprintf( fp, "\n\n");
    return 0;
}

/**
 * One step of a simple Jacobi relaxation.
 */
static double jacobi_relaxation_apply(relaxation_params_t* grp)
{
    struct relaxation_jacobi_hidden_params_s* rp = (struct relaxation_jacobi_hidden_params_s*)grp;
    double diff, sum = 0.0, *n, *o;
    int i, j;

    n = rp->data[(rp->idx + 0) % 2];
    o = rp->data[(rp->idx + 1) % 2];

    // printf("Iteration %d:\n", rp->idx);
    // print_matrix(stdout, o, rp->super.sizex, rp->super.sizey);
    // printf("\n");

    /* The size[x,y] account for the boundary */
    for( i = 1; i < (rp->super.sizey-1); i++ ) {
        for( j = 1; j < (rp->super.sizex-1); j++ ) {
            n[i*rp->super.sizex+j] = 0.25 * (o[ i     * rp->super.sizex + (j-1) ]+  // left
                                       o[ i     * rp->super.sizex + (j+1) ]+  // right
                                       o[ (i-1) * rp->super.sizex + j     ]+  // top
                                       o[ (i+1) * rp->super.sizex + j     ]); // bottom
            diff = n[i*rp->super.sizex+j] - o[i*rp->super.sizex+j];
            sum += diff * diff;
        }
    }

    // printf("Iteration %d\n", rp->idx);
    // print_matrix(stdout, n, rp->super.sizex, rp->super.sizey);

    rp->idx++;

//    printf("Iteration %d: global sum = %.6f\n", rp->idx, sum);
    

    return sum;
}

static double* jacobi_relaxation_get_data(relaxation_params_t* grp)
{
    struct relaxation_jacobi_hidden_params_s* rp = (struct relaxation_jacobi_hidden_params_s*)grp;
    return rp->data[(rp->idx + 1) % 2];
}


const struct relaxation_function_class_s _relaxation_jacobi =
    {
        .type = RELAXATION_JACOBI,
		.method_name = "Sequential Jacobi relaxation",
		._init = jacobi_relaxation_init,
        ._fini = jacobi_relaxation_fini,
        ._coarsen = jacobi_relaxation_coarsen,
        ._print = jacobi_relaxation_print,
        ._relaxation = jacobi_relaxation_apply,
        ._get_data = jacobi_relaxation_get_data,
    };
