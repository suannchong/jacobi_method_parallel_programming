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

/**
 * This is a template for an empty relaxation that should be used as
 * a starting point for the homework.
 */
struct relaxation_template_hidden_params_s {
    struct relaxation_params_s super;
    double* data;
    uint32_t idx;
};
const struct relaxation_function_class_s _relaxation_template;

static struct relaxation_params_s*
template_relaxation_init(struct hw_params_s* hw_params)
{
    struct relaxation_template_hidden_params_s* rp;
    uint32_t np = hw_params->resolution + 2;

    rp = (struct relaxation_template_hidden_params_s*)malloc(sizeof(struct relaxation_template_hidden_params_s));
    if( NULL == rp ) {
        fprintf(stderr, "Cannot allocate memory for the relaxation structure\n");
        return NULL;
    }

    rp->super.sizex = np;
    rp->super.sizey = np;
    rp->super.rel_class = &_relaxation_template;

    rp->idx = 0;
    rp->data = calloc(np*np, sizeof(double));
    if( NULL == rp->data ) {
        fprintf(stderr, "Cannot allocate the memory for the matrices\n");
        goto fail_and_return;
    }

    return (struct relaxation_params_s*)rp;
 fail_and_return:
    if( NULL != rp->data ) free(rp->data);
    free(rp);
    return NULL;
}

static int template_relaxation_fini(relaxation_params_t** prp)
{
    struct relaxation_template_hidden_params_s* rp = (struct relaxation_template_hidden_params_s*)*prp;
    if( NULL != rp->data ) free(rp->data);
    free(rp);
    *prp = NULL;
    return 0;
}

static int template_relaxation_coarsen(relaxation_params_t* grp, double* dst, uint32_t dstx, uint32_t dsty)
{
    struct relaxation_template_hidden_params_s* rp = (struct relaxation_template_hidden_params_s*)grp;
    return coarsen(rp->data, rp->super.sizex, rp->super.sizey,
                   dst, dstx, dsty);
}

static int template_relaxation_print(relaxation_params_t* grp, FILE* fp)
{
    struct relaxation_template_hidden_params_s* rp = (struct relaxation_template_hidden_params_s*)grp;
    fprintf( fp, "\n# Iteration %d\n", rp->idx);
    print_matrix(fp, rp->data, rp->super.sizex, rp->super.sizey);
    fprintf( fp, "\n\n");
    return 0;
}

/**
 * One step of a simple Template relaxation.
 */
static double template_relaxation_apply(relaxation_params_t* grp)
{
    struct relaxation_template_hidden_params_s* rp = (struct relaxation_template_hidden_params_s*)grp;
    fprintf(stdout, "Iteration %d for the template relaxation class. No computations are done!\n",
            rp->idx);
    rp->idx++;
    return 1.0;
}

static double* template_relaxation_get_data(relaxation_params_t* grp)
{
    return ((struct relaxation_template_hidden_params_s*)grp)->data;
}


const struct relaxation_function_class_s _relaxation_template =
    {
        .type = RELAXATION_TEMPLATE,
		.method_name = "Template for Jacobi relaxation (name to be changed)",
        ._init = template_relaxation_init,
        ._fini = template_relaxation_fini,
        ._coarsen = template_relaxation_coarsen,
        ._print = template_relaxation_print,
        ._relaxation = template_relaxation_apply,
        ._get_data = template_relaxation_get_data,
    };
