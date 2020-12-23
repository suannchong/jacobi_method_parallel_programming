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

#ifndef COSC462_HW_HEADER_HAS_BEEN_INCLUDED
#define COSC462_HW_HEADER_HAS_BEEN_INCLUDED

#include <stdio.h>
#include <stdint.h>

#define COSC462_HW_FLAG_GENERATE_HEAT_IMAGE  0x0001
#define COSC462_HW_FLAG_TIME_ALG             0x0002
#define COSC462_HW_FLAG_VERBOSE              0x0010
#define COSC462_HW_FLAG_VERY_VERBOSE         0x0020

typedef struct heat_source_s {
    float x;
    float y;
    double range;
    double temp;
} heat_source_t;

typedef struct hw_params_s {
    uint32_t alg_type;
    uint32_t max_iterations;  /* leave to 0 to only stop on residual */
    uint32_t resolution;
    uint32_t num_threads;
    /* For the OpenSHMEM version */
    uint32_t num_proc_p;
    int32_t flags;
    char* input;
    uint32_t vis_res;  /* visual resolution if the data is dumped in a graphical format */
    int32_t vis_step;  /* number of iteration between generating the graphical heat map.
                        * Set to negative to only get 2 images, one for the initial stage
                        * and one for the final.
                        */
    char* vis_output;
    double *vis_data;
    uint32_t num_sources;
    heat_source_t* heat_sources;
} hw_params_t;

int hw_init(hw_params_t* param, int32_t* argc, char* argv[]);
int hw_fini(hw_params_t* param);

/**
 * Transform, coarsen, a matrix by agragating several point
 * in 2D into a single element. This can be used to prepare
 * the output matrix for generating the graphical visualization.
 */
int coarsen(const double* src, uint32_t srcx, uint32_t srcy,
            double* dst, uint32_t dstx, uint32_t dsty);

#define COSC462_RELAXATION_INIT_HAS_LM  1
/**
 * A function to setup the initial values for the boundaries based on the
 * known heat sources. It is generic for most of the relaxations, for as long
 * as the provided size accounts the boundaries. By providing the correct
 * position of the initial location we can use this function to initialize
 * the matrix in distributed environments.
 */
int relaxation_matrix_set_ex(hw_params_t* hw_params,
                             double* mat,  /* pointer to the local matrix including the ghost regions */
                             uint32_t ln,  /* # local elements including ghost in the x dimension */
                             uint32_t lm,  /* # local elements including ghost in the y dimension */
                             int rowpos,   /* my row position on the global matrix */
                             int colpos);  /* my col position on the global matrix */
/* an initialization routine for square local matrices */
int relaxation_matrix_set(hw_params_t* hw_params,
                          double* mat,  /* pointer to the local matrix including the ghost regions */
                          uint32_t ln,  /* # local elements including ghost */
                          int rowpos,   /* my row position on the global matrix */
                          int colpos);  /* my col position on the global matrix */

/**
 * The number of different values of gray in the output image */
#define MAXVAL 255

/**
 * Dump a matrix into a gray ppm format.
 */
int dump_gray_image(char* template_name, int iter, double* data, uint32_t sizex, uint32_t sizey);

/**
 * Dump a matrix into a RGB ppm format.
 */
int dump_rgb_image(char* template_name, int iter, double* data, uint32_t sizex, uint32_t sizey);

/**
 * Time in microsecond from an initial date. Only the difference between
 * 2 such values is meaningful, as it represent the number of microseconds
 * between the 2 events.
 */
double wtime(void);

int print_matrix( FILE* fp, const double* mat, int sizex, int sizey );

#endif  /* COSC462_HW_HEADER_HAS_BEEN_INCLUDED */
