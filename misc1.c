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
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "header.h"
static double coarsen_block( const double* src, int posx, int posy, int ldsrc, int ncols, int nrows )
{
    double val = 0.0;
    const double* mat;
    int i, j;

    mat = src + posy * ldsrc + posx;
    for( j = 0; j < nrows; j++ ) {
        for( i = 0; i < ncols; i++ ) {
            val += mat[i];
        }
        mat += ldsrc;
    }
    return val / (ncols * nrows);
}

/**
 * Transform, coarsen, a matrix by aggregating several point
 * in 2D into a single element. This can be used to prepare
 * the output matrix for generating the graphical visualization.
 */
int coarsen(const double* src, uint32_t src_rows, uint32_t src_cols,
            double* dst, uint32_t dst_rows, uint32_t dst_cols)
{
    uint32_t i, j, maxx, maxy, pos_src_x, pos_src_y;
    uint32_t stepx, stepy;

    if( (src_rows < dst_rows) || (src_cols < dst_cols) ) {
        fprintf(stderr, "The coarsen functions should only be used to reduce the resolution of the matrix.\n");
        return -1;
    }
    if( (dst_cols <= 2) || (dst_rows <= 2) ) {
        fprintf(stderr, "The destination matrix cannot be a single line in any direction\n");
    }
    stepx = (src_rows - 2 + dst_rows - 2 - 1) / (dst_rows - 2);
    stepy = (src_cols - 2 + dst_cols - 2 - 1) / (dst_cols - 2);

    pos_src_x = 1;
    /* The borders are critical, minimize the coarsening to a single dimension. */
    for( i = 1; i < (dst_cols-1); i++ ) {  /* top and bottom */
        maxx = stepx;
        if( (pos_src_x + maxx) > (src_cols - 1) )
            maxx = src_cols - 1 - pos_src_x;
        dst[i] = coarsen_block( src, pos_src_x, 0, src_cols, maxx, 1);
        dst[(dst_cols * (dst_rows-1)) + i] = coarsen_block( src, pos_src_x, src_rows - 1, src_cols, maxx, 1);
        pos_src_x += maxx;
    }
    pos_src_y = 1;
    for( j = 1; j < (dst_rows-1); j++ ) {  /* left and right */
        maxy = stepy;
        if( (pos_src_y + maxy) > (src_rows - 1) )
            maxx = src_rows - 1 - pos_src_y;
        dst[j * dst_cols] = coarsen_block( src, 0, pos_src_y, src_cols, 1, maxy);
        dst[j * dst_cols + dst_rows - 1] = coarsen_block( src, src_cols - 1, pos_src_y, src_cols, 1, maxy);
        pos_src_y += maxy;
    }
    pos_src_x = 1;
    for( i = 1; i < (dst_cols-1); i++ ) {
        maxx = stepx;
        if( (pos_src_x + maxx) > (src_cols - 1) )
            maxx = src_cols - 1 - pos_src_x;
        pos_src_y = 1;
        for( j = 1; j < (dst_rows-1); j++ ) {
            maxy = stepy;
            if( (pos_src_y + maxy) > (src_rows - 1) )
                maxx = src_rows - 1 - pos_src_y;
            dst[j*dst_cols+i] = coarsen_block( src, pos_src_x, pos_src_y, src_cols, maxx, maxy);
            pos_src_y += maxy; 
        }
        pos_src_x += maxx;
    }
    return 0;
}

int print_matrix( FILE* fp, const double* mat, int sizex, int sizey )
{
    int i, j;

    for( j = 0; j < sizey; j++ ) {
        for( i = 0; i < sizex; i++ ) {
            fprintf(fp, "%10.4lf", mat[j * sizex + i]);
        }
        fprintf(fp, "\n");
    }
    return 0;
}

/**
 * A function to setup the initial values for the boundaries based on the
 * known heat sources. It is generic for most of the relaxations, for as long
 * as the provided size accounts the boundaries. By providing the total number
 * of participants (assuming each participant has it's own matrix), and the
 * identification of the caller on the participants grid (2D) we can set the
 * boundaries only for the processes that are on the border.
 */
int relaxation_matrix_set_ex(hw_params_t* hw_params,
                             double* mat,
                             uint32_t ln,  /* # local elements including ghost in the x dimension */
                             uint32_t lm,  /* # local elements including ghost in the y dimension */
                             int rowpos,   /* my row position on the global matrix */
                             int colpos)   /* my col position on the global matrix */
{
    uint32_t i, j;
    double dist;
    int N = hw_params->resolution;

    for( i = 0; i < hw_params->num_sources; i++ ) {
        /**
         * The heat dissipate linearly in circles around the central
         * point up to the defined range. It only affects the interface
         * between the mediums, so it only has an impact on the boundaries.
         */
        if( 0 == colpos ) {
            for( j = 1; j < ln-1; j++ ) {  /* initialize the top row */
                dist = sqrt( pow((double)(rowpos+j)/(double)(N+1) -
                                 hw_params->heat_sources[i].x, 2) +
                             pow(hw_params->heat_sources[i].y, 2));
                if( dist <= hw_params->heat_sources[i].range ) {
                    mat[j] += ((hw_params->heat_sources[i].range - dist) /
                               hw_params->heat_sources[i].range *
                               hw_params->heat_sources[i].temp);
                }
            }
        }
        if( colpos >= (N - (lm - 2)) ) {
            for( j = 1; j < ln-1; j++ ) {  /* initialize the bottom row */
                dist = sqrt( pow((double)(rowpos+j)/(double)(N+1) -
                                 hw_params->heat_sources[i].x, 2) +
                             pow(1-hw_params->heat_sources[i].y, 2));
                if( dist <= hw_params->heat_sources[i].range ) {
                    mat[(lm-1)*ln+j] += ((hw_params->heat_sources[i].range - dist) /
                                         hw_params->heat_sources[i].range *
                                         hw_params->heat_sources[i].temp);
                }
            }
        }
        if( 0 == rowpos ) {
            for( j = 1; j < lm-1; j++ ) {  /* left-most column */
                dist = sqrt( pow(hw_params->heat_sources[i].x, 2) +
                             pow((double)(colpos+j)/(double)(N+1) -
                                 hw_params->heat_sources[i].y, 2));
                if( dist <= hw_params->heat_sources[i].range ) {
                    mat[j*ln] += ((hw_params->heat_sources[i].range - dist) /
                                  hw_params->heat_sources[i].range *
                                  hw_params->heat_sources[i].temp);
                }
            }
        }
        if( rowpos >= (N - (ln - 2)) ) {
            for( j = 1; j < lm-1; j++ ) {  /* right-most column */
                dist = sqrt( pow(1-hw_params->heat_sources[i].x, 2) +
                             pow((double)(colpos+j)/(double)(N+1) -
                                 hw_params->heat_sources[i].y, 2));
                if( dist <= hw_params->heat_sources[i].range ) {
                    mat[j*ln+(ln-1)] += ((hw_params->heat_sources[i].range - dist) /
                                         hw_params->heat_sources[i].range *
                                         hw_params->heat_sources[i].temp);
                }
            }
        }
    }
    return 0;
}
/**
 * A second version of the initialization where the local part of the distributed
 * matrix is square. It exists only for backward compatibility.
 */
int relaxation_matrix_set(hw_params_t* hw_params,
                          double* mat,
                          uint32_t ln,  /* # local elements including ghost */
                          int rowpos,   /* my row position on the global matrix */
                          int colpos)   /* my col position on the global matrix */
{
    return relaxation_matrix_set_ex(hw_params, mat, ln, ln, rowpos, colpos);
}

#include <float.h>
#include <limits.h>

#include "header.h"

/**
 * Dump a matrix into a gray pgm format.
 * http://netpbm.sourceforge.net/doc/pgm.html
 */
int dump_gray_image(char* template_name, int iter,
                    double* data, uint32_t sizex, uint32_t sizey)
{
    double min = DBL_MAX, max = DBL_MIN, factor;
    int i, j;
    int16_t *val16;
    int8_t *val8;
    FILE* fp;
    char* fname;

    asprintf(&fname, "%s%04d.pgm", template_name, iter);
    fp = fopen(fname, "w");
    if( NULL == fp ) {
        fprintf(stderr, "Can't generate a proper filename for the heatmap image\n");
        free(fname);
        return -1;
    }
    for( i = 0; i < (sizex * sizey); i++ ) {
        if( min > data[i] ) { if( 0.0 != data[i] ) min = data[i]; }
        else if( max < data[i] ) { if ( 0.0 != data[i] ) max = data[i]; }
    }
    factor = (double)MAXVAL / (max - min);
    fprintf(fp, "P5 %d %d %d\n", sizex, sizey, MAXVAL);
    if( MAXVAL >= 256 ) {
        val16 = (int16_t*)malloc(sizex * sizeof(int16_t));
        for( j = 0; j < sizey; j++ ) {
            for( i = 0; i < sizex; i++ ) {
                val16[i] = 0;  /* gracefully handle very small values. */
                if( data[j*sizex+i] > min )
                    val16[i] = (int16_t)ceil(factor * (data[j*sizex+i]-min));
            }
            fwrite(val16, sizeof(int16_t), sizex, fp);
        }
        free(val16);
    } else {
        val8 = (int8_t*)malloc(sizex * sizeof(int8_t));
        for( j = 0; j < sizey; j++ ) {
            for( i = 0; i < sizex; i++ ) {
                val8[i] = 0;  /* gracefully handle very small values. */
                if( data[j*sizex+i] > min )
                    val8[i] = (int8_t)ceil(factor * (data[j*sizex+i]-min));
            }
            fwrite(val8, sizeof(int8_t), sizex, fp);
        }
        free(val8);
    }
    fclose(fp);
    free(fname);
    return 0;
}

/**
 * Dump a matrix into a RGB ppm format.
 */
int dump_rgb_image(char* template_name, int iter,
                   double* data, uint32_t sizex, uint32_t sizey)
{
    return 0;
}

/**
 * Time in microsecond from an initial date. The difference between 2
 * such values is meaningful, as it represent the number of microseconds
 * between the 2 events.
 */
#if HAVE_TIME_H
#include <time.h>
#else
#include <sys/time.h>
#endif  /* HAVE_TIME_H */

double wtime(void)
{
#if HAVE_CLOCK_GETTIME
    struct timespec tp;
    clock_gettime(CLOCK_REALTIME, &tp);
    return (double)(tp.tv_sec * 1000000 + tp.tv_nsec / 1000);
#else
    struct timeval tv;
    gettimeofday( &tv, NULL);
    return (double)(tv.tv_sec * 1000000 + tv.tv_usec);
#endif  /* HAVE_CLOCK_GETTIME */
}
