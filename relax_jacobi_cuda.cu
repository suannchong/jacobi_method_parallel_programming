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

// name mangling 
#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif
    #include "jacobi.h"
    #include "header.h"
#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

struct relaxation_jacobi_cuda_hidden_params_s {
    struct relaxation_params_s super;
    double* data[2];
    uint32_t idx;
    // cuda params
    double* d_data[2];
    // double* d_n, *d_o;
    size_t cuda_malloc_size;
    double* d_sum;
    int* d_nx, *d_ny;
    double* sum_array;
};

    static struct relaxation_params_s* jacobi_cuda_relaxation_init(struct hw_params_s* hw_params);
    static int jacobi_cuda_relaxation_fini(relaxation_params_t** prp);
    static int jacobi_cuda_relaxation_coarsen(relaxation_params_t* grp, double* dst, uint32_t dstx, uint32_t dsty);
    static int jacobi_cuda_relaxation_print(relaxation_params_t* grp, FILE* fp);
    static double jacobi_cuda_relaxation_apply(relaxation_params_t* grp);
    static double* jacobi_cuda_relaxation_get_data(relaxation_params_t* grp);

// name mangling 
#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif
    struct relaxation_function_class_s _relaxation_jacobi_cuda = {
    .type = RELAXATION_JACOBI_CUDA,
    .method_name = NULL ,//"CUDA for Jacobi relaxation",
    ._init = jacobi_cuda_relaxation_init,
    ._fini = jacobi_cuda_relaxation_fini,
    ._coarsen = jacobi_cuda_relaxation_coarsen,
    ._print = jacobi_cuda_relaxation_print, 
    ._relaxation = jacobi_cuda_relaxation_apply,
    ._get_data = jacobi_cuda_relaxation_get_data,
};

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

#define THREADS_PER_BLOCK 32
#define BLOCK_SIZE 16


__global__ void compute_jacobi_kernel(double* d_n, double* d_o, int* d_nx, int* d_ny, double* d_sum){
    int N = (*d_nx) * (*d_ny);

    // for (int index = threadIdx.x + blockIdx.x * blockDim.x; index < N; index += blockDim.x*gridDim.x){
        int index = threadIdx.x + blockIdx.x * blockDim.x;
        int i = index / (*d_nx);
        int j = index % (*d_nx);

        double diff;
        if ((i > 0) && (i < ((*d_ny)-1)) && (j > 0) && (j < ((*d_nx)-1))){
            d_n[i*(*d_nx)+j] = 0.25*(  d_o[i*(*d_nx)+(j-1)] + d_o[i*(*d_nx)+(j+1)] +
                                   d_o[(i-1)*(*d_nx)+j] + d_o[(i+1)*(*d_nx)+j] );
            diff = d_n[i*(*d_nx)+j] - d_o[i*(*d_nx)+j];
            d_sum[i*(*d_nx)+j] =  diff*diff;
        }
    // }
        
}

static struct relaxation_params_s*
jacobi_cuda_relaxation_init(struct hw_params_s* hw_params)
{
    struct relaxation_jacobi_cuda_hidden_params_s* rp;
    uint32_t np = hw_params->resolution + 2;

    rp = (struct relaxation_jacobi_cuda_hidden_params_s*)malloc(sizeof(struct relaxation_jacobi_cuda_hidden_params_s));
    if( NULL == rp ) {
        fprintf(stderr, "Cannot allocate memory for the relaxation structure\n");
        return NULL;
    }

    rp->super.sizex = np;
    rp->super.sizey = np;
    rp->super.rel_class = &_relaxation_jacobi_cuda;
    

    rp->idx = 0;
    rp->data[0] = (double*) calloc(np*np, sizeof(double));
    rp->data[1] = (double*) malloc(np*np*sizeof(double));
    rp->sum_array = (double*) malloc(np*np*sizeof(double));

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

    // Allocate space for device copies 
    rp->cuda_malloc_size = rp->super.sizex*rp->super.sizey*sizeof(double);
    cudaMalloc((void**) &rp->d_data[0], rp->cuda_malloc_size);
    cudaMalloc((void**) &rp->d_data[1], rp->cuda_malloc_size);
    cudaMalloc((void**) &rp->d_sum, rp->cuda_malloc_size);
    cudaMalloc((void**) &rp->d_nx, sizeof(int));
    cudaMalloc((void**) &rp->d_ny, sizeof(int));

    // print_matrix(stdout, o, np, np);
    // printf("\n");

    // copy arrays from host to device
    cudaMemcpy(rp->d_data[1], rp->data[1], rp->cuda_malloc_size, cudaMemcpyHostToDevice);
    cudaMemcpy(rp->d_data[0], rp->data[0], rp->cuda_malloc_size, cudaMemcpyHostToDevice);
    cudaMemcpy(rp->d_nx, &(rp->super.sizex), sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(rp->d_ny, &(rp->super.sizey), sizeof(int), cudaMemcpyHostToDevice);


    return (struct relaxation_params_s*)rp;
 fail_and_return:
    if( NULL != rp->data[0] ) free(rp->data[0]);
    if( NULL != rp->data[1] ) free(rp->data[1]);
    free(rp);
    return NULL;
}

static int jacobi_cuda_relaxation_fini(relaxation_params_t** prp)
{
    struct relaxation_jacobi_cuda_hidden_params_s* rp = (struct relaxation_jacobi_cuda_hidden_params_s*)*prp;
    if( NULL != rp->data[0] ) free(rp->data[0]);
    if( NULL != rp->data[1] ) free(rp->data[1]);

    // free space for device 
    cudaFree(rp->d_data[0]); cudaFree(rp->d_data[1]);
    cudaFree(rp->d_sum); cudaFree(rp->d_nx); cudaFree(rp->d_ny);

    free(rp);
    *prp = NULL;
    return 0;
}

static int jacobi_cuda_relaxation_coarsen(relaxation_params_t* grp, double* dst, uint32_t dstx, uint32_t dsty)
{
    struct relaxation_jacobi_cuda_hidden_params_s* rp = (struct relaxation_jacobi_cuda_hidden_params_s*)grp;

    return coarsen(rp->data[(rp->idx + 1) % 2], rp->super.sizex, rp->super.sizey,
                   dst, dstx, dsty);
}

static int jacobi_cuda_relaxation_print(relaxation_params_t* grp, FILE* fp)
{
    struct relaxation_jacobi_cuda_hidden_params_s* rp = (struct relaxation_jacobi_cuda_hidden_params_s*)grp;
    fprintf( fp, "\n# Iteration %d\n", rp->idx);
    print_matrix(fp, rp->data[(rp->idx + 1) % 2], rp->super.sizex, rp->super.sizey);
    fprintf( fp, "\n\n");
    return 0;
}

/**
 * One step of a simple jacobi_cuda relaxation.
 */
static double jacobi_cuda_relaxation_apply(relaxation_params_t* grp)
{
    struct relaxation_jacobi_cuda_hidden_params_s* rp = (struct relaxation_jacobi_cuda_hidden_params_s*)grp;
    double sum = 0.0, *n, *o;

    // printf("Iteration %d\n", rp->idx);
    // print_matrix(stdout, o, rp->super.sizex, rp->super.sizey);

    // Launch computer_jacobi kernel on GPU
    int num_blocks_x = ceil(rp->super.sizex*rp->super.sizey/THREADS_PER_BLOCK);
    int num_blocks_y = ceil(rp->super.sizey/BLOCK_SIZE);
    dim3 grid(num_blocks_x, 1);
    dim3 block(THREADS_PER_BLOCK, 1);
    // size_t dynamic_size =  rp->super.sizex*rp->super.sizey*sizeof(double);
    // compute_jacobi_kernel<<< 1,1024 >>>(rp->d_data[(rp->idx+0)%2], rp->d_data[(rp->idx+1)%2], rp->d_nx, rp->d_ny, rp->d_sum);
    compute_jacobi_kernel<<< grid,block >>>(rp->d_data[(rp->idx+0)%2], rp->d_data[(rp->idx+1)%2], rp->d_nx, rp->d_ny, rp->d_sum);
    
    // // Copy results back from device to host 
    cudaMemcpy(rp->data[(rp->idx+0)%2], rp->d_data[(rp->idx+0)%2], rp->cuda_malloc_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(rp->data[(rp->idx+1)%2], rp->d_data[(rp->idx+1)%2], rp->cuda_malloc_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(rp->sum_array, rp->d_sum, (rp->super.sizex)*(rp->super.sizey)*sizeof(double), cudaMemcpyDeviceToHost);

    for (int i=0; i<(rp->super.sizex*rp->super.sizey);i++){
        // printf("sum_array[%d] = %.6f\n", i, sum_array[i]);
        sum += rp->sum_array[i];
    }

    // printf("Iteration %d\n", rp->idx);
    // print_matrix(stdout, n, rp->super.sizex, rp->super.sizey);

    rp->idx++;

    // printf("Iteration %d: global sum = %.6f\n", rp->idx, sum);
    

    return sum;
}

static double* jacobi_cuda_relaxation_get_data(relaxation_params_t* grp)
{
    struct relaxation_jacobi_cuda_hidden_params_s* rp = (struct relaxation_jacobi_cuda_hidden_params_s*)grp;
    return rp->data[(rp->idx + 1) % 2];
}