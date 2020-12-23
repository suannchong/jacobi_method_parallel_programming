#include<iostream>
#include<cmath>
#include "gpu_reduce.hpp"

#define BLOCKSIZE 1024

int main()
{
    using my_t = unsigned;

    size_t N0 = 14; // Grid has a side length of 2^N0 + 1
    size_t n = std::pow(2, N0) + 1;
    size_t N = n*n;

    // Populate array
    my_t* A = (my_t*) malloc(sizeof(my_t) * N);
    for (size_t i = 0; i < N; i++)
        A[i] = 1;

    my_t* dA;
    cudaMalloc(&dA, sizeof(my_t)*N); checkCUDAError("Error allocating dA");
    cudaMemcpy(dA, A, sizeof(my_t)*N, cudaMemcpyHostToDevice); checkCUDAError("Error copying A"); 

    my_t tot = 0.;

    size_t nPar1 = (n-1) * (n-1);
    size_t nPar2 = N - nPar1 - 1;

    tot  = GPUReduction<BLOCKSIZE>(dA, nPar1);
    tot += GPUReduction<BLOCKSIZE>(dA + nPar1, nPar2);
    tot += A[nPar1 + nPar2];

    std::cout << "N: " << N << std::endl;
    std::cout << "Result: " << tot << std::endl;
}