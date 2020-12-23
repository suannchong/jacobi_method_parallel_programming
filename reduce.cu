#include <stdio.h>
#include <stdlib.h>

// template <unsigned int blockSize>
__device__ void warpReduce(volatile double* sdata, int tid){
	// if (blockSize >= 64) sdata[tid] += sdata[tid+32];
	// if (blockSize >= 32) sdata[tid] += sdata[tid+16];
	// if (blockSize >= 16) sdata[tid] += sdata[tid+8];
	// if (blockSize >= 8) sdata[tid] += sdata[tid+4];
	// if (blockSize >= 4) sdata[tid] += sdata[tid+2];
	// if (blockSize >= 2) sdata[tid] += sdata[tid+1];
	sdata[tid] += sdata[tid+32];
	sdata[tid] += sdata[tid+16];
	sdata[tid] += sdata[tid+8];
	sdata[tid] += sdata[tid+4];
	sdata[tid] += sdata[tid+2];
 	sdata[tid] += sdata[tid+1];
}

__global__ void reduce(double* g_idata, double*g_odata){
	int N = 16;
	__shared__ double sdata[16];

	// int tid = threadIdx.x;
	// int i = blockIdx.x*(blockDim.x) + threadIdx.x;
	// sdata[tid] = g_idata[i]; //+ g_idata[i+blockDim.x];

	for (int i=blockIdx.x*blockDim.x + threadIdx.x; i<N; i+=blockDim.x * gridDim.x){
		sdata[i] = g_idata[i];
		// printf("g_idata[%d] = %.6f\n", i, g_idata[i] );
		// printf("sdata[%d] = %.6f\n", i, sdata[i] );
		__syncthreads();
	}

	for (int i=blockIdx.x*blockDim.x + threadIdx.x; i<N; i+=blockDim.x * gridDim.x){
		for (unsigned int s=blockDim.x/2; s>0; s*=2){ 
			if (i < s && i+s < N){
				printf("before: sdata[%d] = %.6f \t sdata[%d+%d] = %.6f\n", i, sdata[i], i, s, sdata[i+s]);
				sdata[i] += sdata[i+s];
				printf("reduce loop: sdata[%d] = %.6f\n", i, sdata[i] );
			}

		// 	if ((i+s )< 10){
		// 		printf("before: sdata[%d] = %.6f \t sdata[%d+%d] = %.6f\n", i, sdata[i], i, s, sdata[i+s]);
		// 		sdata[i] += sdata[i+s];
		// 		printf("reduce loop: sdata[%d] = %.6f\n", i, sdata[i] );
		// 	}
			__syncthreads();
		}
	}

	// int i = blockIdx.x*(blockDim.x) + threadIdx.x;
	// // do reduction in shared mem
	// for (unsigned int s=blockDim.x/2; s>0; s>>=1){
		
	// 	if ((i+s )< 10){
	// 		printf("before: sdata[%d] = %.6f \t sdata[%d+%d] = %.6f\n", i, sdata[i], i, s, sdata[i+s]);
	// 		sdata[i] += sdata[i+s];
	// 		printf("reduce loop: sdata[%d] = %.6f\n", i, sdata[i] );
	// 	}
	// 	// printf("sdata[%d] = %.6f\n", i, sdata[i] );
	// 	__syncthreads();
	// } 

	// if (blockSize >= 512) {
	// 	if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }

	// if (blockSize >= 256) {
	// 	if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }

	// if (blockSize >= 128) {
	// 	if (tid < 64) { sdata[tid] += sdata[tid + 64]; } __syncthreads(); }

	// if (tid < 32) warpReduce(sdata, tid);

	// write result for this block to global mem
	// if (tid == 0) {
	// 	g_odata[blockIdx.x] = sdata[0];
	// 	printf("sdata[%d] = %.6f\n", tid, sdata[tid] );
	// }

}

int main(int argc, char** argv){
	int N = 16;

	double* array_in = (double*) calloc(1,N*sizeof(double));
	double* array_out = (double*) calloc(1, N*sizeof(double));
	
	double* d_array_in, *d_array_out;
	int size = N*sizeof(double);

	for (int i=0; i<N;i++){
		array_in[i] = i*0.1;
		// printf("array_in[%d] = %.6f\n", i, array_in[i] );
	}

	cudaMalloc((void**) &d_array_in, size);
	cudaMalloc((void**) &d_array_out, size);

	cudaMemcpy(d_array_in, array_in, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_array_out, array_out, size, cudaMemcpyHostToDevice);

	reduce<<<1,4>>>(d_array_in, d_array_out);

	cudaMemcpy(array_out, d_array_out, size, cudaMemcpyDeviceToHost);

	printf("sum = %.6f\n", array_out[0]);

	free(array_out); free(array_in);
	cudaFree(d_array_in); cudaFree(d_array_out);
	return 0;

}