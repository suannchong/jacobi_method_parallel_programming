#include <stdio.h>
#include <stdlib.h>

__global__ void kernel(double* d_array, double* d_sum){
	// int i = (blockIdx.x*blockDim.x) + threadIdx.x;
	// int x = i % 2560;
	// int y = i % 2560;

	for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < 2560*2560; i += blockDim.x*gridDim.x){
		d_array[i] = i+0.1;
		printf("d_array[%d] = %.2f\n", i ,d_array[i] );
		*d_sum += d_array[i];
	}
}

int main(int argc, char** argv){
	int N = 2560*2560;
	double* array = (double*) calloc(N,sizeof(double));
	double sum = 2.03;
	double* d_array, *d_sum;
	size_t size = N*sizeof(double);

	cudaMalloc((void**) &d_array, size);
	cudaMalloc((void**) &d_sum, sizeof(double));

	cudaMemcpy(d_array,array,size,cudaMemcpyHostToDevice);
	cudaMemcpy(d_sum, &sum, sizeof(double), cudaMemcpyHostToDevice);

	for (int k = 0; k < 100000; k++){
		kernel<<<1,1024>>>(d_array,d_sum);

		cudaMemcpy(array, d_array, size, cudaMemcpyDeviceToHost);
		cudaMemcpy(&sum, d_sum, sizeof(double), cudaMemcpyDeviceToHost);

		// for (int i=0; i<N; i++){
		// 	printf("array[%d] = %.2f\n",i, array[i] );
		// }

		printf("sum = %.2f\n", sum );
	}

	cudaFree(d_array);
	cudaFree(d_sum);
	return 0;
}