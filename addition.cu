#include <stdio.h>
#include <stdlib.h>

#define N 10

__global__ void add(int*a, int*b, int*c){
	c[blockIdx.x] = a[blockIdx.x+1] + b[blockIdx.x];
	
	//*c = *a + *b;
}

int main(void){
	int *a, *b, *c;
	int *d_a, *d_b, *d_c;
	int size = N*sizeof(int);

	// Allocate space for device copies of a, b, c
	cudaMalloc((void**) &d_a, size); 
	cudaMalloc((void**) &d_b, size);
	cudaMalloc((void**) &d_c, size);

	// set up input values
	a = (int*) calloc(1,size); 
	b = (int*) calloc(1,size); 
	c = (int*) calloc(1,size);

	for (int i=0; i<N; i++){
		a[i] = i*1;
		b[i] = i*2;
		printf("a[%d]=%d \t b[%d]=%d\n",i,a[i],i,b[i]);
	
	}

	// copy inputs to device
	cudaMemcpy(d_a, a, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, b, size, cudaMemcpyHostToDevice);

	// Launch add kernel on GPU
	add<<<N,1>>>(d_a, d_b, d_c);

	// copy result back to host
	cudaMemcpy(c, d_c, size, cudaMemcpyDeviceToHost);
	
	for (int i=0; i<N; i++){	
		printf("c[%d] = %d\n", i, c[i]);
	}

	// Clean up
	free(a); free(b); free(c);
	cudaFree(d_a); cudaFree(d_b); cudaFree(d_c);

	return 0;

	
}

