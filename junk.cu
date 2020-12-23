#include <stdlib.h>
#include <stdio.h>

#define SIZE 1365

__global__ void func(int* a, int s)
{
	    int i = (blockIdx.x * 1024) + threadIdx.x;
		    if(i > s)
				        return;
						    a[i] = i;
							    return;
}

int main(int argc, char** argv)
{
int* a;
int i, s, s2;

s = SIZE;
if(s > 1024)
s2 = 1024;
else
s2 = s;

cudaMallocManaged(&a, sizeof(int) * s);

func<<<(s / 1024) + 1, s2>>>(a, s);
cudaDeviceSynchronize();

for(i = 0; i < s; i++)
printf("%d: %d\n", i, a[i]);
return(0);
}
