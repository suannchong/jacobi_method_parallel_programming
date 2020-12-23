#include <stdio.h>
#include <unistd.h>

__global__ void hello( ){
	printf("Hello from device!\n");
}

int main(void){

	hello<<< 1,1 >>>( );
	printf("Hello from Host!\n");

	sleep(1);

	return 0;
}

