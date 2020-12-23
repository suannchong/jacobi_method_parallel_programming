CC = mpicc
LD = mpicc
NVCC = nvcc
CFLAGS += -g3 -std=c11 -D_GNU_SOURCE -fopenmp 
LDFLAGS += -g -lm -lgomp -lcuda
NVFLAGS = -g -std=c++11 -D_GNU_SOURCE
#CFLAGS = -Wall -O3 -DNDEBUG
#LDFLAGS = -O3
AR = ar

.PHONY: all

all: hw_tester

fast: CFLAGS = -std=c11 -D_GNU_SOURCE -O3 
fast: clean hw_tester

lib: libhw.a

hw_tester: % : driver.o libhw.a libcudart.so
	$(LD) -o hw_tester $^ $(LDFLAGS)

libhw.a: relaxation.o relax_jacobi.o relax_template.o misc.o relax_jacobi_mpi_blocking.o relax_jacobi_mpi_nonblocking.o relax_jacobi_mpi_collective.o relax_jacobi_mpi_rma.o relax_jacobi_openmp.o relax_jacobi_cuda.o
	$(AR) cru libhw.a $^

relax_jacobi_cuda.o: relax_jacobi_cuda.cu jacobi.h header.h Makefile
	$(NVCC) -c relax_jacobi_cuda.cu -o $@ $(NVFLAGS)

%.o: %.c jacobi.h header.h Makefile
	$(CC) -c $< -o $@ $(CFLAGS)

clean:
	rm -rf *.o *.a hw_tester
