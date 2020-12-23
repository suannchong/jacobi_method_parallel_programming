Steps to add your algorithm into the provided skeleton. To reiterate on the homework goal: we are looking at providing different communication layers to solve the same Jacobi implementation. The relax_template.c provides an example of the minimal interface you need to implement in order to have your communication scheme used by the driver. This also allows to easily implement multiple relaxation algorithms using the same interface, and the same supporting software infrastructure.

All the following steps assumes you have the most recent copy of the homework and that you are located in the homework directory in your fork.

### Adding a new algorithm

- check the defines in jacobi.h to find the one corresponding to the target communication scheme. As an example for the OpenSHMEM version, the identifier is RELAXATION\_JACOBI_OSHMEM. This is important as it will allow you to correctly add your algorithm.
- make a copy of relax\_template.c as relax\_bstudent\_comm\_flavor.c (or any other name you prefer), and add it as a dependency to the libhw.a line in the Makefile.
- in this renamed template file replace all ```template``` by a unique naming scheme (```oshmem_jacobi``` would make sense for the OSHMEM version but any another unique name of your liking would work).
- at the bottom of the file change the variable \_relaxation\_template (due to the name change in the previous step) to have the field ```.type``` set to the corresponding define identified in the first step (RELAXATION\_JACOBI_OSHMEM for OpenSHMEM).
- in relaxation.c add the prototype of your algorithm and then add it to the \_relaxation\_classes. Taking the OSHMEM version as an example you should have something similar to:

extern const struct relaxation\_function\_class\_s \_relaxation\_oshmem\_jacobi;

const struct relaxation\_function\_class\_s const* \_relaxation\_classes[] =
    {&\_relaxation\_jacobi, &\_relaxation\_template, ..., &\_relaxation\_oshmem\_jacobi, NULL};

With these changes your algorithm should be now fully integrated in the driver software, and will then be selected at runtime via an argument to the hw_tester.

### Reacting to the configuration file

The elements in the configuration file (the .dat files) are mostly fixed. The tester I provide will not understand anything that is not already there, and will therefore not be able to supply them to your algorithm. If you think there is a lack of support for some awesome feature, propose a patch via a pull request and we'll see. You will need however to change the alg (second line of small.dat) to select your implementation (by setting this n umber to the same decimal value as the define in the first step).

In any case the configuration file is loaded by the ```read_input_file``` function into a well defined structure, ```hw_params_t```. This includes the problem size you are executing as well as the number of threads that will be used (in the ```num_threads```). Your algorithm is expected to take this number into account and start the correct number of threads per process and put them all to work in the resolution of the Jacobi relaxation.

If we are looking at a distributed algorithm where you need to implement the communication support, keep in mind that the number of thread is irrelevant, the communications will only be done in a single thread (MPI_THREAD_SERIALIZED in MPI).

### Initialization/execution/finalization

Preparing the environment for communications should be done only once in the initialization (the ```_init``` function). This will basically identify the process placement in the grid of processes, set the boundaries of the data is will handle, and everything else you might consider necessary for the safe setup of your distributed relaxation algorithm.

The ```_fini``` function does all the opposite in order to leave the software infrastructure in a sane and clean state. This function is to be called once per process, so in a multithreaded scenario you will need to turn off your threads.Unlike in the pthread example, in a distributed environment all processes will call all the functions (including ```_init```, ```_fini```, ```_apply``` and ```_coarsen```). The ```_apply``` function is collective, and will therefore be called by all processes together, giving an opportunity to each one to compute the local relaxation, and the global residual. This global residual will be computed by all processes using a reduction operation, and this residual will be returned by each process ```_apply```. The ```_apply``` function executes a single iteration of the Jacobi relaxation.

### Exposing data

Last but not least the ```_get_data``` is a critical function. It is supposed to return the pointer to the local matrix (in a distributed execution each process will return a pointer to it's local part of the entire data) that contains the most up-to-date version of the local data, i.e. the data as updated by the last iteration. This function is what the tester is using to extract the data from your algorithm and to validate it against a known set of data. It is assumed that this pointer points to the entire local data, including the ghost region around the perimeter of the local matrix tile.

### Avoiding conflicts

If you want to simplify your life and prevent your code from generating conflicts with my version, you should avoid as much as possible directly altering any of the provided files (that's why I suggested to make a copy of the template instead of modifying it directly). This way you will be able to easily pull any update on the COSC462 repository and integrate it seamlessly into your copy.

### Visualizing the output

Check the arguments of the driver application, they allow you to specify the verbosity level and to generate images of the heat propagation at different granularity. Once generated hou can merge these gif images together, either using python or anything else you find, to obtain an animation of the heat propagation in your setup.# jacobi_method_parallel_programming
