#include "ga.h"
#include "macdecls.h"
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "del_api.h"
#include "assert.h"
#include <sys/stat.h>

#ifdef MPI
#include <mpi.h>
#else
#include "sndrcv.h"
#endif

int main(int argc, char ** argv) {

    int my_id, nprocs;
#ifdef MPI
    MPI_Init(&argc, &argv);
#else
    PBEGIN_(argc,argv);
#endif

    GA_Initialize();
    int heap = 1000000000, stack = 1000000000;
    MA_init(C_DBL, stack, heap);

    double sub[4];
    double min[4];

    // Tests the 6 cases: nonzero to zero, zero to nonzero, ed > em, and em <= ed, and ed == em,
    // And sub = min
    sub[0] = 0.0; sub[1] = -0.0003125512; 
    sub[2] = 0.00422134; //sub[3] = 0.0000512;
    //sub[3] = 0.00041;
    sub[3] = 0.000123123123;

    min[0] = -.000612311; min[1] = 0.0;
    min[2] = 0.00123412; //min[3] = 0.00003321;
    //min[3] = 0.000961;
    min[3] = 0.000123123123;

    int ga_min, ga_sub;
 
    int dims[1];
    int chunk[1]; chunk[0] = 2;
    dims[0] = 4;   
    ga_min = NGA_Create(C_DBL, 1, dims, "Min GA", chunk);
    ga_sub = GA_Duplicate(ga_min,"Sub GA");

    int lo[1], hi[1];
    lo[0] = 0; hi[0] = 3;
    NGA_Put(ga_sub, lo,hi, sub, dims); 
    NGA_Put(ga_min, lo,hi, min, dims); 

    char * name = calloc(16,strlen("test.delta"));
    strcpy(name,"test.delta");
    computeDeltaFromMemToDisk(ga_sub, ga_min, 32, name, 0, -3.1);
    int ga_del = recoverFromDiskToMem(ga_sub, 32, name);
    double recoveredVals[4];
    int i = 0;
    NGA_Get(ga_del, lo, hi, recoveredVals, dims); 
    for (i = 0; i < dims[0]; i++) {
        if (my_id == 0) printf("val %d: %lf\n", i, recoveredVals[i]);
    }
       
    GA_Destroy(ga_min);
    GA_Destroy(ga_sub);

    GA_Terminate();
#ifdef MPI
    MPI_Finalize();
#else
    PEND_();
#endif

}
