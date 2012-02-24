#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ga.h"
#include "macdecls.h"

#ifdef MPI
#include <mpi.h>
#else
#include "sndrcv.h"
#endif


/* Implementation steps:
- Get, through some means, two streams of data. Store these into memory.
- Figure out the size of the stored arrays. 
- The two streams are a minuend and subtrahend stream, or a base and delta stream.
- Break these up into almost-equal sized contiguous chunks.
- Send thechunks off to local processors.
- On the local processors, either recover the minuend, or compute a delta (using my deltaLE.c program)
- Write this recovered data/delta back into a corresponding spot in the global array.
- On the main thread, do whatever with this computed data - store it somewhere, etc. 
*/

struct gctx_t {
	FILE * minFP;
	FILE * subFP;
	FILE * delFP;
	FILE * outFP;
} ;

struct gctx_t * gctx;

void recoverMinuend(gctx_t * gctx);
void computeDelta(gctx_t *gctx);

int main(int argc, char * argv[]) {
	int me, nprocs;
#ifdef MPI
	MPI_Init(&argc, &argv);
#else 
	PBEGIN_(argc, argv);
#endif

	GA_Initialize();
	int heap = 3000000, stack = 3000000;
	if (!MA_init(C_DBL, stack, heap)) GA_Error("MA_init failed",stack+heap);
	
	me = GA_Nodeid();
	nprocs = GA_Nnodes();
	if (me == 0) {
		printf("\nUsing %d processes\n\n",nprocs); fflush(stdout);
	}

	/* Parse args and determine whether to recover a minuend, or compute a delta. */
	int doRecover = 0;
	if (me == 0) {
	gctx = calloc(1, sizeof(gctx_t));
	char opt;
	int nrFiles = 0;
	while ((opt = getopt(argc, argv, "rmsdo")) != EOF) {
		switch (opt) {
		case 'r':
			doRecover++;
			break;
		case 'm':
			gctx->minFP = fopen(optarg,"r");
			nrFiles++;
			break;
		case 's':
			gctx->subFP = fopen(optarg,"r");
			nrFiles++;
			break;
		case 'd':
			gctx->delFP = fopen(optarg,"r");
			nrFiles++;
			break;
		case 'o':
			gctx->outFP = fopen(optarg,"w");
			nrFiles++;
			break;
		}
	}
	if (nrFiles != 3) {
		printf("Supply either a minuend, subtrahend, and output filename (for delta),\nor a minuend, delta, and output (for recovery)\ne.g. : ./ga-delta.x -s sub -m min -o out\n");
		exit(0);
	}
	}

	if (doRecover) {
		recoverMinuend(gctx);
	} else {
		computeDelta(gctx);
	}
	
	//something();

	if (me == 0) 
		printf("\nTerminating.\n");
	GA_Terminate();

#ifdef MPI
	MPI_Finalize();
#else
	PEND_();
#endif
}

void recoverMinuend(struct gctx_t * gctx) {
	// Read the minuend and subtrahend into buffers of doubles.
	
	// Split these into subarrays for the workers.

	// Do the delta computation, writing deltas out to a local buffer.

	//Combine local buffer into global array.

	//Write global array out or something/
}

void computeDelta(struct gctx_t *gctx) {

}
