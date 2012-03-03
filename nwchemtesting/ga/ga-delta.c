#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "ga.h"
#include "macdecls.h"
#include "zlib.h"
#include "assert.h"

#ifdef MPI
#include <mpi.h>
#else
#include "sndrcv.h"
#endif

#define EXP_MASK 2047 //0b11111111111
#define FRAC_MASK (-1L ^ ((1L << 63) >> 11)) //0b00000000000011....
#define MSB_MASK (1L << 63)
#define BIAS 1023
#define FRAC_WIDTH 52
#define EXP_WIDTH 11



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
	unsigned long int minSize;
	unsigned long int subSize;
	unsigned long int delSize;
	int thresholdExp;
	int doFilter;
};

struct gctx_t * gctx;
int doCompress(double * ds, unsigned long int size, FILE * outFile, int level);
int doDecompress(char * cs, unsigned long int sizeIn, unsigned long int sizeOut, char * outDs);
int recoverMinuend(struct gctx_t * gctx);
int computeDelta(struct gctx_t *gctx);
void printBin(char * msg, long int x);

int main(int argc, char * argv[]) {
	int my_id, nprocs;
	gctx = calloc(1,sizeof(struct gctx_t));
#ifdef MPI
	MPI_Init(&argc, &argv);
#else 
	PBEGIN_(argc, argv);
#endif

	GA_Initialize();
	int heap = 1000000000, stack = 1000000000;
	if (!MA_init(C_DBL, stack, heap)) GA_Error("MA_init failed",stack+heap);
	
	my_id = GA_Nodeid();
	nprocs = GA_Nnodes();
	if (my_id == 0) {
		printf("\n=== Using %d processes\n\n",nprocs); fflush(stdout);
	}

	/* Parse args and determine whether to recover a minuend, or compute a delta. */
	int doRecover = 0;
	char opt;
	int nrFiles = 0;
	while ((opt = getopt(argc, argv, "m:s:d:o:f:r")) != EOF) {
		switch (opt) {
		case 'f':
			gctx->doFilter = 1;
			sscanf(optarg, "%d", &(gctx->thresholdExp));
			if (my_id == 0) printf("=== Using a cutoff of 10^%d\n",gctx->thresholdExp);
            break;
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
	if (nrFiles != 3 && my_id == 0) {
		printf("Supply either a minuend, subtrahend, and output filename (for delta),\nor a minuend, delta, and output (for recovery)\ne.g. : ./ga-delta.x -s sub -m min -o out\n");
		exit(0);
	}

	if (doRecover) {
		recoverMinuend(gctx);
	} else {
	    int res = computeDelta(gctx);
	}

	if (my_id == 0) 
		printf("=== Terminating.\n");
	GA_Terminate();

#ifdef MPI
	MPI_Finalize();
#else
	PEND_();
#endif
}

/************************************
 ************ COMPUTE DELTA *********
 ************************************/
int computeDelta(struct gctx_t * gctx) {
	int my_id = GA_Nodeid();
	int nrProcs = GA_Nnodes();

/* Find file sizes */
	struct stat * minStat = calloc(1,sizeof(struct stat));
	struct stat * subStat = calloc(1,sizeof(struct stat));
	fstat(fileno(gctx->minFP),minStat);
	fstat(fileno(gctx->subFP),subStat);
	if (my_id == 0) {
		printf("Minuend size: %ld\n",minStat->st_size);
		printf("Subtrahend size: %ld\n",subStat->st_size);
	}

	if (minStat->st_size != subStat->st_size) GA_Error("Error: Minuend and Subtrahend not same size", 1);	
	gctx->minSize = minStat->st_size;
	gctx->subSize = subStat->st_size;
	

/* Create the global arrays */
	int dims[1], chunk[1];
	int ga_min, ga_sub, ga_out;

	dims[0] = gctx->minSize/8;
	chunk[0] = dims[0]/nrProcs - 1;
	ga_min = NGA_Create(C_DBL, 1, dims, "Minuend Array", chunk);
	if (!ga_min) GA_Error("Failure to create minuend array", 1);
	if (my_id == 0) printf("=== Created Minuend Array.\n");
	ga_sub = GA_Duplicate(ga_min,"Subtrahend Array");
	ga_out = GA_Duplicate(ga_min,"Output Array");

	if (!ga_sub) GA_Error("duplicating min -> sub failed.",1);
	if (!ga_out) GA_Error("duplicating min -> out failed.",1);
	if (my_id == 0) printf("=== Created Subtrahend and Output Array.\n");
	
	
/* Read the minuend and subtrahend into local buffers, copy this to global array.
 * Maybe split to two processors o_O */
	if (my_id == 0) {
		int globalStride[1];
		globalStride[0] = dims[0];
		double * minDoubles = calloc(dims[0],sizeof(double));
		
		double * subDoubles = calloc(dims[0],sizeof(double));
		unsigned long int i = 0;
		char * buf1 = (char *) malloc(8);
		char * buf2 = (char *) malloc(8);
		printf("=== Filling min/sub array on processor %d\n",my_id);
		while (fscanf(gctx->minFP,"%8c", buf1) != EOF && fscanf(gctx->subFP, "%8c", buf2) != EOF) {
			minDoubles[i] = *((double *) buf1);
			subDoubles[i] = *((double *) buf2);
			i++;
		}	
		if (i != dims[0]) GA_Error("Did not read as many doubles as expected from subtrahend",1);
		printf("=== Filled min/sub arrays on processor %d\n",my_id);
		int lo[1]; lo[0] = 0;
		int hi[1]; hi[0] = dims[0] - 1;
		NGA_Put(ga_sub, lo, hi, subDoubles, globalStride);
		NGA_Put(ga_min, lo, hi, minDoubles, globalStride);
		free(subDoubles);
		free(minDoubles);
	}


	GA_Sync();

	if (my_id == 0) printf("=== Transferring global doubles to local arrays.\n");
/* Move global doubles to local arrays */ 
	int lo[1], hi[1];
	NGA_Distribution(ga_out, my_id, lo,hi);
	int stride[1]; stride[0] = hi[0] - lo[0] + 1;
	
	double * localMinDoubles = calloc(stride[0],sizeof(double));
	NGA_Get(ga_min, lo, hi, localMinDoubles, stride);
	
	double * localSubDoubles = calloc(stride[0],sizeof(double));
	NGA_Get(ga_sub, lo, hi, localSubDoubles, stride); 

	if (my_id == 0) printf("=== Starting delta computation\n");	
/* Do the delta computation, writing deltas out to a local buffer. */
	unsigned long int j;
	double * localOutDoubles = calloc(stride[0],sizeof(double));

	double m;
	double s;
	double o;
	long int exp_m;
	long int exp_d; //delta 1, "diff"
	short sign_d; //Sign bit of delta
	long int frac_d;
	double d;
	short dExpIsGreater = 0; //is the value "SH_AMT" negative? If so, we need
							//to indicate this when encoding the delta.
	long int m_as_bits;
	long int d_as_bits;
	short SH_AMT;
	long NEG_ZERO = 1L << 63;

	printf("=== Processor %d entering delta computation\n",my_id);	
	for (j = 0; j < stride[0]; j++) {
	
		m = localMinDoubles[j];
		s = localSubDoubles[j];
		//Write a 0 if m and s are equal.
		if (m == s) {
			localOutDoubles[j] = 0;
			continue;
		}
	
		//Write a -0 if s is nonzero but m is 0.
		if (m == 0 && s != 0) {
		    long int negzero = 1L << 63;
			localOutDoubles[j] = * (double *) &negzero;
			continue;
		}
		//Compute a normal subtraction difference.
		d = m - s;

		//If |d| < 10^thresholdExp, round to zero and write that out.
		if (gctx->doFilter && 
				(((d < pow(10,gctx->thresholdExp)) && (d >= 0)) || 
				((d > -1 * pow(10,gctx->thresholdExp)) && d < 0))
			) {
			localOutDoubles[j] = 0;
			continue;
		}

		m_as_bits = * (long int *) &m;
	    d_as_bits = * (long int *) &d;
		frac_d = d_as_bits & FRAC_MASK;
		exp_m =  ((m_as_bits >> FRAC_WIDTH) & EXP_MASK) - BIAS;
	    exp_d = ((d_as_bits >> FRAC_WIDTH) & EXP_MASK) - BIAS;
		exp_d = exp_d - exp_m;
		if (exp_d > 0) { //Was the magnitude of d greater than m?
			dExpIsGreater = 1;
			exp_d *= -1;
		} else {
            dExpIsGreater = 0;
        }
	    SH_AMT = (-1) * exp_d + 2;
		frac_d >>= SH_AMT; //Make SH_AMT leading zeros in the fraction
	
		/* Setting this bit indicates to the delta decoder where to "reshift"
		 * the fraction back, also giving the data of how to successfully 
		 * recover the exponent of (m - s). */
		frac_d = frac_d | (1L << (FRAC_WIDTH -  SH_AMT + 1)); //Set new LS zero to 1.
	
		/* Setting this bit indicates that  the difference was larger in magnitude 
	     * than the minuend, i.e., e(d) > e(m). This is needed in computing the 
		 *minuend given the delta and subtrahend. */
	    if (dExpIsGreater) { 
            frac_d = frac_d | (1L << (FRAC_WIDTH - SH_AMT));
        } else {
            frac_d = frac_d & ~(1L << (FRAC_WIDTH - SH_AMT));
        }
		d_as_bits &= (1L << 63); //zero out everything except the sign bit in delta
		d_as_bits |= ((exp_m + BIAS) << 52); //set delta's exponent bits to m's
		d_as_bits |= (frac_d); //set delta's fraction bits.
		d = * (double *) &d_as_bits;  //turn back into double
		localOutDoubles[j] = d;
	}

	printf("=== Processor %d computed %ld deltas.\n",my_id,j);
	//Combine local buffer into global array.
	if (my_id == 0) printf("=== Moving computed deltas to global array.\n");
	NGA_Put(ga_out, lo, hi, localOutDoubles, stride);
	
	GA_Sync();
	if (my_id == 0) {
		double * outDoubles = calloc(dims[0],sizeof(double));
		lo[0] = 0; hi[0] = dims[0] - 1; stride[0] = dims[0];
		NGA_Get(ga_out, lo, hi, outDoubles, stride);
        if (doCompress(outDoubles,dims[0],gctx->outFP, Z_DEFAULT_COMPRESSION) != Z_OK) 
            printf("=== COMPRESSION ERROR???\n");
	}
	if (my_id == 0) printf("Deltas written to specified output file.\n");

    fclose(gctx->minFP);
    fclose(gctx->outFP);
    fclose(gctx->subFP);
	GA_Destroy(ga_sub); 
	GA_Destroy(ga_min);
	GA_Destroy(ga_out);	

    return 0; //Success!
}


/***************************************
 ********** RECOVER MINUEND ***********
 **************************************
 **************************************/
int recoverMinuend(struct gctx_t *gctx) {
	int my_id = GA_Nodeid();
	int nrProcs = GA_Nnodes();

/* Find file sizes */
	struct stat * delStat = calloc(1,sizeof(struct stat));
	struct stat * subStat = calloc(1,sizeof(struct stat));
	fstat(fileno(gctx->delFP),delStat);
	fstat(fileno(gctx->subFP),subStat);
	if (my_id == 0) {
		printf("Delta size: %ld\n",delStat->st_size);
		printf("Subtrahend size: %ld\n",subStat->st_size);
	}
	gctx->delSize = delStat->st_size;
	gctx->subSize = subStat->st_size;
	

/* Create the global arrays */
	int dims[1], chunk[1];
	int ga_del, ga_sub, ga_out;

	dims[0] = gctx->subSize/8;
	chunk[0] = dims[0]/nrProcs - 1;
	ga_del = NGA_Create(C_DBL, 1, dims, "Delta Array", chunk);
	if (!ga_del) GA_Error("Failure to create minuend array", 1);
	if (my_id == 0) printf("=== Created Minuend Array.\n");
	ga_sub = GA_Duplicate(ga_del,"Subtrahend Array");
	ga_out = GA_Duplicate(ga_del,"Output Array");

	if (!ga_sub) GA_Error("duplicating del -> sub failed.",1);
	if (!ga_out) GA_Error("duplicating del -> out failed.",1);
	if (my_id == 0) printf("=== Created Subtrahend and Output Array.\n");
	
	
/* Read the minuend and subtrahend into local buffers, copy this to global array.
 * Maybe split to two processors o_O */
	if (my_id == 0) {
		int globalStride[1];
		globalStride[0] = dims[0];
        char * compDel = calloc(gctx->delSize, sizeof(char));
        unsigned long int read = fread(compDel, 1, gctx->delSize, gctx->delFP);
        printf("=== Read %ld bytes from compressed delta.\n",read);
		double * delDoubles = calloc(dims[0],sizeof(double));
        int ret = doDecompress(compDel, gctx->delSize, gctx->subSize,(char *) delDoubles);
        if (ret != Z_OK) 
            return -1;
		double * subDoubles = calloc(dims[0],sizeof(double));
		unsigned long int i = 0;
		char * subBuf = (char *) malloc(8);
		printf("=== Filling sub array on processor %d\n",my_id);
		while (fscanf(gctx->subFP, "%8c", subBuf) != EOF) {
			subDoubles[i] = *((double *) subBuf);
			i++;
		}	
		if (i != dims[0]) GA_Error("Did not read as many doubles as expected from subtrahend",1);
		printf("=== Filled del/sub arrays on processor %d\n",my_id);
		int lo[1]; lo[0] = 0;
		int hi[1]; hi[0] = dims[0] - 1;
		NGA_Put(ga_sub, lo, hi, subDoubles, globalStride);
		NGA_Put(ga_del, lo, hi, delDoubles, globalStride);
		free(subDoubles);
		free(delDoubles);
	}


	GA_Sync();

	if (my_id == 0) printf("=== Transferring global doubles to local arrays.\n");
/* Move global doubles to local arrays */ 
	int lo[1], hi[1];
	NGA_Distribution(ga_out, my_id, lo,hi);
	int stride[1]; stride[0] = hi[0] - lo[0] + 1;
	
	double * localDelDoubles = calloc(stride[0],sizeof(double));
	NGA_Get(ga_del, lo, hi, localDelDoubles, stride);
	
	double * localSubDoubles = calloc(stride[0],sizeof(double));
	NGA_Get(ga_sub, lo, hi, localSubDoubles, stride); 

	if (my_id == 0) printf("=== Starting delta computation\n");	

	char subBuf[8]; char delBuf[8];
	double dSub;	double dDel;
	double dMin;
	int shamt;
	long int em;
	long int ed;
	int ed_MINUS_em_isPOS;
	int ed_MINUS_em;
	long int fd;
	long int dDelBits;
	double * localOutDoubles = calloc(stride[0], sizeof(double));

	int j = 0;
	printf("=== Processor %d entering recovery for %d vals\n",my_id,stride[0]);
	for (j = 0; j < stride[0]; j++) {
		dDel = localDelDoubles[j];
		dSub = localSubDoubles[j];
//If the delta is 0 or negative zero deal the right way
    	if (dDel == 0 || dDel == (1L << 63)) {
            long int dDelBits = * (long int *) &(localDelDoubles[j]);
			if (dDelBits && MSB_MASK) {
				localOutDoubles[j] = 0;
				continue;
			} else {
				localOutDoubles[j] = dSub;
				continue;
			}
		}

		dDelBits = *((long int *) &dDel);
	

		dDelBits = dDelBits & FRAC_MASK;
		dDelBits <<= 12;
		shamt = 0;
		while ((dDelBits & MSB_MASK) != MSB_MASK) {
			shamt++;
			dDelBits <<= 1;	
		}
        ed_MINUS_em_isPOS = 0;
        if (((dDelBits << 1) & MSB_MASK) < 0) {
            ed_MINUS_em_isPOS = 1;
        } 
		if (ed_MINUS_em_isPOS) {
			ed_MINUS_em = shamt;
		} else { 
			ed_MINUS_em = -1 * shamt;
		}
		dDelBits = *((long int *) &dDel);
		em = ((dDelBits >> FRAC_WIDTH) & EXP_MASK) - BIAS;
		ed = ed_MINUS_em + em;
		ed += BIAS;
		ed <<= FRAC_WIDTH;
		
		dDelBits = *((long int *) &dDel);
		dDelBits &= FRAC_MASK;
		fd = dDelBits << (shamt + 2);
		fd &= FRAC_MASK; //Eliminates possible 1s that pop up after LS
		
		dDelBits = *((long int *) &dDel);
		dDelBits &= MSB_MASK;
		dDelBits |= fd;
		dDelBits |= ed;
		dMin = dSub + *((double *) &dDelBits);	
		localOutDoubles[j] = dMin;
	}
	
	printf("=== Processor %d recovered %d doubles.\n",my_id,j);
	if (my_id == 0) printf("=== Moving recovered minuend values to global array.\n");
	NGA_Put(ga_out, lo, hi, localOutDoubles, stride);
	GA_Sync();
	if (my_id == 0) {
		double * outDoubles = calloc(dims[0],sizeof(double));
		lo[0] = 0; hi[0] = dims[0] - 1; stride[0] = dims[0];
		NGA_Get(ga_out, lo, hi, outDoubles, stride);
		for (j = 0; j < dims[0]; j++) {
			fwrite(&(outDoubles[j]), sizeof(double), 1, gctx->outFP);
		}
	}
	if (my_id == 0) printf("=== Subtrahend doubles written to specified output file.\n");


	GA_Destroy(ga_del);
	GA_Destroy(ga_out);
	GA_Destroy(ga_sub);
    return 0;
}

void printBin(char * msg, long int x) {

	int i = 63;
	while (i >= 0) {
		printf("%ld",(x>>i) & 1);
		if (i == 63) printf("|");
		if (i == 52) printf("|");
		i--;
	}
	printf(": %s ",msg);
	printf("\n");
}



/* Okay for the time being...very large files we'd need to
 * presumably loop and not do it in one shebang 
 * uses zlib to compress our stream of doubles and write it
 * to the specified output file. */
int doCompress(double * ds,unsigned long int size, FILE * outFile, int level) {
    /* Need to compensate for elts being 8 bytes long*/
    size = sizeof(double) * size;
    int ret, flush;
    unsigned long int have;
    z_stream strm;
    unsigned char * in = (char *) ds; 
    unsigned char * out = calloc(size,sizeof(unsigned char)); 
    
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    ret = deflateInit(&strm, level);
    if (ret != Z_OK)
        return ret; 
    strm.avail_in = size;
    strm.next_in = in;
    strm.avail_out = size;
    strm.next_out = out;
    ret = deflate(&strm, Z_FINISH);
    assert (ret != Z_STREAM_ERROR);
    have = size - strm.avail_out; /* Check space in output buffer */
    int res = fwrite(out, 1, have, outFile); /* Write to disk */
    printf("=== %d bytes written to disk\n",res);

    (void) deflateEnd(&strm);
    return Z_OK;
}


/* "sizeOut" should be the size in bytes of the subtrahend array that 
 *  we're using as a base for recovery. 
 * "sizeIn" is the size of the input delta file. */
int doDecompress(char * cs, unsigned long int sizeIn, unsigned long int sizeOut, char * outDs) {
    int ret;
    unsigned have;
    z_stream strm;
    
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    ret = inflateInit(&strm);
    if (ret != Z_OK) 
        return ret;

    strm.avail_in = sizeIn;
    strm.next_in = cs;
    strm.avail_out = sizeOut;
    strm.next_out = outDs;
    ret = inflate(&strm, Z_NO_FLUSH);

    assert(ret != Z_STREAM_ERROR);
    switch (ret) {
    case Z_NEED_DICT:
        ret = Z_DATA_ERROR;
    case Z_DATA_ERROR:
    case Z_MEM_ERROR:
        (void) inflateEnd(&strm);
        return ret; 
    }
    have = sizeOut - strm.avail_out;
    printf("=== Decompressed size: %u\n",have);
    assert(have == sizeOut);

    /* Cast the decompressed byte array back into doubles */
    (void) inflateEnd(&strm);
    return Z_OK;
}

