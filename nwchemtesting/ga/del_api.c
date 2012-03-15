#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "ga.h"
#include "del_api.h"
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

int doCompress(double * ds, unsigned long int size, FILE * outFile, int level);
int doDecompress(char * cs, unsigned long int sizeIn, unsigned long int sizeOut, char * outDs);
int recoverMinuend(struct del_t * gctx);
int computeDelta(struct del_t *gctx);
void printBin(char * msg, long int x);


int computeDeltaFromMemToDisk(int subHandle, int minHandle, unsigned long int size, char * delName, int doFilter, float thresholdExp) {

	int my_id, nprocs;
	struct del_t * gctx = calloc(1,sizeof(struct del_t));
    gctx->subHandle = subHandle;
    gctx->minHandle = minHandle;
    gctx->subSize = size;
    gctx->minSize = size;
    gctx->delName = calloc(1,strlen(delName));
    strcpy(gctx->delName,delName);
    gctx->doFilter = doFilter;
    gctx->thresholdExp = thresholdExp;
    
	
	my_id = GA_Nodeid();
	nprocs = GA_Nnodes();
	if (my_id == 0) {
        printf("====================================================");
		printf("\n=== Using %d processes in computeDeltaFromMemToDisk\n",nprocs); 
        printf("====================================================\n");
	}

	return computeDelta(gctx);
}

int recoverFromDiskToMem(int subHandle, unsigned long int size, char * delName) {
    int my_id = GA_Nodeid();
    int nprocs = GA_Nnodes();
    if (my_id == 0) {
        printf("\n=======================================\n");
        printf("=== Using %d processes in recoverFromDiskToMem\n",nprocs); 
        printf("=======================================\n");
    }
    struct del_t * gctx = calloc(1, sizeof(struct del_t));
    gctx->subHandle = subHandle;
    gctx->subSize = size;
    gctx->delName = calloc(1, strlen(delName));
    strcpy(gctx->delName,delName);
    gctx->delFP = fopen(gctx->delName,"r");
    
    struct stat * delStat = calloc(1, sizeof(struct stat));
    fstat(fileno(gctx->delFP),delStat);
    gctx->delSize =  delStat->st_size;

    int minHandle =  recoverMinuend(gctx);
    if (my_id == 0) {
        printf("=== Recovered min handle is %d\n",minHandle);
    }

    return minHandle;
}

/************************************
 ************ COMPUTE DELTA *********
 ************************************/
int computeDelta(struct del_t * gctx) {
	int my_id = GA_Nodeid();
	int nrProcs = GA_Nnodes();


	int dims[1], chunk[1];
    int ga_min = gctx->minHandle;
    int ga_sub = gctx->subHandle;
    int ga_del;
    
	dims[0] = gctx->minSize/8;
	chunk[0] = dims[0]/nrProcs;
    ga_del = NGA_Create(C_DBL,1, dims, "Delta array", chunk);
	
	GA_Sync();

	if (my_id == 0) printf("=== Transferring global doubles to local arrays.\n");
/* Move global doubles to local arrays */ 
	int lo[1], hi[1];
	NGA_Distribution(ga_del, my_id, lo,hi);
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

	printf("<%d> Processor %d entering delta computation\n",my_id,my_id);	
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

	printf("<%d> Processor %d computed %ld deltas.\n",my_id,my_id,j);
	//Combine local buffer into global array.
	if (my_id == 0) printf("=== Moving computed deltas to global array.\n");
	NGA_Put(ga_del, lo, hi, localOutDoubles, stride);
	
	GA_Sync();
	if (my_id == 0) {
		double * outDoubles = calloc(dims[0],sizeof(double));
		lo[0] = 0; hi[0] = dims[0] - 1; stride[0] = dims[0];
		NGA_Get(ga_del, lo, hi, outDoubles, stride);
        gctx->delFP = fopen(gctx->delName,"w");
        int i;
        if (doCompress(outDoubles,dims[0],gctx->delFP, Z_DEFAULT_COMPRESSION) != Z_OK)  {
            printf("=== COMPRESSION ERROR???\n");
            return -1;
        }
        fclose(gctx->delFP);
	}
	if (my_id == 0) printf("Deltas written to %s.\n",gctx->delName);

	GA_Destroy(ga_del);	

    return 0; //Success!
}


/***************************************
 ********** RECOVER MINUEND ***********
 **************************************
 **************************************/

int recoverMinuend(struct del_t * gctx) {
	int my_id = GA_Nodeid();
	int nrProcs = GA_Nnodes();
    int ga_sub = gctx->subHandle;
    int ga_min, ga_del;
    int dims[1]; dims[0] = gctx->subSize/8;
    int chunk[1]; chunk[0] = dims[0]/nrProcs;
    ga_min = NGA_Create(C_DBL, 1, dims, "Recovered minuend array", chunk);
    if (!ga_min) GA_Error("failed to make output min array",1);
    if (my_id == 0) printf("=== Created output min array\n");
    ga_del = GA_Duplicate(ga_sub,"Delta array");
    if (!ga_del) GA_Error("failed to make input delta array",1);
    if (my_id == 0) printf("=== Created delta array\n");



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
            GA_Error("Trouble decompressing the delta file.\n",1);
		printf("=== Filled del/sub arrays on processor %d\n",my_id);
		int lo[1]; lo[0] = 0;
		int hi[1]; hi[0] = dims[0] - 1;
		NGA_Put(ga_del, lo, hi, delDoubles, globalStride);
		free(delDoubles);
        free(compDel);
	}


	GA_Sync();

	if (my_id == 0) printf("=== Transferring global doubles to local arrays.\n");
/* Move global doubles to local arrays */ 
	int lo[1], hi[1];
	NGA_Distribution(ga_del, my_id, lo,hi);
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
	printf("<%d> Processor %d entering recovery\n",my_id,my_id);
	for (j = 0; j < stride[0]; j++) {
		dDel = localDelDoubles[j];
		dSub = localSubDoubles[j];
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
	
	printf("<%d> Processor %d recovered %d doubles.\n",my_id,my_id,j);
	if (my_id == 0) printf("=== Moving recovered minuend values to global array.\n");
	NGA_Put(ga_min, lo, hi, localOutDoubles, stride);
	GA_Sync();
	GA_Destroy(ga_del);
    return ga_min; //Return handle to the recovered GA
    // If necessary, we could write it to file.
	/*if (my_id == 0) {
		double * outDoubles = calloc(dims[0],sizeof(double));
		lo[0] = 0; hi[0] = dims[0] - 1; stride[0] = dims[0];
		NGA_Get(ga_min, lo, hi, outDoubles, stride);
		for (j = 0; j < dims[0]; j++) {
			fwrite(&(outDoubles[j]), sizeof(double), 1, gctx->outFP);
		}
	}
	if (my_id == 0) printf("=== Subtrahend doubles written to specified output file.\n");*/
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

