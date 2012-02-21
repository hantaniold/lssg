#include <stdio.h>
#include <unistd.h>
#include <pthread.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <math.h>

/** pdeltaLE.c **/
/* Takes in two long-float data files in big endian order,
** and writes their bit-wise deltas to an output file. */
/* Or, takes in a subtrahend, delta file and returns */


/* Constants */

#define EXP_MASK 2047 //0b11111111111
#define FRAC_MASK (-1L ^ ((1L << 63) >> 11)) //0b00000000000011....
#define MSB_MASK (1L << 63)
#define BIAS 1023
#define FRAC_WIDTH 52
#define EXP_WIDTH 11

/* Locks */

pthread_mutex_t statsLock;
pthread_mutex_t gctxLock;
pthread_mutex_t fpMinLock;
pthread_mutex_t fpSubLock;

/* Structs */

struct gctx_t {
	short DO_GN; /* set = use single exponent for all nonzero deltas*/ 
	int GNE; /* group norm exponent */
	short DO_FILTER; /* Should we set a threshold to round to zero?*/
	int THRESHOLD_VALUE;   /* Said threshold */
	int nrThreads;
	int blocksGrabbed; /* How many blocks taken so far..*/
	int blockSize;
};

struct workerArgs {
	int id; /*what thread */
	FILE * fpSub; /* pointer to subtrahend file */
	FILE * fpMin; /* pointer to minuend file */
	FILE * fpDel; /* Pointer to delta file (if exists) */
	FILE * fpOut;
	int blockNr; /* What block we're reading */
};

/* Globals */
struct gctx_t * gctx;
long ZERO = 0;
long NEG_ZERO = (1L << 63);
double POS_SH_AMTS = 0;
double NEG_SH_AMTS = 0;
long int ZERO_SHIFTS = 0;
double POS_COUNT = 0;
double NEG_COUNT = 0;
long int NR_BOTH_ZERO = 0;
long int NR_M_ZERO = 0;
long int NR_S_ZERO = 0;
long int NR_M_POS_WHEN_S_ZERO = 0;
long int NR_M_NEG_WHEN_S_ZERO = 0;
long int NR_M_POS_WHEN_BOTH_NONZERO = 0;
long int NR_M_NEG_WHEN_BOTH_NONZERO = 0;
long int NR_S_POS_WHEN_BOTH_NONZERO = 0;
long int NR_S_NEG_WHEN_BOTH_NONZERO = 0;
long int NR_TIMES_S_AND_M_NONZERO = 0;
int DIFF_BELOW_THRESHOLD = 0;
long int NR_FILTERED = 0;

void * readLoop(void *args);
void * recoverLoop(void * args); 
double gDelta(double s, double m);
void printBin(char * msg, long int x);


int main(int argc, char * argv[]) {
	FILE * fp1;
	FILE * fp2;
	FILE * outputFile;
    if (argc < 4) {
		printf("Delta Usage: delta <SUB> <MIN/DELT> <OUTPUT> [-d] [Flags]\n \
Flags: \n \
\t -gn group-norm-exp \n \
\t -filter filter-threshold-exp \n \
\t -threads numberthreads \n \
\t -b blockSize (KB) \n	If you use -d flag, just provide <SUB> <DELTA> <OUTPUT> flags...");
		exit(0);
    }
	fp1 = fopen(argv[1],"r");
	fp2 = fopen(argv[2],"r");
	outputFile = fopen(argv[3],"w");

	gctx = calloc(1,sizeof(struct gctx_t));
	gctx->DO_FILTER = 0;
	gctx->DO_GN = 0;
	gctx->blockSize = 16;
	gctx->blocksGrabbed = 0;
	gctx->nrThreads = 1;

	int i = 4;
	while (i < argc) {
		if (strcmp(argv[i],"-gn") == 0) {
			gctx->DO_GN = 1;
			sscanf(argv[i+1],"Group norm. enabled - all deltas have exp. %d",&(gctx->GNE));
		}
		if (strcmp(argv[i], "-filter") == 0) {
			gctx->DO_FILTER = 1;
			sscanf(argv[i+1],"%d",&(gctx->THRESHOLD_VALUE));
			printf("Values below 10^%d will be rounded to zero.\n",gctx->THRESHOLD_VALUE);
		}
		if (strcmp(argv[i], "-threads") == 0) {
			sscanf(argv[i+1],"%d",&(gctx->nrThreads));
			printf("Threads: %d",gctx->nrThreads);
		}
		if (strcmp(argv[i], "-b") == 0) {
			sscanf(argv[i+1], "%d", &(gctx->blockSize));
			printf("Block size KB (2^10): %d\n",gctx->blockSize);
		}
		i++;
	}

	pthread_t workerThread;
	for (i = 0; i < gctx->nrThreads; i++) {
		struct workerArgs * wa = calloc(1,sizeof(struct workerArgs));
		wa->id = i;
		wa->fpSub = fp1;	wa->fpMin = fp2;
		wa->fpDel = fp2;	wa->fpOut = outputFile;
		if (strcmp(argv[4],"-d") == 0) { /* Recover? */
			if (pthread_create(&workerThread,NULL,recoverLoop,wa) != 0) {
				perror("Couldn't create worker.");
				exit(0);
			}
		} else {
			if (pthread_create(&workerThread, NULL, readLoop, wa) != 0) {
				perror("Couldn't create worker.");
				free(wa);	free(gctx);	exit(0);
			}
		}
	}
	while(1) {
		sleep(1);
		exit(0);
	}

	/* Below is the infrastructure that should be put into the readLoop */
	char * subBuf = (char *) malloc(8);	char * minBuf = (char *) malloc(8);
	double dSub;						double dMin; 
	double delta;  
	while ((fscanf(fp1, "%8c", subBuf) != EOF) && (fscanf(fp2, "%8c", minBuf) != EOF )) {
		dSub = *((double *) subBuf);	dMin = *((double *) minBuf);
		/* writing 0 to the delta file signifies that dSub = dMin, or 
			both were zero. */
		if ((dSub == 0 && dMin == 0) || (d1 == d2)) { 
			NR_BOTH_ZERO ++;
			fwrite(&ZERO,sizeof(double),1,outputFile); 
		/* writing -0 to the delta file signifies that dSub != 0 but dMin = 0 */
		} else if (dMin == 0) { 
			NR_M_ZERO ++;
			fwrite(&NEG_ZERO,sizeof(double),1,outputFile);
		} else   {
			deltaLF = gDelta(dSub,dMin);
			if (dSub == 0)  { 
				NR_S_ZERO++;
				if (dMin < 0) { NR_M_NEG_WHEN_S_ZERO++; }
				else { NR_M_POS_WHEN_S_ZERO++; }
			} else {
				NR_TIMES_S_AND_M_NONZERO++;
				if (dMin > 0)	{ NR_M_POS_WHEN_BOTH_NONZERO++; }
				else			{ NR_M_NEG_WHEN_BOTH_NONZERO++; }
				if (dSub > 0)	{ NR_S_POS_WHEN_BOTH_NONZERO++; }
				else			{ NR_S_NEG_WHEN_BOTH_NONZERO++; }
			}
			if (DIFF_BELOW_THRESHOLD) {
				DIFF_BELOW_THRESHOLD = 0;
				fwrite(&ZERO,sizeof(double),1,outputFile);
			} else {					
				fwrite(&delta,sizeof(double),1,outputFile);
			}
		}
	}
	printResults(argv);
	fclose(fp1);
	fclose(fp2);
	fclose(outputFile);		
	exit(1);
} 

void * readLoop(void * args) {
	struct workerArgs * wa = (struct workerArgs *) args;
	printf("Hi from %d\n",wa->id);	
	char  h[1];
	while(1) {
		pthread_mutex_lock(&fpMinLock);
		if (fscanf(wa->fpMin,"%c",h) == EOF) {
			pthread_mutex_unlock(&fpMinLock);
			printf("Bye from %d\n",wa->id);
			pthread_exit(NULL);
		}
		printf("%d: %s\n",wa->id,h);
		pthread_mutex_unlock(&fpMinLock);
		sleep(1);
	}
	pthread_exit(NULL);
}

/* "goodDelta" - Given doubles "m"inuend and "s"ubtrahend, compute
the delta which is the sign of (m - s), exponent of m, and fraction of (m - s)'s fraction * 2^(exponent(s) - exponent(m)), divided by two and the e(s) - e(m) + 1'th bit after the decimal set to a 1.  */
double gDelta(double s, double m) {
	long int exp_m;
	long int exp_d; //delta 1, "diff"
	short sign_d; //Sign bit of delta
	long int frac_d;
	double d = m - s;
	if (gctx->DO_FILTER && 
			(((d < pow(10,gctx->THRESHOLD_VALUE)) && (d >= 0)) || 
			((d > -1 * pow(10,gctx->THRESHOLD_VALUE)) && d < 0))
		) {
		DIFF_BELOW_THRESHOLD = 1;
		NR_FILTERED++;
		return d;
	}

	short IS_NEG_SHIFT = 0; //is the value "SH_AMT" negative? If so, we need
							//to indicate this when encoding the delta.
	long int m_as_bits = * (long int *) &m;
    long int d_as_bits = * (long int *) &d;
	frac_d = d_as_bits & FRAC_MASK;
	if (gctx->DO_GN) { 
		exp_m = gctx->GNE;
	} else {
		exp_m =  ((m_as_bits >> FRAC_WIDTH) & EXP_MASK) - BIAS;
	}
    exp_d = ((d_as_bits >> FRAC_WIDTH) & EXP_MASK) - BIAS;
	exp_d = exp_d - exp_m;
	if (exp_d > 0) { //Was the magnitude of d greater than m?
		NEG_SH_AMTS += exp_d;
		exp_d *= -1;
		NEG_COUNT++;	
		IS_NEG_SHIFT = 1;
	} else if (exp_d < 0) {
		POS_SH_AMTS += (-1 * exp_d);
		POS_COUNT++;
	} else {
		ZERO_SHIFTS ++;
	}
	short SH_AMT = (-1) * exp_d + 2;
	frac_d >>= SH_AMT; //Get exp_dMin + 2 leading zeros in the fraction

	/* Setting this bit indicates to the delta decoder where to "reshift"
	 * the fraction back, also giving the data of how to successfully 
	 * recover the exponent of (m - s). */
// .00110010 -> .00000011 (SH_AMT = 4)
// Now FRAC_WIDTH = 8, SH_AMT = 4, so set to: .00100011
	frac_d = frac_d | (1L << (FRAC_WIDTH -  SH_AMT + 1)); //Set new LS zero to 1.

	/* Setting this bit indicates that  the difference was larger in magnitude 
     * than the minuend, i.e., e(d) > e(m). This is needed in computing the 
	 *minuend given the delta and subtrahend. */
// If so, our example becomes .00110011.
	if (IS_NEG_SHIFT) frac_d = frac_d | (1L << (FRAC_WIDTH - SH_AMT));

	d_as_bits &= (1L << 63); //zero out everything except the sign bit in delta
	d_as_bits |= ((exp_m + BIAS) << 52); //set delta's exponent bits to m's
	d_as_bits |= (frac_d); //set delta's fraction bits.
	d = * (double *) &d_as_bits;  //turn back into double
	return d;
}

/*	Do the inverse operation for recovering the minuend from a subtrahend and
	a delta file. If group normalization is enabled, keep that in mind! */
void * recoverLoop(void * args) {
	struct workerArgs * wa = (struct workerArgs *) args;
    FILE * fpSub = wa->fpSub;
	FILE * fpDel = wa->fpDel;
	FILE * fpOut = wa->fpOut;
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
	while (fscanf(fpSub,"%8c",subBuf) != EOF && fscanf(fpDel,"%8c",delBuf) != EOF) {
		dSub = *((double *) subBuf);
		dDel = *((double *) delBuf);
		if (dDel == 0) {
			if (dDel && MSB_MASK) {
				fwrite(&ZERO,sizeof(double),1,fpOut); 
				continue;
			} else {
				fwrite(&dSub,sizeof(double),1,fpOut);
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
		ed_MINUS_em_isPOS = ((dDelBits << 1) & MSB_MASK > 0);
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
		fd = dDelBits << (shamt + 1);
		
		dDelBits = *((long int *) &dDel);
		dDelBits &= MSB_MASK;
		dDelBits |= fd;
		dDelBits |= ed;
		dMin = dSub + *((double *) &dDelBits);	
		fd = fwrite(&dMin,sizeof(double),1,fpOut);
	
	}
	fclose(fpOut);
	pthread_exit(NULL);

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

void printResults(char * argv) { 
	printf("SUBTRAHEND: %s\n",argv[1]);
	printf("MINUEND: %s\n",argv[2]);
	printf("DELTA WRITTEN TO %s\n",argv[3]);
	printf("------------------------\n");
	printf("Times when s and m nonzero: %ld\n",NR_TIMES_S_AND_M_NONZERO);
	printf("    nr m pos when above true: %ld\n",NR_M_POS_WHEN_BOTH_NONZERO);
	printf("    nr m neg when above true: %ld\n",NR_M_NEG_WHEN_BOTH_NONZERO);
	printf("    nr s pos when above true: %ld\n",NR_S_POS_WHEN_BOTH_NONZERO);
	printf("    nr s neg when above true: %ld\n",NR_S_NEG_WHEN_BOTH_NONZERO);
	printf("POS_SH_AMTS/POS_COUNT: %lf\n",POS_SH_AMTS/POS_COUNT);
	printf("POS_COUNT: %ld \nNEG_COUNT: %ld\n",(long int) POS_COUNT,(long int) NEG_COUNT);
	printf("NEG_SH_AMTS/NEG_COUNT: %lf\n",NEG_SH_AMTS/NEG_COUNT);
	printf("AV_SHIFT (neg and pos): %lf\n",(NEG_SH_AMTS + POS_SH_AMTS) / (POS_COUNT + NEG_COUNT));
	printf("ZERO_SHIFTS: %ld\n",ZERO_SHIFTS);
	printf("NR_S_ZERO: %ld\nNR_M_ZERO: %ld\nNR_BOTH_ZERO: %ld\n",NR_S_ZERO,NR_M_ZERO,NR_BOTH_ZERO);
	printf("NR_M_POS_WHEN_S_ZERO: %ld\nNR_M_NEG_WHEN_S_ZERO: %ld \n",NR_M_POS_WHEN_S_ZERO,NR_M_NEG_WHEN_S_ZERO);
	printf("TOTAL DOUBLES PARSED: %ld\n",(long int) (POS_COUNT+NEG_COUNT+ZERO_SHIFTS+NR_M_ZERO+NR_BOTH_ZERO+NR_FILTERED));
	printf("DIFFS ROUNDED TO ZERO: %ld\n",NR_FILTERED);
	return;
}
