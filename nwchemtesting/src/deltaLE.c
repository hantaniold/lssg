#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <math.h>

/** delta.c **/
/* Takes in two long-float data files in big endian order,
** and writes their bit-wise deltas to an output file. */

//add pthreading

#define EXP_MASK 2047 //0b11111111111
#define FRAC_MASK (-1L ^ ((1L << 63) >> 11)) //0b00000000000011....
#define BIAS 1023
#define FRAC_WIDTH 52
#define EXP_WIDTH 11

long ZERO = 0;
long NEG_ZERO = 1L << 63;
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
short DO_GN = 0; /* set = use single exponent for all nonzero deltas*/ 
int GNE = 0; /* group norm exponent */ 
short DO_FILTER = 0; /* Should we set a threshold to round to zero?*/
int THRESHOLD_VALUE = 0;   /* Said threshold */
int DIFF_BELOW_THRESHOLD = 0;
long int NR_FILTERED = 0;

double gDelta(double s, double m);
void printBin(char * msg, long int x);


int main(int argc, char * argv[]) {
	FILE * fp1;
	FILE * fp2;
	FILE * outputFile;
    if (argc < 4) {
		printf("Usage: delta <SUBTRAHEND> <MINUEND> <OUTPUT> <TYPE>\n \
		Types include: -ldf (Perform subtraction as longs but store as double) \n  -lf (perform subtraction as double and store as double) \n -xor (xor the arguments) \n");
		exit(0);
    }
	fp1 = fopen(argv[1],"r");
	fp2 = fopen(argv[2],"r");
	outputFile = fopen(argv[3],"w");

	char * deltaType;
	int i = 4;
	while (i < argc) {
		if (strcmp(argv[i],"-t") == 0)  {
			deltaType = argv[i+1];
		}
		if (strcmp(argv[i],"-gn") == 0) {
			DO_GN = 1;
			sscanf(argv[i+1],"%d",&GNE);
		}
		if (strcmp(argv[i], "-filter") == 0) {
			DO_FILTER = 1;
			sscanf(argv[i+1],"%d",&THRESHOLD_VALUE);
			printf("%d\n",THRESHOLD_VALUE);
	
		}
		i++;
		
	}

	char * buf1 = (char *) malloc(8);  //Holds 8 bytes of input file
	char * buf2 = (char *) malloc(8);
	double d1; //The 8 bytes of buf1 interpreted as a double
    double d2; 
	double deltaLF;  //A delta computed between two numbers, in some fashion
	while ((fscanf(fp1, "%8c", buf1) != EOF) && (fscanf(fp2, "%8c", buf2) != EOF )) {
			d1 = *((double *) buf1);
			d2 = *((double *) buf2);
		/* writing 0 to the delta file signifies that d1 = d2, or 
			both were zero. */
			if ((d1 == 0 && d2 == 0) || (d1 == d2)) { 
				NR_BOTH_ZERO ++;
				fwrite(&ZERO,sizeof(double),1,outputFile); 
		/* writing -0 to the delta file signifies that d1 != 0 but d2 = 0 */
			} else if (d2 == 0) { 
				NR_M_ZERO ++;
				fwrite(&NEG_ZERO,sizeof(double),1,outputFile);
			} else   {
				deltaLF = gDelta(d1,d2);
				if (d1 == 0)  { 
					NR_S_ZERO++;
					if (d2 < 0) { NR_M_NEG_WHEN_S_ZERO++; }
					else { NR_M_POS_WHEN_S_ZERO++; }
				} else {
					NR_TIMES_S_AND_M_NONZERO++;
					if (d2 > 0) {
						NR_M_POS_WHEN_BOTH_NONZERO++; }
					else {
	 					NR_M_NEG_WHEN_BOTH_NONZERO++; }
					if (d1 > 0) {
						NR_S_POS_WHEN_BOTH_NONZERO++; }
					else {
	 					NR_S_NEG_WHEN_BOTH_NONZERO++; }

				}
				if (DIFF_BELOW_THRESHOLD) {
					DIFF_BELOW_THRESHOLD = 0;
					fwrite(&ZERO,sizeof(double),1,outputFile);
				} else {					
					fwrite(&deltaLF,sizeof(double),1,outputFile);
				}
			}
			
		}
	
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
	fclose(fp1);
	fclose(fp2);
	fclose(outputFile);		
	exit(1);
} 

/* "goodDelta" - Given doubles "m"inuend and "s"ubtrahend, compute
the delta which is the sign of (m - s), exponent of m, and fraction of (m - s)'s fraction * 2^(exponent(s) - exponent(m)), divided by two and the e(s) - e(m) + 1'th bit after the decimal set to a 1.  */
double gDelta(double s, double m) {
	long int exp_m;
	long int exp_d; //delta 1, "diff"
	short sign_d; //Sign bit of delta
	long int frac_d;
	double d = m - s;
	if (DO_FILTER && 
			(((d < pow(10,THRESHOLD_VALUE)) && (d >= 0)) || 
			((d > -1 * pow(10,THRESHOLD_VALUE)) && d < 0))
		) {
		DIFF_BELOW_THRESHOLD = 1;
		NR_FILTERED++;
		return d;
	}

	short IS_NEG_SHIFT = 0; //is the value "SH_AMT" negative? If so, we need
							//to indicate this when encoding the delta.
//	printf("m: %lf s: %lf d: %lf \n",m,s,d);
	long int m_as_bits = * (long int *) &m;
    long int d_as_bits = * (long int *) &d;
//	printBin("d_as_bits",d_as_bits);
	frac_d = d_as_bits & FRAC_MASK;
//	printBin("frac_d",frac_d);
	if (DO_GN) { 
		exp_m = GNE;
	} else {
		exp_m =  ((m_as_bits >> FRAC_WIDTH) & EXP_MASK) - BIAS;
	}
    exp_d = ((d_as_bits >> FRAC_WIDTH) & EXP_MASK) - BIAS;
//	printf("exp_m: %ld, exp_d: %ld,",exp_m,exp_d);
	exp_d = exp_d - exp_m;
//	printf(" exp_d - exp_m: %ld,",exp_d);
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
	frac_d >>= SH_AMT; //Get exp_d2 + 1 leading zeros in the fraction
//	printf(" SH_AMT: %d\n",SH_AMT);
//	printBin("frac_d >> SH_AMT",frac_d);

	/* Setting this bit indicates to the delta decoder where to "reshift"
	 * the fraction back, also giving the data of how to successfully 
	 * recover the exponent of (m - s). */
	frac_d = frac_d | (1L << (FRAC_WIDTH -  SH_AMT + 1)); //Set new LS zero to 1.

	/* Setting this bit indicates that  the difference was larger in magnitude 
     * than the minuend, i.e., e(d) > e(m). This is needed in computing the 
	 *minuend given the delta and subtrahend. */
	if (IS_NEG_SHIFT) frac_d = frac_d | (1L << (FRAC_WIDTH - SH_AMT));

//	printBin("frac_d with 1",frac_d);
	d_as_bits &= (1L << 63); //zero out everything except the sign bit in delta
	d_as_bits |= ((exp_m + BIAS) << 52); //set delta's exponent bits to m's
	d_as_bits |= (frac_d); //set delta's fraction bits.
	d = * (double *) &d_as_bits;  //turn back into double
//	printBin("exp_m",(exp_m  + BIAS) << 52);
//	printBin("d_as_bits final",d_as_bits);
	return d;
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
