#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>

/** delta.c **/
/* Takes in two long-float data files in big endian order,
** and writes their bit-wise deltas to an output file. */


#define EXP_MASK 2047 //0b11111111111
#define FRAC_MASK (-1L ^ ((1L << 63) >> 11)) //0b00000000000011....
#define BIAS 1023
#define FRAC_WIDTH 52
#define EXP_WIDTH 11

double gDelta(double s, double m);
int main(int argc, char * argv[]) {
	FILE * fp1;
	FILE * fp2;
	FILE * outputFile;
    if (argc != 5) {
		printf("Usage: delta <SUBTRAHEND> <MINUEND> <OUTPUT> <TYPE>\n \
		Types include: -ldf (Perform subtraction as longs but store as double) \n  -lf (perform subtraction as double and store as double) \n -xor (xor the arguments) \n");
		exit(0);
    }
	fp1 = fopen(argv[1],"r");
	fp2 = fopen(argv[2],"r");
	outputFile = fopen(argv[3],"w");
    char * deltaType = argv[4];
	char * buf1 = (char *) malloc(8);  //Holds 8 bytes of input file
	char * buf2 = (char *) malloc(8);
	double d1; //The 8 bytes of buf1 interpreted as a double
    double d2; 
	double deltaLF;  //A delta computed between two numbers, in some fashion
    long l1;
    long l2;
    long deltaLDF;
    long deltaXOR;
	while ((fscanf(fp1, "%8c", buf1) != EOF) && (fscanf(fp2, "%8c", buf2) != EOF )) {
        if (strcmp(deltaType,"-lf") == 0) {
			d1 = *((double *) buf1);
			d2 = *((double *) buf2);
			deltaLF = d2 - d1;
		    fwrite(&deltaLF,sizeof(double),1,outputFile);
      //    printf("1: %lf \t 2: %lf \n d: %lf \n",a,b,delta);
        } else if (strcmp(deltaType,"-ldf") == 0) {
            l1 = *((long *) buf1);
			l2 = *((long *) buf2);
			deltaLDF = (double) l2 - l1;
			fwrite(&deltaLDF,sizeof(long),1,outputFile);
//			printf("1: %ld \t 2: %ld \n d: %ld \n",l1,l2,deltaLD);
        } else if (strcmp(deltaType,"-xor") == 0) {
            l1 = *((long *) buf1);
			l2 = *((long *) buf2);
			deltaXOR =  l2 ^ l1;
			fwrite(&deltaXOR,sizeof(long),1,outputFile);
		} else {
			d1 = *((double *) buf1);
			d2 = *((double *) buf2);
			if (d1 == 0) { 
				fwrite(&d1,sizeof(double),1,outputFile); 
			} else {
				deltaLF = gDelta(d1,d2);
				fwrite(&deltaLF,sizeof(double),1,outputFile);
			}
			
		}
	}
	fclose(fp1);
	fclose(fp2);
	fclose(outputFile);		
	exit(1);
} 


void printBin(char * msg, long int x) {

	int i = 63;
	while (i > 0) {
		printf("%ld",(x>>i) & 1);
		if (i == 63) printf("|");
		if (i == 52) printf("|");
		i--;
	}
	printf(": %s ",msg);
	printf("\n");
}
/* "goodDelta" - Given doubles "m"inuend and "s"ubtrahend, compute
the delta which is the sign of (m - s), exponent of m, and fraction of (m - s)'s fraction * 2^(exponent(s) - exponent(m)), divided by two and the e(s) - e(m) + 1'th bit after the decimal set to a 1.  */
double gDelta(double s, double m) {
// 0.045 - 0.05
	long int exp_m;
	long int exp_d; //delta 1, "diff"
	short sign_d; //Sign bit of delta
	long int frac_d;
	double d = m - s;
//	printf("m: %lf s: %lf d: %lf \n",m,s,d);
	long int m_as_bits = * (long int *) &m;
    long int d_as_bits = * (long int *) &d;
//	printBin("d_as_bits",d_as_bits);
	frac_d = d_as_bits & FRAC_MASK;
//	printBin("frac_d",frac_d);
	exp_m =  ((m_as_bits >> FRAC_WIDTH) & EXP_MASK) - BIAS;
    exp_d = ((d_as_bits >> FRAC_WIDTH) & EXP_MASK) - BIAS;
//	printf("exp_m: %ld, exp_d: %ld,",exp_m,exp_d);
	exp_d = exp_d - exp_m;
//	printf(" exp_d - exp_m: %ld,",exp_d);
	/* Implement way to deal with exp_d2 >= 0. this only does for < 0 */
    short SH_AMT = (-1) * exp_d + 1;
	frac_d >>= SH_AMT; //Get exp_d2 + 1 leading zeros in the fraction
//	printf(" SH_AMT: %d\n",SH_AMT);
//	printBin("frac_d >> SH_AMT",frac_d);
	frac_d = frac_d | (1L << (FRAC_WIDTH -  SH_AMT)); //Set new LS zero to 1.
//	printBin("frac_d with 1",frac_d);
	d_as_bits &= (1L << 63); //zero out everything except the sign bit in delta
	d_as_bits |= ((exp_m + BIAS) << 52); //set delta's exponent bits to m's
	d_as_bits |= (frac_d); //set delta's fraction bits.
	d = * (double *) &d_as_bits;  //turn back into double
//	printBin("exp_m",(exp_m  + BIAS) << 52);
//	printBin("d_as_bits final",d_as_bits);
	return d;
}

