#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>

/** delta.c **/
/* Takes in two long-float data files in big endian order,
** and writes their bit-wise deltas to an output file. */



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
			deltaLF = gDelta(d1,d2);
			fwrite(&deltaLF,sizeof(double),1,outputFile);
		}
	}
	fclose(fp1);
	fclose(fp2);
	fclose(outputFile);		
	exit(1);
} 

/* "goodDelta" - Given doubles "m"inuend and "s"ubtrahend, compute
the delta which is the sign of (m - s), exponent of m, and fraction of (m - s)'s fraction * 2^(exponent(s) - exponent(m)), divided by two and the e(s) - e(m) + 1'th bit after the decimal set to a 1.  */
double gDelta(double s, double m) {
	long int exp_m;
	long int exp_d1; //delta 1, "diff"
    short exp_d2; //delta 2, "diff'"
	double frac_d1;
	double frac_d2;
	double d1 = m - s;
	short EXP_MASK = 2047; //0b11111111111
    short BIAS = 1023;
	long int m_as_bits = * (long int *) &m;
    long int d1_as_bits = * (long int *) &d1;
	exp_m =  ((m_as_bits >> 52) & EXP_MASK) - BIAS;
    exp_d1 = ((d1_as_bits >> 52) & EXP_MASK) - BIAS;

    printf("m: %ld d1: %ld\n",exp_m,exp_d1);
	return 2.0;
}
