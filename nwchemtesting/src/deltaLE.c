#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

/** delta.c **/
/* Takes in two long-float data files in big endian order,
** and writes their bit-wise deltas to an output file. */

int main(int argc, char * argv[]) {
	FILE * fp1;
	FILE * fp2;
	FILE * outputFile;
	fp1 = fopen(argv[1],"r");
	fp2 = fopen(argv[2],"r");
	outputFile = fopen(argv[3],"w");
	char * buf1 = (char *) malloc(8);  //Holds 8 bytes of input file
	char * buf2 = (char *) malloc(8);
	double a; //The 8 bytes of buf1 interpreted as a long.
    double b; 
	double delta;  //l2 - l1. Written to file.
	while ((fscanf(fp1, "%8c", buf1) != EOF) && (fscanf(fp2, "%8c", buf2) != EOF )) {
		a = *((double *) buf1);
		b = *((double *) buf2);
		delta = b - a;
    //    printf("1: %lf \t 2: %lf \n d: %lf \n",a,b,delta);
		fwrite(&delta,sizeof(double),1,outputFile);
	}
	fclose(fp1);
	fclose(fp2);
	fclose(outputFile);		
	exit(1);
}
