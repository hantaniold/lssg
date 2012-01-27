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
	char * lBuf1 = (char *) malloc(8); //Those bytes flipped
	char * buf2 = (char *) malloc(8);
	char * lBuf2 = (char *) malloc(8);
	double l1; //The 8 bytes of buf1 interpreted as a long.
    double l2; 
	double d;  //l2 - l1. Written to file.
	while ((fscanf(fp1, "%8c", buf1) != EOF) && (fscanf(fp2, "%8c", buf2) != EOF )) {
		for (int i = 0; i < 8; i++) {
			lBuf1[i] = buf1[7-i];
			lBuf2[i] = buf2[7-i];
		//	printf("%c %c\n",lBuf1[i],lBuf2[i]);
		}

		l1 = *((double *) lBuf1);
		l2 = *((double *) lBuf2);
		d = l2 - l1;

        printf("1: %lf \t 2: %lf \n d: %lf \n",l1,l2,d);
	//	printf("%ld\n",l1);
	//	printf("%ld\n",l2);
	//	printf("%ld\n",l2 - l1);
	//	fwrite(&d,sizeof(long),1,outputFile);
		fwrite(&d,sizeof(double),1,outputFile);

	}
		
	fclose(fp1);
	fclose(fp2);
	fclose(outputFile);		
	exit(1);
}
