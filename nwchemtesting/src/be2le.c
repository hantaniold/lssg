#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

/** be2le.c **/
/* Takes in one 64-bit IEEE FP  data file in big endian order,
** and writes it to a specified file in little endian orer. */

int main(int argc, char * argv[]) {
	FILE * fp;
	FILE * outputFile;
	fp = fopen(argv[1],"r");
	outputFile = fopen(argv[2],"w");
	char * buf = (char *) malloc(8);  //Holds 8 bytes of input file
	char * lBuf = (char *) malloc(8); //Those bytes flipped
	long l; //The 8 bytes of buf1 interpreted as a long.
	while (fscanf(fp, "%8c", buf) != EOF) {
		for (int i = 0; i < 8; i++) {
			lBuf[i] = buf[7-i];
		}
		l = *((long *) lBuf);
		fwrite(&l,sizeof(long),1,outputFile);
	}
	fclose(fp);
	fclose(outputFile);		
	exit(1);
}
