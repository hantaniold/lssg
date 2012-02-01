#include <stdio.h>
#include <stdlib.h>
#include <strings.h>


/** filterSmall.c **/
/**Gets rid of all doubles in a file smaller than a certain magnitude **/

int main(int argc, char * argv[]) {
	FILE * fp;
	FILE * outputFile;
	fp = fopen(argv[1],"r");
	outputFile = fopen(argv[2],"w");
    if (argc != 3) {
		printf("Usage: filterSmall <input> <output> \n \
		Filters out doubles from a datafile with value less than 1e-10");
		exit(0);
	} 
    double l;
    double ZERO = 0;
    int retval = 0;
    char * buf = (char *) malloc(8*sizeof(char));
	while ((retval = fscanf(fp, "%8c", buf)) != EOF) {
//		printf("%c%c%c%c%c%c%c%c",buf[0],buf[1],buf[2],buf[3],buf[4],buf[5],buf[6],buf[7]);
        l = * (double *) buf;
		if ((l < 0.0001 && l > 0) || 
			(l > -0.0001 && l < 0)) {
//            printf("a: %lf\n",ZERO);
			fwrite(&ZERO,sizeof(double),1,outputFile);
		} else {
			fwrite(&l,sizeof(double),1,outputFile);
 //           printf("b: %lf\n",l);
		}
	}
	fclose(fp);
	fclose(outputFile);		
	exit(1);
}
