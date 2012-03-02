#include <stdlib.h>
#include <stdio.h>



int main(int argc, char ** argv) {

double * d = calloc(1,sizeof(double));
int * i = calloc(1,sizeof(int));

if (argc != 3) {
	printf("Usage: ./ga2bin input output\n");
}
FILE * fp = fopen(argv[1],"r");
FILE * fpOut = fopen(argv[2],"w");

while (fscanf(fp, "%d %lf",i,d) != EOF) {
	fwrite(d,sizeof(double),1,fpOut);
}
exit(0);
}
