#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <float.h>
#include <getopt.h>

#include "utils.h"


/*****************************************************************************/

void pressEnterToContinue(void)
{
	char c;
	printf("\nPress ENTER to continue ... ");
	fflush(stdin);
	while(1){
		c = getchar();
		if(c == '\n' || c == EOF)
			break;
	}

	return;
}
/*****************************************************************************/


void rTableFile(char filein[], double **table, int nrow, int ncol)
{
	double tmp;
	int col,
		row;
	FILE *mf;

	if((mf = fopen(filein, "r")) == NULL) {
		printf("\n error opening table file\n");
		exit(EXIT_FAILURE);
	}

	for(row = 0; row < nrow; row++) {
		for(col = 0; col < ncol; col++) {
			fscanf(mf, "%lf", &tmp);
			table[row][col] = tmp;
		}
	}
	fclose(mf);

	return;
}
/*****************************************************************************/
