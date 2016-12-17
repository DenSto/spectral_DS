#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "assert.h"

#define LENGTH 2000

int main(int argc, char* argv[]){
	char buf[LENGTH], *token, *fname;
	char ch;
	int nCols = 0, count = 0, lines = 0;
	int i;
	double t,t0, *ave;
	FILE* in;


	if(argc != 3){
		printf("Usage: %s t0\n",argv[0]);
		return 1; 
	}

	fname = argv[1];
	t0 = atof(argv[2]);

	in = fopen(fname,"r");
	assert(in != NULL);
	
	fgets(buf,LENGTH,in);
	//if(buf[0] == '\0') exit(0);

	while(buf[0] != '\0' && !feof(in)){
		if(buf[0] != '#'){
			lines++;
			if(lines == 1){
				token = strtok(buf," 	");
				while( token != NULL){
					nCols++;
					token = strtok(NULL," 	");
				}
				printf("columns: %d\n",nCols);
			}
		}
		fgets(buf,LENGTH,in);
	}

	rewind(in);

	ave = malloc(sizeof(double)*(nCols-1));
	while(buf[0] != '\0' && !feof(in)){
		if(buf[0] != '#'){
			i=0;
			token = strtok(buf," 	");
			t=atof(token);
			if(t >= t0){
				count++;
				token = strtok(NULL," 	");
				while( token != NULL){
					ave[i] += atof(token);
					token = strtok(NULL," 	");
					i++;
				}
			}
		
		}
		fgets(buf,LENGTH,in);
	}


	printf("%d\n",count);
	for(i = 0; i < nCols-1; i++){
		ave[i]/=((double)count);
		printf("%e    ",ave[i]);
	}
	printf("\n");

	fclose(in);

	return 0;
}
