#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#ifdef MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#endif

#include "utils.h"
#include "dconfObj.h"
#include "run.h"

int main(int argc, char* argv[])
{
	
	char iniFile[MAXLENGTH];
	double t1 = 0., t2 = 0.;
	
	if(argc != 2)
  	{
    	printf("Need an iniFile, buddy !\n");
    	exit(EXIT_FAILURE);
  	}
  	else
  	{
    	strcpy(iniFile, argv[1]);
  	}
  			
		
	t1 = time(NULL);
	int thisRank = 0;
	int size = 1;
    
#ifdef MPI
	MPI_Init(&argc, &argv);
 	MPI_Comm_size(MPI_COMM_WORLD, &size);
  	MPI_Comm_rank(MPI_COMM_WORLD, &thisRank);
#endif
        
        /* READING parameter file */
        dconfObj_t simParam = readDconfObj(iniFile);
        simParam->size = size;
        simParam->thisRank = thisRank;

#ifdef MPI
        if(simParam->memoryIntensive == 0)
          fftw_mpi_init();
#endif
        /* RUN analysis */
        analysis(simParam, thisRank, size);
        
#ifdef MPI
        if(simParam->memoryIntensive == 0)
          fftw_mpi_cleanup();
        
	MPI_Finalize();
#endif
	
	t2 = time(NULL);
	printf("Execution took %e s\n", t2-t1);
	dconfObj_del(&simParam);

	return 0;
}
    	
