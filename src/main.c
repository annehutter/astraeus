#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <complex.h>

#ifdef MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "cifog/confObj.h"
#include "cifog/grid.h"

#include "utils.h"
#include "dconfObj.h"
#include "gal_gtree.h"
#include "outgal.h"
#include "nion.h"
#include "domain.h"
#include "run.h"
#include "output.h"

int main(int argc, char *argv[])
{
  int thisRank = 0;
  int size = 1;
  double t1 = 0., t2 = 0.;
  
  char iniFile[MAXLENGTH];
  
#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &thisRank);
  
  t1 = MPI_Wtime();
#else
  t1 = time(NULL);
#endif

  if(argc != 2)
  {
    printf("walking_trees: (C) Use at own risk...\nUSAGE: walking_trees <FILE>\n");
  }
  else
  {
    strcpy(iniFile, argv[1]);
  }
    
  dconfObj_t simParam = readDconfObj(iniFile);
   
  if(thisRank==0)
    copy_iniFile_executable(simParam->outFileName);

  run_astraeus(simParam, thisRank, size);
  
  dconfObj_del(&simParam);
  
#ifdef MPI
  t2 = MPI_Wtime();
  printf("rank %d: Execution took %e s\n", thisRank, t2-t1);
  
  MPI_Finalize();
#else
  t2 = time(NULL);
  printf("Execution took %e s\n", t2-t1);
#endif
}
