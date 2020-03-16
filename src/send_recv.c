#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "utils.h"

#ifdef __MPI
void send_recv_galaxies_int(domain_t *thisDomain, int *numGalOnRanks, int *toSend, int *recvNumGal, int *recvNumGalOnRanks, int **toRecv)
{
    MPI_Status status;
    int size = thisDomain->size;
    int thisRank = thisDomain->thisRank;
    
    int head_offset = 0;
    int tmpValueCounter = 0;
    int *toRecvLoc = NULL;
    
    toRecvLoc = *toRecv;
    
    for(int head=0; head<size; head++)
    {        
        head_offset = 0;
        for(int destRank=0; destRank<size; destRank++)
        {
            if(destRank != head)
            {
                if(thisRank == head)
                {
                    MPI_Send(&numGalOnRanks[destRank], 1, MPI_INT, destRank, 100, MPI_COMM_WORLD);
                    if(numGalOnRanks[destRank] != 0)
                    {
                        MPI_Send(&toSend[head_offset], numGalOnRanks[destRank], MPI_INT, destRank, 101, MPI_COMM_WORLD);
                    }
                }
            }
            
            if(destRank == head && destRank == thisRank)
            {            
                recvNumGalOnRanks[head] = numGalOnRanks[head];
                tmpValueCounter += recvNumGalOnRanks[head];
                toRecvLoc = realloc(toRecvLoc, tmpValueCounter*sizeof(int));
                for(int i=0; i<recvNumGalOnRanks[head]; i++)
                {
                    toRecvLoc[tmpValueCounter - recvNumGalOnRanks[head] + i] = toSend[head_offset+i];
                }
            }
            head_offset += numGalOnRanks[destRank];
        }

        
        if(thisRank != head)
        {
            MPI_Recv(&recvNumGalOnRanks[head], 1, MPI_INT, head, 100, MPI_COMM_WORLD, &status);
            
            tmpValueCounter += recvNumGalOnRanks[head];
            toRecvLoc = realloc(toRecvLoc, tmpValueCounter*sizeof(int));
                
            if(recvNumGalOnRanks[head] != 0)
            {
                MPI_Recv(&toRecvLoc[tmpValueCounter-recvNumGalOnRanks[head]], recvNumGalOnRanks[head], MPI_INT, head, 101, MPI_COMM_WORLD, &status);
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    *recvNumGal = tmpValueCounter;
    *toRecv = toRecvLoc;
}
#endif

/* toSend needs to be arrays in nion, and provide the list to whcih ranks galaxies are sent to as an input */
#ifdef __MPI
void communicate_array_int(int num_gal_on_rank, int *sourceRank, double *zsource, int **thisArray, domain_t *thisDomain, int include_redshift_bins)
{
    int *numGalOnRanks, *recvNumGalOnRanks;
    int *toSend;
    int *toRecv;
    int recvNumGal;

    toSend = NULL;
    toRecv = NULL;
        
    MPI_Barrier(MPI_COMM_WORLD);
            
    /* PACK GALAXIES TO SENd & BROADCAST GALAXIES */
    send_recv_galaxies_int(thisDomain, numGalOnRanks, toSend, &recvNumGal, recvNumGalOnRanks, &toRecv);

    
    /* DELETE CURRENT GALAXIES */
    free(*thisArray);
    *thisArray = NULL;
    *thisArray = toRecv;
    free(toSend);
    
    if(numGalOnRanks != NULL) free(numGalOnRanks);
    if(recvNumGalOnRanks != NULL) free(recvNumGalOnRanks);
}
#endif