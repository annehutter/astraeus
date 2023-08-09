#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "utils.h"
#include "dconfObj.h"
#include "gal_gtree.h"
#include "outgal.h"
#include "output.h"

/*-------------------------------------------------------*/
/* COPYING TO OUTPUT STRUCTS                             */
/*-------------------------------------------------------*/

/* VERTICAL OUPUT */
void copy_tree_to_outgtree(gtree_t *thisTree, outgtree_t *thisOutGTree)
{
    if(thisTree->numGal != thisOutGTree->numGal)
    {
        printf("Tree and OutTree don't have the same numGal, can't do much here !\n");
        exit(EXIT_FAILURE);
    }
    
    for(int i=0; i<thisTree->numGal; i++)
        copy_gal_to_outgaltree(&(thisTree->galaxies[i]), &(thisOutGTree->outgalaxies[i]));
}

void copy_gal_to_outgaltree(gal_t *thisGal, outgal_t *thisOutGal)
{
#if defined WITHROCKSTARID
    thisOutGal->ID = (int64_t) (thisGal->ID);
    thisOutGal->descID = (int64_t) (thisGal->descID);
#endif
  
    thisOutGal->localID = (int32_t) (thisGal->localID);
    thisOutGal->localDescID = (int32_t) (thisGal->localDescID);
    thisOutGal->numProg = (int32_t) (thisGal->numProg);
    thisOutGal->snapnumber = (int32_t) (thisGal->snapnumber);
    
    thisOutGal->scalefactor = (float) (thisGal->scalefactor);
    for(int i=0; i<3; i++)
    {
        thisOutGal->pos[i] = (float) (thisGal->pos[i]);
        thisOutGal->vel[i] = (float) (thisGal->vel[i]);
    }

    thisOutGal->Mvir = (float) (thisGal->Mvir);
    thisOutGal->Mvir_prog = (float) (thisGal->Mvir_prog);
    thisOutGal->Rvir = (float) (thisGal->Rvir);
    thisOutGal->velDisp = (float) (thisGal->velDisp);
    thisOutGal->velMax = (float) (thisGal->velMax);
    thisOutGal->spin = (float) (thisGal->spin);
    thisOutGal->scalefactorLastMajorMerger = (float) (thisGal->scalefactorLastMajorMerger);
    
    thisOutGal->MgasIni = (float) (thisGal->MgasIni);
#ifndef FIRST
    thisOutGal->fracMgasMer = (float) (thisGal->fracMgasMer);
#endif
#if defined WITHMETALS
    thisOutGal->MgasNew = (float) (thisGal->MgasNew);
    thisOutGal->MgasEj = (float) (thisGal->MgasEj);
#endif
    thisOutGal->Mgas = (float) (thisGal->Mgas);
    thisOutGal->Mstar = (float) (thisGal->Mstar);
    
#ifndef FIRST
    thisOutGal->fesc = (float) (thisGal->fesc);
    thisOutGal->Nion = (float) (thisGal->Nion);
    thisOutGal->fej = (float) (thisGal->fej);
#endif
    thisOutGal->feff = (float) (thisGal->feff);
    thisOutGal->fg = (float) (thisGal->fg);
    thisOutGal->zreion = (float) (thisGal->zreion);
    thisOutGal->photHI_bg = (float) (thisGal->photHI_bg);

#if defined WITHMETALS
    for(int i=0; i<3; i++)
    {
        thisOutGal->Mmetal[i] = (float) (thisGal->Mmetal[i]);
        thisOutGal->MmetalIni[i] = (float) (thisGal->MmetalIni[i]);
        thisOutGal->MmetalNew[i] = (float) (thisGal->MmetalNew[i]);
        thisOutGal->fracMmetalMer[i] = (float) (thisGal->fracMmetalMer[i]);
        thisOutGal->MmetalEj[i] = (float) (thisGal->MmetalEj[i]); 
    }
    
    thisOutGal->Mdust = (float) (thisGal->Mdust);
#endif
        
    if(1./thisGal->scalefactor - 1. > thisGal->zreion && thisGal->zreion > 0.)
      printf("SOMETHING WRONG (thisGal): z = %e\t zreion = %e\t numProg = %d\n", 1./thisGal->scalefactor - 1., thisGal->zreion, thisGal->numProg);
    if(1./thisOutGal->scalefactor - 1. > thisOutGal->zreion && thisOutGal->zreion > 0.)
      printf("SOMETHING WRONG (thisOutGal): z = %e\t zreion = %e\n", 1./thisOutGal->scalefactor - 1., thisOutGal->zreion);
}

/* HORIZONTAL OUPUT */
void copy_gal_to_outgalsnap(gal_t *thisGal, outgalsnap_t *thisOutGal)
{
    thisOutGal->scalefactor = (float) (thisGal->scalefactor);
    for(int i=0; i<3; i++)
    {
        thisOutGal->pos[i] = (float) (thisGal->pos[i]);
        thisOutGal->vel[i] = (float) (thisGal->vel[i]);
    }

    thisOutGal->Mvir = (float) (thisGal->Mvir);
    thisOutGal->Mvir_prog = (float) (thisGal->Mvir_prog);
    thisOutGal->Rvir = (float) (thisGal->Rvir);
    thisOutGal->velDisp = (float) (thisGal->velDisp);
    thisOutGal->velMax = (float) (thisGal->velMax);
    thisOutGal->spin = (float) (thisGal->spin);
    thisOutGal->scalefactorLastMajorMerger = (float) (thisGal->scalefactorLastMajorMerger);
    
    thisOutGal->Mgas = (float) (thisGal->Mgas);
    thisOutGal->MgasIni = (float) (thisGal->MgasIni);
#ifndef FIRST
    thisOutGal->fracMgasMer = (float) (thisGal->fracMgasMer);
#endif
#if defined WITHMETALS
    thisOutGal->MgasNew = (float) (thisGal->MgasNew);
    thisOutGal->MgasEj = (float) (thisGal->MgasEj);
#endif
    thisOutGal->Mgas = (float) (thisGal->Mgas);
    thisOutGal->Mstar = (float) (thisGal->Mstar);
    
#ifndef FIRST
    thisOutGal->fesc = (float) (thisGal->fesc);
    thisOutGal->Nion = (float) (thisGal->Nion);
    thisOutGal->fej = (float) (thisGal->fej);
#endif
    thisOutGal->feff = (float) (thisGal->feff);
    thisOutGal->fg = (float) (thisGal->fg);
    thisOutGal->zreion = (float) (thisGal->zreion);
    thisOutGal->photHI_bg = (float) (thisGal->photHI_bg);

    for(int i=0; i<thisGal->snapnumber+1; i++)
        thisOutGal->stellarmasshistory[i] = (float) (thisGal->stellarmasshistory[i]);
    free(thisGal->stellarmasshistory);

#if defined WITHMETALS
    for(int i=0; i<3; i++)
    {
        thisOutGal->Mmetal[i] = (float) (thisGal->Mmetal[i]);
        thisOutGal->MmetalIni[i] = (float) (thisGal->MmetalIni[i]);
        thisOutGal->MmetalNew[i] = (float) (thisGal->MmetalNew[i]);
        thisOutGal->fracMmetalMer[i] = (float) (thisGal->fracMmetalMer[i]);
        thisOutGal->MmetalEj[i] = (float) (thisGal->MmetalEj[i]);
    }

    for(int i=0; i<thisGal->snapnumber+1; i++)
        thisOutGal->metalmasshistory[i] = (float) (thisGal->metalmasshistory[i]);
    free(thisGal->metalmasshistory);
    
    thisOutGal->Mdust = (float) (thisGal->Mdust);
#endif
}

void copy_gal_to_outGalList(gal_t *thisGal, int thisOutGal, outgalsnap_t *thisOutGalList)
{
    thisOutGalList[thisOutGal].scalefactor = (float) (thisGal->scalefactor);
    for(int i=0; i<3; i++)
    {
        thisOutGalList[thisOutGal].pos[i] = (float) (thisGal->pos[i]);
        thisOutGalList[thisOutGal].vel[i] = (float) (thisGal->vel[i]);
    }

    thisOutGalList[thisOutGal].Mvir = (float) (thisGal->Mvir);
    thisOutGalList[thisOutGal].Mvir_prog = (float) (thisGal->Mvir_prog);
    thisOutGalList[thisOutGal].Rvir = (float) (thisGal->Rvir);
    thisOutGalList[thisOutGal].velDisp = (float) (thisGal->velDisp);
    thisOutGalList[thisOutGal].velMax = (float) (thisGal->velMax);
    thisOutGalList[thisOutGal].spin = (float) (thisGal->spin);
    thisOutGalList[thisOutGal].scalefactorLastMajorMerger = (float) (thisGal->scalefactorLastMajorMerger);
    
    thisOutGalList[thisOutGal].MgasIni = (float) (thisGal->MgasIni);
#ifndef FIRST
    thisOutGalList[thisOutGal].fracMgasMer = (float) (thisGal->fracMgasMer);
#endif
#if defined WITHMETALS
    thisOutGalList[thisOutGal].MgasNew = (float) (thisGal->MgasNew);
    thisOutGalList[thisOutGal].MgasEj = (float) (thisGal->MgasEj);
#endif
    thisOutGalList[thisOutGal].Mgas = (float) (thisGal->Mgas);
    thisOutGalList[thisOutGal].Mstar = (float) (thisGal->Mstar);
    
#ifndef FIRST
    thisOutGalList[thisOutGal].fesc = (float) (thisGal->fesc);
    thisOutGalList[thisOutGal].Nion = (float) (thisGal->Nion);
    thisOutGalList[thisOutGal].fej = (float) (thisGal->fej);
#endif
    thisOutGalList[thisOutGal].feff = (float) (thisGal->feff);
    thisOutGalList[thisOutGal].fg = (float) (thisGal->fg);
    thisOutGalList[thisOutGal].zreion = (float) (thisGal->zreion);
    thisOutGalList[thisOutGal].photHI_bg = (float) (thisGal->photHI_bg);

    for(int i=0; i<thisGal->snapnumber+1; i++)
        thisOutGalList[thisOutGal].stellarmasshistory[i] = (float) (thisGal->stellarmasshistory[i]);
    free(thisGal->stellarmasshistory);

#if defined WITHMETALS
    for(int i=0; i<3; i++)
    {
        thisOutGalList[thisOutGal].Mmetal[i] = (float) (thisGal->Mmetal[i]);
        thisOutGalList[thisOutGal].MmetalIni[i] = (float) (thisGal->MmetalIni[i]);
        thisOutGalList[thisOutGal].MmetalNew[i] = (float) (thisGal->MmetalNew[i]);
        thisOutGalList[thisOutGal].fracMmetalMer[i] = (float) (thisGal->fracMmetalMer[i]);
        thisOutGalList[thisOutGal].MmetalEj[i] = (float) (thisGal->MmetalEj[i]);
    }

    for(int i=0; i<thisGal->snapnumber+1; i++)
        thisOutGalList[thisOutGal].metalmasshistory[i] = (float) (thisGal->metalmasshistory[i]);
    free(thisGal->metalmasshistory);
    
   thisOutGalList[thisOutGal].Mdust = (float) (thisGal->Mdust);
#endif
}

/*-------------------------------------------------------*/
/* WRITING OUTPUT FILES                                  */
/*-------------------------------------------------------*/

char *define_outfilename(int snap, const char *genericOutfileName)
{
    char *thisOutfileName = NULL;
    char outfileEnding[MAXLENGTH];
    
    sprintf(outfileEnding, "/Galsnap_%d.dat", snap);
    thisOutfileName = concat(genericOutfileName, outfileEnding);
    
    return thisOutfileName;
}

/* VERTICAL OUTPUT */
void write_treelist(dconfObj_t simParam, int thisRank, int numGtree, gtree_t **thisGtreeList)
{
    char thisOutfileName[MAXLENGTH];
    strcpy(thisOutfileName, simParam->outFileName);
    char char_thisRank[MAXLENGTH];
    outgtree_t *outTree;
    FILE *file;

    sprintf(char_thisRank, "%d", thisRank);
    strcat(thisOutfileName, "tree_");
    strcat(thisOutfileName, char_thisRank);
    strcat(thisOutfileName, ".dat");

    file = fopen(thisOutfileName, "wb");
    if(file == NULL)
      printf("Error opening file to write the trees\n");
    fwrite(&numGtree, sizeof(int32_t), 1, file);

    if(simParam->outputType == 2)
    {
      for(int i=0; i<numGtree; i++)
      {
          outTree = initOutGtree(thisGtreeList[i]->numGal);
          copy_tree_to_outgtree(thisGtreeList[i], outTree);
          fwrite(&(outTree->numGal), sizeof(int32_t), 1, file);
          deallocate_outgtree(outTree);
      }
      
      for(int i=0; i<numGtree; i++)
      {
          outTree = initOutGtree(thisGtreeList[i]->numGal);
          copy_tree_to_outgtree(thisGtreeList[i], outTree);
          fwrite(outTree->outgalaxies, sizeof(outgal_t), outTree->numGal, file);
          deallocate_outgtree(outTree);
      }
    }
    else
    {
      for(int i=0; i<numGtree; i++)
      {
          outTree = initOutGtree(thisGtreeList[i]->numGal);
          copy_tree_to_outgtree(thisGtreeList[i], outTree);
          fwrite(&(outTree->numGal), sizeof(int32_t), 1, file);
          fwrite(outTree->outgalaxies, sizeof(outgal_t), outTree->numGal, file);
          deallocate_outgtree(outTree);
      }
    }

    fclose(file); 
}

/* HORIZONTAL OUTPUT */
void write_galaxies_of_snap_to_file(dconfObj_t simParam, int snap, int numOutGal, outgalsnap_t *outGalList)
{
    char *thisOutfileName = NULL;
#ifdef MPI
    MPI_Status status;
    MPI_File mpiFile;

    /* Create the outgalsnap_t type for MPI */
    MPI_Datatype OUTGALSNAP_T;
#if defined FIRST
    MPI_Type_contiguous(OUTPUTSNAPNUMBER+21, MPI_FLOAT, &OUTGALSNAP_T);
#else
#if defined WITHMETALS
    MPI_Type_contiguous(2*OUTPUTSNAPNUMBER+43, MPI_FLOAT, &OUTGALSNAP_T);
#else
    MPI_Type_contiguous(OUTPUTSNAPNUMBER+25, MPI_FLOAT, &OUTGALSNAP_T);
#endif
#endif
    
    MPI_Type_commit(&OUTGALSNAP_T);
#else
    FILE *file;
#endif

    /* create filename for galaxy output */
    thisOutfileName = define_outfilename(snap, simParam->outFileName); 

    /* delete file if it already existed */
    remove(thisOutfileName);

#ifdef MPI
    /* all processors write galaxies into the file */
    MPI_File_open(MPI_COMM_WORLD, thisOutfileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mpiFile);
    MPI_File_write_shared(mpiFile, outGalList, numOutGal, OUTGALSNAP_T, &status);
    MPI_File_close(&mpiFile);
#else
    file = fopen(thisOutfileName, "wb");
    if(file == NULL)
      printf("Error opening file containing the number of galaxies");
    fwrite(outGalList, sizeof(outgalsnap_t), numOutGal, file);
    fclose(file);
#endif

    /* deallocation */
    free(thisOutfileName);
}

void write_num_galaxies_to_file(dconfObj_t simParam, int *numGalSnap, int numSnaps, int thisRank)
{
    int *numGalList = allocate_array_int(numSnaps, "numGalList");
    char *thisOutfileName = concat(simParam->outFileName, "/numGal.dat");
    
#ifdef MPI
    MPI_Reduce(numGalSnap, numGalList, numSnaps, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
    for(int i=0; i<numSnaps; i++) 
        numGalList[i] = numGalSnap[i];
#endif
    if(thisRank==0)
        for(int i=0; i<numSnaps; i++) 
            printf("snap %d : numgal = %d\n", i, numGalList[i]);
                      
    if(thisRank==0)
    {                 
        FILE *file;
        file = fopen(thisOutfileName, "w");
        if(file == NULL)
            printf("Error opening file containing the number of galaxies");
        fwrite(&numSnaps, sizeof(int), 1, file);
        fwrite(numGalList, sizeof(int), numSnaps, file);
        fclose(file);
    }
    
    free(thisOutfileName);
    free(numGalList);
}

/*-------------------------------------------------------*/
/* DELETING EXISTING FILES                               */
/*-------------------------------------------------------*/

void delete_file(int snap, const char *genericOutfileName)
{
    char *thisOutfileName = NULL;
    
    thisOutfileName = define_outfilename(snap, genericOutfileName);
    remove(thisOutfileName);
    
    free(thisOutfileName);
}

void delete_allfiles(int startSnap, int endSnap, int deltaSnap, const char *genericOutfileName)
{    
    for(int thisSnap=startSnap; thisSnap<=endSnap; thisSnap=thisSnap+deltaSnap)
    {
        delete_file(thisSnap, genericOutfileName);
    }
}

/*-------------------------------------------------------*/
/* COPYINH INIFILE                                       */
/*-------------------------------------------------------*/

void copy_iniFile_executable(char *outputfile)
{
    char *source_ini = "iniFile.ini";
    char *target_ini;
    
    char *source_exe = "astraeus";
    char *target_exe;
    
      
    target_ini = concat_strings(2, outputfile, "/iniFile.ini");
    target_exe = concat_strings(2, outputfile, "/astraeus");
    
    copy_file(source_ini, target_ini);
    copy_file(source_exe, target_exe);
    
    free(target_ini);
    free(target_exe);
}
