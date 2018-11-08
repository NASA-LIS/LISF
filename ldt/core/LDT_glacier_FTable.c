//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center Land Information System (LDT)
//
// See RELEASE_NOTES.txt for more information.
//
// The LDT source code and documentation are not in the public domain
// and may not be freely distributed.  Only qualified entities may receive 
// the source code and documentation. 
//
// Qualified entities must be covered by a Software Usage Agreement. 
// The Software Usage Agreement contains all the terms and conditions
// regarding the release of the LDT software.
//
// NASA GSFC MAKES NO REPRESENTATIONS ABOUT THE SUITABILITY OF THE
// SOFTWARE FOR ANY PURPOSE.  IT IS PROVIDED AS IS WITHOUT EXPRESS OR
// IMPLIED WARRANTY.  NEITHER NASA GSFC NOR THE US GOVERNMENT SHALL BE
// LIABLE FOR ANY DAMAGES SUFFERED BY THE USER OF THIS SOFTWARE.
//
// See the Software Usage Agreement for the full disclaimer of warranty.
//
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
//BOP
//
//  !MODULE: LDT_glacier_FTable
//  
//
// !DESCRIPTION:
//   Function table registries for storing the interface 
//   implementations for managing different sources of 
//   glacier datasets
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"

struct glaciermasknode
{ 
  char *name;
  void (*func)(int*, float*);

  struct glaciermasknode* next;
} ;

struct glaciermasknode* glaciermask_table = NULL; 


//BOP
// !ROUTINE: registerreadglaciermask
// \label{registerreadglaciermask}
// 
// !INTERFACE:
void FTN(registerreadglaciermask)(char *j, void (*func)(int*, float*), int len)
//  
// !DESCRIPTION:
//  Makes an entry in the registry for the routine to read
//  glacier mask data
//
//  The arguments are: 
//  \begin{description}
//   \item[j] 
//    index of the topography source
//  \end{description}
  //EOP
{ 
  struct glaciermasknode* current;
  struct glaciermasknode* pnode; 
  // create node
  
  pnode=(struct glaciermasknode*) malloc(sizeof(struct glaciermasknode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(glaciermask_table == NULL){
    glaciermask_table = pnode;
  }
  else{
    current = glaciermask_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readglaciermask
// \label{readglaciermask}
// 
// !INTERFACE:
void FTN(readglaciermask)(char *j, int *n, float *array, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to read the 
//  elevation data
//
//  The arguments are: 
//  \begin{description}
//   \item[j] 
//    index of the topography source
//   \item[n]
//    index of the nest
//   \item[nt]
//    number of types or bins (bands)
//   \item[fgrd]
//    gridcell fraction of type or values
//   \item[array]
//    pointer to the glacier mask array
//  \end{description}
//EOP
{ 
  struct glaciermasknode* current;
  
  current = glaciermask_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Glacier mask reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,array); 
}



