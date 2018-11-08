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
//  !MODULE: LDT_irrig_FTable
//  
//
// !DESCRIPTION:
//   Function table registries for storing the interface 
//   implementations for managing different sources of 
//   irrigation data
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"

struct irrigtypenode
{ 
  char *name;
  void (*func)(int*, float*, int*);

  struct irrigtypenode* next;
} ;

struct irrigtypenode* irrigtype_table = NULL; 

struct irrigfracnode
{ 
  char *name;
  void (*func)(int*, float*);

  struct irrigfracnode* next;
} ;

struct irrigfracnode* irrigfrac_table = NULL; 


//BOP
// !ROUTINE: registerreadirrigtype
// \label{registerreadirrigtype}
// 
// !INTERFACE:
void FTN(registerreadirrigtype)(char *j, void (*func)(int*, float*, int*), int len)
//  
// !DESCRIPTION:
//  Makes an entry in the registry for the routine to read
//  irrigation data
//
//  The arguments are: 
//  \begin{description}
//   \item[j] 
//    index of the irrigation source
//  \end{description}
  //EOP
{ 
  struct irrigtypenode* current;
  struct irrigtypenode* pnode; 
  // create node
  
  pnode=(struct irrigtypenode*) malloc(sizeof(struct irrigtypenode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(irrigtype_table == NULL){
    irrigtype_table = pnode;
  }
  else{
    current = irrigtype_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readirrigtype
// \label{readirrigtype}
// 
// !INTERFACE:
void FTN(readirrigtype)(char *j, int *n, float *array, int *nt, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to read the 
//  irrigation data
//
//  The arguments are: 
//  \begin{description}
//   \item[j] 
//    index of the irrigation type source
//   \item[n]
//    index of the nest
//   \item[array]
//    pointer to the irrigation type array
//  \end{description}
//EOP
{ 
  struct irrigtypenode* current;
  
  current = irrigtype_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Irrigation reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,array,nt); 
}

//BOP
// !ROUTINE: registerreadirrigfrac
// \label{registerreadirrigfrac}
//
// !INTERFACE:
void FTN(registerreadirrigfrac)(char *j, void (*func)(int*, float*), int len)
//  
// !DESCRIPTION:
//  Makes an entry in the registry for the routine to read
//  irrigfrac data
// 
//  The arguments are: 
//  \begin{description}
//   \item[j] 
//    index of the irrigation source
//   \item[n]
//    index of the nest
//   \item[array]
//    pointer to the irrigation fraction array
//  \end{description}
  //EOP
{ 

  struct irrigfracnode* current;
  struct irrigfracnode* pnode; 
  // create node
  
  pnode=(struct irrigfracnode*) malloc(sizeof(struct irrigfracnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(irrigfrac_table == NULL){
    irrigfrac_table = pnode;
  }
  else{
    current = irrigfrac_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readirrigfrac
// \label{readirrigfrac}
// 
// !INTERFACE:
void FTN(readirrigfrac)(char *j, int *n, float *array, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to read the 
//  irrigfrac data
//
//  The arguments are: 
//  \begin{description}
//   \item[j] 
//    index of the irrigation fraction source
//   \item[n]
//    index of the nest
//   \item[array]
//    pointer to the irrigation fraction array
//  \end{description}
//EOP
{ 
  struct irrigfracnode* current;
  
  current = irrigfrac_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Irrigation fraction reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  
  current->func(n,array); 
}

