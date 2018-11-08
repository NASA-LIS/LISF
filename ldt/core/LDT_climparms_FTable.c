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
// !MODULE: LDT_climparms_FTable
//  
//
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing different sources of 
//  forcing climatology datasets (e.g., forcing downscaling).
//   
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include <string.h>

#include "ftn_drv.h"

struct climpptnode
{ 
  char *name;
  void (*func)(int*, int*, int*, float*, float*);

  struct climpptnode* next;
} ;
struct climpptnode* climppt_table = NULL; 

struct climtminnode
{ 
  char *name;
  void (*func)(int*, float*);

  struct climtminnode* next;
} ;
struct climtminnode* climtmin_table = NULL; 

struct climtmaxnode
{ 
  char *name;
  void (*func)(int*, float*);

  struct climtmaxnode* next;
} ;
struct climtmaxnode* climtmax_table = NULL; 

//BOP
// !ROUTINE: registerreadclimppt
// \label{registerreadclimppt}
// 
// !INTERFACE:
void FTN(registerreadclimppt)(char *j,void (*func)(int*,int*,int*,float*,float*),int len)
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read climppt data
//
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \item[j]
//  index of the climatology data source
//  \end{description}
//EOP
{ 
  struct climpptnode* current;
  struct climpptnode* pnode; 
  // create node
  
  pnode=(struct climpptnode*) malloc(sizeof(struct climpptnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(climppt_table == NULL){
    climppt_table = pnode;
  }
  else{
    current = climppt_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readclimppt
// \label{readclimppt}
//
// !INTERFACE:
void FTN(readclimppt)(char *j, int *n, int *nc, int *nr, float *gridarray, float *array,int len)
//  
// !DESCRIPTION:
// Invokes the routine from the registry to 
// reading climppt data 
// 
// The arguments are: 
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//  index of the climatology data source
//  \item[array]
//  pointer to the climatology data
//  \end{description}
//EOP
{ 
  struct climpptnode* current;
  
  current = climppt_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("CLIMPPT reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,nc,nr,gridarray,array); 
}


//BOP
// !ROUTINE: registerreadclimtmin
// \label{registerreadclimtmin}
// 
// !INTERFACE:
void FTN(registerreadclimtmin)(char *j,void (*func)(int*,float*),int len)
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read climtmin data
//
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \item[j]
//  index of the climatology data source
//  \end{description}
//EOP
{ 
  struct climtminnode* current;
  struct climtminnode* pnode; 
  // create node
  
  pnode=(struct climtminnode*) malloc(sizeof(struct climtminnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(climtmin_table == NULL){
    climtmin_table = pnode;
  }
  else{
    current = climtmin_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readclimtmin
// \label{readclimtmin}
//
// !INTERFACE:
void FTN(readclimtmin)(char *j, int *n,float *array,int len)
//  
// !DESCRIPTION:
// Invokes the routine from the registry to 
// reading climtmin data 
// 
// The arguments are: 
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//  index of the climatology data source
//  \item[array]
//  pointer to the climatology data
//  \end{description}
//EOP
{ 
  struct climtminnode* current;
  
  current = climtmin_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("CLIMTMIN reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,array); 
}


//BOP
// !ROUTINE: registerreadclimtmax
// \label{registerreadclimtmax}
// 
// !INTERFACE:
void FTN(registerreadclimtmax)(char *j,void (*func)(int*,float*),int len)
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read climtmax data
//
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \item[j]
//  index of the climatology ta source
//  \end{description}
//EOP
{ 
  struct climtmaxnode* current;
  struct climtmaxnode* pnode; 
  // create node
  
  pnode=(struct climtmaxnode*) malloc(sizeof(struct climtmaxnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(climtmax_table == NULL){
    climtmax_table = pnode;
  }
  else{
    current = climtmax_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readclimtmax
// \label{readclimtmax}
//
// !INTERFACE:
void FTN(readclimtmax)(char *j, int *n,float *array,int len)
//  
// !DESCRIPTION:
// Invokes the routine from the registry to 
// reading climtmax data 
// 
// The arguments are: 
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//  index of the climatology data source
//  \item[array]
//  pointer to the climatology data
//  \end{description}
//EOP
{ 
  struct climtmaxnode* current;
  
  current = climtmax_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("CLIMTMAX reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,array); 
}




