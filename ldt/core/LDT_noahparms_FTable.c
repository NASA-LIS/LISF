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
// !MODULE: LDT_noahparms_FTable
//  
//
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing different sources of 
//  greenness fraction data
//   
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include <string.h>

#include "ftn_drv.h"

struct tbotnode
{ 
  char *name;
  void (*func)(int*, float*);

  struct tbotnode* next;
} ;
struct tbotnode* tbot_table = NULL; 

//BOP
// !ROUTINE: registerreadtbot
// \label{registerreadtbot}
// 
// !INTERFACE:
void FTN(registerreadtbot)(char *j,void (*func)(int*,float*), int len)
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read tbot data
//
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \item[j]
//  index of the tbot data source
//  \end{description}
//EOP
{ 
  struct tbotnode* current;
  struct tbotnode* pnode; 
  // create node
  
  pnode=(struct tbotnode*) malloc(sizeof(struct tbotnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(tbot_table == NULL){
    tbot_table = pnode;
  }
  else{
    current = tbot_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readtbot
// \label{readtbot}
//
// !INTERFACE:
void FTN(readtbot)(char *j, int *n,float *array, int len)
//  
// !DESCRIPTION:
// Invokes the routine from the registry to 
// reading tbot data 
// 
// The arguments are: 
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//  index of the tbot data source
//  \item[array]
//  pointer to the tbot data
//  \end{description}
//EOP
{ 
  struct tbotnode* current;
  
  current = tbot_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Bottom temp reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,array); 
}

