//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center Land Information System (LIS) v7.2
//
// Copyright (c) 2015 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
//BOP
//
// !MODULE: LIS_emiss_FTable
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

struct emisssetnode
{ 
  char *name;
  void (*func)(int*);

  struct emisssetnode* next;
} ;
struct emisssetnode* emissset_table = NULL; 

struct emissreadnode
{ 
  char *name;
  void (*func)(int*,void*,void*,float*,float*);

  struct emissreadnode* next;
} ;
struct emissreadnode* emissread_table = NULL; 

//BOP
// !ROUTINE: registeremissivitysetup
// \label{registeremissivitysetup}
// 
// !INTERFACE:
void FTN(registeremissivitysetup)(char *j,void (*func)(int*),int len)
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// setup emiss data reading routines
//
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \item[j]
//  index of the greenness data source
//  \end{description}
//EOP
{ 
  struct emisssetnode* current;
  struct emisssetnode* pnode; 
  // create node
  
  pnode=(struct emisssetnode*) malloc(sizeof(struct emisssetnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(emissset_table == NULL){
    emissset_table = pnode;
  }
  else{
    current = emissset_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}

//BOP
// !ROUTINE: emissivitysetup
// \label{emissivitysetup}
//
// !INTERFACE:
void FTN(emissivitysetup)(char *j,int *n, int len)
//  
// !DESCRIPTION:
// Invokes the routine from the registry to 
// setup emiss data reading
// 
// The arguments are: 
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//  index of the greenness data source
//  \end{description}
//EOP
{ 
  struct emisssetnode* current;
  
  current = emissset_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("emissivitysetup routine for runmode %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}

//BOP
// !ROUTINE: registerreademissivity
// \label{registerreademissivity}
// 
// !INTERFACE:
void FTN(registerreademissivity)(char *j,void (*func)(int*, void*, void*, float*, float*),int len)
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read emiss data
//
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \item[j]
//  index of the greenness data source
//  \end{description}
//EOP
{ 
  struct emissreadnode* current;
  struct emissreadnode* pnode; 
  // create node
  
  pnode=(struct emissreadnode*) malloc(sizeof(struct emissreadnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(emissread_table == NULL){
    emissread_table = pnode;
  }
  else{
    current = emissread_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}

//BOP
// !ROUTINE: reademissivity
// \label{reademissivity}
//
// !INTERFACE:
void FTN(reademissivity)(char *j,int *n, void *time1, void *time2, float *array1, float *array2, int len)
//  
// !DESCRIPTION:
// Invokes the routine from the registry to 
// reading emiss data 
// 
// The arguments are: 
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//  index of the greenness data source
//  \item[time]
//  month 
//  \item[array]
//  pointer to the greenness data
//  \end{description}
//EOP
{ 
  struct emissreadnode* current;
  
  current = emissread_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("reademissivity routine for runmode %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,time1, time2, array1, array2); 
}



