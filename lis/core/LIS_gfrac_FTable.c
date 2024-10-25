//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
//BOP
//
// !MODULE: LIS_gfrac_FTable
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

struct gfracsetnode
{ 
  char *name;
  void (*func)(int*);

  struct gfracsetnode* next;
} ;
struct gfracsetnode* gfracset_table = NULL; 

struct gfracreadnode
{ 
  char *name;
  void (*func)(int*,void*, void*,float*,float*);

  struct gfracreadnode* next;
} ;
struct gfracreadnode* gfracread_table = NULL; 

//BOP
// !ROUTINE: registergfracsetup
// \label{registergfracsetup}
// 
// !INTERFACE:
void FTN(registergfracsetup)(char *j,void (*func)(int*),int len)
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// setup gfrac data reading routines
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
  int len1;
  struct gfracsetnode* current;
  struct gfracsetnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct gfracsetnode*) malloc(sizeof(struct gfracsetnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(gfracset_table == NULL){
    gfracset_table = pnode;
  }
  else{
    current = gfracset_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}

//BOP
// !ROUTINE: gfracsetup
// \label{gfracsetup}
//
// !INTERFACE:
void FTN(gfracsetup)(char *j,int *n, int len)
//  
// !DESCRIPTION:
// Invokes the routine from the registry to 
// setup gfrac data reading
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
  struct gfracsetnode* current;
  
  current = gfracset_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("gfracsetup routine for runmode %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}

//BOP
// !ROUTINE: registerreadgfrac
// \label{registerreadgfrac}
// 
// !INTERFACE:
void FTN(registerreadgfrac)(char *j,void (*func)(int*, void*, void*, float*, float*),int len)
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read gfrac data
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
  int len1;
  struct gfracreadnode* current;
  struct gfracreadnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct gfracreadnode*) malloc(sizeof(struct gfracreadnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(gfracread_table == NULL){
    gfracread_table = pnode;
  }
  else{
    current = gfracread_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}

//BOP
// !ROUTINE: readgfrac
// \label{readgfrac}
//
// !INTERFACE:
void FTN(readgfrac)(char *j,int *n, void *time1, void *time2, float *array1, float *array2, int len)
//  
// !DESCRIPTION:
// Invokes the routine from the registry to 
// reading gfrac data 
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
  struct gfracreadnode* current;
  
  current = gfracread_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("readgfrac routine for runmode %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,time1,time2,array1,array2); 
}



