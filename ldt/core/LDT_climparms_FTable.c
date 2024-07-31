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
  int len1;
  struct climpptnode* current;
  struct climpptnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct climpptnode*) malloc(sizeof(struct climpptnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
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
  int len1;
  struct climtminnode* current;
  struct climtminnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct climtminnode*) malloc(sizeof(struct climtminnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
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
  int len1;
  struct climtmaxnode* current;
  struct climtmaxnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct climtmaxnode*) malloc(sizeof(struct climtmaxnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
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




