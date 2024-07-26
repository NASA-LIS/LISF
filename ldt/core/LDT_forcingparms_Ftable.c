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
// !MODULE: LDT_forcingparms_FTable
//  
//
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing different sources of 
//  forcing parameter datasets (e.g., elevation correction).
//   
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include <string.h>

#include "ftn_drv.h"

struct forcelevnode
{ 
  char *name;
  void (*func)(int*, int*, float*, float*);

  struct forcelevnode* next;
} ;
struct forcelevnode* forcelev_table = NULL; 

//BOP
// !ROUTINE: registerreadforcelev
// \label{registerreadforcelev}
// 
// !INTERFACE:
void FTN(registerreadforcelev)(char *j,void (*func)(int*,int*,float*,float*),int len)
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read forcelev data
//
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \item[j]
//  index of the forcing parameter source
//  \end{description}
//EOP
{ 
  int len1;
  struct forcelevnode* current;
  struct forcelevnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct forcelevnode*) malloc(sizeof(struct forcelevnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(forcelev_table == NULL){
    forcelev_table = pnode;
  }
  else{
    current = forcelev_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readforcelev
// \label{readforcelev}
//
// !INTERFACE:
void FTN(readforcelev)(char *j, int *n, int *m, float *elev, float *elevdiff,int len)
//  
// !DESCRIPTION:
// Invokes the routine from the registry to 
// reading forcelev data 
// 
// The arguments are: 
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//  index of the forcing parameter data source
//  \item[array]
//  pointer to the forcing parameter
//  \end{description}
//EOP
{ 
  struct forcelevnode* current;
  
  current = forcelev_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Forcing elev reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,m,elev,elevdiff); 
}

