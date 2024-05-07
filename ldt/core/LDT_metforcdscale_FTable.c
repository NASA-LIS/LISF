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
// !MODULE: LDT_metforcdscale_FTable
//  
//
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing different sources of 
//  forcing parametersets (e.g., forcing temporal downscaling).
//   
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include <string.h>

#include "ftn_drv.h"

struct timedscalenode
{ 
  char *name;
  void (*func)(int*);

  struct timedscalenode* next;
} ;
struct timedscalenode* tmdscale_table = NULL; 

//BOP
// !ROUTINE: registerapplytimedscale
// \label{registerapplytimedscale}
// 
// !INTERFACE:
void FTN(registerapplytimedscale)(char *j,void (*func)(int*),int len)
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
//  read forcing data to be downscaled.
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
  struct timedscalenode* current;
  struct timedscalenode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct timedscalenode*) malloc(sizeof(struct timedscalenode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(tmdscale_table == NULL){
    tmdscale_table = pnode;
  }
  else{
    current = tmdscale_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: applytimedscale
// \label{applytimedscale}
//
// !INTERFACE:
void FTN(applytimedscale)(char *j, int *n, int len)
//  
// !DESCRIPTION:
// Invokes the routine from the registry to 
//  reading forcing to be temporally downscaled.
// 
// The arguments are: 
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//   index of the forcing parameter data source
//  \item[array]
//   pointer to the forcing parameter
//  \end{description}
//EOP
{ 
  struct timedscalenode* current;
  
  current = tmdscale_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Temporal downscaling routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}

