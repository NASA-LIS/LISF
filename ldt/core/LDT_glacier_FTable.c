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

//MN added to support glacier fraction
struct glacierfracnode
{
  char *name;
  void (*func)(int*, float*);

  struct glacierfracnode* next;
} ;

struct glacierfracnode* glacierfrac_table = NULL;

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
  int len1;
  struct glaciermasknode* current;
  struct glaciermasknode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct glaciermasknode*) malloc(sizeof(struct glaciermasknode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
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


//BOP
// !ROUTINE: registerreadglacierfrac
// \label{registerreadglacierfrac}
//
// !INTERFACE:
void FTN(registerreadglacierfrac)(char *j, void (*func)(int*, float*), int len)
//  
// !DESCRIPTION:
//  Makes an entry in the registry for the routine to read
//  glacierfrac data
// 
//  The arguments are: 
//  \begin{description}
//   \item[j] 
//    index of the glacierfrac source
//   \item[n]
//    index of the nest
//   \item[array]
//    pointer to the glacier fraction array
//  \end{description}
  //EOP
{
  int len1;
  struct glacierfracnode* current;
  struct glacierfracnode* pnode;
  // create node

  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct glacierfracnode*) malloc(sizeof(struct glacierfracnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL;

  if(glacierfrac_table == NULL){
    glacierfrac_table = pnode;
  }
  else{
    current = glacierfrac_table;
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode;
  }
}

//BOP
// !ROUTINE: readglacierfrac
// \label{readglacierfrac}
// 
// !INTERFACE:
void FTN(readglacierfrac)(char *j, int *n, float *array, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to read the 
//  glacierfrac data
//
//  The arguments are: 
//  \begin{description}
//   \item[j] 
//    index of the glacier fraction source
//   \item[n]
//    index of the nest
//   \item[array]
//    pointer to the glacier fraction array
//  \end{description}
//EOP
{
  struct glacierfracnode* current;

  current = glacierfrac_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n");
      printf("Glacier fraction reading routine for source %s is not defined\n",j);
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n");
      printf("****************Error****************************\n");
    }
  }

  current->func(n,array);
}
