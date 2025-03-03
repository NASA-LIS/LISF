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
// !MODULE: LIS_runoffdata_FTable
//  
//
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing the operations of different 
//  radiative transfer models
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"

struct runoffdatainitnode
{ 
  char *name;
  void (*func)();

  struct runoffdatainitnode* next;
} ;
struct runoffdatainitnode* runoffdatainit_table = NULL; 


struct runoffdatareadnode
{ 
  char *name;
  void (*func)(int*, float*, float*);

  struct runoffdatareadnode* next;
} ;
struct runoffdatareadnode* runoffdataread_table = NULL; 

//BOP
// !ROUTINE: registerinitrunoffdata
// \label{registerinitrunoffdata}
// 
// !INTERFACE:
void FTN(registerinitrunoffdata)(char *j, void (*func)(), int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  initialize the radiative transfer model
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the runoffdata
//. \end{description}
//EOP
{ 
  int len1;
  struct runoffdatainitnode* current;
  struct runoffdatainitnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct runoffdatainitnode*) malloc(sizeof(struct runoffdatainitnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(runoffdatainit_table == NULL){
    runoffdatainit_table = pnode;
  }
  else{
    current = runoffdatainit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: initrunoffdata
// \label{initrunoffdata}
// 
// !INTERFACE:
void FTN(initrunoffdata)(char *j, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to
//  initialize the runoffdata
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the runoffdata
// \end{description}
//EOP
{ 
  struct runoffdatainitnode* current;
  
  current = runoffdatainit_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("init routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}

//BOP
// !ROUTINE: registerreadrunoffdata
// \label{registerreadrunoffdata}
// 
// !INTERFACE:
void FTN(registerreadrunoffdata)(char *j, void (*func)(int*,float*, float*), int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  specify the forward model integration step of the runoffdata
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the runoffdata
//. \end{description}
//EOP
{ 
  int len1;
  struct runoffdatareadnode* current;
  struct runoffdatareadnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct runoffdatareadnode*) malloc(sizeof(struct runoffdatareadnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(runoffdataread_table == NULL){
    runoffdataread_table = pnode;
  }
  else{
    current = runoffdataread_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: readrunoffdata
// \label{readrunoffdata}
// 
// !INTERFACE:
void FTN(readrunoffdata)(char *j, int *index, float *qs, float *qsb,int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to
//  call the forward model integration step of 
//  the runoffdata
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the runoffdata
// \end{description}
//EOP
{ 
  struct runoffdatareadnode* current;
  
  current = runoffdataread_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("runoffdata run routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(index,qs,qsb); 
}



