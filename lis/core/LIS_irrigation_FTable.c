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
// !MODULE: LIS_irrigation_FTable
//  
//
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing the operations of different 
//  land slide models. 
//
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"

struct irrigationinitnode
{ 
  char *name;
  void (*func)(void*);

  struct irrigationinitnode* next;
} ;
struct irrigationinitnode* irrigationinit_table = NULL; 

struct irrigationapplyupdatenode
{ 
  char *name;
  void (*func)(int*,void*);

  struct irrigationapplyupdatenode* next;
} ;
struct irrigationapplyupdatenode* irrigationapplyupdate_table = NULL; 


struct irrigationfinalnode
{ 
  char *name;
  void (*func)();

  struct irrigationfinalnode* next;
} ;
struct irrigationfinalnode* irrigationfinal_table = NULL; 


//BOP
// !ROUTINE: registerirrigationschemeinit
// \label{registerirrigationschemeinit}
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine to 
//  perform land surface model initialization
// 
// !INTERFACE:
void FTN(registerirrigationschemeinit)(char *j, void (*func)(void*),int len)
//EOP
{ 
  int len1;
  struct irrigationinitnode* current;
  struct irrigationinitnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct irrigationinitnode*) malloc(sizeof(struct irrigationinitnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(irrigationinit_table == NULL){
    irrigationinit_table = pnode;
  }
  else{
    current = irrigationinit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: irrigationschemeinit
// \label{irrigationschemeinit}
//
// !INTERFACE:
void FTN(irrigationschemeinit)(char *j,void *state,int len)
//
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to perform
//  land surface scheme initialization
// 
//  \begin{description}
//  \item[i]
//   index of the LSM
//  \end{description}
//EOP
{ 
  struct irrigationinitnode* current;
  
  current = irrigationinit_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("init routine for irrigation scheme %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(state); 

}
//BOP
// !ROUTINE: registerirrigationupdate
// \label{registerirrigationupdate}
//
// !INTERFACE:
void FTN(registerirrigationupdate)(char *j, void (*func)(int*, void*),int len)
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine
//  to run the land surface model 
// 
//  \begin{description}
//  \item[i]
//   index of the LSM
//  \end{description}
//EOP
{ 
  int len1;
  struct irrigationapplyupdatenode* current;
  struct irrigationapplyupdatenode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct irrigationapplyupdatenode*) malloc(sizeof(struct irrigationapplyupdatenode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(irrigationapplyupdate_table == NULL){
    irrigationapplyupdate_table = pnode;
  }
  else{
    current = irrigationapplyupdate_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: applyirrigationupdates
// \label{applyirrigationupdates}
//
// !INTERFACE:
void FTN(applyirrigationupdates)(char *j,int *n, void *state, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to run the 
//  land surface model 
// 
//  \begin{description}
//  \item[i]
//   index of the LSM
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{
  struct irrigationapplyupdatenode* current;
  
  current = irrigationapplyupdate_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("run routine for irrigation model %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,state); 
}



