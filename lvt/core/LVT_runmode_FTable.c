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
// !MODULE: LVT_runmode_FTable
//  
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing the operation of different 
//  runmode specifications
// 
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"

struct rmodeinitnode
{ 
  char *name;
  void (*func)();

  struct rmodeinitnode* next;
} ;
struct rmodeinitnode* rmodeinit_table = NULL; 

struct rmoderunnode
{ 
  char *name;
  void (*func)();

  struct rmoderunnode* next;
} ;

struct rmoderunnode* rmoderun_table = NULL; 

//BOP
// !ROUTINE: registerlvtinit
// \label{registerlvtinit}
//
// !INTERFACE:
void FTN(registerlvtinit)(char *j,void (*func)(),int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read the runtime domain specifics
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//EOP
{ 
  int len1;
  struct rmodeinitnode* pnode;
  struct rmodeinitnode* current;

  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct rmodeinitnode*) malloc(sizeof(struct rmodeinitnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(rmodeinit_table == NULL){
    rmodeinit_table = pnode;
  }
  else{
    current = rmodeinit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: lvtinit
// \label{lvtinit}
//
// !INTERFACE:
void FTN(lvtinit)(char *j,int len)
// !DESCRIPTION: 
//  Calls the routine from the registry to read the
//  the runtime domain specifics
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//
//EOP
{ 
  struct rmodeinitnode* current;
  
  current = rmodeinit_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("init routine for runmode %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}

//BOP
// !ROUTINE: registerlvtrun
// \label{registerlvtrun}
//
// !INTERFACE:
void FTN(registerlvtrun)(char *j,void (*func)(), int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read the runtime domain specifics
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//EOP
{ 
  int len1;
  struct rmoderunnode* pnode;
  struct rmoderunnode* current;

  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct rmoderunnode*) malloc(sizeof(struct rmoderunnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(rmoderun_table == NULL){
    rmoderun_table = pnode;
  }
  else{
    current = rmoderun_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: lvtrun
// \label{lvtrun}
//
// !INTERFACE:
void FTN(lvtrun)(char *j, int len)
// !DESCRIPTION: 
//  Calls the routine from the registry to read the
//  the runtime domain specifics
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//
//EOP
{ 

  struct rmoderunnode* current;
  
  current = rmoderun_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("run routine for runmode %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}




