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
// !MODULE: LIS_runmode_FTable
//  
//
// !DESCRIPTION:
//   Function table registries for storing the interface 
//   implementations of different running modes in LIS
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

struct rmodefinalnode
{ 
  char *name;
  void (*func)();

  struct rmodefinalnode* next;
} ;
struct rmodefinalnode* rmodefinal_table = NULL; 

//BOP
// !ROUTINE: registerlisinit
// \label{registerlisinit}
//
// !INTERFACE:
void FTN(registerlisinit)(char *j, void (*func)(),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the
//  LIS initialization method for a certain 
//  running mode
// 
// The arguments are:
// \begin{description}
//  \item[j]
//   name of the running mode
// \end{description}
//EOP
{ 
  int len1;
  struct rmodeinitnode* current;
  struct rmodeinitnode* pnode; 
  // create node
  
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
// !ROUTINE: lisinit
// \label{lisinit}
// 
// !INTERFACE:
void FTN(lisinit)(char *j,int len)
//  
// !DESCRIPTION:
//  Invokes the routine from the registry to 
//  perform the LIS initialization for the 
//  specified running mode
// 
// The arguments are:
// \begin{description}
//  \item[j]
//   name of the running mode
// \end{description}
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
// !ROUTINE: registerlisrun
// \label{registerlisrun}
// 
// !INTERFACE:
void FTN(registerlisrun)(char *j, void (*func)(),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the LIS 
//  run method for a certain running mode. 
// 
// The arguments are:
// \begin{description}
//  \item[j]
//   name of the running mode
// \end{description}
//EOP
{ 
  int len1;
  struct rmoderunnode* current;
  struct rmoderunnode* pnode; 
  // create node
  
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
// !ROUTINE: lisrun
// \label{lisrun}
// 
// !INTERFACE:
void FTN(lisrun)(char *j,int len)
//  
// !DESCRIPTION:
//  Invokes the LIS run method from the registry 
//  for the specified running mode
// 
// The arguments are:
// \begin{description}
//  \item[j]
//   name of the running mode
// \end{description}
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

//BOP
// !ROUTINE: registerlisfinalize
// \label{registerlisfinalize}
// 
// !INTERFACE:
void FTN(registerlisfinalize)(char *j, void (*func)(),int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// perform LIS finalization for the specified 
// running mode
// 
// The arguments are:
// \begin{description}
//  \item[j]
//   name of the running mode
// \end{description}
//EOP
{ 
  int len1;
  struct rmodefinalnode* current;
  struct rmodefinalnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct rmodefinalnode*) malloc(sizeof(struct rmodefinalnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(rmodefinal_table == NULL){
    rmodefinal_table = pnode;
  }
  else{
    current = rmodefinal_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: lisfinalize
// \label{lisfinalize}
// 
// !INTERFACE:
void FTN(lisfinalize)(char *j,int len)
//  
// !DESCRIPTION:
//  Invokes the routine from the registry to perform 
//  LIS finalization call for the specified running mode
// 
// The arguments are:
// \begin{description}
//  \item[j]
//   name of the running mode
// \end{description}
//EOP
{ 

  struct rmodefinalnode* current;
  
  current = rmodefinal_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("finalize routine for runmode %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}
