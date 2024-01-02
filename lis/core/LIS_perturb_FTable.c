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
// !MODULE: LIS_perturb_FTable
//  
//
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing the operations of different 
//  perturbation algorithms
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"

struct pertinitnode
{ 
  char *name;
  void (*func)(int*);

  struct pertinitnode* next;
} ;
struct pertinitnode* pertinit_table = NULL; 

struct pertsetnode
{ 
  char *name;
  void (*func)(int*, int*, void*, void*);

  struct pertsetnode* next;
} ;
struct pertsetnode* pertset_table = NULL; 

struct pertmethodnode
{ 
  char *name;
  void (*func)(int*, int*, int*, void*, void*);

  struct pertmethodnode* next;
} ;
struct pertmethodnode* pertmethod_table = NULL; 

struct pertrstnode
{ 
  char *name;
  void (*func)();

  struct pertrstnode* next;
} ;
struct pertrstnode* pertrst_table = NULL; 

struct pertwrtnode
{ 
  char *name;
  void (*func)(int*);

  struct pertwrtnode* next;
} ;
struct pertwrtnode* pertwrt_table = NULL; 

//BOP
// !ROUTINE: registerperturbinit
// \label{registerperturbinit}
// 
// !INTERFACE:
void FTN(registerperturbinit)(char *j, void (*func)(int*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  initialize the pertubation scheme. 
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the perturbation algorithm
//. \end{description}
//EOP
{ 
  int len1;
  struct pertinitnode* current;
  struct pertinitnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct pertinitnode*) malloc(sizeof(struct pertinitnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(pertinit_table == NULL){
    pertinit_table = pnode;
  }
  else{
    current = pertinit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: perturbinit
// \label{perturbinit}
// 
// !INTERFACE:
void FTN(perturbinit)(char *j, int *index,int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to
//  initialize the pertubation scheme
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the perturbation algorithm
// \end{description}
//EOP
{ 
  struct pertinitnode* current;
  
  current = pertinit_table;

  while(strcmp(current->name,j)!=0){
    current = current->next;

    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("init routine for perturbation scheme %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(index); 
}


//BOP
// !ROUTINE: registerperturbsetup
// \label{registerperturbsetup}
// 
// !INTERFACE:
void FTN(registerperturbsetup)(char *j, void (*func)(int*, int*, void*, void*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  initialize the pertubation scheme. 
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the perturbation algorithm
//. \end{description}
//EOP
{ 
  int len1;
  struct pertsetnode* current;
  struct pertsetnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct pertsetnode*) malloc(sizeof(struct pertsetnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(pertset_table == NULL){
    pertset_table = pnode;
  }
  else{
    current = pertset_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: perturbsetup
// \label{perturbsetup}
// 
// !INTERFACE:
void FTN(perturbsetup)(char *j, int *index, int *k,  void *base, void *pert,int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to
//  initialize the pertubation scheme
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the perturbation algorithm
//   \item[index]
//   object of perturbation (forcing, state, obs)
//   \item[k]
//   index of the DA instance
//   \item[base]
//   base state for perturbation
//   \item[pert]
//   output perturbation state 
// \end{description}
//EOP
{ 
  struct pertsetnode* current;
  
  current = pertset_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("setup routine for perturbation scheme %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(index,k,base,pert); 
}


//BOP
// !ROUTINE: registerperturbmethod
// \label{registerperturbmethod}
// 
// !INTERFACE:
void FTN(registerperturbmethod)(char *j, void (*func)(int*, int*, int*, void*, void*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  invoke the pertubation scheme. 
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the perturbation algorithm
//. \end{description}
//EOP
{ 
  int len1;
  struct pertmethodnode* current;
  struct pertmethodnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct pertmethodnode*) malloc(sizeof(struct pertmethodnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(pertmethod_table == NULL){
    pertmethod_table = pnode;
  }
  else{
    current = pertmethod_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: perturbmethod
// \label{perturbmethod}
// 
// !INTERFACE:
void FTN(perturbmethod)(char *j, int *index, int *n, int *k, void *base, void *pert,int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to
//  call the pertubation scheme
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the perturbation algorithm
//   \item[n]
//   index of the nest
//   \item[index]
//   object of perturbation (forcing, state, obs)
//   \item[k]
//   index of the DA instance
//   \item[base]
//   base state for perturbation
//   \item[pert]
//   output perturbation state 
// \end{description}
//EOP
{ 
  struct pertmethodnode* current;
  
  current = pertmethod_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("run routine for perturbation scheme %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(index,n,k,base,pert); 
}

//BOP
// !ROUTINE: registerperturbwriterst
// \label{registerperturbwriterst}
// 
// !INTERFACE:
void FTN(registerperturbwriterst)(char *j, void (*func)(int*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  write restart files for the perturbations scheme. 
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the perturbation algorithm
//. \end{description}
//EOP
{ 
  int len1;
  struct pertwrtnode* current;
  struct pertwrtnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct pertwrtnode*) malloc(sizeof(struct pertwrtnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(pertwrt_table == NULL){
    pertwrt_table = pnode;
  }
  else{
    current = pertwrt_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: writepertrestart
// \label{writepertrestart}
// 
// !INTERFACE:
void FTN(writepertrestart)(char *j, int *nest,int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to
//  write the restart file for the pertubation scheme
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the perturbation algorithm
//   \item[nest]
//    nest index
// \end{description}
//EOP
{ 
  struct pertwrtnode* current;
  
  current = pertwrt_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("write restart routine for perturbation scheme %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(nest); 
}

//BOP
// !ROUTINE: registerperturbreadrst
// \label{registerperturbreadrst}
// 
// !INTERFACE:
void FTN(registerperturbreadrst)(char *j, void (*func)(),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  read restart files for the perturbations scheme. 
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the perturbation algorithm
//. \end{description}
//EOP
{ 
  int len1;
  struct pertrstnode* current;
  struct pertrstnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct pertrstnode*) malloc(sizeof(struct pertrstnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(pertrst_table == NULL){
    pertrst_table = pnode;
  }
  else{
    current = pertrst_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: readpertrestart
// \label{readpertrestart}
// 
// !INTERFACE:
void FTN(readpertrestart)(char *j,int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to
//  read the restart file for the pertubation scheme
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the perturbation algorithm
//   \item[nest]
//    nest index
// \end{description}
//EOP
{ 
  struct pertrstnode* current;
  
  current = pertrst_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("read restart routine for perturbation scheme %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}



