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
// !MODULE: LIS_biasestim_FTable
//  
// !DESCRIPTION:
//  Function table registries for storing the interface
//  implementations for the operation of various onlie bias
//  estimation routines to be used in data assimilation
//
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"
struct biasinitnode
{ 
  char *name;
  void (*func)();

  struct biasinitnode* next;
} ;
struct biasinitnode* biasinit_table = NULL; 

struct biassetupnode
{ 
  char *name;
  void (*func)(int*);

  struct biassetupnode* next;
} ;
struct biassetupnode* biassetup_table = NULL; 

struct biascompnode
{ 
  char *name;
  void (*func)(int*,int*);

  struct biascompnode* next;
} ;
struct biascompnode* biascomp_table = NULL; 

struct biasupdnode
{ 
  char *name;
  void (*func)(int*,int*);

  struct biasupdnode* next;
} ;
struct biasupdnode* biasupd_table = NULL; 

struct biasoutnode
{ 
  char *name;
  void (*func)(int*,int*);

  struct biasoutnode* next;
} ;
struct biasoutnode* biasout_table = NULL; 

struct biasfinalnode
{ 
  char *name;
  void (*func)();

  struct biasfinalnode* next;
} ;
struct biasfinalnode* biasfinal_table = NULL; 

//BOP
// !ROUTINE: registerbiasestimationinit
// \label{registerrbiasestimationinit}
//  
// 
// !INTERFACE:
void FTN(registerbiasestimationinit)(char *j, void (*func)(),int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine that
//  sets up structures for initializing bias estimation algorithm
//  for use in data assimilation
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the bias estimation algorithm. 
//  \end{description}
//EOP
{ 
  int len1;
  struct biasinitnode* current;
  struct biasinitnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct biasinitnode*) malloc(sizeof(struct biasinitnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(biasinit_table == NULL){
    biasinit_table = pnode;
  }
  else{
    current = biasinit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: biasestimationinit
// \label{biasestimationinit}
// 
// !INTERFACE:
void FTN(biasestimationinit)(char *j,int len)
//  
// !DESCRIPTION: 
// Invokes the the selected bias estimation algorithm
// initialization for use in data assimilation
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the bias estimation scheme
//  \end{description}
//
//EOP
{ 
  struct biasinitnode* current;
  
  current = biasinit_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;

    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("init routine for bias estimation scheme %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}

//BOP
// !ROUTINE: registerbiasestimationsetup
// \label{registerrbiasestimationsetup}
//  
// 
// !INTERFACE:
void FTN(registerbiasestimationsetup)(char *j, void (*func)(int*),int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine that
//  sets up structures for setting up bias estimation algorithm
//  for use in data assimilation
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the bias estimation algorithm. 
//  \end{description}
//EOP
{ 
  int len1;
  struct biassetupnode* current;
  struct biassetupnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct biassetupnode*) malloc(sizeof(struct biassetupnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(biassetup_table == NULL){
    biassetup_table = pnode;
  }
  else{
    current = biassetup_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: biasestimationsetup
// \label{biasestimationsetup}
// 
// !INTERFACE:
void FTN(biasestimationsetup)(char *j, int *k,int len)
//  
// !DESCRIPTION: 
// Invokes the the selected bias estimation algorithm
// setup for use in data assimilation
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the bias estimation scheme
//  \item[k]
//   index of the DA instance
//  \end{description}
//
//EOP
{ 
  struct biassetupnode* current;
  
  current = biassetup_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("setup routine for bias estimation scheme %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(k); 
}

//BOP
// !ROUTINE: registerbiasestimationcompute
// \label{registerrbiasestimationcompute}
//  
// 
// !INTERFACE:
void FTN(registerbiasestimationcompute)(char *j, void (*func)(int*, int*),int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine that
//  sets up structures for computing bias correction
//  for use in data assimilation
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the bias estimation algorithm. 
//  \end{description}
//EOP
{ 
  int len1;
  struct biascompnode* current;
  struct biascompnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct biascompnode*) malloc(sizeof(struct biascompnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(biascomp_table == NULL){
    biascomp_table = pnode;
  }
  else{
    current = biascomp_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: computebiascorrection
// \label{computebiascorrection}
// 
// !INTERFACE:
void FTN(computebiascorrection)(char *j, int *n, int *k,int len)
//  
// !DESCRIPTION: 
// Invokes the the selected bias estimation algorithm
// for computing bias correction prior to model propagation
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the bias estimation scheme
//  \item[n]
//   index of the nest
//  \item[k]
//   index of the assimilation instance
//  \end{description}
//
//EOP
{ 
  struct biascompnode* current;
  
  current = biascomp_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("compute routine for bias estimation scheme %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,k); 
}


//BOP
// !ROUTINE: registerbiasestimationupdate
// \label{registerrbiasestimationupdate}
//  
// 
// !INTERFACE:
void FTN(registerbiasestimationupdate)(char *j, void (*func)(int*, int*),int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine that
//  sets up structures for updating bias correction estimation
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the bias estimation algorithm. 
//  \end{description}
//EOP
{ 
  int len1;
  struct biasupdnode* current;
  struct biasupdnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct biasupdnode*) malloc(sizeof(struct biasupdnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(biasupd_table == NULL){
    biasupd_table = pnode;
  }
  else{
    current = biasupd_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: applybiascorrection
// \label{applybiascorrection}
// 
// !INTERFACE:
void FTN(applybiascorrection)(char *j, int *n, int *k,int len)
//  
// !DESCRIPTION: 
// Invokes the the selected bias estimation algorithm
// for updating bias correction estimate
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the bias estimation scheme
//  \item[n]
//   index of the nest
//  \item[k]
//   index of the assimilation instance
//  \end{description}
//
//EOP
{ 
  struct biasupdnode* current;
  
  current = biasupd_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("apply bias outine for bias estimation scheme %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,k); 
}

//BOP
// !ROUTINE: registerbiasestimationrestart
// \label{registerrbiasestimationrestart}
//  
// 
// !INTERFACE:
void FTN(registerbiasestimationrestart)(char *j, void (*func)(int*, int*),int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine that
//  sets up structures for writing bias restart 
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the bias estimation algorithm. 
//  \end{description}
//EOP
{ 
  int len1;
  struct biasoutnode* current;
  struct biasoutnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct biasoutnode*) malloc(sizeof(struct biasoutnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(biasout_table == NULL){
    biasout_table = pnode;
  }
  else{
    current = biasout_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: writebiasrestart
// \label{writebiasrestart}
// 
// !INTERFACE:
void FTN(writebiasrestart)(char *j, int *n, int *k,int len)
//  
// !DESCRIPTION: 
// Invokes the the selected bias estimation algorithm
// for writing bias restart
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the bias estimation scheme
//  \item[n]
//   index of the nest
//  \item[k]
//   index of the assimilation instance
//  \end{description}
//
//EOP
{ 
  struct biasoutnode* current;
  
  current = biasout_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("restart routine for bias estimation scheme %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,k); 
}


//BOP
// !ROUTINE: registerbiasestimationfinalize
// \label{registerrbiasestimationfinalize}
//  
// 
// !INTERFACE:
void FTN(registerbiasestimationfinalize)(char *j, void (*func)(),int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine that
//  finalizes structures related to bias correction algorithm. 
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the bias estimation algorithm. 
//  \end{description}
//EOP
{ 
  int len1;
  struct biasfinalnode* current;
  struct biasfinalnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct biasfinalnode*) malloc(sizeof(struct biasfinalnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(biasfinal_table == NULL){
    biasfinal_table = pnode;
  }
  else{
    current = biasfinal_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: biasestimationfinalize
// \label{biasestimationfinalize}
// 
// !INTERFACE:
void FTN(biasestimationfinalize)(char *j,int len)
//  
// !DESCRIPTION: 
// Finalizes the the selected bias estimation algorithm 
// related structures
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the bias estimation algorithm
//  \end{description}
//
//EOP
{ 
  struct biasfinalnode* current;
  
  current = biasfinal_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("finalize routine for bias estimation scheme %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}
