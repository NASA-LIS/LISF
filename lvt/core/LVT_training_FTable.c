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
// !MODULE: LVT_training_FTable
//  
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing the operation of different 
//  training algorithms
// 
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"

struct trainalginitnode
{ 
  char *name;
  void (*func)();

  struct trainalginitnode* next;
} ;
struct trainalginitnode* trainalginit_table = NULL; 

struct trainalgrunnode
{ 
  char *name;
  void (*func)(int*);

  struct trainalgrunnode* next;
} ;

struct trainalgrunnode* trainalgrun_table = NULL; 

//BOP
// !ROUTINE: registertraininginit
// \label{registertraininginit}
//
// !INTERFACE:
void FTN(registertraininginit)(char *j,void (*func)(),int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// perform the initialization step of the training algorithm
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//EOP
{ 
  int len1;
  struct trainalginitnode* pnode;
  struct trainalginitnode* current;

  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct trainalginitnode*) malloc(sizeof(struct trainalginitnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(trainalginit_table == NULL){
    trainalginit_table = pnode;
  }
  else{
    current = trainalginit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: traininginit
// \label{traininginit}
//
// !INTERFACE:
void FTN(traininginit)(char *j,int len)
// !DESCRIPTION: 
//  Calls the routine from the registry to invoke
//  the initialization step of the training algorithm
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//
//EOP
{ 
  struct trainalginitnode* current;
  
  current = trainalginit_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("init routine for training %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}

//BOP
// !ROUTINE: registertrainingrun
// \label{registertrainingrun}
//
// !INTERFACE:
void FTN(registertrainingrun)(char *j,void (*func)(int*), int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// perform the run step of the training algorithm
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//EOP
{ 
  int len1;
  struct trainalgrunnode* pnode;
  struct trainalgrunnode* current;

  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct trainalgrunnode*) malloc(sizeof(struct trainalgrunnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(trainalgrun_table == NULL){
    trainalgrun_table = pnode;
  }
  else{
    current = trainalgrun_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: trainingrun
// \label{trainingrun}
//
// !INTERFACE:
void FTN(trainingrun)(char *j, int *pass, int len)
// !DESCRIPTION: 
//  Calls the routine from the registry to invoke the 
//  the run step of the training algorithm. 
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//
//EOP
{ 

  struct trainalgrunnode* current;
  
  current = trainalgrun_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("run routine for training %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(pass); 
}




