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
// !MODULE: LDT_anndata_FTable
//  
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing the operation of different 
//  domain specifications
// 
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"

struct anninpsetupnode
{ 
  char *name;
  void (*func)();

  struct anninpsetupnode* next;
} ;
struct anninpsetupnode* anninpsetup_table = NULL; 

struct annoutsetupnode
{ 
  char *name;
  void (*func)();

  struct annoutsetupnode* next;
} ;
struct annoutsetupnode* annoutsetup_table = NULL; 

struct anninputsourcenode
{ 
  char *name;
  void (*func)(int*, int*, int*, int*);

  struct anninputsourcenode* next;
} ;
struct anninputsourcenode* anninputsource_table = NULL; 

struct annoutputsourcenode
{ 
  char *name;
  void (*func)(int*, int*, int*, int*);

  struct annoutputsourcenode* next;
} ;
struct annoutputsourcenode* annoutputsource_table = NULL; 

//BOP
// !ROUTINE: registeranninputsourcesetup
// \label{registeranninputsourcesetup}
//
// !INTERFACE:
void FTN(registeranninputsourcesetup)(char *j,void (*func)(),int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read the runtime domain specifics
// 
// The arguments are: 
// \begin{description}
// \item[j]
//  name of the observation source
//  \end{description}
//EOP
{ 
  int len1;
  struct anninpsetupnode* current;
  struct anninpsetupnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct anninpsetupnode*) malloc(sizeof(struct anninpsetupnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(anninpsetup_table == NULL){
    anninpsetup_table = pnode;
  }
  else{
    current = anninpsetup_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: setupanninputsource
// \label{setupanninputsource}
//
// !INTERFACE:
void FTN(setupanninputsource)(char *j, int len)
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
  struct anninpsetupnode* current;
  
  current = anninpsetup_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Observation setup routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}


//BOP
// !ROUTINE: registerannoutputsourcesetup
// \label{registerannoutputsourcesetup}
//
// !INTERFACE:
void FTN(registerannoutputsourcesetup)(char *j,void (*func)(),int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read the runtime domain specifics
// 
// The arguments are: 
// \begin{description}
// \item[j]
//  name of the observation source
//  \end{description}
//EOP
{ 
  int len1;
  struct annoutsetupnode* current;
  struct annoutsetupnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct annoutsetupnode*) malloc(sizeof(struct annoutsetupnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(annoutsetup_table == NULL){
    annoutsetup_table = pnode;
  }
  else{
    current = annoutsetup_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: setupannoutputsource
// \label{setupannoutputsource}
//
// !INTERFACE:
void FTN(setupannoutputsource)(char *j, int len)
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
  struct annoutsetupnode* current;
  
  current = annoutsetup_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Observation setup routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}


//BOP
// !ROUTINE: registerreadanninputsource
// \label{registerreadanninputsource}
//
// !INTERFACE:
void FTN(registerreadanninputsource)(char *j,void (*func)(int*, int*, int*, int*), int len)
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
  struct anninputsourcenode* current;
  struct anninputsourcenode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct anninputsourcenode*) malloc(sizeof(struct anninputsourcenode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(anninputsource_table == NULL){
    anninputsource_table = pnode;
  }
  else{
    current = anninputsource_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: readanninputsource
// \label{readanninputsource}
//
// !INTERFACE:
void FTN(readanninputsource)(char *j, int *n, int *iomode, int *s, int *e,int len)
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
  struct anninputsourcenode* current;
  
  current = anninputsource_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Read observation source routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,iomode, s,e); 
}


//BOP
// !ROUTINE: registerreadannoutputsource
// \label{registerreadannoutputsource}
//
// !INTERFACE:
void FTN(registerreadannoutputsource)(char *j,void (*func)(int*, int*, int*, int*), int len)
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
  struct annoutputsourcenode* current;
  struct annoutputsourcenode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct annoutputsourcenode*) malloc(sizeof(struct annoutputsourcenode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(annoutputsource_table == NULL){
    annoutputsource_table = pnode;
  }
  else{
    current = annoutputsource_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: readannoutputsource
// \label{readannoutputsource}
//
// !INTERFACE:
void FTN(readannoutputsource)(char *j, int *n, int *iomode, int *s, int *e, int len)
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
  struct annoutputsourcenode* current;
  
  current = annoutputsource_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Read observation source routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,iomode,s, e); 
}





