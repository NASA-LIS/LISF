//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA GSFC Land surface Verification Toolkit (LVT) V1.0
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
//BOP
//
// !MODULE: LVT_obs_FTable
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

struct obsininode
{ 
  char *name;
  void (*func)(int*);

  struct obsininode* next;
} ;
struct obsininode* obsini_table = NULL; 

struct obsreadnode
{ 
  char *name;
  void (*func)(int*);

  struct obsreadnode* next;
} ;

struct obsreadnode* obsread_table = NULL; 

//BOP
// !ROUTINE: registerobssetup
// \label{registerobssetup}
//
// !INTERFACE:
void FTN(registerobssetup)(char *j,void (*func)(int*), int len)
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
  struct obsininode* current;
  struct obsininode* pnode; 
  // create node
  
  pnode=(struct obsininode*) malloc(sizeof(struct obsininode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(obsini_table == NULL){
    obsini_table = pnode;
  }
  else{
    current = obsini_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: observationsetup
// \label{observationsetup}
//
// !INTERFACE:
void FTN(observationsetup)(char *j,  int *source, int len)
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
  struct obsininode* current;
  
  current = obsini_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("observation init routine for source %s is not defined\n",j); 
      printf("please see the configs/lvt.config.master file.....\n"); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(source); 
}


//BOP
// !ROUTINE: registerobsread
// \label{registerobsread}
//
// !INTERFACE:
void FTN(registerobsread)(char *j,void (*func)(int*), int len)
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
  struct obsreadnode* current;
  struct obsreadnode* pnode; 
  // create node
  
  pnode=(struct obsreadnode*) malloc(sizeof(struct obsreadnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(obsread_table == NULL){
    obsread_table = pnode;
  }
  else{
    current = obsread_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: readobservationsource
// \label{readobservationsource}
//
// !INTERFACE:
void FTN(readobservationsource)(char *j, int *source, int len)
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

  struct obsreadnode* current;
  
  current = obsread_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("observation read routine for source %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(source); 
}





