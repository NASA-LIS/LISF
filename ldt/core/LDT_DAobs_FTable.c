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
// !MODULE: LDT_daobs_FTable
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

struct daobsinitnode
{ 
  char *name;
  void (*func)();

  struct daobsinitnode* next;
} ;
struct daobsinitnode* daobsinit_table = NULL; 

struct daobsreadnode
{ 
  char *name;
  void (*func)(int*);

  struct daobsreadnode* next;
} ;
struct daobsreadnode* daobsread_table = NULL; 

//BOP
// !ROUTINE: registerdaobssetup
// \label{registerdaobssetup}
//
// !INTERFACE:
void FTN(registerdaobssetup)(char *j,void (*func)(),int len)
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
  struct daobsinitnode* current;
  struct daobsinitnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct daobsinitnode*) malloc(sizeof(struct daobsinitnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(daobsinit_table == NULL){
    daobsinit_table = pnode;
  }
  else{
    current = daobsinit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: daobservationsetup
// \label{daobservationsetup}
//
// !INTERFACE:
void FTN(daobservationsetup)(char *j, int len)
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
  struct daobsinitnode* current;
  
  current = daobsinit_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Observation setup routine for %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}


//BOP
// !ROUTINE: registerdaobsread
// \label{registerdaobsread}
//
// !INTERFACE:
void FTN(registerdaobsread)(char *j,void (*func)(int*), int len)
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
  struct daobsreadnode* current;
  struct daobsreadnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct daobsreadnode*) malloc(sizeof(struct daobsreadnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(daobsread_table == NULL){
    daobsread_table = pnode;
  }
  else{
    current = daobsread_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: readdaobservationsource
// \label{readdaobservationsource}
//
// !INTERFACE:
void FTN(readdaobservationsource)(char *j, int *n, int len)
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
  struct daobsreadnode* current;
  
  current = daobsread_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Read observation source routine for %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}





