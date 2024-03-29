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
// !MODULE: LVT_domain_FTable
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

struct domaininputnode
{ 
  char *name;
  void (*func)();

  struct domaininputnode* next;
} ;
struct domaininputnode* domaininput_table = NULL; 

struct domainmakenode
{ 
  char *name;
  void (*func)();

  struct domainmakenode* next;
} ;
struct domainmakenode* domainmake_table = NULL; 


//BOP
// !ROUTINE: registerinput
// \label{registerinput}
//
// !INTERFACE:
void FTN(registerinput)(char *j,void (*func)(), int len)
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
  struct domaininputnode* pnode;
  struct domaininputnode* current;

  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct domaininputnode*) malloc(sizeof(struct domaininputnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(domaininput_table == NULL){
    domaininput_table = pnode;
  }
  else{
    current = domaininput_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: readinput
// \label{readinput}
//
// !INTERFACE:
void FTN(readinput)(char *j, int len)
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
  struct domaininputnode* current;
  
  current = domaininput_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("read_input routine for domain %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}

//BOP
// !ROUTINE: registerdomain
// \label{registerdomain}
//  
// 
// !INTERFACE:
void FTN(registerdomain)(char *j,void (*func)(),int len)
// !DESCRIPTION: 
// Creates an entry in the registry for the routine 
// to create grid and tile spaces for a particular 
// domain definition
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//EOP
{ 
  int len1;
  struct domainmakenode* pnode;
  struct domainmakenode* current;

  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct domainmakenode*) malloc(sizeof(struct domainmakenode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(domainmake_table == NULL){
    domainmake_table = pnode;
  }
  else{
    current = domainmake_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: makedomain
// \label{makedomain}
//
// !INTERFACE:
void FTN(makedomain)(char *j,int len)
//  
// !DESCRIPTION: 
// Invokes the routine from the registry 
// to create grid and tile spaces for a particular 
// domain definition
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//EOP
{ 

  struct domainmakenode* current;
  
  current = domainmake_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("make domain routine for domain %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}




