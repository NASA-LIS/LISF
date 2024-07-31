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
//  !MODULE: LDT_irrig_FTable
//  
//
// !DESCRIPTION:
//   Function table registries for storing the interface 
//   implementations for managing different sources of 
//   irrigation data
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"

struct irrigtypenode
{ 
  char *name;
  void (*func)(int*, float*, int*);

  struct irrigtypenode* next;
} ;

struct irrigtypenode* irrigtype_table = NULL; 

struct irrigfracnode
{ 
  char *name;
  void (*func)(int*, float*);

  struct irrigfracnode* next;
} ;

struct irrigfracnode* irrigfrac_table = NULL; 

struct irriggwrationode
{
  char *name;
  void (*func)(int*, float*);

  struct irriggwrationode* next;
} ;
struct irriggwrationode* irriggwratio_table = NULL;

//BOP
// !ROUTINE: registerreadirrigtype
// \label{registerreadirrigtype}
// 
// !INTERFACE:
void FTN(registerreadirrigtype)(char *j, void (*func)(int*, float*, int*), int len)
//  
// !DESCRIPTION:
//  Makes an entry in the registry for the routine to read
//  irrigation data
//
//  The arguments are: 
//  \begin{description}
//   \item[j] 
//    index of the irrigation source
//  \end{description}
  //EOP
{ 
  int len1;
  struct irrigtypenode* current;
  struct irrigtypenode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct irrigtypenode*) malloc(sizeof(struct irrigtypenode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(irrigtype_table == NULL){
    irrigtype_table = pnode;
  }
  else{
    current = irrigtype_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readirrigtype
// \label{readirrigtype}
// 
// !INTERFACE:
void FTN(readirrigtype)(char *j, int *n, float *array, int *nt, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to read the 
//  irrigation data
//
//  The arguments are: 
//  \begin{description}
//   \item[j] 
//    index of the irrigation type source
//   \item[n]
//    index of the nest
//   \item[array]
//    pointer to the irrigation type array
//  \end{description}
//EOP
{ 
  struct irrigtypenode* current;
  
  current = irrigtype_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Irrigation reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,array,nt); 
}

//BOP
// !ROUTINE: registerreadirrigfrac
// \label{registerreadirrigfrac}
//
// !INTERFACE:
void FTN(registerreadirrigfrac)(char *j, void (*func)(int*, float*), int len)
//  
// !DESCRIPTION:
//  Makes an entry in the registry for the routine to read
//  irrigfrac data
// 
//  The arguments are: 
//  \begin{description}
//   \item[j] 
//    index of the irrigation source
//   \item[n]
//    index of the nest
//   \item[array]
//    pointer to the irrigation fraction array
//  \end{description}
  //EOP
{ 

  int len1;
  struct irrigfracnode* current;
  struct irrigfracnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct irrigfracnode*) malloc(sizeof(struct irrigfracnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(irrigfrac_table == NULL){
    irrigfrac_table = pnode;
  }
  else{
    current = irrigfrac_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readirrigfrac
// \label{readirrigfrac}
// 
// !INTERFACE:
void FTN(readirrigfrac)(char *j, int *n, float *array, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to read the 
//  irrigfrac data
//
//  The arguments are: 
//  \begin{description}
//   \item[j] 
//    index of the irrigation fraction source
//   \item[n]
//    index of the nest
//   \item[array]
//    pointer to the irrigation fraction array
//  \end{description}
//EOP
{ 
  struct irrigfracnode* current;
  
  current = irrigfrac_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Irrigation fraction reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  
  current->func(n,array); 
}

//BOP
// !ROUTINE: registerreadirrigGWratio
// \label{registerreadirrigGWratio}
// 
// !INTERFACE:
void FTN(registerreadirriggwratio)(char *j,void (*func)(int*,float*),int len)
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read groundwater irrigation ratio data
//
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \item[j]
//  index of the groundwater irrigation ratio data source
//  \end{description}
//EOP
{
  struct irriggwrationode* current;
  struct irriggwrationode* pnode;
  // create node

  pnode=(struct irriggwrationode*) malloc(sizeof(struct irriggwrationode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL;

  if(irriggwratio_table == NULL){
    irriggwratio_table = pnode;
  }
  else{
    current = irriggwratio_table;
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode;
  }
}

//BOP
// !ROUTINE: readirrigGWratio
// \label{readirrigGWratio}
//
// !INTERFACE:
void FTN(readirriggwratio)(char *j, int *n,float *array,int len)
//  
// !DESCRIPTION:
// Invokes the routine from the registry to 
// reading groundwater irrigation ratio data 
// 
// The arguments are: 
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//  index of the groundwater irrigation ratio data source
//  \item[array]
//  pointer to the groundwater irrigation data
//  \end{description}
//EOP
{
  struct irriggwrationode* current;

  current = irriggwratio_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n");
      printf("irrigGWratio reading routine for source %s is not defined\n",j);
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n");
      printf("****************Error****************************\n");
    }
  }
  current->func(n,array);
}
