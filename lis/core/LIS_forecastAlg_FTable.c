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
// !MODULE: LIS_forecastAlg_FTable
//  
//
// !DESCRIPTION:
//  
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h" 

struct forecastalginitnode
{ 
  char *name;
  void (*func)();

  struct forecastalginitnode* next;
} ;
struct forecastalginitnode* forecastalg_init = NULL; 

struct forecastalgsamplenode
{ 
  char *name;
  void (*func)(int*, int*, int*, int*, int*, int*);

  struct forecastalgsamplenode* next;
} ;
struct forecastalgsamplenode* forecastalg_sample = NULL; 


//BOP
// !ROUTINE: registerforecastalginit
// \label{registerforecastalginit}
// 
// !INTERFACE:
void FTN(registerforecastalginit)(char *j, void (*func)(), int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  initialize the radiative transfer model
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the RTM
//. \end{description}
//EOP
{ 
  int len1;
  struct forecastalginitnode* current;
  struct forecastalginitnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct forecastalginitnode*) malloc(sizeof(struct forecastalginitnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(forecastalg_init == NULL){
    forecastalg_init = pnode;
  }
  else{
    current = forecastalg_init; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: forecastalginit
// \label{forecastalginit}
// 
// !INTERFACE:
void FTN(forecastalginit)(char *j, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to
//  initialize the RTM
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the RTM
// \end{description}
//EOP
{ 
  struct forecastalginitnode* current;
  
  current = forecastalg_init;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("init routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}


//BOP
// !ROUTINE: registerforecastsampledate
// \label{registerforecastsampledate}
// 
// !INTERFACE:
void FTN(registerforecastsampledate)(char *j, void (*func)(int*, int*, int*, int*, int*, int*), int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  initialize the radiative transfer model
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the RTM
//. \end{description}
//EOP
{ 
  int len1;
  struct forecastalgsamplenode* current;
  struct forecastalgsamplenode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct forecastalgsamplenode*) malloc(sizeof(struct forecastalgsamplenode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(forecastalg_sample == NULL){
    forecastalg_sample = pnode;
  }
  else{
    current = forecastalg_sample; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: forecastalgsampledate
// \label{forecastalgsampledate}
// 
// !INTERFACE:
void FTN(forecastalgsampledate)(char *j, int *n, int *kk, int *k, int *yr, int *mo, int* da,int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to
//  sampleialize the RTM
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the RTM
// \end{description}
//EOP
{ 
  struct forecastalgsamplenode* current;
  
  current = forecastalg_sample;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("sample date routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,kk,k,yr,mo,da); 
}

