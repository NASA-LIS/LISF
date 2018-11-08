//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center Land Information System (LIS) v7.2
//
// Copyright (c) 2015 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
//BOP
//
// !MODULE: LIS_roughness_FTable
//  
//
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing different sources of 
//  greenness fraction data
//   
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include <string.h>

#include "ftn_drv.h"

struct roughnesssetnode
{ 
  char *name;
  void (*func)(int*);

  struct roughnesssetnode* next;
} ;
struct roughnesssetnode* roughnessset_table = NULL; 

struct roughnessreadnode
{ 
  char *name;
  void (*func)(int*,void*,void*, float*, float*);

  struct roughnessreadnode* next;
} ;
struct roughnessreadnode* roughnessread_table = NULL; 

//BOP
// !ROUTINE: registerroughnesssetup
// \label{registerroughnesssetup}
// 
// !INTERFACE:
void FTN(registerroughnesssetup)(char *j,void (*func)(int*),int len)
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// setup roughness data reading routines
//
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \item[j]
//  index of the greenness data source
//  \end{description}
//EOP
{ 
  struct roughnesssetnode* current;
  struct roughnesssetnode* pnode; 
  // create node
  
  pnode=(struct roughnesssetnode*) malloc(sizeof(struct roughnesssetnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(roughnessset_table == NULL){
    roughnessset_table = pnode;
  }
  else{
    current = roughnessset_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}

//BOP
// !ROUTINE: roughnesssetup
// \label{roughnesssetup}
//
// !INTERFACE:
void FTN(roughnesssetup)(char *j,int *n, int len)
//  
// !DESCRIPTION:
// Invokes the routine from the registry to 
// setup roughness data reading
// 
// The arguments are: 
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//  index of the greenness data source
//  \end{description}
//EOP
{ 
  struct roughnesssetnode* current;
  
  current = roughnessset_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("roughnesssetup routine for runmode %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}

//BOP
// !ROUTINE: registerreadroughness
// \label{registerreadroughness}
// 
// !INTERFACE:
void FTN(registerreadroughness)(char *j,void (*func)(int*, void*, void*, float*, float*),int len)
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read roughness data
//
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \item[j]
//  index of the greenness data source
//  \end{description}
//EOP
{ 
  struct roughnessreadnode* current;
  struct roughnessreadnode* pnode; 
  // create node
  
  pnode=(struct roughnessreadnode*) malloc(sizeof(struct roughnessreadnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(roughnessread_table == NULL){
    roughnessread_table = pnode;
  }
  else{
    current = roughnessread_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}

//BOP
// !ROUTINE: readroughness
// \label{readroughness}
//
// !INTERFACE:
void FTN(readroughness)(char *j,int *n, void *time1, void *time2, float *array1, float *array2, int len)
//  
// !DESCRIPTION:
// Invokes the routine from the registry to 
// reading roughness data 
// 
// The arguments are: 
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//  index of the greenness data source
//  \item[time]
//  month 
//  \item[array]
//  pointer to the greenness data
//  \end{description}
//EOP
{ 
  struct roughnessreadnode* current;
  
  current = roughnessread_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("readroughness routine for runmode %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,time1,time2, array1,array2); 
}



