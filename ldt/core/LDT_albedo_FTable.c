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
// !MODULE: LDT_albedo_FTable
//  
//
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing different sources of 
//  albedo fraction data
//   
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include <string.h>

#include "ftn_drv.h"


struct albedosetnode
{ 
  char *name;
  void (*func)();

  struct albedosetnode* next;
} ;
struct albedosetnode* albedoset_table = NULL; 

struct albedonode
{ 
  char *name;
  void (*func)(int*, float*);

  struct albedonode* next;
} ;
struct albedonode* albedo_table = NULL; 

struct mxsnoalbnode
{ 
  char *name;
  void (*func)(int*, float*);

  struct mxsnoalbnode* next;
} ;
struct mxsnoalbnode* mxsnoalb_table = NULL; 


//BOP
// !ROUTINE: registersetalbedoattribs
// \label{registersetalbedoattribs}
//  
// 
// !INTERFACE:
void FTN(registersetalbedoattribs)(char *j, void (*func)(),int len)
// !DESCRIPTION:
//  Creates an entry in the registry for the routine to 
//  read albedo data
// 
//  The arguments are: 
//  \begin{description}
//   \item[j]
//    index of the albedo source
//   \end{description}
//EOP
{ 
  int len1;
  struct albedosetnode* current;
  struct albedosetnode* pnode; 
  // create node

  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct albedosetnode*) malloc(sizeof(struct albedosetnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(albedoset_table == NULL){
    albedoset_table = pnode;
  }
  else{
    current = albedoset_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: setalbedoattribs
// \label{setalbedoattribs}
//  
// !DESCRIPTION: 
// Invokes the routine from the registry for 
// reading albedo data. 
//
//  The arguments are: 
//  \begin{description}
//   \item[n]
//    index of the nest
//   \item[j]
//    index of the albedo source
//   \item[array]
//    pointer to the albedo data
//   \item[marray]
//    pointer to the mask data
//  \end{description}
//
// !INTERFACE:
void FTN(setalbedoattribs)(char *j,int len)
//EOP
{ 

  struct albedosetnode* current;
  
  current = albedoset_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;

    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("setAlbedoAttribs routine for source %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func();
}


//BOP
// !ROUTINE: registerreadalbedo
// \label{registerreadalbedo}
// 
// !INTERFACE:
void FTN(registerreadalbedo)(char *j,void (*func)(int*,float*), int len )
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read albedo data
//
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \item[j]
//  index of the albedo data source
//  \end{description}
//EOP
{ 
  int len1;
  struct albedonode* current;
  struct albedonode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct albedonode*) malloc(sizeof(struct albedonode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(albedo_table == NULL){
    albedo_table = pnode;
  }
  else{
    current = albedo_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}

//BOP
// !ROUTINE: readalbedo
// \label{readalbedo}
//
// !INTERFACE:
void FTN(readalbedo)(char *j, int *n,float *array, int len)
//  
// !DESCRIPTION:
// Invokes the routine from the registry to 
// reading albedo data 
// 
// The arguments are: 
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//  index of the albedo data source
//  \item[array]
//  pointer to the albedo data
//  \end{description}
//EOP
{ 
  struct albedonode* current;
  
  current = albedo_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Albedo reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,array); 
}


//BOP
// !ROUTINE: registerreadmxsnoalb
// \label{registerreadmxsnoalb}
// 
// !INTERFACE:
void FTN(registerreadmxsnoalb)(char *j,void (*func)(int*,float*),int len)
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read mxsnoalb data
//
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \item[j]
//  index of the albedo data source
//  \end{description}
//EOP
{ 
  int len1;
  struct mxsnoalbnode* current;
  struct mxsnoalbnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct mxsnoalbnode*) malloc(sizeof(struct mxsnoalbnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(mxsnoalb_table == NULL){
    mxsnoalb_table = pnode;
  }
  else{
    current = mxsnoalb_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readmxsnoalb
// \label{readmxsnoalb}
//
// !INTERFACE:
void FTN(readmxsnoalb)(char *j, int *n,float *array, int len)
//  
// !DESCRIPTION:
// Invokes the routine from the registry to 
// reading mxsnoalb data 
// 
// The arguments are: 
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//  index of the max snow albedo data source
//  \item[array]
//  pointer to the max snow albedo data
//  \end{description}
//EOP
{ 
  struct mxsnoalbnode* current;
  
  current = mxsnoalb_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Max snow albedo reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,array); 
}


