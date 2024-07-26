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
// !MODULE: LDT_landcover_FTable
//  
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing different sources of 
//  landcover data
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"

struct lcsetnode
{ 
  char *name;
  void (*func)();

  struct lcsetnode* next;
} ;
struct lcsetnode* lcset_table = NULL; 

struct lcnode
{ 
  char *name;
  void (*func)(int*, int*, float*, float*);

  struct lcnode* next;
} ;

struct lcnode* lc_table = NULL; 

struct croplcnode
{
  char *name;
  void (*func)(int*, int*, float*);

  struct croplcnode* next;
} ;

struct croplcnode* croplc_table = NULL;


struct drootnode
{
  char *name;
  void (*func)(int*, float*);

  struct drootnode* next;
} ;

struct drootnode* droot_table = NULL;


//BOP
// !ROUTINE: registerrreadlc
// \label{registerreadlc}
//  
// 
// !INTERFACE:
void FTN(registerreadlc)(char *j, void (*func)(int*, int*, float*, float*),int len)
// !DESCRIPTION:
//  Creates an entry in the registry for the routine to 
//  read landcover data
// 
//  The arguments are: 
//  \begin{description}
//   \item[j]
//    index of the landcover source
//   \end{description}
//EOP
{ 
  int len1;
  struct lcnode* current;
  struct lcnode* pnode; 
  // create node

  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lcnode*) malloc(sizeof(struct lcnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lc_table == NULL){
    lc_table = pnode;
  }
  else{
    current = lc_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readlandcover
// \label{readlandcover}
//  
// !DESCRIPTION: 
// Invokes the routine from the registry for 
// reading landcover data. 
//
//  The arguments are: 
//  \begin{description}
//   \item[n]
//    index of the nest
//   \item[j]
//    index of the landcover source
//   \item[array]
//    pointer to the landcover data
//   \item[marray]
//    pointer to the mask data
//  \end{description}
//
// !INTERFACE:
void FTN(readlandcover)(char *j,int *n, int *nt, float *array, float *marray, int len)
//EOP
{ 

  struct lcnode* current;
  
  current = lc_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;

    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Landcover reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,nt,array,marray);
}


//BOP
// !ROUTINE: registerreadcroptype
// \label{registerreadcroptype}
// 
// !INTERFACE:
void FTN(registerreadcroptype)(char *j, void (*func)(int*, int*, float*),int len)
//  
// !DESCRIPTION:
//  Creates an entry in the registry for the routine to 
//  read the crop type data
// 
//  The arguments are: 
//  \begin{description}
//  \item[j]
//   index of the crop type source
//  \end{description}
  //EOP
{
  int len1;
  struct croplcnode* current;
  struct croplcnode* pnode;
  // create node

  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct croplcnode*) malloc(sizeof(struct croplcnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL;

  if(croplc_table == NULL){
    croplc_table = pnode;
  }
  else{
    current = croplc_table;
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode;
  }

}

//BOP
// !ROUTINE: readcroptype
// \label{readcroptype}
// 
// !INTERFACE:
void FTN(readcroptype)(char *j,int *n,int *nct,float *array, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to read the
//  crop type data
// 
//  The arguments are: 
//  \begin{description}
//  \item[n]
//   index of the nest
//  \item[nct]
//   number of crop types
//  \item[array]
//   pointer to the crop type array
//  \end{description}
//EOP
{
  struct croplcnode* current;

  current = croplc_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("**************** Error ****************************\n");
      printf("Crop type reading routine for source %s is not defined\n",j);
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("Program will seg fault.....\n");
      printf("**************** Error ****************************\n");
    }
  }
  current->func(n,nct,array);
}

//BOP
// !ROUTINE: registerreadrootdepth
// \label{registerreadrootdepth}
// 
// !INTERFACE:
void FTN(registerreadrootdepth)(char *j,void (*func)(int*,float*),int len)
//  
// !DESCRIPTION:
//  Creates an entry in the registry for the routine to 
//  read the rootdepth fraction data
// 
//  The arguments are: 
//  \begin{description}
//  \item[j]
//   index of the root depth source
//  \end{description}
  //EOP
{
  int len1;
  struct drootnode* current;
  struct drootnode* pnode;
  // create node

  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct drootnode*) malloc(sizeof(struct drootnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL;

  if(droot_table == NULL){
    droot_table = pnode;
  }
  else{
    current = droot_table;
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode;
  }
}

//BOP
// !ROUTINE: readrootdepth
// \label{readrootdepth}
// 
// !INTERFACE:
void FTN(readrootdepth)(char *j,int *n,float *rtdarray, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to read the
//  rootdepth fraction data
// 
//  The arguments are: 
//  \begin{description}
//  \item[n]
//   index of the nest
//  \item[array]
//   pointer to the rootdepth data
//  \end{description}
//EOP
{

  struct drootnode* current;

  current = droot_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("**************** Error ***************************\n");
      printf("Rootdepth reading routine for source %s is not defined\n",j);
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("Program will seg fault.....\n");
      printf("**************** Error ***************************\n");
    }
  }
  current->func(n,rtdarray);
}

