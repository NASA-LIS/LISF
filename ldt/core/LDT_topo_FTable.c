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
//  !MODULE: LDT_topo_FTable
//  
//
// !DESCRIPTION:
//   Function table registries for storing the interface 
//   implementations for managing different sources of 
//   topography data
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"

struct elevnode
{ 
  char *name;
  void (*func)(int*, int*, float*, float*);

  struct elevnode* next;
} ;

struct elevnode* elev_table = NULL; 

struct slopenode
{ 
  char *name;
  void (*func)(int*, int*, float*, float*);

  struct slopenode* next;
} ;

struct slopenode* slope_table = NULL; 

struct aspectnode
{ 
  char *name;
  void (*func)(int*, int*, float*, float*);

  struct aspectnode* next;
} ;

struct aspectnode* aspect_table = NULL; 

struct curvnode
{ 
  char *name;
  void (*func)(int*, float*);

  struct curvnode* next;
} ;

struct curvnode* curv_table = NULL; 


//BOP
// !ROUTINE: registerreadelev
// \label{registerreadelev}
// 
// !INTERFACE:
void FTN(registerreadelev)(char *j, void (*func)(int*, int*, float*, float*), int len)
//  
// !DESCRIPTION:
//  Makes an entry in the registry for the routine to read
//  elevation data
//
//  The arguments are: 
//  \begin{description}
//   \item[j] 
//    index of the topography source
//  \end{description}
  //EOP
{ 
  int len1;
  struct elevnode* current;
  struct elevnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct elevnode*) malloc(sizeof(struct elevnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(elev_table == NULL){
    elev_table = pnode;
  }
  else{
    current = elev_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readelev
// \label{readelev}
// 
// !INTERFACE:
void FTN(readelev)(char *j, int *n, int *nt, float *fgrd, float *array, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to read the 
//  elevation data
//
//  The arguments are: 
//  \begin{description}
//   \item[j] 
//    index of the topography source
//   \item[n]
//    index of the nest
//   \item[nt]
//    number of types or bins (bands)
//   \item[fgrd]
//    gridcell fraction of type or values
//   \item[array]
//    pointer to the elevation array
//  \end{description}
//EOP
{ 
  struct elevnode* current;
  
  current = elev_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Elevation reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,nt,fgrd,array); 
}

//BOP
// !ROUTINE: registerreadslope
// \label{registerreadslope}
//
// !INTERFACE:
void FTN(registerreadslope)(char *j, void (*func)(int*, int*, float*, float*), int len)
//  
// !DESCRIPTION:
//  Makes an entry in the registry for the routine to read
//  slope data
// 
//  The arguments are: 
//  \begin{description}
//   \item[j] 
//    index of the topography source
//   \item[n]
//    index of the nest
//   \item[nt]
//    number of types or bins (bands)
//   \item[fgrd]
//    gridcell fraction of type or values
//   \item[array]
//    pointer to the slope array
//  \end{description}
  //EOP
{ 
  int len1;
  struct slopenode* current;
  struct slopenode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct slopenode*) malloc(sizeof(struct slopenode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(slope_table == NULL){
    slope_table = pnode;
  }
  else{
    current = slope_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readslope
// \label{readslope}
// 
// !INTERFACE:
void FTN(readslope)(char *j, int *n, int *nt, float *fgrd, float *array, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to read the 
//  slope data
//
//  The arguments are: 
//  \begin{description}
//   \item[j] 
//    index of the topography source
//   \item[n]
//    index of the nest
//   \item[nt]
//    number of types or bins (bands)
//   \item[fgrd]
//    gridcell fraction of type or values
//   \item[array]
//    pointer to the slope array
//  \end{description}
//EOP
{ 
  struct slopenode* current;
  
  current = slope_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Slope reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  
  current->func(n,nt,fgrd,array); 
}


//BOP
// !ROUTINE: registerreadaspect
// \label{registerreadaspect}
// 
// !INTERFACE:
void FTN(registerreadaspect)(char *j, void (*func)(int*, int*, float*, float*), int len)
//  
// !DESCRIPTION:
//  Makes an entry in the registry for the routine to read
//  aspect data
//
//  The arguments are: 
//  \begin{description}
//   \item[j] 
//    index of the topography source
//   \item[n]
//    index of the nest
//   \item[nt]
//    number of types or bins (bands)
//   \item[fgrd]
//    gridcell fraction of type or values
//   \item[array]
//    pointer to the aspect array
//  \end{description}
  //EOP
{ 
  int len1;
  struct aspectnode* current;
  struct aspectnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct aspectnode*) malloc(sizeof(struct aspectnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(aspect_table == NULL){
    aspect_table = pnode;
  }
  else{
    current = aspect_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readaspect
// \label{readaspect}
// 
// !INTERFACE:
void FTN(readaspect)(char *j, int *n, int *nt, float *fgrd, float *array, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to read the 
//  aspect data
//
//  The arguments are: 
//  \begin{description}
//   \item[j] 
//    index of the topography source
//   \item[n]
//    index of the nest
//   \item[nt]
//    number of types or bins (bands)
//   \item[fgrd]
//    gridcell fraction of type or values
//   \item[array]
//    pointer to the aspect array
//  \end{description}
//EOP
{ 
  struct aspectnode* current;
  
  current = aspect_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Aspect reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,nt,fgrd,array); 
}

//BOP
// !ROUTINE: registerreadcurv
// \label{registerreadcurv}
// 
// !INTERFACE:
void FTN(registerreadcurv)(char *j, void (*func)(int*, float*), int len)
//  
// !DESCRIPTION:
//  Makes an entry in the registry for the routine to read
//  curvature data
//
//  The arguments are: 
//  \begin{description}
//   \item[j] 
//    index of the topography source
//  \end{description}
  //EOP
{ 
  int len1;
  struct curvnode* current;
  struct curvnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct curvnode*) malloc(sizeof(struct curvnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(curv_table == NULL){
    curv_table = pnode;
  }
  else{
    current = curv_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}

//BOP
// !ROUTINE: readcurv
// \label{readcurv}
// 
// !INTERFACE:
void FTN(readcurv)(char *j, int *n, float *array, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to read the 
//  curvature data
//
//  The arguments are: 
//  \begin{description}
//   \item[n]
//    index of the nest
//   \item[j] 
//    index of the topography source
//   \item[array]
//    pointer to the curvature array
//  \end{description}
//EOP
{ 
  struct curvnode* current;
  
  current = curv_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Curvature reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,array); 
}




