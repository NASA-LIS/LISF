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
// !MODULE: LIS_landslide_FTable
//  
//
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing the operations of different 
//  land slide models. 
//
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"

struct landslideinitnode
{ 
  char *name;
  void (*func)();

  struct landslideinitnode* next;
} ;
struct landslideinitnode* landslideinit_table = NULL; 

struct landsliderunnode
{ 
  char *name;
  void (*func)(int*);

  struct landsliderunnode* next;
} ;
struct landsliderunnode* landsliderun_table = NULL; 

struct landslideoutnode
{ 
  char *name;
  void (*func)(int*);

  struct landslideoutnode* next;
} ;
struct landslideoutnode* landslideout_table = NULL; 

struct landslidefinalnode
{ 
  char *name;
  void (*func)();

  struct landslidefinalnode* next;
} ;
struct landslidefinalnode* landslidefinal_table = NULL; 

struct landslidesetpedecnode
{ 
  char *name;
  void (*func)(void*);

  struct landslidesetpedecnode* next;
} ;
struct landslidesetpedecnode* landslidesetpedec_table = NULL; 

struct landslideqcpedecnode
{ 
  char *name;
  void (*func)(void*,void*);

  struct landslideqcpedecnode* next;
} ;
struct landslideqcpedecnode* landslideqcpedec_table = NULL; 

struct landslidesetpeprednode
{ 
  char *name;
  void (*func)(void*);

  struct landslidesetpeprednode* next;
} ;
struct landslidesetpeprednode* landslidesetpepred_table = NULL; 

struct landslideresetnode
{ 
  char *name;
  void (*func)();

  struct landslideresetnode* next;
} ;
struct landslideresetnode* landslidereset_table = NULL; 

//BOP
// !ROUTINE: registerlandslidemodelinit
// \label{registerlandslidemodelinit}
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine to 
//  perform land surface model initialization
// 
// !INTERFACE:
void FTN(registerlandslidemodelinit)(char *j, void (*func)(),int len)
//EOP
{ 
  int len1;
  struct landslideinitnode* current;
  struct landslideinitnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct landslideinitnode*) malloc(sizeof(struct landslideinitnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(landslideinit_table == NULL){
    landslideinit_table = pnode;
  }
  else{
    current = landslideinit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: initializelandslidemodel
// \label{initializelandslidemodel}
//
// !INTERFACE:
void FTN(initializelandslidemodel)(char *j,int len)
//
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to perform
//  land surface model initialization
// 
//  \begin{description}
//  \item[i]
//   index of the LSM
//  \end{description}
//EOP
{ 
  struct landslideinitnode* current;
  
  current = landslideinit_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("init routine for landslide model %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 

}
//BOP
// !ROUTINE: registerlandslidemodelrun
// \label{registerlandslidemodelrun}
//
// !INTERFACE:
void FTN(registerlandslidemodelrun)(char *j, void (*func)(int*),int len)
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine
//  to run the land surface model 
// 
//  \begin{description}
//  \item[i]
//   index of the LSM
//  \end{description}
//EOP
{ 
  int len1;
  struct landsliderunnode* current;
  struct landsliderunnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct landsliderunnode*) malloc(sizeof(struct landsliderunnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(landsliderun_table == NULL){
    landsliderun_table = pnode;
  }
  else{
    current = landsliderun_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: runlandslidemodel
// \label{runlandslidemodel}
//
// !INTERFACE:
void FTN(runlandslidemodel)(char *j,int *n,int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to run the 
//  land surface model 
// 
//  \begin{description}
//  \item[i]
//   index of the LSM
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{
  struct landsliderunnode* current;
  
  current = landsliderun_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("run routine for landslide model %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}


//BOP
// !ROUTINE: registerlandslidemodeloutput
// \label{registerlandslidemodeloutput}
//
// !INTERFACE:
void FTN(registerlandslidemodeloutput)(char *j, void (*func)(int*),int len)
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine
//  to run the land surface model 
// 
//  \begin{description}
//  \item[i]
//   index of the LSM
//  \end{description}
//EOP
{ 
  int len1;
  struct landslideoutnode* current;
  struct landslideoutnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct landslideoutnode*) malloc(sizeof(struct landslideoutnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(landslideout_table == NULL){
    landslideout_table = pnode;
  }
  else{
    current = landslideout_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: outputlandslidemodel
// \label{outputlandslidemodel}
//
// !INTERFACE:
void FTN(outputlandslidemodel)(char *j,int *n,int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to run the 
//  land surface model 
// 
//  \begin{description}
//  \item[i]
//   index of the LSM
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{
  struct landslideoutnode* current;
  
  current = landslideout_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("output routine for landslide %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}

//BOP
// !ROUTINE: registerlandslidemodelfinalize
// \label{registerlandslidemodelfinalize}
//
// !INTERFACE:
void FTN(registerlandslidemodelfinalize)(char *j, void (*func)(),int len)
//  
// !DESCRIPTION:
//  Creates an entry in the registry for the routine
//  to cleanup allocated structures specific to the 
//  land surface model
// 
//  \begin{description}
//  \item[i]
//   index of the LANDSLIDEMODEL
//  \end{description}
// 
//EOP
{ 
  int len1;
  struct landslidefinalnode* current;
  struct landslidefinalnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct landslidefinalnode*) malloc(sizeof(struct landslidefinalnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(landslidefinal_table == NULL){
    landslidefinal_table = pnode;
  }
  else{
    current = landslidefinal_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: landslidemodelfinalize
// \label{landslidemodelfinalize}
//
// !INTERFACE:
void FTN(landslidemodelfinalize)(char *j,int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry for cleaning up
//  allocated structures specific to the land surface model
// 
//  \begin{description}
//  \item[i]
//   index of the LANDSLIDEMODEL
//  \end{description}
// 
//EOP
{
  struct landslidefinalnode* current;
  
  current = landslidefinal_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("finalize routine for landslide model %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}

//BOP
// !ROUTINE: registersetuplandslidepeobspred
// \label{registersetuplandslidepeobspred}
// 
// !INTERFACE:
void FTN(registersetuplandslidepeobspred)(char *j, void (*func)(void*),int len)
//  
// !DESCRIPTION: 
//  registers the method to initialize the obs pred for the PE observation
//  (model simulated value of the observation)
//  
//  \begin{description}
//  \item[i]
//   index of the LANDSLIDE
//  \item[j]
//   index of the optimization set
//  \end{description}
//EOP
{ 
  int len1;
  struct landslidesetpeprednode* current;
  struct landslidesetpeprednode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct landslidesetpeprednode*) malloc(sizeof(struct landslidesetpeprednode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(landslidesetpepred_table == NULL){
    landslidesetpepred_table = pnode;
  }
  else{
    current = landslidesetpepred_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: setuplandslidepeobspredspace
// \label{setuplandslidepeobspredspace}
// 
// !INTERFACE:
void FTN(setuplandslidepeobspredspace)(char *j, void *pepred,int len)
//  
// !DESCRIPTION: 
//   invokes the method to compute the obs pred for the PE observation
//  (model simulated value of the observation)
//
//  \begin{description}
//  \item[i]
//   index of the LANDSLIDE
//  \item[j]
//   index of the variable
//  \item[pepred]
//   object containing the PE obspred 
//  \end{description}
//EOP
{ 
  struct landslidesetpeprednode* current;
  
  current = landslidesetpepred_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("setup landslide obspred routine for  %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(pepred); 
}


//BOP
// !ROUTINE: registerqclandslidedecisionspace
// \label{registerqclandslidedecisionspace}
// 
// !INTERFACE:
void FTN(registerqclandslidedecisionspace)(char *j, void (*func)(void*, void*),int len)
//  
// !DESCRIPTION: 
//  Method to registry an interface implementation for 
//   Qc'ing a specified decision space. 
//  \begin{description}
//  \item[i]
//   index of the LSM
//  \item[j]
//   index of the optimization set
//  \end{description}
//EOP
{ 
  int len1;
  struct landslideqcpedecnode* current;
  struct landslideqcpedecnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct landslideqcpedecnode*) malloc(sizeof(struct landslideqcpedecnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(landslideqcpedec_table == NULL){
    landslideqcpedec_table = pnode;
  }
  else{
    current = landslideqcpedec_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: qclandslidedecisionspace
// \label{qclandslidedecisionspace}
// 
// !INTERFACE:
void FTN(qclandslidedecisionspace)(char *j, void *dec, void *feas,int len)
//  
// !DESCRIPTION: 
//  Method to QC a specified decision space.
//
//  \begin{description}
//  \item[i]
//   index of the LSM
//  \item[j]
//   index of the optimization set
//  \item[dec]
//   decision space object
//  \item[feas]
//   feasible space object
//  \end{description}
//EOP
{ 
  struct landslideqcpedecnode* current;
  
  current = landslideqcpedec_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("QC landslide decision space routine for  %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(dec,feas); 
}

//BOP
// !ROUTINE: registersetlandslidedecisionspace
// \label{registersetlandslidedecisionspace}
// 
// !INTERFACE:
void FTN(registersetlandslidedecisionspace)(char *j, void (*func)(void*),int len)
//  Makes an entry in the registry for the routine 
//  to set the decision space for optimization
// !DESCRIPTION: 
// 
//  \begin{description}
//  \item[i]
//   index of the LANDSLIDE
//  \item[j]
//   index of the variable
//  \end{description}
//EOP
{ 
  int len1;
  struct landslidesetpedecnode* current;
  struct landslidesetpedecnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct landslidesetpedecnode*) malloc(sizeof(struct landslidesetpedecnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(landslidesetpedec_table == NULL){
    landslidesetpedec_table = pnode;
  }
  else{
    current = landslidesetpedec_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: setlandslidedecisionspace
// \label{setlandslidedecisionspace}
// 
// !INTERFACE:
void FTN(setlandslidedecisionspace)(char *j, void *state,int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to set the 
//  decision space for optimization
//
//  \begin{description}
//  \item[i]
//   index of the landslide model
//  \item[n]
//   index of the nest
//  \item[statevars]
//   pointer to the prognostic variable array
//  \end{description}
//EOP
{ 
  struct landslidesetpedecnode* current;
  
  current = landslidesetpedec_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("set landslide decision space routine for  %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(state); 
}


//BOP
// !ROUTINE: registerlandslidemodelreset
// \label{registerlandslidemodelreset}
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine to 
//  perform land surface model resetialization
// 
// !INTERFACE:
void FTN(registerlandslidemodelreset)(char *j, void (*func)(),int len)
//EOP
{ 
  int len1;
  struct landslideresetnode* current;
  struct landslideresetnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct landslideresetnode*) malloc(sizeof(struct landslideresetnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(landslidereset_table == NULL){
    landslidereset_table = pnode;
  }
  else{
    current = landslidereset_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: resetlandslidemodel
// \label{resetlandslidemodel}
//
// !INTERFACE:
void FTN(resetlandslidemodel)(char *j,int len)
//
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to perform
//  landslide model reset
// 
//  \begin{description}
//  \item[i]
//   index of the LSM
//  \end{description}
//EOP
{ 
  struct landslideresetnode* current;
  
  current = landslidereset_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("reset routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}

