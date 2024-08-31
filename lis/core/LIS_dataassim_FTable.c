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
// !MODULE: LIS_dataassim_FTable
//  
// !DESCRIPTION:
//  Function table registries for storing the interface
//  implementations for the operation of data assimilation  
//  algorithms and observations.
//
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"
struct dasetupnode
{ 
  char *name;
  void (*func)(int*);

  struct dasetupnode* next;
} ;
struct dasetupnode* dasetup_table = NULL; 

struct dainitnode
{ 
  char *name;
  void (*func)();

  struct dainitnode* next;
} ;
struct dainitnode* dainit_table = NULL; 

struct dacompincrnode
{ 
  char *name;
  void (*func)(int*, int*);

  struct dacompincrnode* next;
} ;
struct dacompincrnode* dacompincr_table = NULL; 

struct daapplyincrnode
{ 
  char *name;
  void (*func)(int*, int*);

  struct daapplyincrnode* next;
} ;
struct daapplyincrnode* daapplyincr_table = NULL; 

struct daoutnode
{ 
  char *name;
  void (*func)(int*, int*);

  struct daoutnode* next;
} ;
struct daoutnode* daout_table = NULL; 

struct dafinalnode
{ 
  char *name;
  void (*func)();

  struct dafinalnode* next;
} ;
struct dafinalnode* dafinal_table = NULL; 

struct daobssetnode
{ 
  char *name;
  void (*func)(int*, void*, void*);

  struct daobssetnode* next;
} ;
struct daobssetnode* daobsset_table = NULL; 

struct daobsreadnode
{ 
  char *name;
  void (*func)(int*, int*, void*, void*);

  struct daobsreadnode* next;
} ;
struct daobsreadnode* daobsread_table = NULL; 

struct daobswritenode
{ 
  char *name;
  void (*func)(int*, int*, void*);

  struct daobswritenode* next;
} ;
struct daobswritenode* daobswrite_table = NULL; 

struct dagetnsonode
{ 
  char *name;
  void (*func)(int*, int*, int*, int*);

  struct dagetnsonode* next;
} ;
struct dagetnsonode* dagetnso_table = NULL; 

struct daobsfinalnode
{ 
  char *name;
  void (*func)();

  struct daobsfinalnode* next;
} ;
struct daobsfinalnode* daobsfinal_table = NULL; 

//BOP
// !ROUTINE: registerdainit
// \label{registerdainit}
//  
// !DESCRIPTION: 
// Creates an entry in the registry for the routine to 
// intialize algorithm specific structures for 
// the data assimilation method used. 
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the assimilation algorithm
//  \end{description}
// 
// !INTERFACE:
void FTN(registerdainit)(char *j, void (*func)(),int len)
//EOP
{ 
  int len1;
  struct dainitnode* current;
  struct dainitnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct dainitnode*) malloc(sizeof(struct dainitnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(dainit_table == NULL){
    dainit_table = pnode;
  }
  else{
    current = dainit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: dataassiminit
// \label{dataassiminit}
//
// !INTERFACE:
void FTN(dataassiminit)(char *j, int len)
//  
// !DESCRIPTION: 
// Invokes the routine from the registry to initialize
// the specific structures for the 
// data assimilation algorithm used. 
// 
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the assimilation algorithm
//  \end{description}
// 
//EOP
{
  struct dainitnode* current;
  
  current = dainit_table;
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
// !ROUTINE: registerdasetup
// \label{registerdasetup}
//  
// !DESCRIPTION: 
// Creates an entry in the registry for the routine to 
// set up algorithm specific structures for 
// the data assimilation method used. 
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the assimilation algorithm
//  \end{description}
// 
// !INTERFACE:
void FTN(registerdasetup)(char *j, void (*func)(int*),int len)
//EOP
{ 
  int len1;
  struct dasetupnode* current;
  struct dasetupnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct dasetupnode*) malloc(sizeof(struct dasetupnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(dasetup_table == NULL){
    dasetup_table = pnode;
  }
  else{
    current = dasetup_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: dataassimsetup
// \label{dataassimsetup}
//
// !INTERFACE:
void FTN(dataassimsetup)(char *j, int *k,int len)
//  
// !DESCRIPTION: 
// Invokes the routine from the registry to setup
// the specific structures for the 
// data assimilation algorithm used. 
// 
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the assimilation algorithm
//  \item[k]
//   index of the DA instance
//  \end{description}
// 
//EOP
{
  struct dasetupnode* current;
  
  current = dasetup_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("setup routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(k); 
}


//BOP
// !ROUTINE: registerapplyincrements
// \label{registerapplyincrements}
//  
// !INTERFACE:
void FTN(registerapplyincrements)(char *j, void (*func)(int*, int*),int len)
// !DESCRIPTION: 
// Creates an entry in the registry for the routine that
// applies the analysis increments
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the assimilation algorithm
//  \end{description}
//  
//EOP
{ 
  int len1;
  struct daapplyincrnode* current;
  struct daapplyincrnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct daapplyincrnode*) malloc(sizeof(struct daapplyincrnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(daapplyincr_table == NULL){
    daapplyincr_table = pnode;
  }
  else{
    current = daapplyincr_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: applyincrements
// \label{applyIncrements}
//  
// !INTERFACE:
void FTN(applyincrements)(char *j, int *n, int *k,int len)
// !DESCRIPTION:
// Invoke the routine from the registry that applies
// analysis increments
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the assimilation algorithm
//  \item[n]
//   index of the nest
//  \item[k]
//   index of the assimilation instance
//  \end{description}
//  
//EOP
{
  struct daapplyincrnode* current;
  
  current = daapplyincr_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("apply increments routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,k); 
}


//BOP
// !ROUTINE: registercomputeincrements
// \label{registercomputeincrements}
//  
// !INTERFACE:
void FTN(registercomputeincrements)(char *j, void (*func)(int*, int*),int len)
// !DESCRIPTION: 
// Creates an entry in the registry for the routine that
// computes the analysis increments
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the assimilation algorithm
//  \end{description}
//  
//EOP
{ 
  int len1;
  struct dacompincrnode* current;
  struct dacompincrnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct dacompincrnode*) malloc(sizeof(struct dacompincrnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(dacompincr_table == NULL){
    dacompincr_table = pnode;
  }
  else{
    current = dacompincr_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: computeincrements
// \label{computeIncrements}
//  
// !INTERFACE:
void FTN(computeincrements)(char *j, int *n, int *k,int len)
// !DESCRIPTION:
// Invoke the routine from the registry that computes 
// analysis increments
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the assimilation algorithm
//  \item[n]
//   index of the nest
//  \item[k]
//   index of assimilation instance
//  \end{description}
//  
//EOP
{
  struct dacompincrnode* current;
  
  current = dacompincr_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("compute increments routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,k); 
}

//BOP
// !ROUTINE: registerdafinalize
// \label{registerdafinalize}
//  
//  
// !INTERFACE:
void FTN(registerdafinalize)(char *j, void (*func)(),int len)
// !DESCRIPTION: 
//  Cretes an entry in the registry for the routine 
//  to cleanup allocated structures
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the assimilation algorithm
//  \end{description}
//EOP
{ 
  int len1;
  struct dafinalnode* current;
  struct dafinalnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct dafinalnode*) malloc(sizeof(struct dafinalnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(dafinal_table == NULL){
    dafinal_table = pnode;
  }
  else{
    current = dafinal_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: dafinalize
// \label{dafinalize}
// 
// !INTERFACE:
void FTN(dafinalize)(char *j,int len)
// !DESCRIPTION: 
// Invokes the routine from the registry that cleans up
// the data assimilation specific structures
// 
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the assimilation algorithm 
//  \end{description}
// 
//EOP
{ 
  struct dafinalnode* current;
  
  current = dafinal_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("finalize routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}


//BOP
// !ROUTINE: registerdaoutput
// \label{registerdaoutput}
//  
//  
// !INTERFACE:
void FTN(registerdaoutput)(char *j, void (*func)(int*, int*),int len)
// !DESCRIPTION: 
//  Cretes an entry in the registry for the 
//  writing diagnostic output from the data
//  assimilation routines
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the assimilation algorithm
//  \end{description}
//EOP
{ 
  int len1;
  struct daoutnode* current;
  struct daoutnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct daoutnode*) malloc(sizeof(struct daoutnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(daout_table == NULL){
    daout_table = pnode;
  }
  else{
    current = daout_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: daoutput
// \label{daoutput}
// 
// !INTERFACE:
void FTN(daoutput)(char *j,int *n, int *k,int len)
// !DESCRIPTION: 
// Invokes the routine from the registry to write 
// data assimilation specific output
// 
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the assimilation algorithm 
//  \item[n]
//   index of the nest
//  \item[k]
//   index of the assimilation instance
//  \end{description}
// 
//EOP
{ 
  struct daoutnode* current;
  
  current = daout_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Output routine for  %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,k); 
}
//BOP
// !ROUTINE: registergetnso
// \label{registergetnso}
//
// !INTERFACE:
void FTN(registergetnso)(char *j, void (*func)(int*, int*, int*, int*),int len)
//  
// !DESCRIPTION: 
// Creates an entry in the registry for the routine that
// computes the number of selected observations for the 
// specified grid point. 
// 
// The arguments are: 
// \begin{description}
//  \item[j]
//   index of the observations variable
//  \end{description}
//EOP
{ 
  int len1;
  struct dagetnsonode* current;
  struct dagetnsonode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct dagetnsonode*) malloc(sizeof(struct dagetnsonode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(dagetnso_table == NULL){
    dagetnso_table = pnode;
  }
  else{
    current = dagetnso_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: getselectedobsnumber 
// \label{getselectedobsnumber}
//
// !INTERFACE: 
void FTN(getselectedobsnumber)(char *j, int *n, int *gid, int *sid, int *eid,int len)
//
// !DESCRIPTION:
// Invokes the routine from the registry that computes
// the number of selected observations for a particular modeling
// point
// 
// The arguments are: 
// \begin{description}
//  \item[j]
//   index of the observations variable
//  \item[n]
//   index of the nest
//  \item[gid]
//   model grid point index
//  \item[sid]
//   starting index of the selected observations
//  \item[eid]
//   ending index of the selected observations
//  \end{description}
//
//EOP
{ 
  struct dagetnsonode* current;
  
  current = dagetnso_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("getNSO (number of selected observations) routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,gid,sid,eid); 
}

//BOP
// !ROUTINE: registerreaddaobssetup
// \label{registerrdaobssetup}
//  
// 
// !INTERFACE:
void FTN(registerdaobssetup)(char *j, void (*func)(int*, void*, void*),int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine that
//  sets up structures for handling observation data for 
//  data assimilation
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   index of the assimilation set
//  \end{description}
//EOP
{ 
  int len1;
  struct daobssetnode* current;
  struct daobssetnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct daobssetnode*) malloc(sizeof(struct daobssetnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(daobsset_table == NULL){
    daobsset_table = pnode;
  }
  else{
    current = daobsset_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: readobsdataconfig
// \label{readobsdataconfig}
// 
// !INTERFACE:
void FTN(readobsdataconfig)(char *j, int *k, void *obs, void *pert,int len)
//  
// !DESCRIPTION: 
// Invokes the routine from the registry to set up 
// structures for handling observation data for data
// assimilation 
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   index of the variable being assimilated
//  \item[k]
//   index of the DA instance
//  \item[obs]
//   ESMF state that contain the observations
//  \item[pert]
//   ESMF state that contain the observation perturbations
//  \end{description}
//
//EOP
{ 
  struct daobssetnode* current;
  
  current = daobsset_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("setup routine for DA obs %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(k,obs,pert); 
}

//BOP
// !ROUTINE: registerreaddaobs
// \label{registerreaddaobs}
//
// !INTERFACE:
void FTN(registerreaddaobs)(char *j, void (*func)(int*, int*, void*, void*),int len)
//  
// !DESCRIPTION: 
// Creates an entry in the registry for the routine that 
// reads observations for data assimilation. 
// 
// The arguments are: 
// \begin{description}
//  \item[j]
//   index of the observations variable
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
// !ROUTINE: readdaobs
// \label{readdaobs}
//
// !INTERFACE: 
void FTN(readdaobs)(char *j,int *n, int *k, void *obsstate,void *obspertstate, int len)
//
// !DESCRIPTION:
// Invokes the routine from the registry to read 
// observations for data assimilation. This routine
// also packages the observations into an ESMF state
// 
// The arguments are: 
// \begin{description}
//  \item[j]
//   index of the observations variable
//  \item[n]
//   index of the nest
//  \item[obsstate]
//   pointer to the ESMF observation state
//  \item[obspertstate]
//   pointer to the ESMF observation perturbation state
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
      printf("read DA obs routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,k, obsstate,obspertstate); 
}
//BOP
// !ROUTINE: registerwritedaobs
// \label{registerwritedaobs}
//
// !INTERFACE:
void FTN(registerwritedaobs)(char *j, void (*func)(int*, int*, void*),int len)
//  
// !DESCRIPTION: 
// Creates an entry in the registry for the routine that 
// writes observations used in data assimilation to disk. 
// 
// The arguments are: 
// \begin{description}
//  \item[j]
//   index of the observations variable
//  \end{description}
//EOP
{ 
  int len1;
  struct daobswritenode* current;
  struct daobswritenode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct daobswritenode*) malloc(sizeof(struct daobswritenode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(daobswrite_table == NULL){
    daobswrite_table = pnode;
  }
  else{
    current = daobswrite_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: writedaobs
// \label{writedaobs}
//
// !INTERFACE: 
void FTN(writedaobs)(char *j,int *n, int *k, void *obsstate,int len)
//
// !DESCRIPTION:
// Invokes the routine from the registry to write
// observations used in data assimilation to disk.
// 
// The arguments are: 
// \begin{description}
//  \item[j]
//   index of the observations variable
//  \item[n]
//   index of the nest
//  \item[obsstate]
//   pointer to the ESMF observation state
//  \end{description}
//
//EOP
{ 
  struct daobswritenode* current;
  
  current = daobswrite_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("write DA obs routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,k,obsstate); 
}



