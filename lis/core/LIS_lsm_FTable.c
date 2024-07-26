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
// !MODULE: LIS_lsm_FTable
//  
//
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing the operations of different 
//  land surface models. The registries also contain 
//  related interface implementations for data assimilation, 
//  WRF/GCE/GFS coupling and parameter estimation
//
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"
struct lsminitnode
{ 
  char *name;
  void (*func)();

  struct lsminitnode* next;
} ;
struct lsminitnode* lsminit_table = NULL; 

struct lsmrunnode
{ 
  char *name;
  void (*func)(int*);

  struct lsmrunnode* next;
} ;
struct lsmrunnode* lsmrun_table = NULL; 

struct lsmfinalnode
{ 
  char *name;
  void (*func)();

  struct lsmfinalnode* next;
} ;
struct lsmfinalnode* lsmfinal_table = NULL; 

struct lsmresetnode
{ 
  char *name;
  void (*func)();

  struct lsmresetnode* next;
} ;
struct lsmresetnode* lsmreset_table = NULL; 

struct lsmsetupnode
{ 
  char *name;
  void (*func)();

  struct lsmsetupnode* next;
} ;
struct lsmsetupnode* lsmsetup_table = NULL;

struct lsmrestartnode
{ 
  char *name;
  void (*func)();

  struct lsmrestartnode* next;
} ;
struct lsmrestartnode* lsmrestart_table = NULL;


struct lsmdynsetnode
{ 
  char *name;
  void (*func)(int*);

  struct lsmdynsetnode* next;
} ;
struct lsmdynsetnode* lsmdynset_table = NULL;

struct lsmf2tnode
{ 
  char *name;
  void (*func)(int*);

  struct lsmf2tnode* next;
} ;
struct lsmf2tnode* lsmf2t_table = NULL;

struct lsmwriterstnode
{ 
  char *name;
  void (*func)(int*);

  struct lsmwriterstnode* next;
} ;
struct lsmwriterstnode* lsmwriterst_table = NULL;

//for DA
struct lsmdainitnode
{ 
  char *name;
  void (*func)(int*);

  struct lsmdainitnode* next;
} ;
struct lsmdainitnode* lsmdainit_table = NULL;

struct lsmdagetvarnode
{ 
  char *name;
  void (*func)(int*, void*);

  struct lsmdagetvarnode* next;
} ;
struct lsmdagetvarnode* lsmdagetvar_table = NULL;

struct lsmdasetvarnode
{ 
  char *name;
  void (*func)(int*, void*);

  struct lsmdasetvarnode* next;
} ;
struct lsmdasetvarnode* lsmdasetvar_table = NULL;

struct lsmdaobstransformnode
{ 
  char *name;
  void (*func)(int*, void*);

  struct lsmdaobstransformnode* next;
} ;
struct lsmdaobstransformnode* lsmdaobstransform_table = NULL;

struct lsmdamapobstolsmnode
{ 
  char *name;
  void (*func)(int*, int*, void*,void*);

  struct lsmdamapobstolsmnode* next;
} ;
struct lsmdamapobstolsmnode* lsmdamapobstolsm_table = NULL;

struct lsmdaobsprednode
{ 
  char *name;
  void (*func)(int*, int*, float*);

  struct lsmdaobsprednode* next;
} ;
struct lsmdaobsprednode* lsmdaobspred_table = NULL;

struct lsmdaqcstatenode
{ 
  char *name;
  void (*func)(int*, void*);

  struct lsmdaqcstatenode* next;
} ;
struct lsmdaqcstatenode* lsmdaqcstate_table = NULL;

struct lsmdaqcobsnode
{ 
  char *name;
  void (*func)(int*, int*, void*);

  struct lsmdaqcobsnode* next;
} ;
struct lsmdaqcobsnode* lsmdaqcobs_table = NULL;

struct lsmdascalenode
{ 
  char *name;
  void (*func)(int*, void*);

  struct lsmdascalenode* next;
} ;
struct lsmdascalenode* lsmdascale_table = NULL;

struct lsmdadescalenode
{ 
  char *name;
  void (*func)(int*, void*, void*);

  struct lsmdadescalenode* next;
} ;
struct lsmdadescalenode* lsmdadescale_table = NULL;

struct lsmdaupdatenode
{ 
  char *name;
  void (*func)(int*, void*,void*);

  struct lsmdaupdatenode* next;
} ;
struct lsmdaupdatenode* lsmdaupdate_table = NULL;

struct lsmdawritevarnode
{ 
  char *name;
  void (*func)(int*, int*,void*);

  struct lsmdawritevarnode* next;
} ;
struct lsmdawritevarnode* lsmdawritevar_table = NULL;

struct lsmdiagfordanode
{ 
  char *name;
  void (*func)(int*);

  struct lsmdiagfordanode* next;
} ;
struct lsmdiagfordanode* lsmdiagforda_table = NULL;

//coupling to atmos. models
struct lsmcplsetexportnode
{ 
  char *name;
  void (*func)(int*);

  struct lsmcplsetexportnode* next;
} ;
struct lsmcplsetexportnode* lsmcplsetexport_table = NULL;

//Parameter estimation
struct lsmpesetdecnode
{ 
  char *name;
  void (*func)(void*, void*);

  struct lsmpesetdecnode* next;
} ;
struct lsmpesetdecnode* lsmpesetdec_table = NULL;

struct lsmpegetdecnode
{ 
  char *name;
  void (*func)(void*);

  struct lsmpegetdecnode* next;
} ;
struct lsmpegetdecnode* lsmpegetdec_table = NULL;

struct lsmpesetupdecnode
{ 
  char *name;
  void (*func)(void*, void*);

  struct lsmpesetupdecnode* next;
} ;
struct lsmpesetupdecnode* lsmpesetupdec_table = NULL;

struct lsmpesetprednode
{ 
  char *name;
  void (*func)(void*);

  struct lsmpesetprednode* next;
} ;
struct lsmpesetprednode* lsmpesetpred_table = NULL;

struct lsmpesetupobsprednode
{ 
  char *name;
  void (*func)(void*);

  struct lsmpesetupobsprednode* next;
} ;
struct lsmpesetupobsprednode* lsmpesetupobspred_table = NULL;

struct lsmpeobsprednode
{ 
  char *name;
  void (*func)(void*);

  struct lsmpeobsprednode* next;
} ;
struct lsmpeobsprednode* lsmpeobspred_table = NULL;

//Routing
struct lsmroutinggetrunoffnode
{ 
  char *name;
  void (*func)(int*);

  struct lsmroutinggetrunoffnode* next;
} ;
struct lsmroutinggetrunoffnode* lsmroutinggetrunoff_table = NULL;

struct lsmroutinggetswsnode
{
  char *name;
  void (*func)(int*);

  struct lsmroutinggetswsnode* next;
} ;
struct lsmroutinggetswsnode* lsmroutinggetsws_table = NULL;

struct lsm2rtmnode
{ 
  char *name;
  void (*func)(int*,void*);

  struct lsm2rtmnode* next;
} ;
struct lsm2rtmnode* lsm2rtm_table = NULL; 

struct lsmirriggetnode
{ 
  char *name;
  void (*func)(int*,void*);

  struct lsmirriggetnode* next;
} ;
struct lsmirriggetnode* lsmirrigget_table = NULL; 


//BOP
// !ROUTINE: registerlsmini
// \label{registerlsmini}
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine to 
//  perform land surface model initialization
// 
// !INTERFACE:
void FTN(registerlsminit)(char *j, void (*func)(),int len)
//EOP
{ 
  int len1;
  struct lsminitnode* current;
  struct lsminitnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsminitnode*) malloc(sizeof(struct lsminitnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsminit_table == NULL){
    lsminit_table = pnode;
  }
  else{
    current = lsminit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: lsminit
// \label{lsminit}
//
// !INTERFACE:
void FTN(lsminit)(char *j,int len)
//
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to perform
//  land surface model initialization
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \end{description}
//EOP
{ 

  struct lsminitnode* current;
  int found ; 

  current = lsminit_table;

  while(strcmp(current->name,j)!=0){
    current = current->next;

    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("init routine for LSM %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 

}
//BOP
// !ROUTINE: registerlsmrun
// \label{registerlsmrun}
//
// !INTERFACE:
void FTN(registerlsmrun)(char *j, void (*func)(int*),int len)
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine
//  to run the land surface model 
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmrunnode* current;
  struct lsmrunnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmrunnode*) malloc(sizeof(struct lsmrunnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmrun_table == NULL){
    lsmrun_table = pnode;
  }
  else{
    current = lsmrun_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lsmrun
// \label{lsmrun}
//
// !INTERFACE:
void FTN(lsmrun)(char *j,int *n,int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to run the 
//  land surface model 
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{
  struct lsmrunnode* current;
  
  current = lsmrun_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;  
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("run routine for LSM %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}

//BOP
// !ROUTINE: registerlsmfinalize
// \label{registerlsmfinalize}
//
// !INTERFACE:
void FTN(registerlsmfinalize)(char *j, void (*func)(),int len)
//  
// !DESCRIPTION:
//  Creates an entry in the registry for the routine
//  to cleanup allocated structures specific to the 
//  land surface model
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \end{description}
// 
//EOP
{ 
  int len1;
  struct lsmfinalnode* current;
  struct lsmfinalnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmfinalnode*) malloc(sizeof(struct lsmfinalnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmfinal_table == NULL){
    lsmfinal_table = pnode;
  }
  else{
    current = lsmfinal_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: lsmfinalize
// \label{lsmfinalize}
//
// !INTERFACE:
void FTN(lsmfinalize)(char *j,int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry for cleaning up
//  allocated structures specific to the land surface model
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \end{description}
// 
//EOP
{
  struct lsmfinalnode* current;
  
  current = lsmfinal_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("finalize routine for LSM %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}

//BOP
// !ROUTINE: registerlsmreset
// \label{registerlsmreset}
//
// !INTERFACE:
void FTN(registerlsmreset)(char *j, void (*func)(),int len)
//  
// !DESCRIPTION:
//  Creates an entry in the registry for the routine
//  to cleanup allocated structures specific to the 
//  land surface model
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \end{description}
// 
//EOP
{ 
  int len1;
  struct lsmresetnode* current;
  struct lsmresetnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmresetnode*) malloc(sizeof(struct lsmresetnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmreset_table == NULL){
    lsmreset_table = pnode;
  }
  else{
    current = lsmreset_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: lsmreset
// \label{lsmreset}
//
// !INTERFACE:
void FTN(lsmreset)(char *j,int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry for cleaning up
//  allocated structures specific to the land surface model
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \end{description}
// 
//EOP
{
  struct lsmresetnode* current;
  
  current = lsmreset_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("reset routine for LSM %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}

//BOP
// !ROUTINE: registerlsmsetup
// \label{registerlsmsetup}
//
// !INTERFACE:
void FTN(registerlsmsetup)(char *j, void (*func)(),int len)
//  
// !DESCRIPTION:
//  Makes an entry in the registry for the routine
//  to set up land surface model parameters 
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmsetupnode* current;
  struct lsmsetupnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmsetupnode*) malloc(sizeof(struct lsmsetupnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmsetup_table == NULL){
    lsmsetup_table = pnode;
  }
  else{
    current = lsmsetup_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lsmsetup
// \label{lsmsetup}
//
// !INTERFACE:
void FTN(lsmsetup)(char *j, int len)
//  
// !DESCRIPTION:  
//  Invokes the routine in the registry to set up 
//  land surface model parameters  
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \end{description}
//EOP
{ 
  struct lsmsetupnode* current;
  
  current = lsmsetup_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("setup routine for LSM %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}

//BOP
// !ROUTINE: registerlsmrestart
// \label{registerlsmrestart}
// 
// !INTERFACE:
void FTN(registerlsmrestart)(char *j, void (*func)(),int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// restart the land surface model from a 
// previously saved state
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmrestartnode* current;
  struct lsmrestartnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmrestartnode*) malloc(sizeof(struct lsmrestartnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmrestart_table == NULL){
    lsmrestart_table = pnode;
  }
  else{
    current = lsmrestart_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lsmrestart
// \label{lsmrestart}
//
// !INTERFACE:
void FTN(lsmrestart)(char *j, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to 
//  restart the land surface model from a previously
//  saved state
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \end{description}
//EOP
{ 
  struct lsmrestartnode* current;
  
  current = lsmrestart_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("read restart routine for LSM %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}

//BOP
// !ROUTINE: registerlsmdynsetup
// \label{registerlsmdynsetup}
// 
// !INTERFACE:
void FTN(registerlsmdynsetup)(char *j, void (*func)(int*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  set the time dependent land surface parameters
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmdynsetnode* current;
  struct lsmdynsetnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmdynsetnode*) malloc(sizeof(struct lsmdynsetnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmdynset_table == NULL){
    lsmdynset_table = pnode;
  }
  else{
    current = lsmdynset_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lsmdynsetup
// \label{lsmdynsetup}
// 
// !INTERFACE:
void FTN(lsmdynsetup)(char *j, int *n, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to set the 
//  time dependent land surface parameters
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{ 
  struct lsmdynsetnode* current;
  
  current = lsmdynset_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;

    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("dynamic setup routine for LSM %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}


//BOP
// !ROUTINE: registerlsmf2t
// \label{registerlsmf2t}
// 
// !INTERFACE:
void FTN(registerlsmf2t)(char *j, void (*func)(int*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  transfer forcing to model tiles
// 
//  \begin{description}
//  \item[j]
//   name of the LSM + name of the running mode
//  \item[j]
//   index of the runmode
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmf2tnode* current;
  struct lsmf2tnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmf2tnode*) malloc(sizeof(struct lsmf2tnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmf2t_table == NULL){
    lsmf2t_table = pnode;
  }
  else{
    current = lsmf2t_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lsmf2t
// \label{lsmf2t}
// 
// !INTERFACE:
void FTN(lsmf2t)(char *j, int *n, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to 
//  transfer forcing to model tiles
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \item[j]
//   index of the runmode
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{ 
  struct lsmf2tnode* current;
  
  current = lsmf2t_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;

    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("f2t writing routine for LSM and running mode %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}

//BOP
// !ROUTINE: registerlsmwrst
// \label{registerlsmwrst}
// 
// !INTERFACE:
void FTN(registerlsmwrst)(char *j, void (*func)(int*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine 
//  to write restart files for a land surface model 
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmwriterstnode* current;
  struct lsmwriterstnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmwriterstnode*) malloc(sizeof(struct lsmwriterstnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmwriterst_table == NULL){
    lsmwriterst_table = pnode;
  }
  else{
    current = lsmwriterst_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: lsmwrst
// \label{lsmwrst}
// 
// !INTERFACE:
void FTN(lsmwrst)(char *j, int *n, int len)
//  
// !DESCRIPTION:
//  Invokes the routine from the registry to write
//  restart files for a land surface model
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{ 

  struct lsmwriterstnode* current;
  
  current = lsmwriterst_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("write restart writing routine for LSM %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}

//BOP
// !ROUTINE: registerlsmdainit
// \label{registerlsmdainit}
// 
// !INTERFACE:
void FTN(registerlsmdainit)(char *j, void (*func)(int*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine 
//  for initializing DA related LSM settings
// 
//  \begin{description}
//  \item[j]
//   name of the LSM + DA instance
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmdainitnode* current;
  struct lsmdainitnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmdainitnode*) malloc(sizeof(struct lsmdainitnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmdainit_table == NULL){
    lsmdainit_table = pnode;
  }
  else{
    current = lsmdainit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: lsmdainit
// \label{lsmdainit}
// 
// !INTERFACE:
void FTN(lsmdainit)(char *j, int *k, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry for initializing
//  DA related LSM settings. 
//  \begin{description}
//  \item[j]
//   name of the LSM + DA instance
//  \end{description}
//EOP
{ 
  struct lsmdainitnode* current;
  
  current = lsmdainit_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("init routine for LSM + DA instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(k); 
}

//BOP
// !ROUTINE: registerlsmdagetstatevar
// \label{registerlsmdagetstatevar}
// 
// !INTERFACE:
void FTN(registerlsmdagetstatevar)(char *j, void (*func)(int*, void*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine 
//  for obtaining the specified prognostic variables from the 
//  land surface model (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the LSM + DA instance
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmdagetvarnode* current;
  struct lsmdagetvarnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmdagetvarnode*) malloc(sizeof(struct lsmdagetvarnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmdagetvar_table == NULL){
    lsmdagetvar_table = pnode;
  }
  else{
    current = lsmdagetvar_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: lsmdagetstatevar
// \label{lsmdagetstatevar}
// 
// !INTERFACE:
void FTN(lsmdagetstatevar)(char *j, int *n, void *state, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry for obtaining
//  the specified prognostic variables from the land surface model
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the LSM + DA instance
//  \item[n]
//   index of the nest
//  \item[state]
//   pointer to the prognostic variable state
//  \end{description}
//EOP
{ 
  struct lsmdagetvarnode* current;
  
  current = lsmdagetvar_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("get LSM variable routine for LSM + DA instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,state); 
}

//BOP
// !ROUTINE: registerlsmdasetstatevar
// \label{registerlsmdasetstatevar}
// 
// !INTERFACE:
void FTN(registerlsmdasetstatevar)(char *j, void (*func)(int*, void*),int len)
//  
// !DESCRIPTION:
//  Makes an entry in the registry for updating the specified
//  state variable in a land surface model 
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the LSM + DA instance
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmdasetvarnode* current;
  struct lsmdasetvarnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmdasetvarnode*) malloc(sizeof(struct lsmdasetvarnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmdasetvar_table == NULL){
    lsmdasetvar_table = pnode;
  }
  else{
    current = lsmdasetvar_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: lsmdasetstatevar
// \label{lsmdasetstatevar}
// 
// !INTERFACE:
void FTN(lsmdasetstatevar)(char *j,int *n, void *statevar, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry for updating
//  the specified state variable in a land surface model 
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the LSM + DA instance
//  \item[n]
//   index of the nest
//  \item[statevars]
//   pointer to the prognostic variable state
//  \end{description}
//EOP
{ 
  struct lsmdasetvarnode* current;
  
  current = lsmdasetvar_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;

    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("set LSM variable routine for LSM + DA instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,statevar); 
}


//BOP
// !ROUTINE: registerlsmdaobstransform
// \label{registerlsmdaobstransform}
// 
// !INTERFACE:
void FTN(registerlsmdaobstransform)(char *j, void (*func)(int*, void*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry to perform the  
//  translation of observations to state variables 
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the LSM + DA instance
//  \item[k]
//   index of the observation data
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmdaobstransformnode* current;
  struct lsmdaobstransformnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmdaobstransformnode*) malloc(sizeof(struct lsmdaobstransformnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmdaobstransform_table == NULL){
    lsmdaobstransform_table = pnode;
  }
  else{
    current = lsmdaobstransform_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}

//BOP
// !ROUTINE: lsmdaobstransform
// \label{lsmdaobstransform}
//
// !INTERFACE:
void FTN(lsmdaobstransform)(char *j, int *n, void *obs, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to translate
//  the observations to the prognostic variable space
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the LSM + DA instance
//  \item[n]
//   index of the nest
//   \item[obs]
//   transformed variable
//  \end{description}
//EOP
{ 
  struct lsmdaobstransformnode* current;
  
  current = lsmdaobstransform_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("obs transform routine for LSM + DA instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,obs); 
}

//BOP
// !ROUTINE: registerlsmdagetobspred
// \label{registerlsmdagetobspred}
// 
// !INTERFACE:
void FTN(registerlsmdagetobspred)(char *j, void (*func)(int*,int*,float*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine 
//  that provides an LSM's estimate of the observations.
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the LSM + DA instance
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmdaobsprednode* current;
  struct lsmdaobsprednode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmdaobsprednode*) malloc(sizeof(struct lsmdaobsprednode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmdaobspred_table == NULL){
    lsmdaobspred_table = pnode;
  }
  else{
    current = lsmdaobspred_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}

//BOP
// !ROUTINE: lsmdagetobspred
// \label{lsmdagetobspred}
//
// !INTERFACE:
void FTN(lsmdagetobspred)(char *j, int *n, int *k,float *pred, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to translate
//  the observations to state variables
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \item[n]
//   index of the nest
//  \item[pred]
//   model's estimated observation prediction
//  \end{description}
//EOP
{ 
  struct lsmdaobsprednode* current;
  
  current = lsmdaobspred_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;

    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("obspred routine for LSM + DA instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,k,pred); 
}

//BOP
// !ROUTINE: registerlsmdadiagnosevars
// \label{registerlsmdadiagnosevars}
// 
// !INTERFACE:
void FTN(registerlsmdadiagnosevars)(char *j, void (*func)(int*),int len)
//  
// !DESCRIPTION:
//  Makes an entry in the registry for the routine to 
//  perform land surface model output
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmdiagfordanode* current;
  struct lsmdiagfordanode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmdiagfordanode*) malloc(sizeof(struct lsmdiagfordanode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmdiagforda_table == NULL){
    lsmdiagforda_table = pnode;
  }
  else{
    current = lsmdiagforda_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lsmdadiagnosevars
// \label{lsmdadiagnosevars}
//
// !INTERFACE:
void FTN(lsmdadiagnosevars)(char *j, int *n, int len)
//  
// !DESCRIPTION:
//  Invokes the routine from the registry to perform
//  land surface model output
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{ 

  struct lsmdiagfordanode* current;
  
  current = lsmdiagforda_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("diagnose for DA routine for LSM +DAset %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}


//BOP
// !ROUTINE: registerlsmdamapobstolsm
// \label{registerlsmdamapobstolsm}
// 
// !INTERFACE:
void FTN(registerlsmdamapobstolsm)(char *j, void (*func)(int*, int*, void*, void*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry to perform the  
//  translation of observations to state variables 
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the LSM + DA instance
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmdamapobstolsmnode* current;
  struct lsmdamapobstolsmnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmdamapobstolsmnode*) malloc(sizeof(struct lsmdamapobstolsmnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmdamapobstolsm_table == NULL){
    lsmdamapobstolsm_table = pnode;
  }
  else{
    current = lsmdamapobstolsm_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}

//BOP
// !ROUTINE: lsmdamapobstolsm
// \label{lsmdamapobstolsm}
//
// !INTERFACE:
void FTN(lsmdamapobstolsm)(char *j, int *n, int *k, void *obs, void *lsm, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to translate
//  the observations to state variables
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the LSM + DA instance
//  \item[n]
//   index of the nest
//  \item[obs]
//   observations to be mapped
//  \item[lsm]
//   updated lsm states  
//  \end{description}
//EOP
{ 
  struct lsmdamapobstolsmnode* current;
  
  current = lsmdamapobstolsm_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("map obs to LSM routine for LSM + DA instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,k,obs,lsm); 
}

//BOP
// !ROUTINE: registerlsmdaqcstate
// \label{registerlsmdaqcstate}
// 
// !INTERFACE:
void FTN(registerlsmdaqcstate)(char *j, void (*func)(int*, void*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  QC the updated LSM state
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the LSM + DA instance
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmdaqcstatenode* current;
  struct lsmdaqcstatenode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmdaqcstatenode*) malloc(sizeof(struct lsmdaqcstatenode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmdaqcstate_table == NULL){
    lsmdaqcstate_table = pnode;
  }
  else{
    current = lsmdaqcstate_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lsmdaqcstate
// \label{lsmdaqcstate}
// 
// !INTERFACE:
void FTN(lsmdaqcstate)(char *j,int *n, void *LSM_State, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to set the 
//  QC the updated LSM variables
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \item[n]
//   index of the nest
//  \item[LSM\_State]
//   The LSM state being qc'd
//  \end{description}
//EOP
{ 
  struct lsmdaqcstatenode* current;
  
  current = lsmdaqcstate_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("map obs to LSM routine for LSM + DA instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,LSM_State); 
}

//BOP
// !ROUTINE: registerlsmdaqcobsstate
// \label{registerlsmdaqcobsstate}
// 
// !INTERFACE:
void FTN(registerlsmdaqcobsstate)(char *j, void (*func)(int*, int*, void*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  QC the OBS state based on LSM variables and states
//  (for data assimilation).
// 
//  \begin{description}
//  \item[j]
//   name of the LSM + DA instance
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmdaqcobsnode* current;
  struct lsmdaqcobsnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmdaqcobsnode*) malloc(sizeof(struct lsmdaqcobsnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmdaqcobs_table == NULL){
    lsmdaqcobs_table = pnode;
  }
  else{
    current = lsmdaqcobs_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lsmdaqcobsstate
// \label{lsmdaqcobsstate}
// 
// !INTERFACE:
void FTN(lsmdaqcobsstate)(char *j, int *n, int *k, void *LSM_State, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to set the 
//  QC the observation state based on LSM variables and states
//  (for data assimilation).
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \item[n]
//   index of the nest
//  \item[LSM\_State]
//   The LSM state being qc'd
//  \end{description}
//EOP
{ 
  struct lsmdaqcobsnode* current;
  
  current = lsmdaqcobs_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("qc obs state routine for LSM + DA instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,k,LSM_State); 
}


//BOP
// !ROUTINE: registerlsmdascalestatevar
// \label{registerlsmdascalestatevar}
// 
// !INTERFACE:
void FTN(registerlsmdascalestatevar)(char *j, void (*func)(int*, void*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  scale the LSM state variables
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the LSM + DA instance
//  \item[j]
//   index of the assimilated variable
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmdascalenode* current;
  struct lsmdascalenode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmdascalenode*) malloc(sizeof(struct lsmdascalenode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmdascale_table == NULL){
    lsmdascale_table = pnode;
  }
  else{
    current = lsmdascale_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lsmdascalestatevar
// \label{lsmdascalestatevar}
// 
// !INTERFACE:
void FTN(lsmdascalestatevar)(char *j, int *n, void *LSM_State, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to scale the LSM 
//  state variables
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the LSM + DA instance
//  \item[n]
//   index of the nest
//  \item[LSM\_State]
//   The LSM state being scaled
//  \end{description}
//EOP
{ 
  struct lsmdascalenode* current;
  
  current = lsmdascale_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("scale variable routine for LSM + DA instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,LSM_State); 
}


//BOP
// !ROUTINE: registerlsmdadescalestatevar
// \label{registerlsmdadescalestatevar}
// 
// !INTERFACE:
void FTN(registerlsmdadescalestatevar)(char *j, void (*func)(int*, void*, void*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  descale the LSM state variables
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \item[j]
//   index of the assimilated variable
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmdadescalenode* current;
  struct lsmdadescalenode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmdadescalenode*) malloc(sizeof(struct lsmdadescalenode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmdadescale_table == NULL){
    lsmdadescale_table = pnode;
  }
  else{
    current = lsmdadescale_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lsmdadescalestatevar
// \label{lsmdadescalestatevar}
// 
// !INTERFACE:
void FTN(lsmdadescalestatevar)(char *j,int *n, void *LSM_State, void *LSM_Incr_State, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry descale the
//  LSM state variables 
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \item[n]
//   index of the nest
//  \item[LSM\_State]
//   The LSM state being descaled
//  \end{description}
//EOP
{ 
  struct lsmdadescalenode* current;
  
  current = lsmdadescale_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("descale variables routine for LSM + DA instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,LSM_State,LSM_Incr_State); 
}

//BOP
// !ROUTINE: registerlsmdaupdatestate
// \label{registerlsmdaupdatestate}
// 
// !INTERFACE:
void FTN(registerlsmdaupdatestate)(char *j, void (*func)(int*, void*, void*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  apply the LSM state increments to LSM state
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the LSM + DA instance
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmdaupdatenode* current;
  struct lsmdaupdatenode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmdaupdatenode*) malloc(sizeof(struct lsmdaupdatenode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmdaupdate_table == NULL){
    lsmdaupdate_table = pnode;
  }
  else{
    current = lsmdaupdate_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lsmdaupdatestate
// \label{lsmdaupdatestate}
// 
// !INTERFACE:
void FTN(lsmdaupdatestate)(char *j, int *n, void *LSM_State, void *LSM_Incr_State, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to apply the
//  LSM state increments to LSM State
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the LSM + DA instance
//  \item[n]
//   index of the nest
//  \item[LSM\_State]
//   The LSM state being updated
//  \item[LSM\_Incr\_State]
//   The LSM incr state
//  \end{description}
//EOP
{ 
  struct lsmdaupdatenode* current;
  
  current = lsmdaupdate_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("update state routine for LSM + DA instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,LSM_State, LSM_Incr_State);
}

//BOP
// !ROUTINE: registerwritestatevar
// \label{registerwritestatevar}
// 
// !INTERFACE:
//void FTN(registerwritestatevar)(char *j, int *j, void (*func)(int *, int*, void*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine 
//  for obtaining the specified prognostic variables from the 
//  land surface model 
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \item[j]
//   index of the variable
//  \end{description}
//EOP
//{ 
//  ft_check_index(*i, FT_MAX_LSM, "registerwritestatevar");
//  ft_check_index(*j, FT_MAX_DAVAROBS, "registerwritestatevar");
//  wrtvar[*i][*j].func = func; 
//}
//BOP
// !ROUTINE: writelsmstatevar
// \label{writelsmstatevar}
// 
// !INTERFACE:
//void FTN(writelsmstatevar)(char *j, int *j, int *ftn, int *n, void *state)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry for obtaining
//  the specified prognostic variables from the land surface model
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \item[j]
//   index of the variable
//  \item[n]
//   index of the nest
//  \item[statevars]
//   pointer to the prognostic variable array
//  \end{description}
//EOP
//{ 
//  if(wrtvar[*i][*j].func==NULL) {
//    printf("****************Error****************************\n"); 
//    printf("routine that writes state prognostic variables for lsm %d and variable %d is not defined \n",*i,*j); 
//    printf("program will segfault.....\n"); 
//    printf("****************Error****************************\n"); 
//    exit;
//  }
//  wrtvar[*i][*j].func(ftn, n, state); 
//}

//BOP
// !ROUTINE: registerlsmcplsetexport
// \label{registerlsmcplsetexport}
// 
// !INTERFACE:
void FTN(registerlsmcplsetexport)(char *j, void (*func)(int*),int len)
//  
// !DESCRIPTION:
//  Makes an entry in the registry for the routine 
//  to set export states for a LSM
// 
//  \begin{description}
//  \item[j]
//   name of the LSM + instance of run mode
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmcplsetexportnode* current;
  struct lsmcplsetexportnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmcplsetexportnode*) malloc(sizeof(struct lsmcplsetexportnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmcplsetexport_table == NULL){
    lsmcplsetexport_table = pnode;
  }
  else{
    current = lsmcplsetexport_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }


}
//BOP
// !ROUTINE: lsmcplsetexport
// \label{lsmcplsetexport}
// 
// !INTERFACE:
void FTN(lsmcplsetexport)(char *j, int *n, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine in the registry to set the 
//  export states for a LSM
// 
//  \begin{description}
//  \item[j]
//   name of the LSM + instance of run mode
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{ 
  struct lsmcplsetexportnode* current;
  
  current = lsmcplsetexport_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("set export LSM routine for LSM + CPL instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}


//BOP
// !ROUTINE: registerlsmpesetdecisionspace
// \label{registerlsmpesetdecisionspace}
// 
// !INTERFACE:
void FTN(registerlsmpesetdecisionspace)(char *j, void (*func)(void*,void*),int len)
//
//  Makes an entry in the registry for the routine 
//  to set the decision space for parameter estimation
// !DESCRIPTION: 
// 
//  \begin{description}
//  \item[j]
//   name of the LSM + instance of the runmode
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmpesetdecnode* current;
  struct lsmpesetdecnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmpesetdecnode*) malloc(sizeof(struct lsmpesetdecnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmpesetdec_table == NULL){
    lsmpesetdec_table = pnode;
  }
  else{
    current = lsmpesetdec_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lsmpesetdecisionspace
// \label{lsmpesetdecisionspace}
// 
// !INTERFACE:
void FTN(lsmpesetdecisionspace)(char *j, void *dec, void *feas, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to set the 
//  decision space for parameter estimation
//
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \item[dec]
//   decision space object
//  \item[feas]
//   feasible space object
//  \end{description}
//EOP
{ 
  struct lsmpesetdecnode* current;
  
  current = lsmpesetdec_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("set decision space routine for LSM %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(dec,feas);
}

//BOP
// !ROUTINE: registerlsmpegetdecisionspace
// \label{registerlsmpegetdecisionspace}
// 
// !INTERFACE:
void FTN(registerlsmpegetdecisionspace)(char *j, void (*func)(void*),int len)
//  Makes an entry in the registry for the routine 
//  to get the decision space for parameter estimation
// !DESCRIPTION: 
// 
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmpegetdecnode* current;
  struct lsmpegetdecnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmpegetdecnode*) malloc(sizeof(struct lsmpegetdecnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmpegetdec_table == NULL){
    lsmpegetdec_table = pnode;
  }
  else{
    current = lsmpegetdec_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lsmpegetdecisionspace
// \label{lsmpegetdecisionspace}
// 
// !INTERFACE:
void FTN(lsmpegetdecisionspace)(char *j, void *state, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to get the 
//  decision space for parameter estimation
//
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \item[statevars]
//   pointer to the prognostic variable array
//  \end{description}
//EOP
{ 
  struct lsmpegetdecnode* current;
  
  current = lsmpegetdec_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("get decision space routine for LSM %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(state);
}

//BOP
// !ROUTINE: registerlsmpesetupdecisionspace
// \label{registerlsmpesetupdecisionspace}
// 
// !INTERFACE:
void FTN(registerlsmpesetupdecisionspace)(char *j, void (*func)(void*, void*),int len)
//  
// !DESCRIPTION: 
//  Method to registry an interface implementation for 
//   setting up the LSM decision space. 
//
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmpesetupdecnode* current;
  struct lsmpesetupdecnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmpesetupdecnode*) malloc(sizeof(struct lsmpesetupdecnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmpesetupdec_table == NULL){
    lsmpesetupdec_table = pnode;
  }
  else{
    current = lsmpesetupdec_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lsmpesetupdecisionspace
// \label{lsmpesetupdecisionspace}
// 
// !INTERFACE:
void FTN(lsmpesetupdecisionspace)(char *j, void *dec, void *feas, int len)
//  
// !DESCRIPTION: 
//  Method to setup  the LSM decision space
//
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \item[dec]
//   decision space object
//  \item[feas]
//   feasible space object
//  \end{description}
//EOP
{ 
  struct lsmpesetupdecnode* current;
  
  current = lsmpesetupdec_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("setup decision space outine for LSM %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(dec,feas); 
}

//BOP
// !ROUTINE: registerlsmpesetupobspred
// \label{registerlsmpesetupobspred}
// 
// !INTERFACE:
void FTN(registerlsmpesetupobspred)(char *j, void (*func)(void*),int len)
//  
// !DESCRIPTION: 
//  registers the method to initialize the obs pred for the PE observation
//  (model simulated value of the observation)
//  
//  \begin{description}
//  \item[j]
//   name of the LSM + PE instance
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmpesetupobsprednode* current;
  struct lsmpesetupobsprednode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmpesetupobsprednode*) malloc(sizeof(struct lsmpesetupobsprednode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmpesetupobspred_table == NULL){
    lsmpesetupobspred_table = pnode;
  }
  else{
    current = lsmpesetupobspred_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lsmpesetupobspredspace
// \label{lsmpesetupobspredspace}
// 
// !INTERFACE:
void FTN(lsmpesetupobspredspace)(char *j, void *pepred, int len)
//  
// !DESCRIPTION: 
//   invokes the method to compute the obs pred for the PE observation
//  (model simulated value of the observation)
//
//  \begin{description}
//  \item[j]
//   name of the LSM + PE instance
//  \item[pepred]
//   object containing the PE obspred 
//  \end{description}
//EOP
{ 
  struct lsmpesetupobsprednode* current;
  
  current = lsmpesetupobspred_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("setup obspred routine for LSM + PE instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(pepred); 
}

//BOP
// !ROUTINE: registerlsmpegetobspred
// \label{registerlsmpegetobspred}
// 
// !INTERFACE:
void FTN(registerlsmpegetobspred)(char *j, void (*func)(void*),int len)
//  
// !DESCRIPTION: 
//  registers the method to compute the obs pred for the PE observation
//  (model simulated value of the observation)
//  
//  \begin{description}
//  \item[j]
//   name of the LSM
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmpeobsprednode* current;
  struct lsmpeobsprednode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmpeobsprednode*) malloc(sizeof(struct lsmpeobsprednode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmpeobspred_table == NULL){
    lsmpeobspred_table = pnode;
  }
  else{
    current = lsmpeobspred_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lsmpegetobspred
// \label{lsmpegetobspred}
// 
// !INTERFACE:
void FTN(lsmpegetobspred)(char *j, void *pepred, int len)
//  
// !DESCRIPTION: 
//   invokes the method to compute the obs pred for the PE observation
//  (model simulated value of the observation)
//
//  \begin{description}
//  \item[j]
//   name of the LSM + PE instance
//  \item[pepred]
//   object containing the PE obspred 
//  \end{description}
//EOP
{ 
  struct lsmpeobsprednode* current;
  
  current = lsmpeobspred_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("PE obspred routine for LSM + PE instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(pepred); 
}


//BOP
// !ROUTINE: registerlsmroutinggetrunoff
// \label{registerlsmroutinggetrunoff}
// 
// !INTERFACE:
void FTN(registerlsmroutinggetrunoff)(char *j, void (*func)(int*),int len)
//  
// !DESCRIPTION: 
//  creates an entry in the registry for the routine to 
//  return runoff fields from a LSM.   
//
//  \begin{description}
//  \item[j]
//   name of the LSM + routing instance
//  \end{description}
//EOP
{ 
  int len1;
  struct lsmroutinggetrunoffnode* current;
  struct lsmroutinggetrunoffnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmroutinggetrunoffnode*) malloc(sizeof(struct lsmroutinggetrunoffnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmroutinggetrunoff_table == NULL){
    lsmroutinggetrunoff_table = pnode;
  }
  else{
    current = lsmroutinggetrunoff_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}

//BOP
// !ROUTINE: lsmroutinggetrunoff
// \label{lsmroutinggetrunoff}
// 
// !INTERFACE:
void FTN(lsmroutinggetrunoff)(char *j, int *n, int len)
//  
// !DESCRIPTION: 
//  Invokes the registered routine that returns the
//  runoff fields from a LSM. 
//
//  \begin{description}
//  \item[j]
//   name of the LSM + routing instance
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{ 
  struct lsmroutinggetrunoffnode* current;
  
  current = lsmroutinggetrunoff_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("get runoff routine for LSM + routing instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}

//BOP
// !ROUTINE: registerlsmroutinggetsws
// \label{registerlsmroutinggetsws}
//
// !INTERFACE:
void FTN(registerlsmroutinggetsws)(char *j, void (*func)(int*),int len)
//
// !DESCRIPTION:
//  creates an entry in the registry for the routine to
//  set the surface water storage fields from the routing
//  model within the LSM
//
//  \begin{description}
//  \item[j]
//   name of the LSM + routing instance
//  \end{description}
//EOP
{
  int len1;
  struct lsmroutinggetswsnode* current;
  struct lsmroutinggetswsnode* pnode;
  // create node

  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmroutinggetswsnode*) malloc(sizeof(struct lsmroutinggetswsnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL;

  if(lsmroutinggetsws_table == NULL){
    lsmroutinggetsws_table = pnode;
  }
  else{
    current = lsmroutinggetsws_table;
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode;
  }
}

//BOP
// !ROUTINE: lsmroutinggetsws
// \label{lsmroutinggetsws}
//
// !INTERFACE:
void FTN(lsmroutinggetsws)(char *j, int *n, int len)
//
// !DESCRIPTION:
//  Invokes the registered routine that sets the
//  surface water storage fields from the routing model
//  within the LSM
//
//  \begin{description}
//  \item[j]
//   name of the LSM + routing instance
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{
  struct lsmroutinggetswsnode* current;

  current = lsmroutinggetsws_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n");
      printf("set sws routine for LSM + routing instance %s is not defined\n",j);
      printf("program will seg fault.....\n");
      printf("****************Error****************************\n");
    }
  }
  current->func(n);
}

//BOP
// !ROUTINE: registerlsm2rtm
// \label{registerlsm2rtm}
// 
// !INTERFACE:
void FTN(registerlsm2rtm)(char *j,  void (*func)(int*, void*), int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  transfer surface properties to the RTM data structure 
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the RTM
//. \end{description}
//EOP
{ 
  int len1;
  struct lsm2rtmnode* current;
  struct lsm2rtmnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsm2rtmnode*) malloc(sizeof(struct lsm2rtmnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsm2rtm_table == NULL){
    lsm2rtm_table = pnode;
  }
  else{
    current = lsm2rtm_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: lsm2rtm
// \label{lsm2rtm}
// 
// !INTERFACE:
void FTN(lsm2rtm)(char *j, int *index, void *sfc, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to
//  map surface properties to the RTM data structure
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the RTM
// \end{description}
//EOP
{ 
  struct lsm2rtmnode* current;
  
  current = lsm2rtm_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("lsm2rtm routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(index,sfc); 
}

//-------------------- Irrigation updates -------------------------


//BOP
// !ROUTINE: registerlsmirrigationgetstates
// \label{registerlsmirrigationgetstates}
// 
// !INTERFACE:
void FTN(registerlsmirrigationgetstates)(char *j,  void (*func)(int*, void*), int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  initialize the LSM irrigation states
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the irrigation scheme
//. \end{description}
//EOP
{ 
  int len1;
  struct lsmirriggetnode* current;
  struct lsmirriggetnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lsmirriggetnode*) malloc(sizeof(struct lsmirriggetnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmirrigget_table == NULL){
    lsmirrigget_table = pnode;
  }
  else{
    current = lsmirrigget_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: getirrigationlsmstates
// \label{getirrigationlsmstates}
// 
// !INTERFACE:
void FTN(getirrigationlsmstates)(char *j, int *n, void *irrigstate, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to
//  getialize the LSM irrigation states
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the irrigation model
// \end{description}
//EOP
{ 
  struct lsmirriggetnode* current;
  
  current = lsmirrigget_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("get routine for lsm and irrigation scheme %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n, irrigstate); 
}




