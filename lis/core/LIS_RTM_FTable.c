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
// !MODULE: LIS_RTM_FTable
//  
//
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing the operations of different 
//  radiative transfer models
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"

struct rtminitnode
{ 
  char *name;
  void (*func)(int*);

  struct rtminitnode* next;
} ;
struct rtminitnode* rtminit_table = NULL; 

struct rtmf2tnode
{ 
  char *name;
  void (*func)(int*);

  struct rtmf2tnode* next;
} ;
struct rtmf2tnode* rtmf2t_table = NULL; 

struct geo2rtmnode
{ 
  char *name;
  void (*func)(int*);

  struct geo2rtmnode* next;
} ;
struct geo2rtmnode* geo2rtm_table = NULL; 

struct rtmrunnode
{ 
  char *name;
  void (*func)(int*);

  struct rtmrunnode* next;
} ;
struct rtmrunnode* rtmrun_table = NULL; 

//Parameter estimation
struct rtmpesetdecnode
{ 
  char *name;
  void (*func)(void*, void*);

  struct rtmpesetdecnode* next;
} ;
struct rtmpesetdecnode* rtmpesetdec_table = NULL;

struct rtmpegetdecnode
{ 
  char *name;
  void (*func)(void*);

  struct rtmpegetdecnode* next;
} ;
struct rtmpegetdecnode* rtmpegetdec_table = NULL;

struct rtmpesetupdecnode
{ 
  char *name;
  void (*func)(void*, void*);

  struct rtmpesetupdecnode* next;
} ;
struct rtmpesetupdecnode* rtmpesetupdec_table = NULL;

struct rtmpesetprednode
{ 
  char *name;
  void (*func)(void*);

  struct rtmpesetprednode* next;
} ;
struct rtmpesetprednode* rtmpesetpred_table = NULL;

struct rtmpesetupobsprednode
{ 
  char *name;
  void (*func)(void*);

  struct rtmpesetupobsprednode* next;
} ;
struct rtmpesetupobsprednode* rtmpesetupobspred_table = NULL;

struct rtmpeobsprednode
{ 
  char *name;
  void (*func)(void*);

  struct rtmpeobsprednode* next;
} ;
struct rtmpeobsprednode* rtmpeobspred_table = NULL;

struct rtmresetnode
{ 
  char *name;
  void (*func)();

  struct rtmresetnode* next;
} ;
struct rtmresetnode* rtmreset_table = NULL; 

//BOP
// !ROUTINE: registerrtminit
// \label{registerrtminit}
// 
// !INTERFACE:
void FTN(registerrtminit)(char *j, void (*func)(int*), int len)
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
  struct rtminitnode* current;
  struct rtminitnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct rtminitnode*) malloc(sizeof(struct rtminitnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(rtminit_table == NULL){
    rtminit_table = pnode;
  }
  else{
    current = rtminit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: rtminit
// \label{rtminit}
// 
// !INTERFACE:
void FTN(rtminitialize)(char *j, int *index, int len)
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
  struct rtminitnode* current;
  
  current = rtminit_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("init routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(index); 
}

//BOP
// !ROUTINE: registerrtmf2t
// \label{registerrtmf2t}
// 
// !INTERFACE:
void FTN(registerrtmf2t)(char *j, void (*func)(int*), int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  transfer atmospheric profiles to the RTM data structure 
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the RTM
//. \end{description}
//EOP
{ 
  int len1;
  struct rtmf2tnode* current;
  struct rtmf2tnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct rtmf2tnode*) malloc(sizeof(struct rtmf2tnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(rtmf2t_table == NULL){
    rtmf2t_table = pnode;
  }
  else{
    current = rtmf2t_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: rtmf2t
// \label{rtmf2t}
// 
// !INTERFACE:
void FTN(rtmf2t)(char *j, int *index, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to
//  transfer atmospheric profiles to the RTM 
//  data structure
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the RTM
// \end{description}
//EOP
{ 
  struct rtmf2tnode* current;
  
  current = rtmf2t_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("rtmf2t routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(index); 
}


//BOP
// !ROUTINE: registergeometry2rtm
// \label{registergeometry2rtm}
// 
// !INTERFACE:
void FTN(registergeometry2rtm)(char *j, void (*func)(int*), int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  specify the sensor geometries in the RTM
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the RTM
//. \end{description}
//EOP
{ 
  int len1;
  struct geo2rtmnode* current;
  struct geo2rtmnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct geo2rtmnode*) malloc(sizeof(struct geo2rtmnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(geo2rtm_table == NULL){
    geo2rtm_table = pnode;
  }
  else{
    current = geo2rtm_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: geometry2rtm
// \label{geometry2rtm}
// 
// !INTERFACE:
void FTN(geometry2rtm)(char *j, int *index, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to
//  map sensor geometry information to the RTM 
//  data structure
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the RTM
// \end{description}
//EOP
{ 
  struct geo2rtmnode* current;
  
  current = geo2rtm_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("geometry2rtm routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(index); 
}


//BOP
// !ROUTINE: registerrtmrun
// \label{registerrtmrun}
// 
// !INTERFACE:
void FTN(registerrtmrun)(char *j, void (*func)(int*), int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  specify the forward model integration step of the RTM
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the RTM
//. \end{description}
//EOP
{ 
  int len1;
  struct rtmrunnode* current;
  struct rtmrunnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct rtmrunnode*) malloc(sizeof(struct rtmrunnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(rtmrun_table == NULL){
    rtmrun_table = pnode;
  }
  else{
    current = rtmrun_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: rtmrun
// \label{rtmrun}
// 
// !INTERFACE:
void FTN(rtmrun)(char *j, int *index, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to
//  call the forward model integration step of 
//  the RTM
//  
//  The arguments are: 
//  \begin{description}
//   \item[j]
//   name of the RTM
// \end{description}
//EOP
{ 
  struct rtmrunnode* current;
  
  current = rtmrun_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("RTM run routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(index); 
}


//BOP
// !ROUTINE: registerrtmpesetdecisionspace
// \label{registerrtmpesetdecisionspace}
// 
// !INTERFACE:
void FTN(registerrtmpesetdecisionspace)(char *j, void (*func)(void*,void*),int len)
//
//  Makes an entry in the registry for the routine 
//  to set the decision space for parameter estimation
// !DESCRIPTION: 
// 
//  \begin{description}
//  \item[j]
//   name of the RTM + instance of the runmode
//  \end{description}
//EOP
{ 
  int len1;
  struct rtmpesetdecnode* current;
  struct rtmpesetdecnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct rtmpesetdecnode*) malloc(sizeof(struct rtmpesetdecnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(rtmpesetdec_table == NULL){
    rtmpesetdec_table = pnode;
  }
  else{
    current = rtmpesetdec_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: rtmpesetdecisionspace
// \label{rtmpesetdecisionspace}
// 
// !INTERFACE:
void FTN(rtmpesetdecisionspace)(char *j, void *dec, void *feas, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to set the 
//  decision space for parameter estimation
//
//  \begin{description}
//  \item[j]
//   name of the RTM
//  \item[dec]
//   decision space object
//  \item[feas]
//   feasible space object
//  \end{description}
//EOP
{ 
  struct rtmpesetdecnode* current;
  
  current = rtmpesetdec_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("set decision space routine for RTM %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(dec,feas);
}

//BOP
// !ROUTINE: registerrtmpegetdecisionspace
// \label{registerrtmpegetdecisionspace}
// 
// !INTERFACE:
void FTN(registerrtmpegetdecisionspace)(char *j, void (*func)(void*),int len)
//  Makes an entry in the registry for the routine 
//  to get the decision space for parameter estimation
// !DESCRIPTION: 
// 
//  \begin{description}
//  \item[j]
//   name of the RTM
//  \end{description}
//EOP
{ 
  int len1;
  struct rtmpegetdecnode* current;
  struct rtmpegetdecnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct rtmpegetdecnode*) malloc(sizeof(struct rtmpegetdecnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(rtmpegetdec_table == NULL){
    rtmpegetdec_table = pnode;
  }
  else{
    current = rtmpegetdec_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: rtmpegetdecisionspace
// \label{rtmpegetdecisionspace}
// 
// !INTERFACE:
void FTN(rtmpegetdecisionspace)(char *j, void *state, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to get the 
//  decision space for parameter estimation
//
//  \begin{description}
//  \item[j]
//   name of the RTM
//  \item[statevars]
//   pointer to the prognostic variable array
//  \end{description}
//EOP
{ 
  struct rtmpegetdecnode* current;
  
  current = rtmpegetdec_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("get decision space routine for RTM %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(state);
}

//BOP
// !ROUTINE: registerrtmpesetupdecisionspace
// \label{registerrtmpesetupdecisionspace}
// 
// !INTERFACE:
void FTN(registerrtmpesetupdecisionspace)(char *j, void (*func)(void*, void*),int len)
//  
// !DESCRIPTION: 
//  Method to registry an interface implementation for 
//   setting up the RTM decision space. 
//
//  \begin{description}
//  \item[j]
//   name of the RTM
//  \end{description}
//EOP
{ 
  int len1;
  struct rtmpesetupdecnode* current;
  struct rtmpesetupdecnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct rtmpesetupdecnode*) malloc(sizeof(struct rtmpesetupdecnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(rtmpesetupdec_table == NULL){
    rtmpesetupdec_table = pnode;
  }
  else{
    current = rtmpesetupdec_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: rtmpesetupdecisionspace
// \label{rtmpesetupdecisionspace}
// 
// !INTERFACE:
void FTN(rtmpesetupdecisionspace)(char *j, void *dec, void *feas, int len)
//  
// !DESCRIPTION: 
//  Method to setup  the RTM decision space
//
//  \begin{description}
//  \item[j]
//   name of the RTM
//  \item[dec]
//   decision space object
//  \item[feas]
//   feasible space object
//  \end{description}
//EOP
{ 
  struct rtmpesetupdecnode* current;
  
  current = rtmpesetupdec_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("setup decision space outine for RTM %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(dec,feas); 
}

//BOP
// !ROUTINE: registerrtmpesetupobspred
// \label{registerrtmpesetupobspred}
// 
// !INTERFACE:
void FTN(registerrtmpesetupobspred)(char *j, void (*func)(void*),int len)
//  
// !DESCRIPTION: 
//  registers the method to initialize the obs pred for the PE observation
//  (model simulated value of the observation)
//  
//  \begin{description}
//  \item[j]
//   name of the RTM + PE instance
//  \end{description}
//EOP
{ 
  int len1;
  struct rtmpesetupobsprednode* current;
  struct rtmpesetupobsprednode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct rtmpesetupobsprednode*) malloc(sizeof(struct rtmpesetupobsprednode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(rtmpesetupobspred_table == NULL){
    rtmpesetupobspred_table = pnode;
  }
  else{
    current = rtmpesetupobspred_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: rtmpesetupobspredspace
// \label{rtmpesetupobspredspace}
// 
// !INTERFACE:
void FTN(rtmpesetupobspredspace)(char *j, void *pepred, int len)
//  
// !DESCRIPTION: 
//   invokes the method to compute the obs pred for the PE observation
//  (model simulated value of the observation)
//
//  \begin{description}
//  \item[j]
//   name of the RTM + PE instance
//  \item[pepred]
//   object containing the PE obspred 
//  \end{description}
//EOP
{ 
  struct rtmpesetupobsprednode* current;
  
  current = rtmpesetupobspred_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("setup obspred routine for RTM + PE instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(pepred); 
}

//BOP
// !ROUTINE: registerrtmpegetobspred
// \label{registerrtmpegetobspred}
// 
// !INTERFACE:
void FTN(registerrtmpegetobspred)(char *j, void (*func)(void*),int len)
//  
// !DESCRIPTION: 
//  registers the method to compute the obs pred for the PE observation
//  (model simulated value of the observation)
//  
//  \begin{description}
//  \item[j]
//   name of the RTM
//  \end{description}
//EOP
{ 
  int len1;
  struct rtmpeobsprednode* current;
  struct rtmpeobsprednode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct rtmpeobsprednode*) malloc(sizeof(struct rtmpeobsprednode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(rtmpeobspred_table == NULL){
    rtmpeobspred_table = pnode;
  }
  else{
    current = rtmpeobspred_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: rtmpegetobspred
// \label{rtmpegetobspred}
// 
// !INTERFACE:
void FTN(rtmpegetobspred)(char *j, void *pepred, int len)
//  
// !DESCRIPTION: 
//   invokes the method to compute the obs pred for the PE observation
//  (model simulated value of the observation)
//
//  \begin{description}
//  \item[j]
//   name of the RTM + PE instance
//  \item[pepred]
//   object containing the PE obspred 
//  \end{description}
//EOP
{ 
  struct rtmpeobsprednode* current;
  
  current = rtmpeobspred_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("PE obspred routine for RTM + PE instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(pepred); 
}


//BOP
// !ROUTINE: registerrtmreset
// \label{registerrtmreset}
//  
// !INTERFACE:
void FTN(registerrtmreset)(char *j, void (*func)(), int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine to 
//  perform reset of RTM variables. 
// 
//EOP
{ 
  int len1;
  struct rtmresetnode* current;
  struct rtmresetnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct rtmresetnode*) malloc(sizeof(struct rtmresetnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(rtmreset_table == NULL){
    rtmreset_table = pnode;
  }
  else{
    current = rtmreset_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: resetrtm
// \label{resetrtm}
//
// !INTERFACE:
void FTN(resetrtm)(char *j, int len)
//
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to perform
//  rtm model reset
// 
//  \begin{description}
//  \item[j]
//   name of the RTM
//  \end{description}
//EOP
{ 
  struct rtmresetnode* current;
  
  current = rtmreset_table;
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





