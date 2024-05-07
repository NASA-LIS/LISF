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
// !MODULE: LIS_optUE_FTable
//  
// !DESCRIPTION:
//  Function table registries for storing the interface
//  implementations for the operation of optimization   
//  uncertainty estimation algorithms and observations.
//
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"
struct optueinitnode
{ 
  char *name;
  void (*func)();

  struct optueinitnode* next;
} ;
struct optueinitnode* optueinit_table = NULL; 

struct optuesetnode
{ 
  char *name;
  void (*func)();

  struct optuesetnode* next;
} ;
struct optuesetnode* optueset_table = NULL; 

struct optuecchknode
{ 
  char *name;
  void (*func)(void*);

  struct optuecchknode* next;
} ;
struct optuecchknode* optuecchk_table = NULL; 

struct optuerunnode
{ 
  char *name;
  void (*func)();

  struct optuerunnode* next;
} ;
struct optuerunnode* optuerun_table = NULL; 

struct optuegetdecnode
{ 
  char *name;
  void (*func)(int*);

  struct optuegetdecnode* next;
} ;
struct optuegetdecnode* optuegetdec_table = NULL; 

struct optuesetdecnode
{ 
  char *name;
  void (*func)(int*,void*);

  struct optuesetdecnode* next;
} ;
struct optuesetdecnode* optuesetdec_table = NULL; 

struct optuegetnpnode
{ 
  char *name;
  void (*func)(void*);

  struct optuegetnpnode* next;
} ;
struct optuegetnpnode* optuegetnp_table = NULL; 

struct optuerstnode
{ 
  char *name;
  void (*func)();

  struct optuerstnode* next;
} ;
struct optuerstnode* optuerst_table = NULL; 

struct peobsinitnode
{ 
  char *name;
  void (*func)(void*);

  struct peobsinitnode* next;
} ;
struct peobsinitnode* peobsinit_table = NULL; 

struct peobsgetnode
{ 
  char *name;
  void (*func)(void*);

  struct peobsgetnode* next;
} ;
struct peobsgetnode* peobsget_table = NULL; 

struct peobswritenode
{ 
  char *name;
  void (*func)(void*);

  struct peobswritenode* next;
} ;
struct peobswritenode* peobswrite_table = NULL; 

struct peobsresetnode
{ 
  char *name;
  void (*func)(void*);

  struct peobsresetnode* next;
} ;
struct peobsresetnode* peobsreset_table = NULL; 

struct objfuncinitnode
{ 
  char *name;
  void (*func)();

  struct objfuncinitnode* next;
} ;
struct objfuncinitnode* objfuncinit_table = NULL; 

struct objfunccomputenode
{ 
  char *name;
  void (*func)();

  struct objfunccomputenode* next;
} ;
struct objfunccomputenode* objfunccompute_table = NULL; 

struct objfuncupdatenode
{ 
  char *name;
  void (*func)();

  struct objfuncupdatenode* next;
} ;
struct objfuncupdatenode* objfuncupdate_table = NULL; 

struct objfuncresetnode
{ 
  char *name;
  void (*func)();

  struct objfuncresetnode* next;
} ;
struct objfuncresetnode* objfuncreset_table = NULL; 

struct objfuncevalnode
{ 
  char *name;
  void (*func)();

  struct objfuncevalnode* next;
} ;
struct objfuncevalnode* objfunceval_table = NULL; 
//BOP
// !ROUTINE: registeroptuealginit
// \label{registeroptuealginit}
//  
// 
// !INTERFACE:
void FTN(registeroptuealginit)(char *j, void (*func)(),int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine that
//  sets up structures for initializing optimization 
//  / uncertainty estimation algorithms
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   index of the optmization/UE algorithm
//  \end{description}
//EOP
{ 
  int len1;
  struct optueinitnode* current;
  struct optueinitnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct optueinitnode*) malloc(sizeof(struct optueinitnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(optueinit_table == NULL){
    optueinit_table = pnode;
  }
  else{
    current = optueinit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }


}
//BOP
// !ROUTINE: optuealginit
// \label{optuealginit}
// 
// !INTERFACE:
void FTN(optuealginit)(char *j,int len)
//  
// !DESCRIPTION: 
// Invokes the routine from the registry to set up 
// structures for initializing optimization or
// uncertainty estimation algorithm
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the optimization/UE algorithm
//  \end{description}
//
//EOP
{ 
  struct optueinitnode* current;
  
  current = optueinit_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("init routine for  %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}

//BOP
// !ROUTINE: registeroptuealgsetup
// \label{registeroptuealgsetup}
//  
// 
// !INTERFACE:
void FTN(registeroptuealgsetup)(char *j, void (*func)(),int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine that
//  sets up structures for optimization/ uncertainty estimation  
//  algorithms
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   index of the optmization/UE algorithm
//  \end{description}
//EOP
{ 
  int len1;
  struct optuesetnode* current;
  struct optuesetnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct optuesetnode*) malloc(sizeof(struct optuesetnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(optueset_table == NULL){
    optueset_table = pnode;
  }
  else{
    current = optueset_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: optuealgsetup
// \label{optuealgsetup}
// 
// !INTERFACE:
void FTN(optuealgsetup)(char *j,int len)
//  
// !DESCRIPTION: 
// Invokes the routine from the registry to set up 
// structures for optimization/uncertainty estimation
// algorithm
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the optimization/uncertainty estimation
//   algorithm
//  \end{description}
//
//EOP
{ 
  struct optuesetnode* current;
  
  current = optueset_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("setup routine for OPT/UE %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}


//BOP
// !ROUTINE: registeroptueconvergencecheck
// \label{registeroptueconvergencecheck}
//  
// 
// !INTERFACE:
void FTN(registeroptueconvergencecheck)(char *j, void (*func)(void*),int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine that
//  checks for convergence for an optimization algorithm
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the optimization algorithm
//  \end{description}
//EOP
{ 
  int len1;
  struct optuecchknode* current;
  struct optuecchknode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct optuecchknode*) malloc(sizeof(struct optuecchknode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(optuecchk_table == NULL){
    optuecchk_table = pnode;
  }
  else{
    current = optuecchk_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: checkconvergence
// \label{checkconvergence}
// 
// !INTERFACE:
void FTN(checkconvergence)(char *j, void* check,int len)
//  
// !DESCRIPTION: 
// Invokes the routine from the registry to check for 
// convergence of an optimization algorithm. 
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the optimization algorithm
//  \end{description}
//
//EOP
{ 
  struct optuecchknode* current;
  
  current = optuecchk_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("check convergence routine for OPT/UE %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(check); 
}

//BOP
// !ROUTINE: registeroptuealgrun
// \label{registeroptuealgrun}
//  
// 
// !INTERFACE:
void FTN(registeroptuealgrun)(char *j, void (*func)(),int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine that
//  executes the optimization algorithm run method
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the optimization algorithm. 
//  \end{description}
//EOP
{ 
  int len1;
  struct optuerunnode* current;
  struct optuerunnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct optuerunnode*) malloc(sizeof(struct optuerunnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(optuerun_table == NULL){
    optuerun_table = pnode;
  }
  else{
    current = optuerun_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: runoptue
// \label{runoptue}
// 
// !INTERFACE:
void FTN(runoptue)(char *j,int len)
//  
// !DESCRIPTION: 
// Invokes the routine from the registry to execute
// the optimization algorithm run method
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the optimization algorithm
//  \end{description}
//
//EOP
{ 
  struct optuerunnode* current;
  
  current = optuerun_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("run routine for OPT/UE  %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}

//BOP
// !ROUTINE: registergetdecisionspace
// \label{registergetdecisionspace}
//  
// 
// !INTERFACE:
void FTN(registeroptuegetdecisionspace)(char *j, void (*func)(int*),int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine that
//  retrieves decision space values from the optimization
//  algorithm datastructures
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the optimization algorithm
//  \end{description}
//EOP
{ 
  int len1;
  struct optuegetdecnode* current;
  struct optuegetdecnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct optuegetdecnode*) malloc(sizeof(struct optuegetdecnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(optuegetdec_table == NULL){
    optuegetdec_table = pnode;
  }
  else{
    current = optuegetdec_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: 
// \label{getoptuealgdecspace}
// 
// !INTERFACE:
void FTN(getoptuealgdecspace)(char *j, int *n,int len)
//  
// !DESCRIPTION: 
// Invokes the routine from the registry to invoke
// an objective function method
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the optimization algorithm
//  \item[n]
//   index of the nest
//  \end{description}
//
//EOP
{ 
  struct optuegetdecnode* current;
  
  current = optuegetdec_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("get decision space routine for OPT/UE %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}

//BOP
// !ROUTINE: registeroptuesetdecisionspace
// \label{registeroptuesetdecisionspace}
//  
// 
// !INTERFACE:
void FTN(registeroptuesetdecisionspace)(char *j, void (*func)(int*, void*),int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine that
//  assigns decision space values to the optimization
//  algorithm datastructures
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the optimization algorithm
//  \end{description}
//EOP
{ 
  int len1;
  struct optuesetdecnode* current;
  struct optuesetdecnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct optuesetdecnode*) malloc(sizeof(struct optuesetdecnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(optuesetdec_table == NULL){
    optuesetdec_table = pnode;
  }
  else{
    current = optuesetdec_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: 
// \label{setoptuealgdecspace}
// 
// !INTERFACE:
void FTN(setoptuealgdecspace)(char *j, int *n, void *decspace,int len)
//  
// !DESCRIPTION: 
// Invokes the routine from the registry to set
// a given decision space to the optimization algorithm
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the optimization algorithm
//  \item[n]
//   index of the nest
//  \item[decspace]
//   decision space object
//  \end{description}
//
//EOP
{ 
  struct optuesetdecnode* current;
  
  current = optuesetdec_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("routine for setting the decision space of OPT/UE algorithm %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,decspace); 
}


//BOP
// !ROUTINE: registeroptuegetnparam
// \label{registeroptuegetnparam}
//  
// 
// !INTERFACE:
void FTN(registeroptuegetnparam)(char *j, void (*func)(void*),int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine that
//  retrieves the size of the decision space 
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the optimization algorithm
//  \end{description}
//EOP
{ 
  int len1;
  struct optuegetnpnode* current;
  struct optuegetnpnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct optuegetnpnode*) malloc(sizeof(struct optuegetnpnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(optuegetnp_table == NULL){
    optuegetnp_table = pnode;
  }
  else{
    current = optuegetnp_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: 
// \label{getoptuealgnparam}
// 
// !INTERFACE:
void FTN(getoptuealgnparam)(char *j, void *nparam,int len)
//  
// !DESCRIPTION: 
// Invokes the routine from the registry to retrieve 
// the size of the decision space
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the optimization algorithm
//  \item[nparam]
//   size of the decision space
//  \end{description}
//
//EOP
{ 
  struct optuegetnpnode* current;
  
  current = optuegetnp_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("routine for getdecisionspace size for OPT/UE %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(nparam); 
}

//BOP
// !ROUTINE: registeroptuereadrestart
// \label{registeroptuereadrestart}
//  
// 
// !INTERFACE:
void FTN(registeroptuereadrestart)(char *j, void (*func)(void*),int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine that
//  reads a restart file for the OPT/UE algorithm. 
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the optimization algorithm
//  \end{description}
//EOP
{ 
  int len1;
  struct optuerstnode* current;
  struct optuerstnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct optuerstnode*) malloc(sizeof(struct optuerstnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(optuerst_table == NULL){
    optuerst_table = pnode;
  }
  else{
    current = optuerst_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: 
// \label{optuereadrestart}
// 
// !INTERFACE:
void FTN(optuereadrestart)(char *j,int len)
//  
// !DESCRIPTION: 
// Invokes the routine from the registry to retrieve 
// the size of the decision space
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the optimization algorithm
//  \item[rstflag]
//   flag indicating the type of optue start (cold/restart)
//  \end{description}
//
//EOP
{ 
  struct optuerstnode* current;
  
  current = optuerst_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("readrestart routine for OPT/UE algorithm %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}

//BOP
// !ROUTINE: registerpeobssetup
// \label{registerpeobssetup}
//  
// 
// !INTERFACE:
void FTN(registerpeobssetup)(char *j, void (*func)(void*),int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine that
//  sets up structures for handling observation data for 
//  optimization.
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the PE obs source
//  \end{description}
//EOP
{ 
  int len1;
  struct peobsinitnode* current;
  struct peobsinitnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct peobsinitnode*) malloc(sizeof(struct peobsinitnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(peobsinit_table == NULL){
    peobsinit_table = pnode;
  }
  else{
    current = peobsinit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: setuppeobsspace
// \label{setuppeobsspace}
// 
// !INTERFACE:
void FTN(setuppeobsspace)(char *j, void *obs,int len)
//  
// !DESCRIPTION: 
// Invokes the routine from the registry to set up 
// structures for handling observation data for 
// optimization
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the PE obs source 
//  \item[obs]
//   ESMF state that contain the observations
//  \end{description}
//
//EOP
{ 
  struct peobsinitnode* current;
  
  current = peobsinit_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("setup PE obs routine for  %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(obs); 
}

//BOP
// !ROUTINE: registerpeobsreset
// \label{registerpeobsreset}
//  
// 
// !INTERFACE:
void FTN(registerpeobsreset)(char *j, void (*func)(void*),int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine that
//  sets up structures for handling observation data for 
//  optimization.
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the PE obs source
//  \end{description}
//EOP
{ 
  int len1;
  struct peobsresetnode* current;
  struct peobsresetnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct peobsresetnode*) malloc(sizeof(struct peobsresetnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(peobsreset_table == NULL){
    peobsreset_table = pnode;
  }
  else{
    current = peobsreset_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: resetpeobsspace
// \label{resetpeobsspace}
// 
// !INTERFACE:
void FTN(resetpeobsspace)(char *j, void *obs,int len)
//  
// !DESCRIPTION: 
// Invokes the routine from the registry to reset
// structures for each observation data for 
// optimization
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the PE obs source
//  \item[obs]
//   ESMF state that contain the observations
//  \end{description}
//
//EOP
{ 
  struct peobsresetnode* current;
  
  current = peobsreset_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("reset PE obs routine for  %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(obs); 
}

//BOP
// !ROUTINE: registergetpeobs
// \label{registergetpeobs}
//  
// 
// !INTERFACE:
void FTN(registergetpeobs)(char *j, void (*func)(void*),int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine that
//  sets up structures for reading observation data for 
//  optimization
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the PE obs source
//  \end{description}
//EOP
{ 
  int len1;
  struct peobsgetnode* current;
  struct peobsgetnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct peobsgetnode*) malloc(sizeof(struct peobsgetnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(peobsget_table == NULL){
    peobsget_table = pnode;
  }
  else{
    current = peobsget_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: getpeobs
// \label{getpeobs}
// 
// !INTERFACE:
void FTN(getpeobs)(char *j, void *peobs,int len)
//  
// !DESCRIPTION: 
// Invokes the routine from the registry to set up 
// structures for reading observation data for 
// optimization 
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the PE observation source
//  \item[peobs]
//   ESMF state that contain the observations
//  \end{description}
//
//EOP
{ 
  struct peobsgetnode* current;
  
  current = peobsget_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("get PE obs routine for  %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(peobs); 
}


//BOP
// !ROUTINE: registerwritepeobs
// \label{registerwritepeobs}
//  
// 
// !INTERFACE:
void FTN(registerwritepeobs)(char *j, void (*func)(void*),int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine that
//  writes the processed observations (objective space) 
//  to disk
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the PE obs source
//  \end{description}
//EOP
{ 
  int len1;
  struct peobswritenode* current;
  struct peobswritenode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct peobswritenode*) malloc(sizeof(struct peobswritenode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(peobswrite_table == NULL){
    peobswrite_table = pnode;
  }
  else{
    current = peobswrite_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: writepeobs
// \label{writepeobs}
// 
// !INTERFACE:
void FTN(writepeobs)(char *j, void *peobs,int len)
//  
// !DESCRIPTION: 
// Invokes the routine from the registry to set up 
// structures for writing observation data for 
// optimization 
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the PE obs source
//  \item[peobs]
//   ESMF state that contain the observations
//  \end{description}
//
//EOP
{ 
  struct peobswritenode* current;
  
  current = peobswrite_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("write PE obs routine for  %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(peobs); 
}


//BOP
// !ROUTINE: registerinitobjfunctype
// \label{registerinitobjfunctype}
//  
// 
// !INTERFACE:
void FTN(registerinitobjfunctype)(char *j, void (*func)(),int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine that
//  sets up structures for initializing the method of 
//  objective function evaluation
// 
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the objective function method
//  \end{description}
//EOP
{ 
  int len1;
  struct objfuncinitnode* current;
  struct objfuncinitnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct objfuncinitnode*) malloc(sizeof(struct objfuncinitnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(objfuncinit_table == NULL){
    objfuncinit_table = pnode;
  }
  else{
    current = objfuncinit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: objectivefunctypeinit
// \label{objectivefunctypeinit}
// 
// !INTERFACE:
void FTN(objectivefunctypeinit)(char *j,int len)
//  
// !DESCRIPTION: 
// Invokes the routine from the registry to set up 
// structures for initializing the method of 
//  objective function evaluation
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the objective function method
//  \end{description}
//
//EOP
{ 
  struct objfuncinitnode* current;
  
  current = objfuncinit_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("init objective function routine for  %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}


//BOP
// !ROUTINE: registercomputeobjfunctype
// \label{registercomputeobjfunctype}
//  
// 
// !INTERFACE:
void FTN(registercomputeobjfunctype)(char *j, void (*func)(),int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine that
//  computes objective function based on the desired method
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the objective function method
//  \end{description}
//EOP
{ 
  int len1;
  struct objfunccomputenode* current;
  struct objfunccomputenode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct objfunccomputenode*) malloc(sizeof(struct objfunccomputenode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(objfunccompute_table == NULL){
    objfunccompute_table = pnode;
  }
  else{
    current = objfunccompute_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: computeobjectivefunctype
// \label{computeobjectivefunctype}
// 
// !INTERFACE:
void FTN(computeobjectivefunctype)(char *j,int len)
//  
// !DESCRIPTION: 
// Invokes the routine from the registry to invoke
// an objective function method
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the objective function method
//  \end{description}
//
//EOP
{ 
  struct objfunccomputenode* current;
  
  current = objfunccompute_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("compute objective function routine for  %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}


//BOP
// !ROUTINE: registerupdateobjfunctype
// \label{registerupdateobjfunctype}
//  
// 
// !INTERFACE:
void FTN(registerupdateobjfunctype)(char *j, void (*func)(),int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine that
//  updates objective function based on the desired method
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the objective function method
//  \end{description}
//EOP
{ 
  int len1;
  struct objfuncupdatenode* current;
  struct objfuncupdatenode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct objfuncupdatenode*) malloc(sizeof(struct objfuncupdatenode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(objfuncupdate_table == NULL){
    objfuncupdate_table = pnode;
  }
  else{
    current = objfuncupdate_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: updateobjectivefunctype
// \label{updateobjectivefunctype}
// 
// !INTERFACE:
void FTN(updateobjectivefunc)(char *j,int len)
//  
// !DESCRIPTION: 
// Invokes the routine from the registry to invoke
// an objective function method
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the objective function method
//  \end{description}
//
//EOP
{ 
  struct objfuncupdatenode* current;
  
  current = objfuncupdate_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("update objective function routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}


//BOP
// !ROUTINE: registerresetobjfunctype
// \label{registerresetobjfunctype}
//  
// 
// !INTERFACE:
void FTN(registerresetobjfunctype)(char *j, void (*func)(),int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine that
//  resets objective function 
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the objective function method
//  \end{description}
//EOP
{ 
  int len1;
  struct objfuncresetnode* current;
  struct objfuncresetnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct objfuncresetnode*) malloc(sizeof(struct objfuncresetnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(objfuncreset_table == NULL){
    objfuncreset_table = pnode;
  }
  else{
    current = objfuncreset_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: resetobjectivefunctype
// \label{resetobjectivefunctype}
// 
// !INTERFACE:
void FTN(resetobjectivefunctype)(char *j,int len)
//  
// !DESCRIPTION: 
// Invokes the routine from the registry to invoke
// an objective function method
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the objective function method
//  \end{description}
//
//EOP
{ 
  struct objfuncresetnode* current;
  
  current = objfuncreset_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("reset objective function routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}
//NOTE: ****This is deprecated - check and remove******
//BOP
// !ROUTINE: registeroptuetypeobjfunceval
// \label{registeroptuetypeobjfunceval}
//  
// 
// !INTERFACE:
void FTN(registeroptuetypeobjfunceval)(char *j, void (*func)(),int len)
//
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine that
//  evaluates objective function
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   index of the runmode (param estimation, for e.g)
//  \end{description}
//EOP
{ 
  int len1;
  struct objfuncevalnode* current;
  struct objfuncevalnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct objfuncevalnode*) malloc(sizeof(struct objfuncevalnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(objfunceval_table == NULL){
    objfunceval_table = pnode;
  }
  else{
    current = objfunceval_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: evaluateobjfunction
// \label{evaluateobjfunction}
// 
// !INTERFACE:
void FTN(evaluateobjfunction)(char *j,int len)
//  
// !DESCRIPTION: 
// Invokes the routine from the registry to invoke
// an objective function evaluation
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the objective function evaluation
//  \end{description}
//
//EOP
{ 
  struct objfuncevalnode* current;
  
  current = objfunceval_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("evaluate objective function routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}


