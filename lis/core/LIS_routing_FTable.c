//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center Land Information System (LIS) v7.2
//
// Copyright (c) 2015 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
//BOP
//
// !MODULE: LIS_routing_FTable
//  
//
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing the operations of different 
//  routing models. 
//
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"
struct routinginitnode
{ 
  char *name;
  void (*func)();

  struct routinginitnode* next;
} ;
struct routinginitnode* routinginit_table = NULL; 

struct routingrrstnode
{ 
  char *name;
  void (*func)();

  struct routingrrstnode* next;
} ;
struct routingrrstnode* routingrrst_table = NULL; 

struct routingwrstnode
{ 
  char *name;
  void (*func)(int*);

  struct routingwrstnode* next;
} ;
struct routingwrstnode* routingwrst_table = NULL; 

struct routingrunnode
{ 
  char *name;
  void (*func)(int*);

  struct routingrunnode* next;
} ;
struct routingrunnode* routingrun_table = NULL; 

struct routingoutnode
{ 
  char *name;
  void (*func)(int*);

  struct routingoutnode* next;
} ;
struct routingoutnode* routingout_table = NULL; 

//for DA
struct routingdainitnode
{ 
  char *name;
  void (*func)(int*);

  struct routingdainitnode* next;
} ;
struct routingdainitnode* routingdainit_table = NULL;

struct routingdagetvarnode
{ 
  char *name;
  void (*func)(int*, void*);

  struct routingdagetvarnode* next;
} ;
struct routingdagetvarnode* routingdagetvar_table = NULL;

struct routingdasetvarnode
{ 
  char *name;
  void (*func)(int*, void*);

  struct routingdasetvarnode* next;
} ;
struct routingdasetvarnode* routingdasetvar_table = NULL;

struct routingdaobstransformnode
{ 
  char *name;
  void (*func)(int*, void*);

  struct routingdaobstransformnode* next;
} ;
struct routingdaobstransformnode* routingdaobstransform_table = NULL;

struct routingdamapobstoroutingnode
{ 
  char *name;
  void (*func)(int*, int*, void*,void*);

  struct routingdamapobstoroutingnode* next;
} ;
struct routingdamapobstoroutingnode* routingdamapobstorouting_table = NULL;

struct routingdaobsprednode
{ 
  char *name;
  void (*func)(int*, int*, float*);

  struct routingdaobsprednode* next;
} ;
struct routingdaobsprednode* routingdaobspred_table = NULL;

struct routingdaqcstatenode
{ 
  char *name;
  void (*func)(int*, void*);

  struct routingdaqcstatenode* next;
} ;
struct routingdaqcstatenode* routingdaqcstate_table = NULL;

struct routingdaqcobsnode
{ 
  char *name;
  void (*func)(int*, int*, void*);

  struct routingdaqcobsnode* next;
} ;
struct routingdaqcobsnode* routingdaqcobs_table = NULL;

struct routingdascalenode
{ 
  char *name;
  void (*func)(int*, void*);

  struct routingdascalenode* next;
} ;
struct routingdascalenode* routingdascale_table = NULL;

struct routingdadescalenode
{ 
  char *name;
  void (*func)(int*, void*, void*);

  struct routingdadescalenode* next;
} ;
struct routingdadescalenode* routingdadescale_table = NULL;

struct routingdaupdatenode
{ 
  char *name;
  void (*func)(int*, void*,void*);

  struct routingdaupdatenode* next;
} ;
struct routingdaupdatenode* routingdaupdate_table = NULL;

struct routingdawritevarnode
{ 
  char *name;
  void (*func)(int*, int*,void*);

  struct routingdawritevarnode* next;
} ;
struct routingdawritevarnode* routingdawritevar_table = NULL;

struct routingdiagfordanode
{ 
  char *name;
  void (*func)(int*);

  struct routingdiagfordanode* next;
} ;
struct routingdiagfordanode* routingdiagforda_table = NULL;

struct routingdastatesizenode
{ 
  char *name;
  void (*func)(int*,int*);

  struct routingdastatesizenode* next;
} ;
struct routingdastatesizenode* routingdastatesize_table = NULL;


struct routingdasetpertnode
{ 
  char *name;
  void (*func)(int*,int*, void*,void*);

  struct routingdasetpertnode* next;
} ;
struct routingdasetpertnode* routingdasetpert_table = NULL;


//BOP
// !ROUTINE: registerroutinginit
// \label{registerroutinginit}
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine to 
//  perform routing model initialization
// 
// !INTERFACE:
void FTN(registerroutinginit)(char *j, void (*func)(),int len)
//EOP
{ 
  int len1;
  struct routinginitnode* current;
  struct routinginitnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct routinginitnode*) malloc(sizeof(struct routinginitnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(routinginit_table == NULL){
    routinginit_table = pnode;
  }
  else{
    current = routinginit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: routinginit
// \label{routinginit}
//
// !INTERFACE:
void FTN(routinginit)(char *j,int len)
//
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to perform
//  routing model initialization
// 
//  \begin{description}
//  \item[i]
//   index of the routing model
//  \end{description}
//EOP
{ 
  struct routinginitnode* current;
  
  current = routinginit_table;
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
// !ROUTINE: registerroutingreadrestart
// \label{registerroutingreadrestart}
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine to 
//  perform routing model read restart. 
// 
// !INTERFACE:
void FTN(registerroutingreadrestart)(char *j, void (*func)(),int len)
//EOP
{ 
  int len1;
  struct routingrrstnode* current;
  struct routingrrstnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct routingrrstnode*) malloc(sizeof(struct routingrrstnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingrrst_table == NULL){
    routingrrst_table = pnode;
  }
  else{
    current = routingrrst_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: routingreadrestart
// \label{routingreadrestart}
//
// !INTERFACE:
void FTN(routingreadrestart)(char *j,int len)
//
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to perform
//  routing model initialization
// 
//  \begin{description}
//  \item[i]
//   index of the routing model
//  \end{description}
//EOP
{ 
  struct routingrrstnode* current;
  
  current = routingrrst_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("read restart routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}

//BOP
// !ROUTINE: registerroutingrun
// \label{registerroutingrun}
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine to 
//  perform routing model run
// 
// !INTERFACE:
void FTN(registerroutingrun)(char *j, void (*func)(int*),int len)
//EOP
{ 
  int len1;
  struct routingrunnode* current;
  struct routingrunnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct routingrunnode*) malloc(sizeof(struct routingrunnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingrun_table == NULL){
    routingrun_table = pnode;
  }
  else{
    current = routingrun_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: routingrun
// \label{routingrun}
//
// !INTERFACE:
void FTN(routingrun)(char *j, int *n,int len)
//
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to perform
//  routing model initialization
// 
//  \begin{description}
//  \item[i]
//   index of the routing model
//  \end{description}
//EOP
{ 
  struct routingrunnode* current;
  
  current = routingrun_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("run routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}

//BOP
// !ROUTINE: registerroutingoutput
// \label{registerroutingoutput}
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine to 
//  perform routing model output
// 
// !INTERFACE:
void FTN(registerroutingoutput)(char *j, void (*func)(int*),int len)
//EOP
{ 
  int len1;
  struct routingoutnode* current;
  struct routingoutnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct routingoutnode*) malloc(sizeof(struct routingoutnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingout_table == NULL){
    routingout_table = pnode;
  }
  else{
    current = routingout_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: routingoutput
// \label{routingoutput}
//
// !INTERFACE:
void FTN(routingoutput)(char *j, int *n,int len)
//
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to perform
//  routing model initialization
// 
//  \begin{description}
//  \item[i]
//   index of the routing model
//  \end{description}
//EOP
{ 
  struct routingoutnode* current;
  
  current = routingout_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("output routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}

//BOP
// !ROUTINE: registerroutingwriterestart
// \label{registerroutingwriterestart}
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine to 
//  perform routing model writerestart
// 
// !INTERFACE:
void FTN(registerroutingwriterestart)(char *j, void (*func)(int*),int len)
//EOP
{ 
  int len1;
  struct routingwrstnode* current;
  struct routingwrstnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct routingwrstnode*) malloc(sizeof(struct routingwrstnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingwrst_table == NULL){
    routingwrst_table = pnode;
  }
  else{
    current = routingwrst_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: routingwriterestart
// \label{routingwriterestart}
//
// !INTERFACE:
void FTN(routingwriterestart)(char *j, int *n,int len)
//
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to perform
//  routing model initialization
// 
//  \begin{description}
//  \item[i]
//   index of the routing model
//  \end{description}
//EOP
{ 
  struct routingwrstnode* current;
  
  current = routingwrst_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("write restart routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}

//BOP
// !ROUTINE: registerroutingdainit
// \label{registerroutingdainit}
// 
// !INTERFACE:
void FTN(registerroutingdainit)(char *j, void (*func)(int*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine 
//  for initializing DA related ROUTING settings
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING + DA instance
//  \end{description}
//EOP
{ 
  struct routingdainitnode* current;
  struct routingdainitnode* pnode; 
  // create node
  
  pnode=(struct routingdainitnode*) malloc(sizeof(struct routingdainitnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingdainit_table == NULL){
    routingdainit_table = pnode;
  }
  else{
    current = routingdainit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: routingdainit
// \label{routingdainit}
// 
// !INTERFACE:
void FTN(routingdainit)(char *j, int *k, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry for initializing
//  DA related ROUTING settings. 
//  \begin{description}
//  \item[j]
//   name of the ROUTING + DA instance
//  \end{description}
//EOP
{ 
  struct routingdainitnode* current;
  
  current = routingdainit_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("init routine for ROUTING + DA instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(k); 
}

//BOP
// !ROUTINE: registerroutingdagetstatevar
// \label{registerroutingdagetstatevar}
// 
// !INTERFACE:
void FTN(registerroutingdagetstatevar)(char *j, void (*func)(int*, void*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine 
//  for obtaining the specified prognostic variables from the 
//  land surface model (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING + DA instance
//  \end{description}
//EOP
{ 
  struct routingdagetvarnode* current;
  struct routingdagetvarnode* pnode; 
  // create node
  
  pnode=(struct routingdagetvarnode*) malloc(sizeof(struct routingdagetvarnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingdagetvar_table == NULL){
    routingdagetvar_table = pnode;
  }
  else{
    current = routingdagetvar_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: routingdagetstatevar
// \label{routingdagetstatevar}
// 
// !INTERFACE:
void FTN(routingdagetstatevar)(char *j, int *n, void *state, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry for obtaining
//  the specified prognostic variables from the land surface model
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING + DA instance
//  \item[n]
//   index of the nest
//  \item[state]
//   pointer to the prognostic variable state
//  \end{description}
//EOP
{ 
  struct routingdagetvarnode* current;
  
  current = routingdagetvar_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("get ROUTING variable routine for ROUTING + DA instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,state); 
}

//BOP
// !ROUTINE: registerroutingdasetstatevar
// \label{registerroutingdasetstatevar}
// 
// !INTERFACE:
void FTN(registerroutingdasetstatevar)(char *j, void (*func)(int*, void*),int len)
//  
// !DESCRIPTION:
//  Makes an entry in the registry for updating the specified
//  state variable in a land surface model 
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING + DA instance
//  \end{description}
//EOP
{ 
  struct routingdasetvarnode* current;
  struct routingdasetvarnode* pnode; 
  // create node
  
  pnode=(struct routingdasetvarnode*) malloc(sizeof(struct routingdasetvarnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingdasetvar_table == NULL){
    routingdasetvar_table = pnode;
  }
  else{
    current = routingdasetvar_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: routingdasetstatevar
// \label{routingdasetstatevar}
// 
// !INTERFACE:
void FTN(routingdasetstatevar)(char *j,int *n, void *statevar, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry for updating
//  the specified state variable in a land surface model 
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING + DA instance
//  \item[n]
//   index of the nest
//  \item[statevars]
//   pointer to the prognostic variable state
//  \end{description}
//EOP
{ 
  struct routingdasetvarnode* current;
  
  current = routingdasetvar_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;

    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("set ROUTING variable routine for ROUTING + DA instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,statevar); 
}


//BOP
// !ROUTINE: registerroutingdaobstransform
// \label{registerroutingdaobstransform}
// 
// !INTERFACE:
void FTN(registerroutingdaobstransform)(char *j, void (*func)(int*, void*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry to perform the  
//  translation of observations to state variables 
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING + DA instance
//  \item[k]
//   index of the observation data
//  \end{description}
//EOP
{ 
  struct routingdaobstransformnode* current;
  struct routingdaobstransformnode* pnode; 
  // create node
  
  pnode=(struct routingdaobstransformnode*) malloc(sizeof(struct routingdaobstransformnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingdaobstransform_table == NULL){
    routingdaobstransform_table = pnode;
  }
  else{
    current = routingdaobstransform_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}

//BOP
// !ROUTINE: routingdaobstransform
// \label{routingdaobstransform}
//
// !INTERFACE:
void FTN(routingdaobstransform)(char *j, int *n, void *obs, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to translate
//  the observations to the prognostic variable space
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING + DA instance
//  \item[n]
//   index of the nest
//   \item[obs]
//   transformed variable
//  \end{description}
//EOP
{ 
  struct routingdaobstransformnode* current;
  
  current = routingdaobstransform_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("obs transform routine for ROUTING + DA instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,obs); 
}

//BOP
// !ROUTINE: registerroutingdagetobspred
// \label{registerroutingdagetobspred}
// 
// !INTERFACE:
void FTN(registerroutingdagetobspred)(char *j, void (*func)(int*,int*,float*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine 
//  that provides an ROUTING's estimate of the observations.
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING + DA instance
//  \end{description}
//EOP
{ 
  struct routingdaobsprednode* current;
  struct routingdaobsprednode* pnode; 
  // create node
  
  pnode=(struct routingdaobsprednode*) malloc(sizeof(struct routingdaobsprednode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingdaobspred_table == NULL){
    routingdaobspred_table = pnode;
  }
  else{
    current = routingdaobspred_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}

//BOP
// !ROUTINE: routingdagetobspred
// \label{routingdagetobspred}
//
// !INTERFACE:
void FTN(routingdagetobspred)(char *j, int *n, int *k,float *pred, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to translate
//  the observations to state variables
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING
//  \item[n]
//   index of the nest
//  \item[pred]
//   model's estimated observation prediction
//  \end{description}
//EOP
{ 
  struct routingdaobsprednode* current;
  
  current = routingdaobspred_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;

    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("obspred routine for ROUTING + DA instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,k,pred); 
}

//BOP
// !ROUTINE: registerroutingdadiagnosevars
// \label{registerroutingdadiagnosevars}
// 
// !INTERFACE:
void FTN(registerroutingdadiagnosevars)(char *j, void (*func)(int*),int len)
//  
// !DESCRIPTION:
//  Makes an entry in the registry for the routine to 
//  perform land surface model output
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING
//  \end{description}
//EOP
{ 

  struct routingdiagfordanode* current;
  struct routingdiagfordanode* pnode; 
  // create node
  
  pnode=(struct routingdiagfordanode*) malloc(sizeof(struct routingdiagfordanode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingdiagforda_table == NULL){
    routingdiagforda_table = pnode;
  }
  else{
    current = routingdiagforda_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: routingdadiagnosevars
// \label{routingdadiagnosevars}
//
// !INTERFACE:
void FTN(routingdadiagnosevars)(char *j, int *n, int len)
//  
// !DESCRIPTION:
//  Invokes the routine from the registry to perform
//  land surface model output
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{ 

  struct routingdiagfordanode* current;
  
  current = routingdiagforda_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("diagnose for DA routine for ROUTING +DAset %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}


//BOP
// !ROUTINE: registerroutingdamapobstorouting
// \label{registerroutingdamapobstorouting}
// 
// !INTERFACE:
void FTN(registerroutingdamapobstorouting)(char *j, void (*func)(int*, int*, void*, void*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry to perform the  
//  translation of observations to state variables 
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING + DA instance
//  \end{description}
//EOP
{ 
  struct routingdamapobstoroutingnode* current;
  struct routingdamapobstoroutingnode* pnode; 
  // create node
  
  pnode=(struct routingdamapobstoroutingnode*) malloc(sizeof(struct routingdamapobstoroutingnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingdamapobstorouting_table == NULL){
    routingdamapobstorouting_table = pnode;
  }
  else{
    current = routingdamapobstorouting_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}

//BOP
// !ROUTINE: routingdamapobstorouting
// \label{routingdamapobstorouting}
//
// !INTERFACE:
void FTN(routingdamapobstorouting)(char *j, int *n, int *k, void *obs, void *routing, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to translate
//  the observations to state variables
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING + DA instance
//  \item[n]
//   index of the nest
//  \item[obs]
//   observations to be mapped
//  \item[routing]
//   updated routing states  
//  \end{description}
//EOP
{ 
  struct routingdamapobstoroutingnode* current;
  
  current = routingdamapobstorouting_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("map obs to ROUTING routine for ROUTING + DA instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,k,obs,routing); 
}

//BOP
// !ROUTINE: registerroutingdaqcstate
// \label{registerroutingdaqcstate}
// 
// !INTERFACE:
void FTN(registerroutingdaqcstate)(char *j, void (*func)(int*, void*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  QC the updated ROUTING state
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING + DA instance
//  \end{description}
//EOP
{ 
  struct routingdaqcstatenode* current;
  struct routingdaqcstatenode* pnode; 
  // create node
  
  pnode=(struct routingdaqcstatenode*) malloc(sizeof(struct routingdaqcstatenode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingdaqcstate_table == NULL){
    routingdaqcstate_table = pnode;
  }
  else{
    current = routingdaqcstate_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: routingdaqcstate
// \label{routingdaqcstate}
// 
// !INTERFACE:
void FTN(routingdaqcstate)(char *j,int *n, void *ROUTING_State, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to set the 
//  QC the updated ROUTING variables
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING
//  \item[n]
//   index of the nest
//  \item[ROUTING\_State]
//   The ROUTING state being qc'd
//  \end{description}
//EOP
{ 
  struct routingdaqcstatenode* current;
  
  current = routingdaqcstate_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("map obs to ROUTING routine for ROUTING + DA instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,ROUTING_State); 
}

//BOP
// !ROUTINE: registerroutingdascalestatevar
// \label{registerroutingdascalestatevar}
// 
// !INTERFACE:
void FTN(registerroutingdascalestatevar)(char *j, void (*func)(int*, void*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  scale the ROUTING state variables
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING + DA instance
//  \item[j]
//   index of the assimilated variable
//  \end{description}
//EOP
{ 
  struct routingdascalenode* current;
  struct routingdascalenode* pnode; 
  // create node
  
  pnode=(struct routingdascalenode*) malloc(sizeof(struct routingdascalenode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingdascale_table == NULL){
    routingdascale_table = pnode;
  }
  else{
    current = routingdascale_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: routingdascalestatevar
// \label{routingdascalestatevar}
// 
// !INTERFACE:
void FTN(routingdascalestatevar)(char *j, int *n, void *ROUTING_State, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to scale the ROUTING 
//  state variables
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING + DA instance
//  \item[n]
//   index of the nest
//  \item[ROUTING\_State]
//   The ROUTING state being scaled
//  \end{description}
//EOP
{ 
  struct routingdascalenode* current;
  
  current = routingdascale_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("scale variable routine for ROUTING + DA instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,ROUTING_State); 
}


//BOP
// !ROUTINE: registerroutingdadescalestatevar
// \label{registerroutingdadescalestatevar}
// 
// !INTERFACE:
void FTN(registerroutingdadescalestatevar)(char *j, void (*func)(int*, void*, void*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  descale the ROUTING state variables
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING
//  \item[j]
//   index of the assimilated variable
//  \end{description}
//EOP
{ 
  struct routingdadescalenode* current;
  struct routingdadescalenode* pnode; 
  // create node
  
  pnode=(struct routingdadescalenode*) malloc(sizeof(struct routingdadescalenode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingdadescale_table == NULL){
    routingdadescale_table = pnode;
  }
  else{
    current = routingdadescale_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: routingdadescalestatevar
// \label{routingdadescalestatevar}
// 
// !INTERFACE:
void FTN(routingdadescalestatevar)(char *j,int *n, void *ROUTING_State, void *ROUTING_Incr_State, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry descale the
//  ROUTING state variables 
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING
//  \item[n]
//   index of the nest
//  \item[ROUTING\_State]
//   The ROUTING state being descaled
//  \end{description}
//EOP
{ 
  struct routingdadescalenode* current;
  
  current = routingdadescale_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("descale variables routine for ROUTING + DA instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,ROUTING_State,ROUTING_Incr_State); 
}

//BOP
// !ROUTINE: registerroutingdaupdatestate
// \label{registerroutingdaupdatestate}
// 
// !INTERFACE:
void FTN(registerroutingdaupdatestate)(char *j, void (*func)(int*, void*, void*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  apply the ROUTING state increments to ROUTING state
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING + DA instance
//  \end{description}
//EOP
{ 
  struct routingdaupdatenode* current;
  struct routingdaupdatenode* pnode; 
  // create node
  
  pnode=(struct routingdaupdatenode*) malloc(sizeof(struct routingdaupdatenode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingdaupdate_table == NULL){
    routingdaupdate_table = pnode;
  }
  else{
    current = routingdaupdate_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: routingdaupdatestate
// \label{routingdaupdatestate}
// 
// !INTERFACE:
void FTN(routingdaupdatestate)(char *j, int *n, void *ROUTING_State, void *ROUTING_Incr_State, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to apply the
//  ROUTING state increments to ROUTING State
//  (for data assimilation)
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING + DA instance
//  \item[n]
//   index of the nest
//  \item[ROUTING\_State]
//   The ROUTING state being updated
//  \item[ROUTING\_Incr\_State]
//   The ROUTING incr state
//  \end{description}
//EOP
{ 
  struct routingdaupdatenode* current;
  
  current = routingdaupdate_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("update state routine for ROUTING + DA instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,ROUTING_State, ROUTING_Incr_State);
}

//BOP
// !ROUTINE: registerroutingdaqcobsstate
// \label{registerroutingdaqcobsstate}
// 
// !INTERFACE:
void FTN(registerroutingdaqcobsstate)(char *j, void (*func)(int*, int*, void*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  QC the OBS state based on ROUTING variables and states
//  (for data assimilation).
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING + DA instance
//  \end{description}
//EOP
{ 
  struct routingdaqcobsnode* current;
  struct routingdaqcobsnode* pnode; 
  // create node
  
  pnode=(struct routingdaqcobsnode*) malloc(sizeof(struct routingdaqcobsnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingdaqcobs_table == NULL){
    routingdaqcobs_table = pnode;
  }
  else{
    current = routingdaqcobs_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: routingdaqcobsstate
// \label{routingdaqcobsstate}
// 
// !INTERFACE:
void FTN(routingdaqcobsstate)(char *j, int *n, int *k, void *ROUTING_State, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to set the 
//  QC the observation state based on ROUTING variables and states
//  (for data assimilation).
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING
//  \item[n]
//   index of the nest
//  \item[ROUTING\_State]
//   The ROUTING state being qc'd
//  \end{description}
//EOP
{ 
  struct routingdaqcobsnode* current;
  
  current = routingdaqcobs_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("qc obs state routine for ROUTING + DA instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,k,ROUTING_State); 
}

//BOP
// !ROUTINE: registerroutingdagetstatespacesize
// \label{registerroutingdagetstatespacesize}
// 
// !INTERFACE:
void FTN(registerroutingdagetstatespacesize)(char *j, void (*func)(int*, int*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  QC the OBS state based on ROUTING variables and states
//  (for data assimilation).
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING + DA instance
//  \end{description}
//EOP
{ 

  struct routingdastatesizenode* current;
  struct routingdastatesizenode* pnode; 
  // create node
  
  pnode=(struct routingdastatesizenode*) malloc(sizeof(struct routingdastatesizenode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingdastatesize_table == NULL){
    routingdastatesize_table = pnode;
  }
  else{
    current = routingdastatesize_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: routingdagetstatespacesize
// \label{routingdagetstatespacesize}
// 
// !INTERFACE:
void FTN(routingdagetstatespacesize)(char *j, int *n, int *size, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to set the 
//  QC the observation state based on ROUTING variables and states
//  (for data assimilation).
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING
//  \item[n]
//   index of the nest
//  \item[ROUTING\_State]
//   The ROUTING state being qc'd
//  \end{description}
//EOP
{ 
  struct routingdastatesizenode* current;
  
  current = routingdastatesize_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("getstatespacesize for ROUTING + DA instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,size); 
}


//BOP
// !ROUTINE: registerroutingdasetpertstates
// \label{registerroutingdasetpertstates}
// 
// !INTERFACE:
void FTN(registerroutingdasetpertstates)(char *j, void (*func)(int*, int*, void*, void*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  QC the OBS state based on ROUTING variables and states
//  (for data assimilation).
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING + DA instance
//  \end{description}
//EOP
{ 

  struct routingdasetpertnode* current;
  struct routingdasetpertnode* pnode; 
  // create node
  
  pnode=(struct routingdasetpertnode*) malloc(sizeof(struct routingdasetpertnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingdasetpert_table == NULL){
    routingdasetpert_table = pnode;
  }
  else{
    current = routingdasetpert_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: routingdasetpertstates
// \label{routingdasetpertstates}
// 
// !INTERFACE:
void FTN(routingdasetpertstates)(char *j, int *n, int *Nstate, void *pstate, void *progpert, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to set the 
//  QC the observation state based on ROUTING variables and states
//  (for data assimilation).
// 
//  \begin{description}
//  \item[j]
//   name of the ROUTING
//  \item[n]
//   index of the nest
//  \item[ROUTING\_State]
//   The ROUTING state being qc'd
//  \end{description}
//EOP
{ 
  struct routingdasetpertnode* current;
  
  current = routingdasetpert_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("setpertstates for ROUTING + DA instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,Nstate, pstate, progpert); 
}


