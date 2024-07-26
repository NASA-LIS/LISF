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
// !MODULE: LIS_glaciermodel_FTable
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
struct glaciermodelinitnode
{ 
  char *name;
  void (*func)(int*);

  struct glaciermodelinitnode* next;
} ;
struct glaciermodelinitnode* glaciermodelinit_table = NULL; 

struct glaciermodelrunnode
{ 
  char *name;
  void (*func)(int*,int*);

  struct glaciermodelrunnode* next;
} ;
struct glaciermodelrunnode* glaciermodelrun_table = NULL; 

struct glaciermodelfinalnode
{ 
  char *name;
  void (*func)();

  struct glaciermodelfinalnode* next;
} ;
struct glaciermodelfinalnode* glaciermodelfinal_table = NULL; 

struct glaciermodelsetupnode
{ 
  char *name;
  void (*func)(int*);

  struct glaciermodelsetupnode* next;
} ;
struct glaciermodelsetupnode* glaciermodelsetup_table = NULL;

struct glaciermodelrestartnode
{ 
  char *name;
  void (*func)(int*);

  struct glaciermodelrestartnode* next;
} ;
struct glaciermodelrestartnode* glaciermodelrestart_table = NULL;

struct glaciermodeloutputnode
{ 
  char *name;
  void (*func)(int*, int*);

  struct glaciermodeloutputnode* next;
} ;
struct glaciermodeloutputnode* glaciermodeloutput_table = NULL;

struct glaciermodeldynsetnode
{ 
  char *name;
  void (*func)(int*,int*);

  struct glaciermodeldynsetnode* next;
} ;
struct glaciermodeldynsetnode* glaciermodeldynset_table = NULL;

struct glaciermodelf2tnode
{ 
  char *name;
  void (*func)(int*,int*);

  struct glaciermodelf2tnode* next;
} ;
struct glaciermodelf2tnode* glaciermodelf2t_table = NULL;

struct glaciermodelwriterstnode
{ 
  char *name;
  void (*func)(int*);

  struct glaciermodelwriterstnode* next;
} ;
struct glaciermodelwriterstnode* glaciermodelwriterst_table = NULL;

//Routing
struct glacierroutinggetrunoffnode
{ 
  char *name;
  void (*func)(int*);

  struct glacierroutinggetrunoffnode* next;
} ;
struct glacierroutinggetrunoffnode* glacierroutinggetrunoff_table = NULL;


//BOP
// !ROUTINE: registerglaciermodelini
// \label{registerglaciermodelini}
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine to 
//  perform land surface model initialization
// 
// !INTERFACE:
void FTN(registerglaciermodelinit)(char *j, void (*func)(int*),int len)
//EOP
{ 
  int len1;
  struct glaciermodelinitnode* current;
  struct glaciermodelinitnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct glaciermodelinitnode*) malloc(sizeof(struct glaciermodelinitnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(glaciermodelinit_table == NULL){
    glaciermodelinit_table = pnode;
  }
  else{
    current = glaciermodelinit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: glaciermodelinit
// \label{glaciermodelinit}
//
// !INTERFACE:
void FTN(glaciermodelinit)(char *j,int *mtype, int len)
//
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to perform
//  land surface model initialization
// 
//  \begin{description}
//  \item[j]
//   name of the GLACIERMODEL
//  \end{description}
//EOP
{ 

  struct glaciermodelinitnode* current;
  
  current = glaciermodelinit_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("init routine for GLACIERMODEL %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(mtype); 
}
//BOP
// !ROUTINE: registerglaciermodelrun
// \label{registerglaciermodelrun}
//
// !INTERFACE:
void FTN(registerglaciermodelrun)(char *j, void (*func)(int*, int*),int len)
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine
//  to run the land surface model 
// 
//  \begin{description}
//  \item[j]
//   name of the GLACIERMODEL
//  \end{description}
//EOP
{ 
  int len1;
  struct glaciermodelrunnode* current;
  struct glaciermodelrunnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct glaciermodelrunnode*) malloc(sizeof(struct glaciermodelrunnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(glaciermodelrun_table == NULL){
    glaciermodelrun_table = pnode;
  }
  else{
    current = glaciermodelrun_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: glaciermodelrun
// \label{glaciermodelrun}
//
// !INTERFACE:
void FTN(glaciermodelrun)(char *j,int *n,int *mtype,int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to run the 
//  land surface model 
// 
//  \begin{description}
//  \item[j]
//   name of the GLACIERMODEL
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{
  struct glaciermodelrunnode* current;
  
  current = glaciermodelrun_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("run routine for GLACIERMODEL %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,mtype); 
}

//BOP
// !ROUTINE: registerglaciermodelfinalize
// \label{registerglaciermodelfinalize}
//
// !INTERFACE:
void FTN(registerglaciermodelfinalize)(char *j, void (*func)(),int len)
//  
// !DESCRIPTION:
//  Creates an entry in the registry for the routine
//  to cleanup allocated structures specific to the 
//  land surface model
// 
//  \begin{description}
//  \item[j]
//   name of the GLACIERMODEL
//  \end{description}
// 
//EOP
{ 
  int len1;
  struct glaciermodelfinalnode* current;
  struct glaciermodelfinalnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct glaciermodelfinalnode*) malloc(sizeof(struct glaciermodelfinalnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(glaciermodelfinal_table == NULL){
    glaciermodelfinal_table = pnode;
  }
  else{
    current = glaciermodelfinal_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: glaciermodelfinalize
// \label{glaciermodelfinalize}
//
// !INTERFACE:
void FTN(glaciermodelfinalize)(char *j,int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry for cleaning up
//  allocated structures specific to the land surface model
// 
//  \begin{description}
//  \item[j]
//   name of the GLACIERMODEL
//  \end{description}
// 
//EOP
{
  struct glaciermodelfinalnode* current;
  
  current = glaciermodelfinal_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("finalize routine for GLACIERMODEL %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}

//BOP
// !ROUTINE: registerglaciermodelsetup
// \label{registerglaciermodelsetup}
//
// !INTERFACE:
void FTN(registerglaciermodelsetup)(char *j, void (*func)(int*),int len)
//  
// !DESCRIPTION:
//  Makes an entry in the registry for the routine
//  to set up land surface model parameters 
// 
//  \begin{description}
//  \item[j]
//   name of the GLACIERMODEL
//  \end{description}
//EOP
{ 
  int len1;
  struct glaciermodelsetupnode* current;
  struct glaciermodelsetupnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct glaciermodelsetupnode*) malloc(sizeof(struct glaciermodelsetupnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(glaciermodelsetup_table == NULL){
    glaciermodelsetup_table = pnode;
  }
  else{
    current = glaciermodelsetup_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: glaciermodelsetup
// \label{glaciermodelsetup}
//
// !INTERFACE:
void FTN(glaciermodelsetup)(char *j,int *mtype, int len)
//  
// !DESCRIPTION:  
//  Invokes the routine in the registry to set up 
//  land surface model parameters  
// 
//  \begin{description}
//  \item[j]
//   name of the GLACIERMODEL
//  \end{description}
//EOP
{ 
  struct glaciermodelsetupnode* current;
  
  current = glaciermodelsetup_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("setup routine for GLACIERMODEL %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(mtype); 
}

//BOP
// !ROUTINE: registerglaciermodelrestart
// \label{registerglaciermodelrestart}
// 
// !INTERFACE:
void FTN(registerglaciermodelrestart)(char *j, void (*func)(int*),int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// restart the land surface model from a 
// previously saved state
// 
//  \begin{description}
//  \item[j]
//   name of the GLACIERMODEL
//  \end{description}
//EOP
{ 
  int len1;
  struct glaciermodelrestartnode* current;
  struct glaciermodelrestartnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct glaciermodelrestartnode*) malloc(sizeof(struct glaciermodelrestartnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(glaciermodelrestart_table == NULL){
    glaciermodelrestart_table = pnode;
  }
  else{
    current = glaciermodelrestart_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: glaciermodelrestart
// \label{glaciermodelrestart}
//
// !INTERFACE:
void FTN(glaciermodelrestart)(char *j,int *mtype, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to 
//  restart the land surface model from a previously
//  saved state
// 
//  \begin{description}
//  \item[j]
//   name of the GLACIERMODEL
//  \end{description}
//EOP
{ 
  struct glaciermodelrestartnode* current;
  
  current = glaciermodelrestart_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("read restart routine for GLACIERMODEL %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(mtype); 
}

//BOP
// !ROUTINE: registerglaciermodeldynsetup
// \label{registerglaciermodeldynsetup}
// 
// !INTERFACE:
void FTN(registerglaciermodeldynsetup)(char *j, void (*func)(int*, int*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  set the time dependent land surface parameters
// 
//  \begin{description}
//  \item[j]
//   name of the GLACIERMODEL
//  \end{description}
//EOP
{ 
  int len1;
  struct glaciermodeldynsetnode* current;
  struct glaciermodeldynsetnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct glaciermodeldynsetnode*) malloc(sizeof(struct glaciermodeldynsetnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(glaciermodeldynset_table == NULL){
    glaciermodeldynset_table = pnode;
  }
  else{
    current = glaciermodeldynset_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: glaciermodeldynsetup
// \label{glaciermodeldynsetup}
// 
// !INTERFACE:
void FTN(glaciermodeldynsetup)(char *j, int *n, int *mtype,int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to set the 
//  time dependent land surface parameters
// 
//  \begin{description}
//  \item[j]
//   name of the GLACIERMODEL
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{ 
  struct glaciermodeldynsetnode* current;
  
  current = glaciermodeldynset_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("dynamic setup routine for GLACIERMODEL %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,mtype); 
}

//BOP
// !ROUTINE: registerglaciermodeloutput
// \label{registerglaciermodeloutput}
// 
// !INTERFACE:
void FTN(registerglaciermodeloutput)(char *j, void (*func)(int*, int*),int len)
//  
// !DESCRIPTION:
//  Makes an entry in the registry for the routine to 
//  perform land surface model output
// 
//  \begin{description}
//  \item[j]
//   name of the GLACIERMODEL
//  \end{description}
//EOP
{ 
  int len1;
  struct glaciermodeloutputnode* current;
  struct glaciermodeloutputnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct glaciermodeloutputnode*) malloc(sizeof(struct glaciermodeloutputnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(glaciermodeloutput_table == NULL){
    glaciermodeloutput_table = pnode;
  }
  else{
    current = glaciermodeloutput_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: glaciermodeloutput
// \label{glaciermodeloutput}
//
// !INTERFACE:
void FTN(glaciermodeloutput)(char *j, int *n, int *mtype, int len)
//  
// !DESCRIPTION:
//  Invokes the routine from the registry to perform
//  land surface model output
// 
//  \begin{description}
//  \item[j]
//   name of the GLACIERMODEL
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{ 

  struct glaciermodeloutputnode* current;
  
  current = glaciermodeloutput_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("output writing routine for GLACIERMODEL %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n, mtype); 
}

//BOP
// !ROUTINE: registerglaciermodelf2t
// \label{registerglaciermodelf2t}
// 
// !INTERFACE:
void FTN(registerglaciermodelf2t)(char *j, void (*func)(int*, int*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  transfer forcing to model tiles
// 
//  \begin{description}
//  \item[j]
//   name of the GLACIERMODEL + name of the running mode
//  \item[j]
//   index of the runmode
//  \end{description}
//EOP
{ 
  int len1;
  struct glaciermodelf2tnode* current;
  struct glaciermodelf2tnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct glaciermodelf2tnode*) malloc(sizeof(struct glaciermodelf2tnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(glaciermodelf2t_table == NULL){
    glaciermodelf2t_table = pnode;
  }
  else{
    current = glaciermodelf2t_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: glaciermodelf2t
// \label{glaciermodelf2t}
// 
// !INTERFACE:
void FTN(glaciermodelf2t)(char *j, int *n, int *mtype, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to 
//  transfer forcing to model tiles
// 
//  \begin{description}
//  \item[j]
//   name of the GLACIERMODEL
//  \item[j]
//   index of the runmode
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{ 
  struct glaciermodelf2tnode* current;
  
  current = glaciermodelf2t_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("f2t writing routine for GLACIERMODEL and running mode %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,mtype); 
}

//BOP
// !ROUTINE: registerglaciermodelwrst
// \label{registerglaciermodelwrst}
// 
// !INTERFACE:
void FTN(registerglaciermodelwrst)(char *j, void (*func)(int*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine 
//  to write restart files for a land surface model 
// 
//  \begin{description}
//  \item[j]
//   name of the GLACIERMODEL
//  \end{description}
//EOP
{ 
  int len1;
  struct glaciermodelwriterstnode* current;
  struct glaciermodelwriterstnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct glaciermodelwriterstnode*) malloc(sizeof(struct glaciermodelwriterstnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(glaciermodelwriterst_table == NULL){
    glaciermodelwriterst_table = pnode;
  }
  else{
    current = glaciermodelwriterst_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: glaciermodelwrst
// \label{glaciermodelwrst}
// 
// !INTERFACE:
void FTN(glaciermodelwrst)(char *j, int *n, int len)
//  
// !DESCRIPTION:
//  Invokes the routine from the registry to write
//  restart files for a land surface model
// 
//  \begin{description}
//  \item[j]
//   name of the GLACIERMODEL
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{ 

  struct glaciermodelwriterstnode* current;
  
  current = glaciermodelwriterst_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("write restart writing routine for GLACIERMODEL %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}

//BOP
// !ROUTINE: registerglacierroutinggetrunoff
// \label{registerglacierroutinggetrunoff}
// 
// !INTERFACE:
void FTN(registerglacierroutinggetrunoff)(char *j, void (*func)(int*),int len)
//  
// !DESCRIPTION: 
//  creates an entry in the registry for the routine to 
//  return runoff fields from a GLACIER.   
//
//  \begin{description}
//  \item[j]
//   name of the GLACIER + routing instance
//  \end{description}
//EOP
{ 
  int len1;
  struct glacierroutinggetrunoffnode* current;
  struct glacierroutinggetrunoffnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct glacierroutinggetrunoffnode*) malloc(sizeof(struct glacierroutinggetrunoffnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(glacierroutinggetrunoff_table == NULL){
    glacierroutinggetrunoff_table = pnode;
  }
  else{
    current = glacierroutinggetrunoff_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}

//BOP
// !ROUTINE: glacierroutinggetrunoff
// \label{glacierroutinggetrunoff}
// 
// !INTERFACE:
void FTN(glacierroutinggetrunoff)(char *j, int *n, int len)
//  
// !DESCRIPTION: 
//  Invokes the registered routine that returns the
//  runoff fields from a GLACIER. 
//
//  \begin{description}
//  \item[j]
//   name of the GLACIER + routing instance
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{ 
  struct glacierroutinggetrunoffnode* current;
  
  current = glacierroutinggetrunoff_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("get runoff routine for GLACIER + routing instance %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}





