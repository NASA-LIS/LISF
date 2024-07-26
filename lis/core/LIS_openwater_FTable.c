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
// !MODULE: LIS_openwater_FTable
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
struct openwaterinitnode
{ 
  char *name;
  void (*func)(int*);

  struct openwaterinitnode* next;
} ;
struct openwaterinitnode* openwaterinit_table = NULL; 

struct openwaterrunnode
{ 
  char *name;
  void (*func)(int*,int*);

  struct openwaterrunnode* next;
} ;
struct openwaterrunnode* openwaterrun_table = NULL; 

struct openwaterfinalnode
{ 
  char *name;
  void (*func)();

  struct openwaterfinalnode* next;
} ;
struct openwaterfinalnode* openwaterfinal_table = NULL; 

struct openwatersetupnode
{ 
  char *name;
  void (*func)(int*);

  struct openwatersetupnode* next;
} ;
struct openwatersetupnode* openwatersetup_table = NULL;

struct openwaterrestartnode
{ 
  char *name;
  void (*func)(int*);

  struct openwaterrestartnode* next;
} ;
struct openwaterrestartnode* openwaterrestart_table = NULL;

struct openwateroutputnode
{ 
  char *name;
  void (*func)(int*, int*);

  struct openwateroutputnode* next;
} ;
struct openwateroutputnode* openwateroutput_table = NULL;

struct openwaterdynsetnode
{ 
  char *name;
  void (*func)(int*,int*);

  struct openwaterdynsetnode* next;
} ;
struct openwaterdynsetnode* openwaterdynset_table = NULL;

struct openwaterf2tnode
{ 
  char *name;
  void (*func)(int*,int*);

  struct openwaterf2tnode* next;
} ;
struct openwaterf2tnode* openwaterf2t_table = NULL;

struct openwaterwriterstnode
{ 
  char *name;
  void (*func)(int*);

  struct openwaterwriterstnode* next;
} ;
struct openwaterwriterstnode* openwaterwriterst_table = NULL;


//BOP
// !ROUTINE: registeropenwaterini
// \label{registeropenwaterini}
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine to 
//  perform land surface model initialization
// 
// !INTERFACE:
void FTN(registeropenwaterinit)(char *j, void (*func)(int*),int len)
//EOP
{ 
  int len1;
  struct openwaterinitnode* current;
  struct openwaterinitnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct openwaterinitnode*) malloc(sizeof(struct openwaterinitnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(openwaterinit_table == NULL){
    openwaterinit_table = pnode;
  }
  else{
    current = openwaterinit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: openwaterinit
// \label{openwaterinit}
//
// !INTERFACE:
void FTN(openwaterinit)(char *j,int *mtype, int len)
//
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to perform
//  land surface model initialization
// 
//  \begin{description}
//  \item[j]
//   name of the OPENWATER
//  \end{description}
//EOP
{ 

  struct openwaterinitnode* current;
  
  current = openwaterinit_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("init routine for OPENWATER %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(mtype); 
}
//BOP
// !ROUTINE: registeropenwaterrun
// \label{registeropenwaterrun}
//
// !INTERFACE:
void FTN(registeropenwaterrun)(char *j, void (*func)(int*, int*),int len)
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine
//  to run the land surface model 
// 
//  \begin{description}
//  \item[j]
//   name of the OPENWATER
//  \end{description}
//EOP
{ 
  int len1;
  struct openwaterrunnode* current;
  struct openwaterrunnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct openwaterrunnode*) malloc(sizeof(struct openwaterrunnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(openwaterrun_table == NULL){
    openwaterrun_table = pnode;
  }
  else{
    current = openwaterrun_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: openwaterrun
// \label{openwaterrun}
//
// !INTERFACE:
void FTN(openwaterrun)(char *j,int *n,int *mtype,int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to run the 
//  land surface model 
// 
//  \begin{description}
//  \item[j]
//   name of the OPENWATER
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{
  struct openwaterrunnode* current;
  
  current = openwaterrun_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("run routine for OPENWATER %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,mtype); 
}

//BOP
// !ROUTINE: registeropenwaterfinalize
// \label{registeropenwaterfinalize}
//
// !INTERFACE:
void FTN(registeropenwaterfinalize)(char *j, void (*func)(),int len)
//  
// !DESCRIPTION:
//  Creates an entry in the registry for the routine
//  to cleanup allocated structures specific to the 
//  land surface model
// 
//  \begin{description}
//  \item[j]
//   name of the OPENWATER
//  \end{description}
// 
//EOP
{ 
  int len1;
  struct openwaterfinalnode* current;
  struct openwaterfinalnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct openwaterfinalnode*) malloc(sizeof(struct openwaterfinalnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(openwaterfinal_table == NULL){
    openwaterfinal_table = pnode;
  }
  else{
    current = openwaterfinal_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: openwaterfinalize
// \label{openwaterfinalize}
//
// !INTERFACE:
void FTN(openwaterfinalize)(char *j,int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry for cleaning up
//  allocated structures specific to the land surface model
// 
//  \begin{description}
//  \item[j]
//   name of the OPENWATER
//  \end{description}
// 
//EOP
{
  struct openwaterfinalnode* current;
  
  current = openwaterfinal_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("finalize routine for OPENWATER %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}

//BOP
// !ROUTINE: registeropenwatersetup
// \label{registeropenwatersetup}
//
// !INTERFACE:
void FTN(registeropenwatersetup)(char *j, void (*func)(int*),int len)
//  
// !DESCRIPTION:
//  Makes an entry in the registry for the routine
//  to set up land surface model parameters 
// 
//  \begin{description}
//  \item[j]
//   name of the OPENWATER
//  \end{description}
//EOP
{ 
  int len1;
  struct openwatersetupnode* current;
  struct openwatersetupnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct openwatersetupnode*) malloc(sizeof(struct openwatersetupnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(openwatersetup_table == NULL){
    openwatersetup_table = pnode;
  }
  else{
    current = openwatersetup_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: openwatersetup
// \label{openwatersetup}
//
// !INTERFACE:
void FTN(openwatersetup)(char *j,int *mtype, int len)
//  
// !DESCRIPTION:  
//  Invokes the routine in the registry to set up 
//  land surface model parameters  
// 
//  \begin{description}
//  \item[j]
//   name of the OPENWATER
//  \end{description}
//EOP
{ 
  struct openwatersetupnode* current;
  
  current = openwatersetup_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("setup routine for OPENWATER %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(mtype); 
}

//BOP
// !ROUTINE: registeropenwaterrestart
// \label{registeropenwaterrestart}
// 
// !INTERFACE:
void FTN(registeropenwaterrestart)(char *j, void (*func)(int*),int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// restart the land surface model from a 
// previously saved state
// 
//  \begin{description}
//  \item[j]
//   name of the OPENWATER
//  \end{description}
//EOP
{ 
  int len1;
  struct openwaterrestartnode* current;
  struct openwaterrestartnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct openwaterrestartnode*) malloc(sizeof(struct openwaterrestartnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(openwaterrestart_table == NULL){
    openwaterrestart_table = pnode;
  }
  else{
    current = openwaterrestart_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: openwaterrestart
// \label{openwaterrestart}
//
// !INTERFACE:
void FTN(openwaterrestart)(char *j,int *mtype, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to 
//  restart the land surface model from a previously
//  saved state
// 
//  \begin{description}
//  \item[j]
//   name of the OPENWATER
//  \end{description}
//EOP
{ 
  struct openwaterrestartnode* current;
  
  current = openwaterrestart_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("read restart routine for OPENWATER %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(mtype); 
}

//BOP
// !ROUTINE: registeropenwaterdynsetup
// \label{registeropenwaterdynsetup}
// 
// !INTERFACE:
void FTN(registeropenwaterdynsetup)(char *j, void (*func)(int*, int*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  set the time dependent land surface parameters
// 
//  \begin{description}
//  \item[j]
//   name of the OPENWATER
//  \end{description}
//EOP
{ 
  int len1;
  struct openwaterdynsetnode* current;
  struct openwaterdynsetnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct openwaterdynsetnode*) malloc(sizeof(struct openwaterdynsetnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(openwaterdynset_table == NULL){
    openwaterdynset_table = pnode;
  }
  else{
    current = openwaterdynset_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: openwaterdynsetup
// \label{openwaterdynsetup}
// 
// !INTERFACE:
void FTN(openwaterdynsetup)(char *j, int *n, int *mtype,int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to set the 
//  time dependent land surface parameters
// 
//  \begin{description}
//  \item[j]
//   name of the OPENWATER
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{ 
  struct openwaterdynsetnode* current;
  
  current = openwaterdynset_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("dynamic setup routine for OPENWATER %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,mtype); 
}

//BOP
// !ROUTINE: registeropenwateroutput
// \label{registeropenwateroutput}
// 
// !INTERFACE:
void FTN(registeropenwateroutput)(char *j, void (*func)(int*, int*),int len)
//  
// !DESCRIPTION:
//  Makes an entry in the registry for the routine to 
//  perform land surface model output
// 
//  \begin{description}
//  \item[j]
//   name of the OPENWATER
//  \end{description}
//EOP
{ 
  int len1;
  struct openwateroutputnode* current;
  struct openwateroutputnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct openwateroutputnode*) malloc(sizeof(struct openwateroutputnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(openwateroutput_table == NULL){
    openwateroutput_table = pnode;
  }
  else{
    current = openwateroutput_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: openwateroutput
// \label{openwateroutput}
//
// !INTERFACE:
void FTN(openwateroutput)(char *j, int *n, int *mtype, int len)
//  
// !DESCRIPTION:
//  Invokes the routine from the registry to perform
//  land surface model output
// 
//  \begin{description}
//  \item[j]
//   name of the OPENWATER
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{ 

  struct openwateroutputnode* current;
  
  current = openwateroutput_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("output writing routine for OPENWATER %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n, mtype); 
}

//BOP
// !ROUTINE: registeropenwaterf2t
// \label{registeropenwaterf2t}
// 
// !INTERFACE:
void FTN(registeropenwaterf2t)(char *j, void (*func)(int*, int*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  transfer forcing to model tiles
// 
//  \begin{description}
//  \item[j]
//   name of the OPENWATER + name of the running mode
//  \item[j]
//   index of the runmode
//  \end{description}
//EOP
{ 
  int len1;
  struct openwaterf2tnode* current;
  struct openwaterf2tnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct openwaterf2tnode*) malloc(sizeof(struct openwaterf2tnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(openwaterf2t_table == NULL){
    openwaterf2t_table = pnode;
  }
  else{
    current = openwaterf2t_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: openwaterf2t
// \label{openwaterf2t}
// 
// !INTERFACE:
void FTN(openwaterf2t)(char *j, int *n, int *mtype, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to 
//  transfer forcing to model tiles
// 
//  \begin{description}
//  \item[j]
//   name of the OPENWATER
//  \item[j]
//   index of the runmode
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{ 
  struct openwaterf2tnode* current;
  
  current = openwaterf2t_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("f2t writing routine for OPENWATER and running mode %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,mtype); 
}

//BOP
// !ROUTINE: registeropenwaterwrst
// \label{registeropenwaterwrst}
// 
// !INTERFACE:
void FTN(registeropenwaterwrst)(char *j, void (*func)(int*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine 
//  to write restart files for a land surface model 
// 
//  \begin{description}
//  \item[j]
//   name of the OPENWATER
//  \end{description}
//EOP
{ 
  int len1;
  struct openwaterwriterstnode* current;
  struct openwaterwriterstnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct openwaterwriterstnode*) malloc(sizeof(struct openwaterwriterstnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(openwaterwriterst_table == NULL){
    openwaterwriterst_table = pnode;
  }
  else{
    current = openwaterwriterst_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: openwaterwrst
// \label{openwaterwrst}
// 
// !INTERFACE:
void FTN(openwaterwrst)(char *j, int *n, int len)
//  
// !DESCRIPTION:
//  Invokes the routine from the registry to write
//  restart files for a land surface model
// 
//  \begin{description}
//  \item[j]
//   name of the OPENWATER
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{ 

  struct openwaterwriterstnode* current;
  
  current = openwaterwriterst_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("write restart writing routine for OPENWATER %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}






