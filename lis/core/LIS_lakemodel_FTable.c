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
// !MODULE: LIS_lakemodel_FTable
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
struct lakemodelinitnode
{ 
  char *name;
  void (*func)(int*);

  struct lakemodelinitnode* next;
} ;
struct lakemodelinitnode* lakemodelinit_table = NULL; 

struct lakemodelrunnode
{ 
  char *name;
  void (*func)(int*,int*);

  struct lakemodelrunnode* next;
} ;
struct lakemodelrunnode* lakemodelrun_table = NULL; 

struct lakemodelfinalnode
{ 
  char *name;
  void (*func)();

  struct lakemodelfinalnode* next;
} ;
struct lakemodelfinalnode* lakemodelfinal_table = NULL; 

struct lakemodelsetupnode
{ 
  char *name;
  void (*func)(int*);

  struct lakemodelsetupnode* next;
} ;
struct lakemodelsetupnode* lakemodelsetup_table = NULL;

struct lakemodelrestartnode
{ 
  char *name;
  void (*func)(int*);

  struct lakemodelrestartnode* next;
} ;
struct lakemodelrestartnode* lakemodelrestart_table = NULL;

struct lakemodeloutputnode
{ 
  char *name;
  void (*func)(int*, int*);

  struct lakemodeloutputnode* next;
} ;
struct lakemodeloutputnode* lakemodeloutput_table = NULL;

struct lakemodeldynsetnode
{ 
  char *name;
  void (*func)(int*,int*);

  struct lakemodeldynsetnode* next;
} ;
struct lakemodeldynsetnode* lakemodeldynset_table = NULL;

struct lakemodelf2tnode
{ 
  char *name;
  void (*func)(int*,int*);

  struct lakemodelf2tnode* next;
} ;
struct lakemodelf2tnode* lakemodelf2t_table = NULL;

struct lakemodelwriterstnode
{ 
  char *name;
  void (*func)(int*);

  struct lakemodelwriterstnode* next;
} ;
struct lakemodelwriterstnode* lakemodelwriterst_table = NULL;


//BOP
// !ROUTINE: registerlakemodelini
// \label{registerlakemodelini}
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine to 
//  perform land surface model initialization
// 
// !INTERFACE:
void FTN(registerlakemodelinit)(char *j, void (*func)(int*),int len)
//EOP
{ 
  int len1;
  struct lakemodelinitnode* current;
  struct lakemodelinitnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lakemodelinitnode*) malloc(sizeof(struct lakemodelinitnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lakemodelinit_table == NULL){
    lakemodelinit_table = pnode;
  }
  else{
    current = lakemodelinit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: lakemodelinit
// \label{lakemodelinit}
//
// !INTERFACE:
void FTN(lakemodelinit)(char *j,int *mtype, int len)
//
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to perform
//  land surface model initialization
// 
//  \begin{description}
//  \item[j]
//   name of the LAKEMODEL
//  \end{description}
//EOP
{ 

  struct lakemodelinitnode* current;
  
  current = lakemodelinit_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("init routine for LAKEMODEL %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(mtype); 
}
//BOP
// !ROUTINE: registerlakemodelrun
// \label{registerlakemodelrun}
//
// !INTERFACE:
void FTN(registerlakemodelrun)(char *j, void (*func)(int*, int*),int len)
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine
//  to run the land surface model 
// 
//  \begin{description}
//  \item[j]
//   name of the LAKEMODEL
//  \end{description}
//EOP
{ 
  int len1;
  struct lakemodelrunnode* current;
  struct lakemodelrunnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lakemodelrunnode*) malloc(sizeof(struct lakemodelrunnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lakemodelrun_table == NULL){
    lakemodelrun_table = pnode;
  }
  else{
    current = lakemodelrun_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lakemodelrun
// \label{lakemodelrun}
//
// !INTERFACE:
void FTN(lakemodelrun)(char *j,int *n,int *mtype,int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to run the 
//  land surface model 
// 
//  \begin{description}
//  \item[j]
//   name of the LAKEMODEL
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{
  struct lakemodelrunnode* current;
  
  current = lakemodelrun_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("run routine for LAKEMODEL %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,mtype); 
}

//BOP
// !ROUTINE: registerlakemodelfinalize
// \label{registerlakemodelfinalize}
//
// !INTERFACE:
void FTN(registerlakemodelfinalize)(char *j, void (*func)(),int len)
//  
// !DESCRIPTION:
//  Creates an entry in the registry for the routine
//  to cleanup allocated structures specific to the 
//  land surface model
// 
//  \begin{description}
//  \item[j]
//   name of the LAKEMODEL
//  \end{description}
// 
//EOP
{ 
  int len1;
  struct lakemodelfinalnode* current;
  struct lakemodelfinalnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lakemodelfinalnode*) malloc(sizeof(struct lakemodelfinalnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lakemodelfinal_table == NULL){
    lakemodelfinal_table = pnode;
  }
  else{
    current = lakemodelfinal_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: lakemodelfinalize
// \label{lakemodelfinalize}
//
// !INTERFACE:
void FTN(lakemodelfinalize)(char *j,int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry for cleaning up
//  allocated structures specific to the land surface model
// 
//  \begin{description}
//  \item[j]
//   name of the LAKEMODEL
//  \end{description}
// 
//EOP
{
  struct lakemodelfinalnode* current;
  
  current = lakemodelfinal_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("finalize routine for LAKEMODEL %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}

//BOP
// !ROUTINE: registerlakemodelsetup
// \label{registerlakemodelsetup}
//
// !INTERFACE:
void FTN(registerlakemodelsetup)(char *j, void (*func)(int*),int len)
//  
// !DESCRIPTION:
//  Makes an entry in the registry for the routine
//  to set up land surface model parameters 
// 
//  \begin{description}
//  \item[j]
//   name of the LAKEMODEL
//  \end{description}
//EOP
{ 
  int len1;
  struct lakemodelsetupnode* current;
  struct lakemodelsetupnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lakemodelsetupnode*) malloc(sizeof(struct lakemodelsetupnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lakemodelsetup_table == NULL){
    lakemodelsetup_table = pnode;
  }
  else{
    current = lakemodelsetup_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lakemodelsetup
// \label{lakemodelsetup}
//
// !INTERFACE:
void FTN(lakemodelsetup)(char *j,int *mtype, int len)
//  
// !DESCRIPTION:  
//  Invokes the routine in the registry to set up 
//  land surface model parameters  
// 
//  \begin{description}
//  \item[j]
//   name of the LAKEMODEL
//  \end{description}
//EOP
{ 
  struct lakemodelsetupnode* current;
  
  current = lakemodelsetup_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("setup routine for LAKEMODEL %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(mtype); 
}

//BOP
// !ROUTINE: registerlakemodelrestart
// \label{registerlakemodelrestart}
// 
// !INTERFACE:
void FTN(registerlakemodelrestart)(char *j, void (*func)(int*),int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// restart the land surface model from a 
// previously saved state
// 
//  \begin{description}
//  \item[j]
//   name of the LAKEMODEL
//  \end{description}
//EOP
{ 
  int len1;
  struct lakemodelrestartnode* current;
  struct lakemodelrestartnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lakemodelrestartnode*) malloc(sizeof(struct lakemodelrestartnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lakemodelrestart_table == NULL){
    lakemodelrestart_table = pnode;
  }
  else{
    current = lakemodelrestart_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lakemodelrestart
// \label{lakemodelrestart}
//
// !INTERFACE:
void FTN(lakemodelrestart)(char *j,int *mtype, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to 
//  restart the land surface model from a previously
//  saved state
// 
//  \begin{description}
//  \item[j]
//   name of the LAKEMODEL
//  \end{description}
//EOP
{ 
  struct lakemodelrestartnode* current;
  
  current = lakemodelrestart_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("read restart routine for LAKEMODEL %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(mtype); 
}

//BOP
// !ROUTINE: registerlakemodeldynsetup
// \label{registerlakemodeldynsetup}
// 
// !INTERFACE:
void FTN(registerlakemodeldynsetup)(char *j, void (*func)(int*, int*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  set the time dependent land surface parameters
// 
//  \begin{description}
//  \item[j]
//   name of the LAKEMODEL
//  \end{description}
//EOP
{ 
  int len1;
  struct lakemodeldynsetnode* current;
  struct lakemodeldynsetnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lakemodeldynsetnode*) malloc(sizeof(struct lakemodeldynsetnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lakemodeldynset_table == NULL){
    lakemodeldynset_table = pnode;
  }
  else{
    current = lakemodeldynset_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lakemodeldynsetup
// \label{lakemodeldynsetup}
// 
// !INTERFACE:
void FTN(lakemodeldynsetup)(char *j, int *n, int *mtype,int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to set the 
//  time dependent land surface parameters
// 
//  \begin{description}
//  \item[j]
//   name of the LAKEMODEL
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{ 
  struct lakemodeldynsetnode* current;
  
  current = lakemodeldynset_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("dynamic setup routine for LAKEMODEL %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,mtype); 
}

//BOP
// !ROUTINE: registerlakemodeloutput
// \label{registerlakemodeloutput}
// 
// !INTERFACE:
void FTN(registerlakemodeloutput)(char *j, void (*func)(int*, int*),int len)
//  
// !DESCRIPTION:
//  Makes an entry in the registry for the routine to 
//  perform land surface model output
// 
//  \begin{description}
//  \item[j]
//   name of the LAKEMODEL
//  \end{description}
//EOP
{ 
  int len1;
  struct lakemodeloutputnode* current;
  struct lakemodeloutputnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lakemodeloutputnode*) malloc(sizeof(struct lakemodeloutputnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lakemodeloutput_table == NULL){
    lakemodeloutput_table = pnode;
  }
  else{
    current = lakemodeloutput_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lakemodeloutput
// \label{lakemodeloutput}
//
// !INTERFACE:
void FTN(lakemodeloutput)(char *j, int *n, int *mtype, int len)
//  
// !DESCRIPTION:
//  Invokes the routine from the registry to perform
//  land surface model output
// 
//  \begin{description}
//  \item[j]
//   name of the LAKEMODEL
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{ 

  struct lakemodeloutputnode* current;
  
  current = lakemodeloutput_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("output writing routine for LAKEMODEL %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n, mtype); 
}

//BOP
// !ROUTINE: registerlakemodelf2t
// \label{registerlakemodelf2t}
// 
// !INTERFACE:
void FTN(registerlakemodelf2t)(char *j, void (*func)(int*, int*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine to 
//  transfer forcing to model tiles
// 
//  \begin{description}
//  \item[j]
//   name of the LAKEMODEL + name of the running mode
//  \item[j]
//   index of the runmode
//  \end{description}
//EOP
{ 
  int len1;
  struct lakemodelf2tnode* current;
  struct lakemodelf2tnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lakemodelf2tnode*) malloc(sizeof(struct lakemodelf2tnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lakemodelf2t_table == NULL){
    lakemodelf2t_table = pnode;
  }
  else{
    current = lakemodelf2t_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lakemodelf2t
// \label{lakemodelf2t}
// 
// !INTERFACE:
void FTN(lakemodelf2t)(char *j, int *n, int *mtype, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to 
//  transfer forcing to model tiles
// 
//  \begin{description}
//  \item[j]
//   name of the LAKEMODEL
//  \item[j]
//   index of the runmode
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{ 
  struct lakemodelf2tnode* current;
  
  current = lakemodelf2t_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("f2t writing routine for LAKEMODEL and running mode %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,mtype); 
}

//BOP
// !ROUTINE: registerlakemodelwrst
// \label{registerlakemodelwrst}
// 
// !INTERFACE:
void FTN(registerlakemodelwrst)(char *j, void (*func)(int*),int len)
//  
// !DESCRIPTION: 
//  Makes an entry in the registry for the routine 
//  to write restart files for a land surface model 
// 
//  \begin{description}
//  \item[j]
//   name of the LAKEMODEL
//  \end{description}
//EOP
{ 
  int len1;
  struct lakemodelwriterstnode* current;
  struct lakemodelwriterstnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lakemodelwriterstnode*) malloc(sizeof(struct lakemodelwriterstnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lakemodelwriterst_table == NULL){
    lakemodelwriterst_table = pnode;
  }
  else{
    current = lakemodelwriterst_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: lakemodelwrst
// \label{lakemodelwrst}
// 
// !INTERFACE:
void FTN(lakemodelwrst)(char *j, int *n, int len)
//  
// !DESCRIPTION:
//  Invokes the routine from the registry to write
//  restart files for a land surface model
// 
//  \begin{description}
//  \item[j]
//   name of the LAKEMODEL
//  \item[n]
//   index of the nest
//  \end{description}
//EOP
{ 

  struct lakemodelwriterstnode* current;
  
  current = lakemodelwriterst_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("write restart writing routine for LAKEMODEL %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}






