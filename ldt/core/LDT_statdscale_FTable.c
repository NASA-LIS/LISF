//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center Land Data Toolkit (LDT) V7.0 BETA
// Released January 2008
//
// See SOFTWARE DISTRIBUTION POLICY for software distribution policies
//
// The LDT source code and documentation are in the public domain,
// available without fee for educational, research, non-commercial and
// commercial purposes.  Users may distribute the binary or source
// code to third parties provided this statement appears on all copies and
// that no charge is made for such copies.
//
// NASA GSFC MAKES NO REPRESENTATIONS ABOUT THE SUITABILITY OF THE
// SOFTWARE FOR ANY PURPOSE.  IT IS PROVIDED AS IS WITHOUT EXPRESS OR
// IMPLIED WARRANTY.  NEITHER NASA GSFC NOR THE US GOVERNMENT SHALL BE
// LIABLE FOR ANY DAMAGES SUFFERED BY THE USER OF THIS SOFTWARE.
//
// See COPYRIGHT.TXT for copyright details.
//
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
//BOP
//
// !MODULE: LDT_statdscale_FTable
//  
// !DESCRIPTION:
// 
// 
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"

struct statdscaleinitnode
{ 
  char *name;
  void (*func)();

  struct statdscaleinitnode* next;
} ;
struct statdscaleinitnode* statdscaleinit_table = NULL; 

struct statdscalediagnode
{ 
  char *name;
  void (*func)(int*, int*);

  struct statdscalediagnode* next;
} ;
struct statdscalediagnode* statdscalediag_table = NULL; 

struct statdscalecomputenode
{ 
  char *name;
  void (*func)(int*);

  struct statdscalecomputenode* next;
} ;
struct statdscalecomputenode* statdscalecompute_table = NULL; 

struct statdscaleoutputnode
{ 
  char *name;
  void (*func)(int*);

  struct statdscaleoutputnode* next;
} ;
struct statdscaleoutputnode* statdscaleoutput_table = NULL; 

//BOP
// !ROUTINE: registerinitstatdscale
// \label{registerinitstatdscale}
//
// !INTERFACE:
void FTN(registerinitstatdscale)(char *j,void (*func)(), int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read the runtime domain specifics
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//EOP
{ 
  struct statdscaleinitnode* current;
  struct statdscaleinitnode* pnode; 
  // create node
  
  pnode=(struct statdscaleinitnode*) malloc(sizeof(struct statdscaleinitnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(statdscaleinit_table == NULL){
    statdscaleinit_table = pnode;
  }
  else{
    current = statdscaleinit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: initstatdscale
// \label{initstatdscale}
//
// !INTERFACE:
void FTN(initstatdscale)(char *j, int len)
// !DESCRIPTION: 
//  Calls the routine from the registry to read the
//  the runtime domain specifics
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//
//EOP
{ 
  struct statdscaleinitnode* current;
  
  current = statdscaleinit_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("init routine for downscaling method %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}

//BOP
// !ROUTINE: registerdiagnosestatdscale
// \label{registerdiagnosestatdscale}
//
// !INTERFACE:
void FTN(registerdiagnosestatdscale)(char *j,void (*func)(int*, int*), int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read the runtime domain specifics
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//EOP
{ 
  struct statdscalediagnode* current;
  struct statdscalediagnode* pnode; 
  // create node
  
  pnode=(struct statdscalediagnode*) malloc(sizeof(struct statdscalediagnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(statdscalediag_table == NULL){
    statdscalediag_table = pnode;
  }
  else{
    current = statdscalediag_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: diagnosestatdscale
// \label{diagnosestatdscale}
//
// !INTERFACE:
void FTN(diagnosestatdscale)(char *j, int *n, int *pass,int len)
// !DESCRIPTION: 
//  Calls the routine from the registry to read the
//  the runtime domain specifics
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//
//EOP
{ 
  struct statdscalediagnode* current;
  
  current = statdscalediag_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Diagnose routine for downscaling method %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,pass); 
}


//BOP
// !ROUTINE: registercomputestatdscale
// \label{registercomputestatdscale}
//
// !INTERFACE:
void FTN(registercomputestatdscale)(char *j,void (*func)(int*), int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read the runtime domain specifics
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//EOP
{ 
  struct statdscalecomputenode* current;
  struct statdscalecomputenode* pnode; 
  // create node
  
  pnode=(struct statdscalecomputenode*) malloc(sizeof(struct statdscalecomputenode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(statdscalecompute_table == NULL){
    statdscalecompute_table = pnode;
  }
  else{
    current = statdscalecompute_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: computestatdscale
// \label{computestatdscale}
//
// !INTERFACE:
void FTN(computestatdscale)(char *j, int *pass,int len)
// !DESCRIPTION: 
//  Calls the routine from the registry to read the
//  the runtime domain specifics
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//
//EOP
{ 
  struct statdscalecomputenode* current;
  
  current = statdscalecompute_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Compute routine for downscaling method %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(pass); 
}


//BOP
// !ROUTINE: registeroutputstatdscale
// \label{registeroutputstatdscale}
//
// !INTERFACE:
void FTN(registeroutputstatdscale)(char *j,void (*func)(int*), int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read the runtime domain specifics
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//EOP
{ 
  struct statdscaleoutputnode* current;
  struct statdscaleoutputnode* pnode; 
  // create node
  
  pnode=(struct statdscaleoutputnode*) malloc(sizeof(struct statdscaleoutputnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(statdscaleoutput_table == NULL){
    statdscaleoutput_table = pnode;
  }
  else{
    current = statdscaleoutput_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: outputstatdscale
// \label{outputstatdscale}
//
// !INTERFACE:
void FTN(outputstatdscale)(char *j, int *pass,int len)
// !DESCRIPTION: 
//  Calls the routine from the registry to read the
//  the runtime domain specifics
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//
//EOP
{ 
  struct statdscaleoutputnode* current;
  
  current = statdscaleoutput_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Output routine for downscaling method %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(pass); 
}

