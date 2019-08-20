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
// !MODULE: LDT_anndata_FTable
//  
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing the operation of different 
//  domain specifications
// 
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"

struct anninpsetupnode
{ 
  char *name;
  void (*func)();

  struct anninpsetupnode* next;
} ;
struct anninpsetupnode* anninpsetup_table = NULL; 

struct annoutsetupnode
{ 
  char *name;
  void (*func)();

  struct annoutsetupnode* next;
} ;
struct annoutsetupnode* annoutsetup_table = NULL; 

struct anninputsourcenode
{ 
  char *name;
  void (*func)(int*, int*, int*, int*);

  struct anninputsourcenode* next;
} ;
struct anninputsourcenode* anninputsource_table = NULL; 

struct annoutputsourcenode
{ 
  char *name;
  void (*func)(int*, int*, int*, int*);

  struct annoutputsourcenode* next;
} ;
struct annoutputsourcenode* annoutputsource_table = NULL; 

//BOP
// !ROUTINE: registeranninputsourcesetup
// \label{registeranninputsourcesetup}
//
// !INTERFACE:
void FTN(registeranninputsourcesetup)(char *j,void (*func)(),int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read the runtime domain specifics
// 
// The arguments are: 
// \begin{description}
// \item[j]
//  name of the observation source
//  \end{description}
//EOP
{ 
  struct anninpsetupnode* current;
  struct anninpsetupnode* pnode; 
  // create node
  
  pnode=(struct anninpsetupnode*) malloc(sizeof(struct anninpsetupnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(anninpsetup_table == NULL){
    anninpsetup_table = pnode;
  }
  else{
    current = anninpsetup_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: setupanninputsource
// \label{setupanninputsource}
//
// !INTERFACE:
void FTN(setupanninputsource)(char *j, int len)
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
  struct anninpsetupnode* current;
  
  current = anninpsetup_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Observation setup routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}


//BOP
// !ROUTINE: registerannoutputsourcesetup
// \label{registerannoutputsourcesetup}
//
// !INTERFACE:
void FTN(registerannoutputsourcesetup)(char *j,void (*func)(),int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read the runtime domain specifics
// 
// The arguments are: 
// \begin{description}
// \item[j]
//  name of the observation source
//  \end{description}
//EOP
{ 
  struct annoutsetupnode* current;
  struct annoutsetupnode* pnode; 
  // create node
  
  pnode=(struct annoutsetupnode*) malloc(sizeof(struct annoutsetupnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(annoutsetup_table == NULL){
    annoutsetup_table = pnode;
  }
  else{
    current = annoutsetup_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: setupannoutputsource
// \label{setupannoutputsource}
//
// !INTERFACE:
void FTN(setupannoutputsource)(char *j, int len)
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
  struct annoutsetupnode* current;
  
  current = annoutsetup_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Observation setup routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}


//BOP
// !ROUTINE: registerreadanninputsource
// \label{registerreadanninputsource}
//
// !INTERFACE:
void FTN(registerreadanninputsource)(char *j,void (*func)(int*, int*, int*, int*), int len)
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
  struct anninputsourcenode* current;
  struct anninputsourcenode* pnode; 
  // create node
  
  pnode=(struct anninputsourcenode*) malloc(sizeof(struct anninputsourcenode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(anninputsource_table == NULL){
    anninputsource_table = pnode;
  }
  else{
    current = anninputsource_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: readanninputsource
// \label{readanninputsource}
//
// !INTERFACE:
void FTN(readanninputsource)(char *j, int *n, int *iomode, int *s, int *e,int len)
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
  struct anninputsourcenode* current;
  
  current = anninputsource_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Read observation source routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,iomode, s,e); 
}


//BOP
// !ROUTINE: registerreadannoutputsource
// \label{registerreadannoutputsource}
//
// !INTERFACE:
void FTN(registerreadannoutputsource)(char *j,void (*func)(int*, int*, int*, int*), int len)
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
  struct annoutputsourcenode* current;
  struct annoutputsourcenode* pnode; 
  // create node
  
  pnode=(struct annoutputsourcenode*) malloc(sizeof(struct annoutputsourcenode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(annoutputsource_table == NULL){
    annoutputsource_table = pnode;
  }
  else{
    current = annoutputsource_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: readannoutputsource
// \label{readannoutputsource}
//
// !INTERFACE:
void FTN(readannoutputsource)(char *j, int *n, int *iomode, int *s, int *e, int len)
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
  struct annoutputsourcenode* current;
  
  current = annoutputsource_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Read observation source routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,iomode,s, e); 
}





