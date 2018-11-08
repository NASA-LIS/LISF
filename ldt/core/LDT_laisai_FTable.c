//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center Land Information System (LDT)
//
// See RELEASE_NOTES.txt for more information.
//
// The LDT source code and documentation are not in the public domain
// and may not be freely distributed.  Only qualified entities may receive 
// the source code and documentation. 
//
// Qualified entities must be covered by a Software Usage Agreement. 
// The Software Usage Agreement contains all the terms and conditions
// regarding the release of the LDT software.
//
// NASA GSFC MAKES NO REPRESENTATIONS ABOUT THE SUITABILITY OF THE
// SOFTWARE FOR ANY PURPOSE.  IT IS PROVIDED AS IS WITHOUT EXPRESS OR
// IMPLIED WARRANTY.  NEITHER NASA GSFC NOR THE US GOVERNMENT SHALL BE
// LIABLE FOR ANY DAMAGES SUFFERED BY THE USER OF THIS SOFTWARE.
//
// See the Software Usage Agreement for the full disclaimer of warranty.
//
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
//BOP
//
// !MODULE: LDT_laisai_FTable
//  
//
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing different sources of 
//  LAI/SAI data
//   
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include <string.h>

#include "ftn_drv.h"

struct laisetnode
{ 
  char *name;
  void (*func)();

  struct laisetnode* next;
} ;
struct laisetnode* laiset_table = NULL; 

struct lainode
{ 
  char *name;
  void (*func)(int*, float*, float*);

  struct lainode* next;
} ;

struct lainode* lai_table = NULL; 

struct sainode
{ 
  char *name;
  void (*func)(int*, float*);

  struct sainode* next;
} ;

struct sainode* sai_table = NULL; 

struct laiminnode
{ 
  char *name;
  void (*func)(int*, float*);

  struct laiminnode* next;
} ;
struct laiminnode* laimin_table = NULL; 

struct laimaxnode
{ 
  char *name;
  void (*func)(int*, float*);

  struct laimaxnode* next;
} ;
struct laimaxnode* laimax_table = NULL; 

//BOP
// !ROUTINE: registersetlaiattribs
// \label{registersetlaiattribs}
//  
// 
// !INTERFACE:
void FTN(registersetlaiattribs)(char *j, void (*func)(),int len)
// !DESCRIPTION:
//  Creates an entry in the registry for the routine to 
//  read lai data
// 
//  The arguments are: 
//  \begin{description}
//   \item[j]
//    index of the lai source
//   \end{description}
//EOP
{ 
  struct laisetnode* current;
  struct laisetnode* pnode; 
  // create node

  pnode=(struct laisetnode*) malloc(sizeof(struct laisetnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(laiset_table == NULL){
    laiset_table = pnode;
  }
  else{
    current = laiset_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: setlaiattribs
// \label{setlaiattribs}
//  
// !DESCRIPTION: 
// Invokes the routine from the registry for 
// reading lai data. 
//
//  The arguments are: 
//  \begin{description}
//   \item[n]
//    index of the nest
//   \item[j]
//    index of the lai source
//   \item[array]
//    pointer to the lai data
//   \item[marray]
//    pointer to the mask data
//  \end{description}
//
// !INTERFACE:
void FTN(setlaiattribs)(char *j,int len)
//EOP
{ 

  struct laisetnode* current;
  
  current = laiset_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;

    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("setLAIAttribs routine for source %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func();
}


//BOP
// !ROUTINE: registerreadlai
// \label{registerreadlai}
// 
// !INTERFACE:
void FTN(registerreadlai)(char *j,void (*func)(int*,float*,float*), int len)
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read lai data
//
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \item[j]
//  index of the LAI data source
//  \end{description}
//EOP
{ 

  struct lainode* current;
  struct lainode* pnode; 
  // create node
  
  pnode=(struct lainode*) malloc(sizeof(struct lainode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(lai_table == NULL){
    lai_table = pnode;
  }
  else{
    current = lai_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readlai
// \label{readlai}
//
// !INTERFACE:
void FTN(readlai)(char *j, int *n,float *array,float *marray,int len)
//  
// !DESCRIPTION:
// Invokes the routine from the registry to 
// reading lai data 
// 
// The arguments are: 
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//  index of the LAI data source
//  \item[array]
//  pointer to the LAI data
//  \item[marray]
//  pointer to the mask data
//  \end{description}
//EOP
{ 
  struct lainode* current;
  
  current = lai_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;

    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("LAI reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,array,marray); 
} 


//BOP
// !ROUTINE: registerreadsai
// \label{registerreadsai}
// 
// !INTERFACE:
void FTN(registerreadsai)(char *j,void (*func)(int*,float*), int len)
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read sai data
//
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \item[j]
//  index of the SAI data source
//  \end{description}
//EOP
{ 

  struct sainode* current;
  struct sainode* pnode; 
  // create node
  
  pnode=(struct sainode*) malloc(sizeof(struct sainode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(sai_table == NULL){
    sai_table = pnode;
  }
  else{
    current = sai_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readsai
// \label{readsai}
//
// !INTERFACE:
void FTN(readsai)(char *j, int *n,float *array, int len)
//  
// !DESCRIPTION:
// Invokes the routine from the registry to 
// reading sai data 
// 
// The arguments are: 
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//  index of the SAI data source
//  \item[array]
//  pointer to the SAI data
//  \end{description}
//EOP
{ 

  struct sainode* current;
  
  current = sai_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("SAI reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,array); 
}


//BOP
// !ROUTINE: registerreadlaimin
// \label{registerreadlaimin}
// 
// !INTERFACE:
void FTN(registerreadlaimin)(char *j,void (*func)(int*,float*),int len)
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read laimin data
//
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \item[j]
//  index of the greenness data source
//  \end{description}
//EOP
{ 
  struct laiminnode* current;
  struct laiminnode* pnode; 
  // create node
  
  pnode=(struct laiminnode*) malloc(sizeof(struct laiminnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(laimin_table == NULL){
    laimin_table = pnode;
  }
  else{
    current = laimin_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readlaimin
// \label{readlaimin}
//
// !INTERFACE:
void FTN(readlaimin)(char *j, int *n,float *array,int len)
//  
// !DESCRIPTION:
// Invokes the routine from the registry to 
// reading laimin data 
// 
// The arguments are: 
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//  index of the greenness data source
//  \item[array]
//  pointer to the greenness data
//  \end{description}
//EOP
{ 
  struct laiminnode* current;
  
  current = laimin_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Min LAI reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,array); 
}


//BOP
// !ROUTINE: registerreadlaimax
// \label{registerreadlaimax}
// 
// !INTERFACE:
void FTN(registerreadlaimax)(char *j,void (*func)(int*,float*),int len)
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read laimax data
//
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \item[j]
//  index of the greenness data source
//  \end{description}
//EOP
{ 
  struct laimaxnode* current;
  struct laimaxnode* pnode; 
  // create node
  
  pnode=(struct laimaxnode*) malloc(sizeof(struct laimaxnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(laimax_table == NULL){
    laimax_table = pnode;
  }
  else{
    current = laimax_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readlaimax
// \label{readlaimax}
//
// !INTERFACE:
void FTN(readlaimax)(char *j, int *n,float *array,int len)
//  
// !DESCRIPTION:
// Invokes the routine from the registry to 
// reading laimax data 
// 
// The arguments are: 
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//  index of the greenness data source
//  \item[array]
//  pointer to the greenness data
//  \end{description}
//EOP
{ 
  struct laimaxnode* current;
  
  current = laimax_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Max LAI reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,array); 
}






