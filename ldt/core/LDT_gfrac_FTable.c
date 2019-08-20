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
// !MODULE: LDT_gfrac_FTable
//  
//
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing different sources of 
//  greenness fraction data
//   
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include <string.h>

#include "ftn_drv.h"

struct gfracnode
{ 
  char *name;
  void (*func)(int*, float*, float*);

  struct gfracnode* next;
} ;
struct gfracnode* gfrac_table = NULL; 

struct gfracsetnode
{ 
  char *name;
  void (*func)();

  struct gfracsetnode* next;
} ;
struct gfracsetnode* gfracset_table = NULL; 


struct shdminnode
{ 
  char *name;
  void (*func)(int*, float*);

  struct shdminnode* next;
} ;
struct shdminnode* shdmin_table = NULL; 

struct shdmaxnode
{ 
  char *name;
  void (*func)(int*, float*);

  struct shdmaxnode* next;
} ;
struct shdmaxnode* shdmax_table = NULL; 

//BOP
// !ROUTINE: registersetgfracattribs
// \label{registersetgfracattribs}
//  
// 
// !INTERFACE:
void FTN(registersetgfracattribs)(char *j, void (*func)(),int len)
// !DESCRIPTION:
//  Creates an entry in the registry for the routine to 
//  read gfrac data
// 
//  The arguments are: 
//  \begin{description}
//   \item[j]
//    index of the gfrac source
//   \end{description}
//EOP
{ 
  struct gfracsetnode* current;
  struct gfracsetnode* pnode; 
  // create node

  pnode=(struct gfracsetnode*) malloc(sizeof(struct gfracsetnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(gfracset_table == NULL){
    gfracset_table = pnode;
  }
  else{
    current = gfracset_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: setgfracattribs
// \label{setgfracattribs}
//  
// !DESCRIPTION: 
// Invokes the routine from the registry for 
// reading gfrac data. 
//
//  The arguments are: 
//  \begin{description}
//   \item[n]
//    index of the nest
//   \item[j]
//    index of the gfrac source
//   \item[array]
//    pointer to the gfrac data
//   \item[marray]
//    pointer to the mask data
//  \end{description}
//
// !INTERFACE:
void FTN(setgfracattribs)(char *j,int len)
//EOP
{ 

  struct gfracsetnode* current;
  
  current = gfracset_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;

    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("setGfracAttribs routine for source %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func();
}

//BOP
// !ROUTINE: registerreadgfrac
// \label{registerreadgfrac}
// 
// !INTERFACE:
void FTN(registerreadgfrac)(char *j,void (*func)(int*,float*,float*),int len)
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read gfrac data
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
  struct gfracnode* current;
  struct gfracnode* pnode; 
  // create node
  
  pnode=(struct gfracnode*) malloc(sizeof(struct gfracnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(gfrac_table == NULL){
    gfrac_table = pnode;
  }
  else{
    current = gfrac_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readgfrac
// \label{readgfrac}
//
// !INTERFACE:
void FTN(readgfrac)(char *j, int *n,float *array,float *marray,int len)
//  
// !DESCRIPTION:
// Invokes the routine from the registry to 
// reading gfrac data 
// 
// The arguments are: 
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//  index of the greenness data source
//  \item[array]
//  pointer to the greenness data
//  \item[marray]
//  pointer to the mask array data (optional)
//  \end{description}
//EOP
{ 
  struct gfracnode* current;
  
  current = gfrac_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Greenness fraction reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,array,marray); 
}


//BOP
// !ROUTINE: registerreadshdmin
// \label{registerreadshdmin}
// 
// !INTERFACE:
void FTN(registerreadshdmin)(char *j,void (*func)(int*,float*),int len)
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read shdmin data
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
  struct shdminnode* current;
  struct shdminnode* pnode; 
  // create node
  
  pnode=(struct shdminnode*) malloc(sizeof(struct shdminnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(shdmin_table == NULL){
    shdmin_table = pnode;
  }
  else{
    current = shdmin_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readshdmin
// \label{readshdmin}
//
// !INTERFACE:
void FTN(readshdmin)(char *j, int *n,float *array,int len)
//  
// !DESCRIPTION:
// Invokes the routine from the registry to 
// reading shdmin data 
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
  struct shdminnode* current;
  
  current = shdmin_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Min GVF reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,array); 
}


//BOP
// !ROUTINE: registerreadshdmax
// \label{registerreadshdmax}
// 
// !INTERFACE:
void FTN(registerreadshdmax)(char *j,void (*func)(int*,float*),int len)
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read shdmax data
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
  struct shdmaxnode* current;
  struct shdmaxnode* pnode; 
  // create node
  
  pnode=(struct shdmaxnode*) malloc(sizeof(struct shdmaxnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(shdmax_table == NULL){
    shdmax_table = pnode;
  }
  else{
    current = shdmax_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readshdmax
// \label{readshdmax}
//
// !INTERFACE:
void FTN(readshdmax)(char *j, int *n,float *array,int len)
//  
// !DESCRIPTION:
// Invokes the routine from the registry to 
// reading shdmax data 
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
  struct shdmaxnode* current;
  
  current = shdmax_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Max GVF reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,array); 
}




