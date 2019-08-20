//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center Land Information System (LIS)
//
// See RELEASE.NOTES for more information.
//
// See SOFTWARE DISTRIBUTION POLICY for software distribution policies.
//
// The LIS source code and documentation are in the public domain,
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
// !MODULE: LDT_mask_FTable
//  
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing different sources of 
//  mask information
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"

struct mnode
{ 
  char *name;
  void (*func)(int*, float*);

  struct mnode* next;
} ;
struct mnode* mask_table = NULL; 

struct rgmnode
{
  char *name;
  void (*func)(int*, float*);

  struct rgmnode* next;
} ;
struct rgmnode* regmask_table = NULL;


//BOP
// !ROUTINE: registerrreadmask
// \label{registerreadmask}
//
// !INTERFACE: 
void FTN(registerreadmask)(char *j, void (*func)(int*, float*),int len)
//  
// !DESCRIPTION:
//  Creates an entry in the registry for the routine to 
//  read landmask data
// 
//  The arguments are: 
//  \begin{description}
//   \item[j]
//    index of the landcover source
//   \end{description}
//EOP
{ 
  struct mnode* current;
  struct mnode* pnode; 
  // create node
  
  pnode=(struct mnode*) malloc(sizeof(struct mnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(mask_table == NULL){
    mask_table = pnode;
  }
  else{
    current = mask_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readlandmask
// \label{readlandmask}
//
// !INTERFACE:
void FTN(readlandmask)(char *j, int *n, float *array, int len)
//  
// !DESCRIPTION: 
// Invokes the routine from the registry for 
// reading landmask data. 
//
//  The arguments are: 
//  \begin{description}
//   \item[n]
//    index of the nest
//   \item[j]
//    index of the landcover source
//   \item[array]
//    pointer to the landmask data
//   \end{description}
//EOP
{ 
  struct mnode* current;
  
  current = mask_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Landmask reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,array); 
}


//BOP
// !ROUTINE: registerrreadregmask
// \label{registerreadregmask}
//
// !INTERFACE:
void FTN(registerreadregmask)(char *j, void (*func)(int*, float*),int len)
//
// !DESCRIPTION:
//  Creates an entry in the registry for the routine to
//  read regional mask data
//
//  The arguments are:
//  \begin{description}
//   \item[j]
//    index of the regional mask source
//   \end{description}
//EOP
{
  struct rgmnode* current;
  struct rgmnode* pnode;
  // create node

  pnode=(struct rgmnode*) malloc(sizeof(struct rgmnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL;

  if(regmask_table == NULL){
    regmask_table = pnode;
  }
  else{
    current = regmask_table;
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode;
  }
}

//BOP
// !ROUTINE: readregmask
// \label{readregmask}
//
// !INTERFACE:
void FTN(readregmask)(char *j, int *n, float *array, int len)
//
// !DESCRIPTION:
// Invokes the routine from the registry for
// reading regional mask data.
//
//  The arguments are:
//  \begin{description}
//   \item[n]
//    index of the nest
//   \item[j]
//    index of the regional mask source
//   \item[array]
//    pointer to the regional mask data
//   \end{description}
//EOP
{
  struct rgmnode* current;

  current = regmask_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;

    if(current==NULL) {
      printf("****************Error****************************\n");
      printf("Regional mask reading routine for source %s is not defined\n",j);
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n");
      printf("****************Error****************************\n");
    }
  }
  current->func(n,array);
}

