//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center Land Information System (LDT) V5.0 BETA
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
// !MODULE: LDT_domain_FTable
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

struct domaininputnode
{ 
  char *name;
  void (*func)();

  struct domaininputnode* next;
} ;
struct domaininputnode* domaininput_table = NULL; 

struct domainmakenode
{ 
  char *name;
  void (*func)();

  struct domainmakenode* next;
} ;
struct domainmakenode* domainmake_table = NULL; 

//BOP
// !ROUTINE: registerinput
// \label{registerinput}
//
// !INTERFACE:
void FTN(registerinput)(char *j,void (*func)(), int len)
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
  struct domaininputnode* current;
  struct domaininputnode* pnode; 
  // create node
  
  pnode=(struct domaininputnode*) malloc(sizeof(struct domaininputnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(domaininput_table == NULL){
    domaininput_table = pnode;
  }
  else{
    current = domaininput_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: readinput
// \label{readinput}
//
// !INTERFACE:
void FTN(readinput)(char *j, int len)
// !DESCRIPTION: 
//  Calls the routine from the registry to read the
//  the runtime domain specifics
// 
// The arguments are: 
// \begin{description}
//  \item[j]
//   Character-based domain option
//  \end{description}
//
//EOP
{ 
  struct domaininputnode* current;
  
  current = domaininput_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Read domain input routine for %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}

//BOP
// !ROUTINE: registerdomain
// \label{registerdomain}
//  
// 
// !INTERFACE:
void FTN(registerdomain)(char *j,void (*func)(), int len)
// !DESCRIPTION: 
// Creates an entry in the registry for the routine 
// to create grid and tile spaces for a particular 
// domain definition
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//EOP
{ 
  struct domainmakenode* current;
  struct domainmakenode* pnode; 
  // create node
  
  pnode=(struct domainmakenode*) malloc(sizeof(struct domainmakenode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(domainmake_table == NULL){
    domainmake_table = pnode;
  }
  else{
    current = domainmake_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}
//BOP
// !ROUTINE: makedomain
// \label{makedomain}
//
// !INTERFACE:
void FTN(makedomain)(char *j, int len)
//  
// !DESCRIPTION: 
// Invokes the routine from the registry 
// to create grid and tile spaces for a particular 
// domain definition
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//EOP
{ 
  struct domainmakenode* current;
  
  current = domainmake_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Make domain routine for %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}




