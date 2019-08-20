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
// !MODULE: LDT_routingparams_FTable
//  
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing the operation of different 
//  routingparams specifications
// 
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"

struct routingparamprocinitnode
{ 
  char *name;
  void (*func)();

  struct routingparamprocinitnode* next;
} ;
struct routingparamprocinitnode* routingparamprocinit_table = NULL; 

struct routingparamprocwheadernode
{ 
  char *name;
  void (*func)(int*, int*, int*, int*);

  struct routingparamprocwheadernode* next;
} ;
struct routingparamprocwheadernode* routingparamprocwheader_table = NULL; 

struct routingparamprocwdatanode
{ 
  char *name;
  void (*func)(int*, int*);

  struct routingparamprocwdatanode* next;
} ;
struct routingparamprocwdatanode* routingparamprocwdata_table = NULL; 



//BOP
// !ROUTINE: registerroutingparamprocinit
// \label{registerroutingparamprocinit}
//
// !INTERFACE:
void FTN(registerroutingparamprocinit)(char *j,void (*func)(), int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read the initialize the specific routing parameter processing
// 
// The arguments are: 
// \begin{description}
// \item[j]
//  name of the routing
//  \end{description}
//EOP
{ 
  struct routingparamprocinitnode* current;
  struct routingparamprocinitnode* pnode; 
  // create node
  
  pnode=(struct routingparamprocinitnode*) malloc(sizeof(struct routingparamprocinitnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingparamprocinit_table == NULL){
    routingparamprocinit_table = pnode;
  }
  else{
    current = routingparamprocinit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: routingparamprocinit
// \label{routingparamprocinit}
//
// !INTERFACE:
void FTN(routingparamprocinit)(char *j, int len)
// !DESCRIPTION: 
//  Calls the routine from the registry to initialize the 
//  the routing parameter processing init routine
// 
// The arguments are: 
// \begin{description}
//  \item[j]
//   Character-based routing name
//  \end{description}
//
//EOP
{ 
  struct routingparamprocinitnode* current;
  
  current = routingparamprocinit_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Read routingparam input routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}


//BOP
// !ROUTINE: registerroutingparamprocwriteheader
// \label{registerroutingparamprocwriteheader}
//
// !INTERFACE:
void FTN(registerroutingparamprocwriteheader)(char *j,void (*func)(int*, int*, int*, int*), int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine that
// writes the header information for a particular set of routing 
// parameters. 
// 
// The arguments are: 
// \begin{description}
// \item[j]
//  name of the routing
//  \end{description}
//EOP
{ 
  struct routingparamprocwheadernode* current;
  struct routingparamprocwheadernode* pnode; 
  // create node
  
  pnode=(struct routingparamprocwheadernode*) malloc(sizeof(struct routingparamprocwheadernode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingparamprocwheader_table == NULL){
    routingparamprocwheader_table = pnode;
  }
  else{
    current = routingparamprocwheader_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: routingparamprocwriteheader
// \label{routingparamprocwriteheader}
//
// !INTERFACE:
void FTN(routingparamprocwriteheader)(char *j, int *n, int *ftn, int *dimId, int *monthId, int len)
// !DESCRIPTION: 
//  Calls the routine from the registry to write the 
//  header information for a particular set of routing 
//  parameters. 
// 
// The arguments are: 
// \begin{description}
//  \item[j]
//   Character-based routing name
//  \end{description}
//
//EOP
{ 
  struct routingparamprocwheadernode* current;
  
  current = routingparamprocwheader_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("write header routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,ftn,dimId,monthId); 
}



//BOP
// !ROUTINE: registerroutingparamprocwritedata
// \label{registerroutingparamprocwritedata}
//
// !INTERFACE:
void FTN(registerroutingparamprocwritedata)(char *j,void (*func)(int*, int*), int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine that
// writes the data information for a particular set of routing 
// parameters. 
// 
// The arguments are: 
// \begin{description}
// \item[j]
//  name of the routing
//  \end{description}
//EOP
{ 
  struct routingparamprocwdatanode* current;
  struct routingparamprocwdatanode* pnode; 
  // create node
  
  pnode=(struct routingparamprocwdatanode*) malloc(sizeof(struct routingparamprocwdatanode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingparamprocwdata_table == NULL){
    routingparamprocwdata_table = pnode;
  }
  else{
    current = routingparamprocwdata_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: routingparamprocwritedata
// \label{routingparamprocwritedata}
//
// !INTERFACE:
void FTN(routingparamprocwritedata)(char *j, int *n, int *ftn, int len)
// !DESCRIPTION: 
//  Calls the routine from the registry to write the 
//  data information for a particular set of routing 
//  parameters. 
// 
// The arguments are: 
// \begin{description}
//  \item[j]
//   Character-based routing name
//  \end{description}
//
//EOP
{ 
  struct routingparamprocwdatanode* current;
  
  current = routingparamprocwdata_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("write data routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,ftn); 
}




