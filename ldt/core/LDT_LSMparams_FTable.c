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
// !MODULE: LDT_LSMparams_FTable
//  
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing the operation of different 
//  LSMparams specifications
// 
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"

struct lsmparamprocinitnode
{ 
  char *name;
  void (*func)();

  struct lsmparamprocinitnode* next;
} ;
struct lsmparamprocinitnode* lsmparamprocinit_table = NULL; 

struct lsmparamprocwheadernode
{ 
  char *name;
  void (*func)(int*, int*, int*, int*);

  struct lsmparamprocwheadernode* next;
} ;
struct lsmparamprocwheadernode* lsmparamprocwheader_table = NULL; 

struct lsmparamprocwdatanode
{ 
  char *name;
  void (*func)(int*, int*);

  struct lsmparamprocwdatanode* next;
} ;
struct lsmparamprocwdatanode* lsmparamprocwdata_table = NULL; 



//BOP
// !ROUTINE: registerlsmparamprocinit
// \label{registerlsmparamprocinit}
//
// !INTERFACE:
void FTN(registerlsmparamprocinit)(char *j,void (*func)(), int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read the initialize the specific LSM parameter processing
// 
// The arguments are: 
// \begin{description}
// \item[j]
//  name of the LSM
//  \end{description}
//EOP
{ 
  struct lsmparamprocinitnode* current;
  struct lsmparamprocinitnode* pnode; 
  // create node
  
  pnode=(struct lsmparamprocinitnode*) malloc(sizeof(struct lsmparamprocinitnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmparamprocinit_table == NULL){
    lsmparamprocinit_table = pnode;
  }
  else{
    current = lsmparamprocinit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lsmparamprocinit
// \label{lsmparamprocinit}
//
// !INTERFACE:
void FTN(lsmparamprocinit)(char *j, int len)
// !DESCRIPTION: 
//  Calls the routine from the registry to initialize the 
//  the LSM parameter processing init routine
// 
// The arguments are: 
// \begin{description}
//  \item[j]
//   Character-based LSM name
//  \end{description}
//
//EOP
{ 
  struct lsmparamprocinitnode* current;
  
  current = lsmparamprocinit_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Read LSMparam input routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}


//BOP
// !ROUTINE: registerlsmparamprocwriteheader
// \label{registerlsmparamprocwriteheader}
//
// !INTERFACE:
void FTN(registerlsmparamprocwriteheader)(char *j,void (*func)(int*, int*, int*, int*), int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine that
// writes the header information for a particular set of LSM 
// parameters. 
// 
// The arguments are: 
// \begin{description}
// \item[j]
//  name of the LSM
//  \end{description}
//EOP
{ 
  struct lsmparamprocwheadernode* current;
  struct lsmparamprocwheadernode* pnode; 
  // create node
  
  pnode=(struct lsmparamprocwheadernode*) malloc(sizeof(struct lsmparamprocwheadernode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmparamprocwheader_table == NULL){
    lsmparamprocwheader_table = pnode;
  }
  else{
    current = lsmparamprocwheader_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lsmparamprocwriteheader
// \label{lsmparamprocwriteheader}
//
// !INTERFACE:
void FTN(lsmparamprocwriteheader)(char *j, int *n, int *ftn, int *dimId, int *monthId, int len)
// !DESCRIPTION: 
//  Calls the routine from the registry to write the 
//  header information for a particular set of LSM 
//  parameters. 
// 
// The arguments are: 
// \begin{description}
//  \item[j]
//   Character-based LSM name
//  \end{description}
//
//EOP
{ 
  struct lsmparamprocwheadernode* current;
  
  current = lsmparamprocwheader_table;
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
// !ROUTINE: registerlsmparamprocwritedata
// \label{registerlsmparamprocwritedata}
//
// !INTERFACE:
void FTN(registerlsmparamprocwritedata)(char *j,void (*func)(int*, int*), int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine that
// writes the data information for a particular set of LSM 
// parameters. 
// 
// The arguments are: 
// \begin{description}
// \item[j]
//  name of the LSM
//  \end{description}
//EOP
{ 
  struct lsmparamprocwdatanode* current;
  struct lsmparamprocwdatanode* pnode; 
  // create node
  
  pnode=(struct lsmparamprocwdatanode*) malloc(sizeof(struct lsmparamprocwdatanode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(lsmparamprocwdata_table == NULL){
    lsmparamprocwdata_table = pnode;
  }
  else{
    current = lsmparamprocwdata_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lsmparamprocwritedata
// \label{lsmparamprocwritedata}
//
// !INTERFACE:
void FTN(lsmparamprocwritedata)(char *j, int *n, int *ftn, int len)
// !DESCRIPTION: 
//  Calls the routine from the registry to write the 
//  data information for a particular set of LSM 
//  parameters. 
// 
// The arguments are: 
// \begin{description}
//  \item[j]
//   Character-based LSM name
//  \end{description}
//
//EOP
{ 
  struct lsmparamprocwdatanode* current;
  
  current = lsmparamprocwdata_table;
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




