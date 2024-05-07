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
// !MODULE: LDT_lakeparams_FTable
//  
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing the operation of different 
//  lakeparams specifications
// 
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"

struct lakeparamprocinitnode
{ 
  char *name;
  void (*func)();

  struct lakeparamprocinitnode* next;
} ;
struct lakeparamprocinitnode* lakeparamprocinit_table = NULL; 

struct lakeparamprocwheadernode
{ 
  char *name;
  void (*func)(int*, int*, int*, int*);

  struct lakeparamprocwheadernode* next;
} ;
struct lakeparamprocwheadernode* lakeparamprocwheader_table = NULL; 

struct lakeparamprocwdatanode
{ 
  char *name;
  void (*func)(int*, int*);

  struct lakeparamprocwdatanode* next;
} ;
struct lakeparamprocwdatanode* lakeparamprocwdata_table = NULL; 



//BOP
// !ROUTINE: registerlakeparamprocinit
// \label{registerlakeparamprocinit}
//
// !INTERFACE:
void FTN(registerlakeparamprocinit)(char *j,void (*func)(), int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine to 
// read the initialize the specific lake parameter processing
// 
// The arguments are: 
// \begin{description}
// \item[j]
//  name of the lake
//  \end{description}
//EOP
{ 
  int len1;
  struct lakeparamprocinitnode* current;
  struct lakeparamprocinitnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lakeparamprocinitnode*) malloc(sizeof(struct lakeparamprocinitnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lakeparamprocinit_table == NULL){
    lakeparamprocinit_table = pnode;
  }
  else{
    current = lakeparamprocinit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lakeparamprocinit
// \label{lakeparamprocinit}
//
// !INTERFACE:
void FTN(lakeparamprocinit)(char *j, int len)
// !DESCRIPTION: 
//  Calls the routine from the registry to initialize the 
//  the lake parameter processing init routine
// 
// The arguments are: 
// \begin{description}
//  \item[j]
//   Character-based lake name
//  \end{description}
//
//EOP
{ 
  struct lakeparamprocinitnode* current;
  
  current = lakeparamprocinit_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Read lakeparam input routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}


//BOP
// !ROUTINE: registerlakeparamprocwriteheader
// \label{registerlakeparamprocwriteheader}
//
// !INTERFACE:
void FTN(registerlakeparamprocwriteheader)(char *j,void (*func)(int*, int*, int*, int*), int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine that
// writes the header information for a particular set of lake 
// parameters. 
// 
// The arguments are: 
// \begin{description}
// \item[j]
//  name of the lake
//  \end{description}
//EOP
{ 
  int len1;
  struct lakeparamprocwheadernode* current;
  struct lakeparamprocwheadernode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lakeparamprocwheadernode*) malloc(sizeof(struct lakeparamprocwheadernode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lakeparamprocwheader_table == NULL){
    lakeparamprocwheader_table = pnode;
  }
  else{
    current = lakeparamprocwheader_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lakeparamprocwriteheader
// \label{lakeparamprocwriteheader}
//
// !INTERFACE:
void FTN(lakeparamprocwriteheader)(char *j, int *n, int *ftn, int *dimId, int *monthId, int len)
// !DESCRIPTION: 
//  Calls the routine from the registry to write the 
//  header information for a particular set of lake 
//  parameters. 
// 
// The arguments are: 
// \begin{description}
//  \item[j]
//   Character-based lake name
//  \end{description}
//
//EOP
{ 
  struct lakeparamprocwheadernode* current;
  
  current = lakeparamprocwheader_table;
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
// !ROUTINE: registerlakeparamprocwritedata
// \label{registerlakeparamprocwritedata}
//
// !INTERFACE:
void FTN(registerlakeparamprocwritedata)(char *j,void (*func)(int*, int*), int len)
//  
// !DESCRIPTION: 
// Makes an entry in the registry for the routine that
// writes the data information for a particular set of lake 
// parameters. 
// 
// The arguments are: 
// \begin{description}
// \item[j]
//  name of the lake
//  \end{description}
//EOP
{ 
  int len1;
  struct lakeparamprocwdatanode* current;
  struct lakeparamprocwdatanode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct lakeparamprocwdatanode*) malloc(sizeof(struct lakeparamprocwdatanode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(lakeparamprocwdata_table == NULL){
    lakeparamprocwdata_table = pnode;
  }
  else{
    current = lakeparamprocwdata_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: lakeparamprocwritedata
// \label{lakeparamprocwritedata}
//
// !INTERFACE:
void FTN(lakeparamprocwritedata)(char *j, int *n, int *ftn, int len)
// !DESCRIPTION: 
//  Calls the routine from the registry to write the 
//  data information for a particular set of lake 
//  parameters. 
// 
// The arguments are: 
// \begin{description}
//  \item[j]
//   Character-based lake name
//  \end{description}
//
//EOP
{ 
  struct lakeparamprocwdatanode* current;
  
  current = lakeparamprocwdata_table;
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




