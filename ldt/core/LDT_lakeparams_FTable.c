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
  struct lakeparamprocinitnode* current;
  struct lakeparamprocinitnode* pnode; 
  // create node
  
  pnode=(struct lakeparamprocinitnode*) malloc(sizeof(struct lakeparamprocinitnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
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
  struct lakeparamprocwheadernode* current;
  struct lakeparamprocwheadernode* pnode; 
  // create node
  
  pnode=(struct lakeparamprocwheadernode*) malloc(sizeof(struct lakeparamprocwheadernode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
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
  struct lakeparamprocwdatanode* current;
  struct lakeparamprocwdatanode* pnode; 
  // create node
  
  pnode=(struct lakeparamprocwdatanode*) malloc(sizeof(struct lakeparamprocwdatanode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
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




