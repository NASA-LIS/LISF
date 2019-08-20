//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center Land Information System (LIS) v7.2
//
// Copyright (c) 2015 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
//BOP
//
// !MODULE: LIS_routing_FTable
//  
//
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing the operations of different 
//  routing models. 
//
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"
struct routinginitnode
{ 
  char *name;
  void (*func)();

  struct routinginitnode* next;
} ;
struct routinginitnode* routinginit_table = NULL; 

struct routingrrstnode
{ 
  char *name;
  void (*func)();

  struct routingrrstnode* next;
} ;
struct routingrrstnode* routingrrst_table = NULL; 

struct routingwrstnode
{ 
  char *name;
  void (*func)(int*);

  struct routingwrstnode* next;
} ;
struct routingwrstnode* routingwrst_table = NULL; 

struct routingrunnode
{ 
  char *name;
  void (*func)(int*);

  struct routingrunnode* next;
} ;
struct routingrunnode* routingrun_table = NULL; 

struct routingoutnode
{ 
  char *name;
  void (*func)(int*);

  struct routingoutnode* next;
} ;
struct routingoutnode* routingout_table = NULL; 

//BOP
// !ROUTINE: registerroutinginit
// \label{registerroutinginit}
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine to 
//  perform routing model initialization
// 
// !INTERFACE:
void FTN(registerroutinginit)(char *j, void (*func)(),int len)
//EOP
{ 
  struct routinginitnode* current;
  struct routinginitnode* pnode; 
  // create node
  
  pnode=(struct routinginitnode*) malloc(sizeof(struct routinginitnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(routinginit_table == NULL){
    routinginit_table = pnode;
  }
  else{
    current = routinginit_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: routinginit
// \label{routinginit}
//
// !INTERFACE:
void FTN(routinginit)(char *j,int len)
//
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to perform
//  routing model initialization
// 
//  \begin{description}
//  \item[i]
//   index of the routing model
//  \end{description}
//EOP
{ 
  struct routinginitnode* current;
  
  current = routinginit_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("init routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}


//BOP
// !ROUTINE: registerroutingreadrestart
// \label{registerroutingreadrestart}
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine to 
//  perform routing model read restart. 
// 
// !INTERFACE:
void FTN(registerroutingreadrestart)(char *j, void (*func)(),int len)
//EOP
{ 
  struct routingrrstnode* current;
  struct routingrrstnode* pnode; 
  // create node
  
  pnode=(struct routingrrstnode*) malloc(sizeof(struct routingrrstnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingrrst_table == NULL){
    routingrrst_table = pnode;
  }
  else{
    current = routingrrst_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: routingreadrestart
// \label{routingreadrestart}
//
// !INTERFACE:
void FTN(routingreadrestart)(char *j,int len)
//
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to perform
//  routing model initialization
// 
//  \begin{description}
//  \item[i]
//   index of the routing model
//  \end{description}
//EOP
{ 
  struct routingrrstnode* current;
  
  current = routingrrst_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("read restart routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(); 
}

//BOP
// !ROUTINE: registerroutingrun
// \label{registerroutingrun}
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine to 
//  perform routing model run
// 
// !INTERFACE:
void FTN(registerroutingrun)(char *j, void (*func)(int*),int len)
//EOP
{ 
  struct routingrunnode* current;
  struct routingrunnode* pnode; 
  // create node
  
  pnode=(struct routingrunnode*) malloc(sizeof(struct routingrunnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingrun_table == NULL){
    routingrun_table = pnode;
  }
  else{
    current = routingrun_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: routingrun
// \label{routingrun}
//
// !INTERFACE:
void FTN(routingrun)(char *j, int *n,int len)
//
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to perform
//  routing model initialization
// 
//  \begin{description}
//  \item[i]
//   index of the routing model
//  \end{description}
//EOP
{ 
  struct routingrunnode* current;
  
  current = routingrun_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("run routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}

//BOP
// !ROUTINE: registerroutingoutput
// \label{registerroutingoutput}
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine to 
//  perform routing model output
// 
// !INTERFACE:
void FTN(registerroutingoutput)(char *j, void (*func)(int*),int len)
//EOP
{ 
  struct routingoutnode* current;
  struct routingoutnode* pnode; 
  // create node
  
  pnode=(struct routingoutnode*) malloc(sizeof(struct routingoutnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingout_table == NULL){
    routingout_table = pnode;
  }
  else{
    current = routingout_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: routingoutput
// \label{routingoutput}
//
// !INTERFACE:
void FTN(routingoutput)(char *j, int *n,int len)
//
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to perform
//  routing model initialization
// 
//  \begin{description}
//  \item[i]
//   index of the routing model
//  \end{description}
//EOP
{ 
  struct routingoutnode* current;
  
  current = routingout_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("output routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}

//BOP
// !ROUTINE: registerroutingwriterestart
// \label{registerroutingwriterestart}
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine to 
//  perform routing model writerestart
// 
// !INTERFACE:
void FTN(registerroutingwriterestart)(char *j, void (*func)(int*),int len)
//EOP
{ 
  struct routingwrstnode* current;
  struct routingwrstnode* pnode; 
  // create node
  
  pnode=(struct routingwrstnode*) malloc(sizeof(struct routingwrstnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(routingwrst_table == NULL){
    routingwrst_table = pnode;
  }
  else{
    current = routingwrst_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: routingwriterestart
// \label{routingwriterestart}
//
// !INTERFACE:
void FTN(routingwriterestart)(char *j, int *n,int len)
//
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to perform
//  routing model initialization
// 
//  \begin{description}
//  \item[i]
//   index of the routing model
//  \end{description}
//EOP
{ 
  struct routingwrstnode* current;
  
  current = routingwrst_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("write restart routine for %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n); 
}

