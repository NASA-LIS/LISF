//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center Land Data Toolkit (LDT)
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
// !MODULE: LDT_metforcing_FTable
//  
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing different meteorological 
//  analyses implmented as "metforcings"
// 
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"

struct forcinginputnode
{ 
  char *name;
  void (*func)(int*);

  struct forcinginputnode* next;
} ;
struct forcinginputnode* forcinginput_table = NULL; 

struct forcinggetnode
{ 
  char *name;
  void (*func)(int*,int*);

  struct forcinggetnode* next;
} ;
struct forcinggetnode* forcingget_table = NULL; 

struct forcingtinterpnode
{ 
  char *name;
  void (*func)(int*,int*);

  struct forcingtinterpnode* next;
} ;
struct forcingtinterpnode* forcingtinterp_table = NULL; 

struct forcingresetnode
{ 
  char *name;
  void (*func)(int*);

  struct forcingresetnode* next;
} ;
struct forcingresetnode* forcingreset_table = NULL; 

struct forcingfinalnode
{ 
  char *name;
  void (*func)(int*);

  struct forcingfinalnode* next;
} ;
struct forcingfinalnode* forcingfinal_table = NULL; 

//BOP
// !ROUTINE: registerinitmetforc
// \label{registerinitmetforc}
// 
// !INTERFACE:
void FTN(registerinitmetforc)(char *j, void (*func)(int*), int len)
//  
// !DESCRIPTION: 
// Creates an entry in the registry for the routine
// that define the native domain of the met 
// forcing scheme. 
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the met forcing scheme
// \end{description}
//EOP
{ 
  struct forcinginputnode* current;
  struct forcinginputnode* pnode; 
  // create node
  
  pnode=(struct forcinginputnode*) malloc(sizeof(struct forcinginputnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(forcinginput_table == NULL){
    forcinginput_table = pnode;
  }
  else{
    current = forcinginput_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: initmetforc
// \label{initmetforc}
// 
// !INTERFACE:
void FTN(initmetforc)(char *j,int *findex, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry that 
//  defines the native domain of the met
//  forcine scheme. 
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the met forcing scheme
//  \item[findex]
//   order of the met forcing scheme among all selected
//   forcings in a LDT instance
// \end{description}
//EOP
{   
  struct forcinginputnode* current;
  
  current = forcinginput_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("**************** Error ****************************\n"); 
      printf("Init routine for forcing %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("**************** Error ****************************\n"); 
    }
  }
  current->func(findex); 

}

//BOP
// !ROUTINE: registerretrievemetforc
// \label{registerretrievemetforc}
//
// !INTERFACE:
void FTN(registerretrievemetforc)(char *j,void (*func)(int*, int*), int len)
//  
// !DESCRIPTION: 
// Creates an entry in the registry for the routines to open and 
// read met forcing
// 
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the met forcing scheme
// \end{description}
//EOP
{ 

  struct forcinggetnode* current;
  struct forcinggetnode* pnode; 
  // create node
  
  pnode=(struct forcinggetnode*) malloc(sizeof(struct forcinggetnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(forcingget_table == NULL){
    forcingget_table = pnode;
  }
  else{
    current = forcingget_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: retrievemetforc
// \label{retrievemetforc}
//  
// !INTERFACE:
void FTN(retrievemetforc)(char *j, int *n, int *findex, int len)
// !DESCRIPTION: 
//  Invokes the routine from the registry to open and read
//  the met forcing
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the met forcing scheme
//  \item[n]
//   index of the nest
//  \item[findex]
//   order of the met forcing scheme among all selected
//   forcings in a LDT instance
// \end{description}
//EOP
{ 
  struct forcinggetnode* current;
  
  current = forcingget_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("**************** Error ****************************\n"); 
      printf("Get routine for forcing %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("**************** Error ****************************\n"); 
    }
  }
  current->func(n,findex); 
}
//BOP
// !ROUTINE: registertimeinterpmetforc
// \label{registertimeinterpmetforc}
// 
// !INTERFACE:
void FTN(registertimeinterpmetforc)(char *j,void (*func)(int*, int*), int len)
//  
// !DESCRIPTION: 
//  Creats an entry in the registry for the routine to 
//  perform temporal interpolation
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the met forcing scheme
// \end{description}
//EOP
{ 
  struct forcingtinterpnode* current;
  struct forcingtinterpnode* pnode; 
  // create node
  
  pnode=(struct forcingtinterpnode*) malloc(sizeof(struct forcingtinterpnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(forcingtinterp_table == NULL){
    forcingtinterp_table = pnode;
  }
  else{
    current = forcingtinterp_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: timeinterpmetforc
// \label{timeinterpmetforc}
// 
// !INTERFACE:
void FTN(timeinterpmetforc)(char *j, int *n, int *findex, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to 
//  perform temporal interpolation
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the met forcing scheme
//  \item[n]
//   index of the nest
//  \item[findex]
//   order of the met forcing scheme among all selected
//   forcings in a LDT instance
// \end{description}
//EOP
{ 
  struct forcingtinterpnode* current;
  
  current = forcingtinterp_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("**************** Error ****************************\n"); 
      printf("Time interp routine for forcing %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("**************** Error ****************************\n"); 
    }
  }
  current->func(n,findex); 
}

//BOP
// !ROUTINE: registerfinalmetforc
// \label{registerfinalmetforc}
// 
// !INTERFACE:
void FTN(registerfinalmetforc)(char *j,void (*func)(int*), int len)
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine
//  to cleanup allocated structures
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the met forcing scheme
// \end{description}
//EOP
{ 
  struct forcingfinalnode* current;
  struct forcingfinalnode* pnode; 
  // create node
  
  pnode=(struct forcingfinalnode*) malloc(sizeof(struct forcingfinalnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(forcingfinal_table == NULL){
    forcingfinal_table = pnode;
  }
  else{
    current = forcingfinal_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: finalmetforc
// \label{finalmetforc}
// 
// !INTERFACE:
void FTN(finalmetforc)(char *j, int *findex, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to cleanup the 
//  allocated structures associated with the met forcing 
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the met forcing scheme
//  \item[findex]
//   order of the met forcing scheme among all selected
//   forcings in a LDT instance
// \end{description}
//EOP
{ 
  struct forcingfinalnode* current;
  
  current = forcingfinal_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("**************** Error ****************************\n"); 
      printf("Finalize routine for forcing %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("**************** Error ****************************\n"); 
    }
  }
  current->func(findex); 
}

//BOP
// !ROUTINE: registerresetmetforc
// \label{registerresetmetforc}
// 
// !INTERFACE:
void FTN(registerresetmetforc)(char *j,void (*func)(int*), int len)
//  
// !DESCRIPTION: 
//  Creates an entry in the registry for the routine
//  to reset the data structures
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the met forcing scheme
// \end{description}
//EOP
{ 
  struct forcingresetnode* current;
  struct forcingresetnode* pnode; 
  // create node
  
  pnode=(struct forcingresetnode*) malloc(sizeof(struct forcingresetnode));
  pnode->name=(char*) malloc(len*sizeof(char));
  strcpy(pnode->name,j);
  pnode->func = func;
  pnode->next = NULL; 

  if(forcingreset_table == NULL){
    forcingreset_table = pnode;
  }
  else{
    current = forcingreset_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}
//BOP
// !ROUTINE: resetmetforc
// \label{resetmetforc}
// 
// !INTERFACE:
void FTN(resetmetforc)(char *j, int *findex, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to reset the 
//  data structures associated with the met forcing 
//
// The arguments are: 
// \begin{description}
//  \item[j]
//   name of the met forcing scheme
//  \item[findex]
//   order of the met forcing scheme among all selected
//   forcings in a LDT instance
// \end{description}
//EOP
{ 
  struct forcingresetnode* current;
  
  current = forcingreset_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("**************** Error ****************************\n"); 
      printf("Reset routine for forcing %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("**************** Error ****************************\n"); 
    }
  }
  current->func(findex); 
}

