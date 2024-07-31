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
// !MODULE: LDT_soils_FTable
//  
//
// !DESCRIPTION:
//   Function table registries for storing the interface implementations
//   for managing different sources of soil parameters data
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<string.h>

#include "ftn_drv.h"

struct soilfracsetnode
{ 
  char *name;
  void (*func)();

  struct soilfracsetnode* next;
} ;
struct soilfracsetnode* soilfracset_table = NULL; 

struct soilfracnode
{ 
  char *name;
  void (*func)(int*, int*, float*, float*, float*, float*);

  struct soilfracnode* next;
} ;

struct soilfracnode* soilfrac_table = NULL; 

struct colornode
{ 
  char *name;
  void (*func)(int*, float*);

  struct colornode* next;
} ;

struct colornode* color_table = NULL; 

struct txtsetnode
{ 
  char *name;
  void (*func)();

  struct txtsetnode* next;
} ;
struct txtsetnode* txtset_table = NULL; 

struct txtnode
{ 
  char *name;
  void (*func)(int*, int*, float*, float*);

  struct txtnode* next;
} ;

struct txtnode* txt_table = NULL; 

struct porosnode
{ 
  char *name;
  void (*func)(int*, float*, float*);

  struct porosnode* next;
} ;

struct porosnode* poros_table = NULL; 


//struct drootnode
//{ 
//  char *name;
//  void (*func)(int*, float*);
//
//  struct drootnode* next;
//} ;

//struct drootnode* droot_table = NULL; 

struct hsgsetnode
{
  char *name;
  void (*func)();

  struct hsgsetnode* next;
} ;
struct hsgsetnode* hsgset_table = NULL;

struct hsgnode
{
  char *name;
  void (*func)(int*, float*);

  struct hsgnode* next;
} ;

struct hsgnode* hsg_table = NULL;

struct dsoilnode
{ 
  char *name;
  void (*func)(int*, float*);

  struct dsoilnode* next;
} ;

struct dsoilnode* dsoil_table = NULL; 


//BOP
// !ROUTINE: registersetsoilfractionattribs
// \label{registersetsoilfractionattribs}
//  
// 
// !INTERFACE:
void FTN(registersetsoilfractionattribs)(char *j, void (*func)(),int len)
// !DESCRIPTION:
//  Creates an entry in the registry for the routine to 
//  read soilfraction data
// 
//  The arguments are: 
//  \begin{description}
//   \item[j]
//    index of the soilfraction source
//   \end{description}
//EOP
{ 
  int len1;
  struct soilfracsetnode* current;
  struct soilfracsetnode* pnode; 
  // create node

  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct soilfracsetnode*) malloc(sizeof(struct soilfracsetnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(soilfracset_table == NULL){
    soilfracset_table = pnode;
  }
  else{
    current = soilfracset_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: setsoilfractionattribs
// \label{setsoilfractionattribs}
//  
// !DESCRIPTION: 
// Invokes the routine from the registry for 
// reading soilfraction data. 
//
//  The arguments are: 
//  \begin{description}
//   \item[n]
//    index of the nest
//   \item[j]
//    index of the soilfraction source
//   \item[array]
//    pointer to the soilfraction data
//   \item[marray]
//    pointer to the mask data
//  \end{description}
//
// !INTERFACE:
void FTN(setsoilfractionattribs)(char *j,int len)
//EOP
{ 

  struct soilfracsetnode* current;
  
  current = soilfracset_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;

    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("setSoilfractionAttribs routine for source %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func();
}


//BOP
// !ROUTINE: registerreadsoilfrac
// \label{registerreadsoilfrac}
// 
// !INTERFACE:
void FTN(registerreadsoilfrac)(char *j, void (*func)(int*, int*, 
                float*, float*, float*, float*), int len)
//  
// !DESCRIPTION:
//  Creates an entry in the registry for the routine to 
//  read the soilfrac, clay and silt fraction data
// 
//  The arguments are: 
//  \begin{description}
//  \item[j]
//   index of the soils source
//  \end{description}
  //EOP
{ 
  int len1;
  struct soilfracnode* current;
  struct soilfracnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct soilfracnode*) malloc(sizeof(struct soilfracnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(soilfrac_table == NULL){
    soilfrac_table = pnode;
  }
  else{
    current = soilfrac_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readsoilfrac
// \label{readsoilfrac}
// 
// !INTERFACE:
void FTN(readsoilfrac)(char *j,int *n, int *nt, float *soilsfgrd,  
                       float *sand, float *clay, float *silt,int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to read the soilfrac
//  fraction data
// 
//  The arguments are: 
//  \begin{description}
//  \item[j]
//   index of the soils source
//  \item[n]
//   index of the nest
//   \item[nt]
//    number of types or bins (bands)
//   \item[fgrd]
//    gridcell fraction of type or values
//  \item[array]
//   pointer to the soilfrac array
//  \end{description}
//EOP
{ 
  struct soilfracnode* current;
  
  current = soilfrac_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Silt fraction reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,nt,soilsfgrd,sand,clay,silt); 
}

//BOP
// !ROUTINE: registersettextureattribs
// \label{registersettextureattribs}
//  
// 
// !INTERFACE:
void FTN(registersettextureattribs)(char *j, void (*func)(),int len)
// !DESCRIPTION:
//  Creates an entry in the registry for the routine to 
//  read texture data
// 
//  The arguments are: 
//  \begin{description}
//   \item[j]
//    index of the texture source
//   \end{description}
//EOP
{ 
  int len1;
  struct txtsetnode* current;
  struct txtsetnode* pnode; 
  // create node

  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct txtsetnode*) malloc(sizeof(struct txtsetnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(txtset_table == NULL){
    txtset_table = pnode;
  }
  else{
    current = txtset_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: settextureattribs
// \label{settextureattribs}
//  
// !DESCRIPTION: 
// Invokes the routine from the registry for 
// reading texture data. 
//
//  The arguments are: 
//  \begin{description}
//   \item[n]
//    index of the nest
//   \item[j]
//    index of the texture source
//   \item[array]
//    pointer to the texture data
//   \item[marray]
//    pointer to the mask data
//  \end{description}
//
// !INTERFACE:
void FTN(settextureattribs)(char *j,int len)
//EOP
{ 

  struct txtsetnode* current;
  
  current = txtset_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;

    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("setTextureAttribs routine for source %s is not defined\n",j); 
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func();
}

//BOP
// !ROUTINE: registerreadtexture
// \label{registerreadtexture}
//  
// !DESCRIPTION:
//  Creates an entry in the registry for the routine to 
//  read the soil texture data
// 
//  The arguments are: 
//  \begin{description}
//  \item[j]
//   index of the soils source
//  \end{description}
// 
// !INTERFACE:
void FTN(registerreadsoiltexture)(char *j,void (*func)(int*, int*, float*, float*),int len)
//EOP 
{ 
  int len1;
  struct txtnode* current;
  struct txtnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct txtnode*) malloc(sizeof(struct txtnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(txt_table == NULL){
    txt_table = pnode;
  }
  else{
    current = txt_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}

//BOP
// !ROUTINE: readsoiltexture
// \label{readsoiltexture}
// 
// !INTERFACE:
void FTN(readsoiltexture)(char *j,int *n,int *nt,float *array, float *arrayopt, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to read the
//  soil texture data
// 
//  The arguments are: 
//  \begin{description}
//  \item[n]
//   index of the nest
//  \item[nt]
//   number of types or bins
//  \item[j]
//   index of the soils source
//  \item[array]
//   pointer to the texture array
//  \item[arrayopt]
//   pointer to optional texture array with soil layers
//  \end{description}
//EOP
{ 
  struct txtnode* current;
  
  current = txt_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Texture reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,nt,array,arrayopt); 
}


//BOP
// !ROUTINE: registerreadporosity
// \label{registerreadporosity}
// 
// !INTERFACE:
void FTN(registerreadporosity)(char *j,void (*func)(int*, float*, float*),int len)
//  
// !DESCRIPTION:
//  Creates an entry in the registry for the routine to 
//  read the porosity data
// 
//  The arguments are: 
//  \begin{description}
//  \item[j]
//   index of the soils source
//  \end{description}
  //EOP
{
  int len1;
  struct porosnode* current;
  struct porosnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct porosnode*) malloc(sizeof(struct porosnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(poros_table == NULL){
    poros_table = pnode;
  }
  else{
    current = poros_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readporosity
// \label{readporosity}
// 
// !INTERFACE:
void FTN(readporosity)(char *j,int *n, float *array, float *marray, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to read the
//  porosity data
// 
//  The arguments are: 
//  \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//   index of the soils source
//  \item[array]
//   pointer to the porosity array
//  \end{description}
//EOP
{ 

  struct porosnode* current;
  
  current = poros_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Porosity reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,array,marray); 
}

//BOP
// !ROUTINE: registerreadrootdepth
// \label{registerreadrootdepth}
// 
// !INTERFACE:
//void FTN(registerreadrootdepth)(char *j,void (*func)(int*, float*),int len)
//  
// !DESCRIPTION:
//  Creates an entry in the registry for the routine to 
//  read the rootdepth fraction data
// 
//  The arguments are: 
//  \begin{description}
//  \item[j]
//   index of the soils source
//  \end{description}
  //EOP
//{ 

//  int len1;
//  struct drootnode* current;
// struct drootnode* pnode; 
  // create node
  
//  len1 = len + 1; // ensure that there is space for terminating null
//  pnode=(struct drootnode*) malloc(sizeof(struct drootnode));
//  pnode->name=(char*) calloc(len1,sizeof(char));
//  strncpy(pnode->name,j,len);
//  pnode->func = func;
//  pnode->next = NULL; 

//  if(droot_table == NULL){
//    droot_table = pnode;
//  }
//  else{
//    current = droot_table; 
//    while(current->next!=NULL){
//      current = current->next;
//    }
//    current->next = pnode; 
//  }
//}

//BOP
// !ROUTINE: readrootdepth
// \label{readrootdepth}
// 
// !INTERFACE:
//void FTN(readrootdepth)(char *j,int *n,float *array, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to read the
//  rootdepth fraction data
// 
//  The arguments are: 
//  \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//   index of the soils source
//  \item[array]
//   pointer to the rootdepth data
//  \end{description}
//EOP
//{ 

//  struct drootnode* current;
  
//  current = droot_table;
//  while(strcmp(current->name,j)!=0){
//    current = current->next;
//    if(current==NULL) {
//      printf("****************Error****************************\n"); 
//      printf("Rootdepth reading routine for source %s is not defined\n",j); 
//      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
//      printf("program will seg fault.....\n"); 
//      printf("****************Error****************************\n"); 
//    }
//  }
//  current->func(n,array); 
//}

//BOP
// !ROUTINE: registerreadcolor
// \label{registerreadcolor}
// 
// !INTERFACE:
void FTN(registerreadcolor)(char *j, void (*func)(int*, float*),int len)
//  
// !DESCRIPTION:
//  Creates an entry in the registry for the routine to 
//  read the soil color data
// 
//  The arguments are: 
//  \begin{description}
//  \item[j]
//   index of the soils source
//  \end{description}
  //EOP
{ 
  int len1;
  struct colornode* current;
  struct colornode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct colornode*) malloc(sizeof(struct colornode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(color_table == NULL){
    color_table = pnode;
  }
  else{
    current = color_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }

}

//BOP
// !ROUTINE: readcolor
// \label{readcolor}
// 
// !INTERFACE:
void FTN(readcolor)(char *j,int *n,float *array, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to read the
//  soil color data
// 
//  The arguments are: 
//  \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//   index of the soils source
//  \item[array]
//   pointer to the color array
//  \end{description}
//EOP
{ 
  struct colornode* current;
  
  current = color_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Soil color reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,array); 
}


//BOP
// !ROUTINE: registersethsgattribs
// \label{registersethsgattribs}
//  
// 
// !INTERFACE:
void FTN(registersethsgattribs)(char *j, void (*func)(),int len)
// !DESCRIPTION:
//  Creates an entry in the registry for the routine to 
//  read HSG data
// 
//  The arguments are: 
//  \begin{description}
//   \item[j]
//    index of the HSG source
//   \end{description}
//EOP
{
  int len1;
  struct hsgsetnode* current;
  struct hsgsetnode* pnode;
  // create node

  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct hsgsetnode*) malloc(sizeof(struct hsgsetnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL;

  if(hsgset_table == NULL){
    hsgset_table = pnode;
  }
  else{
    current = hsgset_table;
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode;
  }
}

//BOP
// !ROUTINE: sethsgattribs
// \label{sethsgattribs}
//  
// !DESCRIPTION: 
// Invokes the routine from the registry for 
// reading HSG data. 
//
//  The arguments are: 
//  \begin{description}
//  \end{description}
//
// !INTERFACE:
void FTN(sethsgattribs)(char *j,int len)
//EOP
{

  struct hsgsetnode* current;

  current = hsgset_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;

    if(current==NULL) {
      printf("****************Error****************************\n");
      printf("setTextureAttribs routine for source %s is not defined\n",j);
      printf("program will seg fault.....\n");
      printf("****************Error****************************\n");
    }
  }
  current->func();
}

//BOP
// !ROUTINE: registerreadhsg
// \label{registerreadhsg}
//
// !INTERFACE:
void FTN(registerreadhsg)(char *j, void (*func)(int*, float*), int len)
//
// !DESCRIPTION:
//  Makes an entry in the registry for the routine to read
//  hydrological soil group data
//
//  The arguments are:
//  \begin{description}
//   \item[j]
//    index of the topography source
//  \end{description}
  //EOP
{
  int len1;
  struct hsgnode* current;
  struct hsgnode* pnode;
  // create node

  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct hsgnode*) malloc(sizeof(struct hsgnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL;

  if(hsg_table == NULL){
    hsg_table = pnode;
  }
  else{
    current = hsg_table;
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode;
  }
}

//BOP
// !ROUTINE: readhsg
// \label{readhsg}
//
// !INTERFACE:
void FTN(readhsg)(char *j, int *n, float *array, int len)
//
// !DESCRIPTION:
//  Invokes the routine from the registry to read the
//  hydrological soils group data
//
//  The arguments are:
//  \begin{description}
//   \item[n]
//    index of the nest
//   \item[j]
//    index of the topography source
//   \item[array]
//    pointer to the hsg array
//  \end{description}
//EOP
{
  struct hsgnode* current;

  current = hsg_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
    if(current==NULL) {
      printf("****************Error****************************\n");
      printf("hsg reading routine for source %s is not defined\n",j);
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n");
      printf("****************Error****************************\n");
    }
  }
  current->func(n,array);
}

//BOP
// !ROUTINE: registerreadsoildepth
// \label{registerreadsoildepth}
// 
// !INTERFACE:
void FTN(registerreadsoildepth)(char *j,void (*func)(int*, float*),int len)
//  
// !DESCRIPTION:
//  Creates an entry in the registry for the routine to 
//  read the soildepth fraction data
// 
//  The arguments are: 
//  \begin{description}
//  \item[j]
//   index of the soils source
//  \end{description}
  //EOP
{ 
  int len1;
  struct dsoilnode* current;
  struct dsoilnode* pnode; 
  // create node
  
  len1 = len + 1; // ensure that there is space for terminating null
  pnode=(struct dsoilnode*) malloc(sizeof(struct dsoilnode));
  pnode->name=(char*) calloc(len1,sizeof(char));
  strncpy(pnode->name,j,len);
  pnode->func = func;
  pnode->next = NULL; 

  if(dsoil_table == NULL){
    dsoil_table = pnode;
  }
  else{
    current = dsoil_table; 
    while(current->next!=NULL){
      current = current->next;
    }
    current->next = pnode; 
  }
}

//BOP
// !ROUTINE: readsoildepth
// \label{readsoildepth}
// 
// !INTERFACE:
void FTN(readsoildepth)(char *j,int *n,float *array, int len)
//  
// !DESCRIPTION: 
//  Invokes the routine from the registry to read the
//  soildepth fraction data
// 
//  The arguments are: 
//  \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//   index of the soils source
//  \item[array]
//   pointer to the soildepth data
//  \end{description}
//EOP
{ 

  struct dsoilnode* current;
  
  current = dsoil_table;
  while(strcmp(current->name,j)!=0){
    current = current->next;
  
    if(current==NULL) {
      printf("****************Error****************************\n"); 
      printf("Soildepth reading routine for source %s is not defined\n",j); 
      printf("Please refer to configs/ldt.config_master or LDT User's Guide for options.\n");
      printf("program will seg fault.....\n"); 
      printf("****************Error****************************\n"); 
    }
  }
  current->func(n,array); 
}
