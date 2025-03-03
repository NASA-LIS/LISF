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
// !MODULE: LVT_metric_FTable
//  
// !DESCRIPTION:
//  Function table registries for storing the interface 
//  implementations for managing the operation of different 
//  metrics
// 
//EOP
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>

#include "ftn_drv.h"
#include "FTable.h"
typedef struct
{ 
  void (*func)(void*, void*, void*);
} metric_ini_TABLE; 
metric_ini_TABLE metric_ini[FT_MAX_METRIC];

typedef struct
{ 
  void (*func)(int*); 
} metric_diag_TABLE; 
metric_diag_TABLE metric_diag[FT_MAX_METRIC];

typedef struct
{ 
  void (*func)(int*,void*); 
} metric_comp_TABLE; 
metric_comp_TABLE metric_comp[FT_MAX_METRIC];


typedef struct
{ 
  void (*func)(int*,int*, int*, void*,void*); 
} metric_write_TABLE; 
metric_write_TABLE metric_write[FT_MAX_METRIC];

typedef struct
{ 
  void (*func)(int*,int*); 
} metric_wrst_TABLE; 
metric_wrst_TABLE metric_wrst[FT_MAX_METRIC];

typedef struct
{ 
  void (*func)(int*); 
} metric_rrst_TABLE; 
metric_rrst_TABLE metric_rrst[FT_MAX_METRIC];


typedef struct
{ 
  void (*func)(void*); 
} metric_reset_TABLE; 
metric_reset_TABLE metric_reset[FT_MAX_METRIC];


//BOP
// !ROUTINE: registermetricinit
// \label{registermetricinit}
//
// !INTERFACE:
void FTN(registermetricinit)(int *i,void (*func)(void*, void*, void*))
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
  ft_check_index(*i, FT_MAX_METRIC, "registermetricinit");
  metric_ini[*i].func = func; 
}
//BOP
// !ROUTINE: initmetric
// \label{initmetric}
//
// !INTERFACE:
void FTN(initmetric)(int *i, void *nlevs, void *stats, void *metric)
// !DESCRIPTION: 
//  Calls the routine from the registry to read the
//  the runtime domain specifics
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//
//EOP
{ 
  if(metric_ini[*i].func==NULL) {
    printf("****************Error****************************\n"); 
    printf("init routine for metric %d is not defined\n",*i); 
    printf("program will seg fault.....\n"); 
    printf("****************Error****************************\n"); 
  }	
  metric_ini[*i].func(nlevs,stats,metric); 
}

//BOP
// !ROUTINE: registermetricdiagnose
// \label{registermetricdiagnose}
//
// !INTERFACE:
void FTN(registermetricdiagnose)(int *i,void (*func)(int*))
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
  ft_check_index(*i, FT_MAX_METRIC, "registermetricdiagnose");
  metric_diag[*i].func = func; 
}
//BOP
// !ROUTINE: diagnosemetric
// \label{diagnosemetric}
//
// !INTERFACE:
void FTN(diagnosemetric)(int *i, int *pass)
// !DESCRIPTION: 
//  Calls the routine from the registry to read the
//  the runtime domain specifics
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//
//EOP
{ 
  if(metric_diag[*i].func==NULL) {
    printf("****************Error****************************\n"); 
    printf("diagnose routine for metric %d is not defined\n",*i); 
    printf("program will seg fault.....\n"); 
    printf("****************Error****************************\n"); 
  }	
  metric_diag[*i].func(pass); 
}


//BOP
// !ROUTINE: registermetriccompute
// \label{registermetriccompute}
//
// !INTERFACE:
void FTN(registermetriccompute)(int *i,void (*func)(int*,void*))
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
  ft_check_index(*i, FT_MAX_METRIC, "registermetriccompute");
  metric_comp[*i].func = func; 
}
//BOP
// !ROUTINE: computemetric
// \label{computemetric}
//
// !INTERFACE:
void FTN(computemetric)(int *i, int *pass, void *alarm)
// !DESCRIPTION: 
//  Calls the routine from the registry to read the
//  the runtime domain specifics
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//
//EOP
{ 
  if(metric_comp[*i].func==NULL) {
    printf("****************Error****************************\n"); 
    printf("compute routine for metric %d is not defined\n",*i); 
    printf("program will seg fault.....\n"); 
    printf("****************Error****************************\n"); 
  }	
  metric_comp[*i].func(pass,alarm); 
}

//BOP
// !ROUTINE: registermetricwriterestart
// \label{registermetricwriterestart}
//
// !INTERFACE:
void FTN(registermetricwriterestart)(int *i,void (*func)(int*, int*))
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
  ft_check_index(*i, FT_MAX_METRIC, "registermetricwriterestart");
  metric_wrst[*i].func = func; 
}
//BOP
// !ROUTINE: writemetricrestart
// \label{writemetricrestart}
//
// !INTERFACE:
void FTN(writemetricrestart)(int *i, int *ftn, int *pass)
// !DESCRIPTION: 
//  Calls the routine from the registry to read the
//  the runtime domain specifics
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//
//EOP
{ 
  if(metric_wrst[*i].func==NULL) {
    printf("****************Error****************************\n"); 
    printf("write restart routine for metric %d is not defined\n",*i); 
    printf("program will seg fault.....\n"); 
    printf("****************Error****************************\n"); 
  }	
  metric_wrst[*i].func(ftn,pass); 
}

//BOP
// !ROUTINE: registermetricreadrestart
// \label{registermetricreadrestart}
//
// !INTERFACE:
void FTN(registermetricreadrestart)(int *i,void (*func)(int*))
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
  ft_check_index(*i, FT_MAX_METRIC, "registermetricreadrestart");
  metric_rrst[*i].func = func; 
}
//BOP
// !ROUTINE: readmetricrestart
// \label{readmetricrestart}
//
// !INTERFACE:
void FTN(readmetricrestart)(int *i, int *ftn)
// !DESCRIPTION: 
//  Calls the routine from the registry to read the
//  the runtime domain specifics
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//
//EOP
{ 
  if(metric_rrst[*i].func==NULL) {
    printf("****************Error****************************\n"); 
    printf("read restart routine for metric %d is not defined\n",*i); 
    printf("program will seg fault.....\n"); 
    printf("****************Error****************************\n"); 
  }	
  metric_rrst[*i].func(ftn); 
}

//BOP
// !ROUTINE: registermetricwriteentry
// \label{registermetricwriteentry}
//
// !INTERFACE:
void FTN(registermetricwriteentry)(int *i,void (*func)(int*,int*, int*, void*,void*))
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
  ft_check_index(*i, FT_MAX_METRIC, "registermetricwriteentry");
  metric_write[*i].func = func; 
}
//BOP
// !ROUTINE: writemetricentry
// \label{writemetricentry}
//
// !INTERFACE:
void FTN(writemetricentry)(int *i, int *pass, int *final, int *vlevels,void *stats, void *obs)
// !DESCRIPTION: 
//  Calls the routine from the registry to read the
//  the runtime domain specifics
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//
//EOP
{ 
  if(metric_write[*i].func==NULL) {
    printf("****************Error****************************\n"); 
    printf("writeentry routine for metric %d is not defined\n",*i); 
    printf("program will seg fault.....\n"); 
    printf("****************Error****************************\n"); 
  }	
  metric_write[*i].func(pass,final,vlevels,stats,obs); 
}


//BOP
// !ROUTINE: registermetricreset
// \label{registermetricreset}
//
// !INTERFACE:
void FTN(registermetricreset)(int *i,void (*func)(void*))
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
  ft_check_index(*i, FT_MAX_METRIC, "registermetricreset");
  metric_reset[*i].func = func; 
}
//BOP
// !ROUTINE: resetmetric
// \label{resetmetric}
//
// !INTERFACE:
void FTN(resetmetric)(int *i, void *flag)
// !DESCRIPTION: 
//  Calls the routine from the registry to read the
//  the runtime domain specifics
// 
// The arguments are: 
// \begin{description}
// \item[i]
//  index of the domain
//  \end{description}
//
//EOP
{ 
  if(metric_reset[*i].func==NULL) {
    printf("****************Error****************************\n"); 
    printf("reset routine for metric %d is not defined\n",*i); 
    printf("program will seg fault.....\n"); 
    printf("****************Error****************************\n"); 
  }	
  metric_reset[*i].func(flag); 
}






