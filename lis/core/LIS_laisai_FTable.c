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
// !MODULE: LIS_laisai_FTable
//
// !DESCRIPTION:
//  Function table registries for storing the interface
//  implementations for managing different sources of
//  LAI and SAI data
//EOP
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include<string.h>

#include "ftn_drv.h"
struct laisetnode {
    char *name;
    void (*func)(int*);
    struct laisetnode* next;
};
struct laisetnode* laiset_table = NULL;

struct saisetnode {
    char *name;
    void (*func)(int*);
    struct saisetnode* next;
};
struct saisetnode* saiset_table = NULL;

struct laireadnode {
    char *name;
    void (*func)(int*, float*, float*, float*, float*);
    struct laireadnode* next;
};
struct laireadnode* lairead_table = NULL;

struct saireadnode {
    char *name;
    void (*func)(int*, void*, float*);
    struct saireadnode* next;
};
struct saireadnode* sairead_table = NULL;
//BOP
// !ROUTINE: registerlaisetup
// \label{registerlaisetup}
//
// !INTERFACE:
void FTN(registerlaisetup)(char *j, void (*func)(int*), int len)
// !DESCRIPTION:
// Makes an entry in the registry for the routine to
// setup lai data reading routines
//
// The arguments are:
// \begin{description}
// \item[i]
//  index of the domain
//  \item[j]
//  index of the lai data source
//  \end{description}
//EOP
{
    int len1;
    struct laisetnode* current;
    struct laisetnode* pnode;
    // create node

    len1 = len + 1; // ensure that there is space for terminating null
    pnode = (struct laisetnode*) malloc(sizeof(struct laisetnode));
    pnode->name = (char*) calloc(len1, sizeof(char));
    strncpy(pnode->name, j, len);
    pnode->func = func;
    pnode->next = NULL;

    if (laiset_table == NULL) {
        laiset_table = pnode;
    } else {
        current = laiset_table;
        while (current->next != NULL) {
            current = current->next;
        }
        current->next = pnode;
    }
}

//BOP
// !ROUTINE: laisetup
// \label{laisetup}
//
// !INTERFACE:
void FTN(laisetup)(char *j, int *n, int len)
//
// !DESCRIPTION:
// Invokes the routine from the registry to
// setup lai data reading
//
// The arguments are:
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//  index of the lai data source
//  \end{description}
//EOP
{
    struct laisetnode* current;

    current = laiset_table;
    while (strcmp(current->name, j) != 0) {
        current = current->next;
        if(current == NULL) {
            printf("****************Error****************************\n");
            printf("laisetup routine for  %s is not defined\n",j);
            printf("program will seg fault.....\n");
            printf("****************Error****************************\n");
        }
    }
    current->func(n);
}

//BOP
// !ROUTINE: registersaisetup
// \label{registersaisetup}
//
// !INTERFACE:
void FTN(registersaisetup)(char *j, void (*func)(int*), int len)
// !DESCRIPTION:
// Makes an entry in the registry for the routine to
// setup sai data reading routines
//
// The arguments are:
// \begin{description}
// \item[i]
//  index of the domain
//  \item[j]
//  index of the sai data source
//  \end{description}
//EOP
{
    int len1;
    struct saisetnode* current;
    struct saisetnode* pnode;
    // create node

    len1 = len + 1; // ensure that there is space for terminating null
    pnode = (struct saisetnode*) malloc(sizeof(struct saisetnode));
    pnode->name = (char*) calloc(len1, sizeof(char));
    strncpy(pnode->name, j, len);
    pnode->func = func;
    pnode->next = NULL;

    if (saiset_table == NULL) {
        saiset_table = pnode;
    } else {
        current = saiset_table;
        while (current->next != NULL) {
            current = current->next;
        }
        current->next = pnode;
    }
}

//BOP
// !ROUTINE: saisetup
// \label{saisetup}
//
// !INTERFACE:
void FTN(saisetup)(char *j, int *n, int len)
//
// !DESCRIPTION:
// Invokes the routine from the registry to
// setup sai data reading
//
// The arguments are:
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//  index of the sai data source
//  \end{description}
//EOP
{
    struct saisetnode* current;

    current = saiset_table;
    while (strcmp(current->name, j) != 0) {
        current = current->next;
        if (current == NULL) {
            printf("****************Error****************************\n");
            printf("saisetup routine for  %s is not defined\n",j);
            printf("program will seg fault.....\n");
            printf("****************Error****************************\n");
        }
    }
    current->func(n);
}


//BOP
// !ROUTINE: registerreadlai
// \label{registerreadlai}
//
// !INTERFACE:
void FTN(registerreadlai)(char *j,
                          void (*func)(int*, float*, float*, float*, float*),
                          int len)
//
// !DESCRIPTION:
//  Creates an entry in the registry for the routines to
//  read LAI data
//
//  The arguments are:
//  \begin{description}
//   \item[j]
//    index of the LAI source
//   \end{description}
//
//EOP
{
    int len1;
    struct laireadnode* current;
    struct laireadnode* pnode;
    // create node

    len1 = len + 1; // ensure that there is space for terminating null
    pnode = (struct laireadnode*) malloc(sizeof(struct laireadnode));
    pnode->name = (char*) calloc(len1, sizeof(char));
    strncpy(pnode->name, j, len);
    pnode->func = func;
    pnode->next = NULL;

    if (lairead_table == NULL) {
        lairead_table = pnode;
    } else {
        current = lairead_table;
        while (current->next != NULL) {
            current = current->next;
        }
        current->next = pnode;
    }
}

//BOP
// !ROUTINE: readlai
// \label{readlai}
//
// !INTERFACE:
void FTN(readlai)(char *j, int *n, float *wt1, float *wt2,
                  float *array1, float *array2, int len)
//
// !DESCRIPTION:
// Invokes the routine from the registry to read LAI data.
//
//  The arguments are:
//  \begin{description}
//   \item[n]
//    index of the nest
//   \item[j]
//    index of the LAI source
//  \item[time]
//   time file search begins
//  \item[array]
//  pointer to the LAI data
//  \item[time]
//   time of found file
//   \end{description}
//EOP
{
    struct laireadnode* current;

    current = lairead_table;
    while (strcmp(current->name, j) != 0) {
        current = current->next;
        if (current == NULL) {
            printf("****************Error****************************\n");
            printf("readlai routine for  %s is not defined\n",j);
            printf("program will seg fault.....\n");
            printf("****************Error****************************\n");
        }
    }
    current->func(n, wt1, wt2, array1, array2);
}
//BOP
// !ROUTINE: registerreadsai
// \label{registerreadsai}
//
// !INTERFACE:
void FTN(registerreadsai)(char *j, void (*func)(int*, void*, float*), int len)
// !DESCRIPTION:
//  Creates an entry in the registry for the routines to
//  read SAI data
//
//  The arguments are:
//  \begin{description}
//   \item[j]
//    index of the SAI source
//   \end{description}
//EOP
{
    int len1;
    struct saireadnode* current;
    struct saireadnode* pnode;
    // create node

    len1 = len + 1; // ensure that there is space for terminating null
    pnode = (struct saireadnode*) malloc(sizeof(struct saireadnode));
    pnode->name = (char*) calloc(len1, sizeof(char));
    strncpy(pnode->name, j, len);
    pnode->func = func;
    pnode->next = NULL;

    if (sairead_table == NULL){
        sairead_table = pnode;
    } else {
        current = sairead_table;
        while (current->next != NULL ){
            current = current->next;
        }
        current->next = pnode;
    }
}
//BOP
// !ROUTINE: readsai
// \label{readsai}
//
// !INTERFACE:
void FTN(readsai)(char *j, int *n, void *time, float* array, int len)
//
// !DESCRIPTION:
// Invokes the routine from the registry to read SAI data.
//
//  The arguments are:
//  \begin{description}
//   \item[n]
//    index of nest
//   \item[j]
//    index of the SAI source
//  \item[time]
//    time
//  \item[array]
//  pointer to the SAI data
//   \end{description}
//EOP
{
    struct saireadnode* current;

    current = sairead_table;
    while (strcmp(current->name, j) != 0) {
        current = current->next;
        if (current == NULL) {
            printf("****************Error****************************\n");
            printf("readsai routing for %s is not defined\n",j);
            printf("program will seg fault.....\n");
            printf("****************Error****************************\n");
        }
    }
  current->func(n, time, array);
}


