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
// !MODULE: LIS_alb_FTable
//
//
// !DESCRIPTION:
//  Function table registries for storing the interface
//  implementations for managing different sources of
//  albedo parameter data
//EOP
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ftn_drv.h"

struct albsetnode {
    char *name;
    void (*func)(int*);
    struct albsetnode* next;
};
struct albsetnode* albset_table = NULL;

struct albreadnode {
    char *name;
    void (*func)(int*, float*, float*, float*, float*);
    struct albreadnode* next;
};
struct albreadnode* albread_table = NULL;

//BOP
// !ROUTINE: registeralbedosetup
// \label{registeralbedosetup}
//
// !INTERFACE:
void FTN(registeralbedosetup)(char *j,
                              void (*func)(int*),
                              int len)
// !DESCRIPTION:
// Makes an entry in the registry for the routine to
// setup albedo data reading routines
//
// The arguments are:
// \begin{description}
// \item[i]
//  index of the domain
//  \item[j]
//  index of the albedo data source
//  \end{description}
//EOP
{
    int len1;
    struct albsetnode* current;
    struct albsetnode* pnode;
    // create node

    len1 = len + 1; // ensure that there is space for terminating null
    pnode = (struct albsetnode*) malloc(sizeof(struct albsetnode));
    pnode->name = (char*) calloc(len1, sizeof(char));
    strncpy(pnode->name, j, len);
    pnode->func = func;
    pnode->next = NULL;

    if (albset_table == NULL) {
        albset_table = pnode;
    } else {
        current = albset_table;
        while (current->next != NULL) {
            current = current->next;
        }
        current->next = pnode;
    }
}

//BOP
// !ROUTINE: albedosetup
// \label{albedosetup}
//
// !INTERFACE:
void FTN(albedosetup)(char *j, int *n, int len)
//
// !DESCRIPTION:
// Invokes the routine from the registry to
// setup albedo data reading
//
// The arguments are:
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//  index of the albedo data source
//  \end{description}
//EOP
{
    struct albsetnode* current;

    current = albset_table;
    while (strcmp(current->name, j) != 0) {
        current = current->next;

        if (current == NULL) {
            printf("****************Error****************************\n");
            printf("albedo setup routine for %s is not defined\n",j);
            printf("program will seg fault.....\n");
            printf("****************Error****************************\n");
        }
    }
    current->func(n);
}

//BOP
// !ROUTINE: registerreadalbedo
// \label{registerreadalbedo}
//
// !INTERFACE:
void FTN(registerreadalbedo)(char *j,
                            void (*func)(int*, float*, float*, float*, float*),
                            int len)
// !DESCRIPTION:
// Creates an entry in the registry for the routine to
// read albedo climatology data
//
// The arguments are:
// \begin{description}
//  \item[j]
//   index of the albedo source
//  \end{description}
//EOP
{
    int len1;
    struct albreadnode* current;
    struct albreadnode* pnode;
    // create node

    len1 = len + 1; // ensure that there is space for terminating null
    pnode = (struct albreadnode*) malloc(sizeof(struct albreadnode));
    pnode->name = (char*) calloc(len1, sizeof(char));
    strncpy(pnode->name, j, len);
    pnode->func = func;
    pnode->next = NULL;

    if (albread_table == NULL) {
        albread_table = pnode;
    } else {
        current = albread_table;
        while (current->next != NULL) {
            current = current->next;
        }
        current->next = pnode;
    }
}

//BOP
// !ROUTINE: readalbedo
// \label{readalbedo}
//
// !INTERFACE:
//EMK: Changed time1 and time2 arguments
//to floats with relative weights, consistant with underlying
//read_ALMIPIIalbedo Fortran routine.
void FTN(readalbedo)(char *j, int *n, float *wt1, float *wt2,
                     float *array1, float *array2, int len)
//
// !DESCRIPTION:
// Invokes the routines from the registry for
// reading climatology albedo files
//
// The arguments are:
// \begin{description}
//  \item[j]
//   index of the albedo source
//  \item[n]
//   index of the nest
//  \item[wt1]
//   weight of first time level
//  \item[wt2]
//   weight of second time level
//  \item[array1]
//   pointer to first array containing the albedo data.
//  \item[array2]
//   pointer to second array containing the albedo data.

//  \end{description}
//EOP
{
    struct albreadnode* current;

    current = albread_table;
    while (strcmp(current->name, j) != 0) {
        current = current->next;

        if (current == NULL) {
            printf("****************Error****************************\n");
            printf("readalbedo routine for %s is not defined\n",j);
            printf("program will seg fault.....\n");
            printf("****************Error****************************\n");
        }
    }
    current->func(n, wt1, wt2, array1, array2);
}
