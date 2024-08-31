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
// !MODULE: LDT_OSSEmaskdata_FTable
//
// !DESCRIPTION:
//  Function table registries for storing the interface
//  implementations for managing the operation of different
//  nature run data sources.
//
//EOP
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ftn_drv.h"

struct ossemaskdatasetupnode {
    char *name;
    void (*func)();
    struct ossemaskdatasetupnode* next;
};
struct ossemaskdatasetupnode* ossemaskdatasetup_table = NULL;

struct ossemasksourcenode {
    char *name;
    void (*func)(int*);
    struct ossemasksourcenode* next;
};
struct ossemasksourcenode* ossemasksourcenode_table = NULL;


//BOP
// !ROUTINE: registerossemasksourcesetup
// \label{registerossemasksourcesetup}
//
// !INTERFACE:
void FTN(registerossemasksourcesetup)(char *j, void (*func)(), int len)
//
// !DESCRIPTION:
// Makes an entry in the registry for the routine to
// read the runtime domain specifics
//
// The arguments are:
// \begin{description}
// \item[j]
//  name of the observation source
//  \end{description}
//EOP
{
    struct ossemaskdatasetupnode* current;
    struct ossemaskdatasetupnode* pnode;

    // create node
    pnode = (struct ossemaskdatasetupnode*)
        malloc(sizeof(struct ossemaskdatasetupnode));
    pnode->name = (char*) malloc(len*sizeof(char));
    strcpy(pnode->name, j);
    pnode->func = func;
    pnode->next = NULL;

    if (ossemaskdatasetup_table == NULL) {
        ossemaskdatasetup_table = pnode;
    } else {
        current = ossemaskdatasetup_table;
        while (current->next != NULL) {
            current = current->next;
        }
        current->next = pnode;
    }
}

//BOP
// !ROUTINE: setupossemasksource
// \label{setupossemasksource}
//
// !INTERFACE:
void FTN(setupossemasksource)(char *j, int len)
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
    struct ossemaskdatasetupnode* current;

    current = ossemaskdatasetup_table;
    while (strcmp(current->name, j) != 0) {
        current = current->next;
        if (current == NULL) {
            printf("****************Error****************************\n");
            printf("OSSE mask source setup routine for %s is not defined\n",j);
            printf("program will seg fault.....\n");
            printf("****************Error****************************\n");
        }
    }
  current->func();
}

//BOP
// !ROUTINE: registerreadossemasksource
// \label{registerreadossemasksource}
//
// !INTERFACE:
void FTN(registerreadossemasksource)(char *j, void (*func)(int*), int len)
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
    struct ossemasksourcenode* current;
    struct ossemasksourcenode* pnode;

    // create node
    pnode = (struct ossemasksourcenode*)
        malloc(sizeof(struct ossemasksourcenode));
    pnode->name = (char*) malloc(len*sizeof(char));
    strcpy(pnode->name, j);
    pnode->func = func;
    pnode->next = NULL;

    if (ossemasksourcenode_table == NULL) {
        ossemasksourcenode_table = pnode;
    } else {
        current = ossemasksourcenode_table;
        while (current->next != NULL) {
            current = current->next;
        }
        current->next = pnode;
    }
}

//BOP
// !ROUTINE: readossemasksource
// \label{readossemasksource}
//
// !INTERFACE:
void FTN(readossemasksource)(char *j, int *n, int len)
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
    struct ossemasksourcenode* current;

    current = ossemasksourcenode_table;
    while (strcmp(current->name, j) != 0) {
        current = current->next;
        if (current == NULL) {
            printf("****************Error****************************\n");
            printf("Read OSSE mask source routine for %s is not defined\n", j);
            printf("program will seg fault.....\n");
            printf("****************Error****************************\n");
        }
    }
  current->func(n);
}
