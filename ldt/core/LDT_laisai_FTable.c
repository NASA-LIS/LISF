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
// !MODULE: LDT_laisai_FTable
//
//
// !DESCRIPTION:
//  Function table registries for storing the interface
//  implementations for managing different sources of
//  LAI/SAI data
//
//EOP
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ftn_drv.h"

struct laisetnode {
    char *name;
    void (*func)();
    struct laisetnode* next;
};
struct laisetnode* laiset_table = NULL;

struct lainode {
    char *name;
    //EMK Removed third argument, since not used by any Fortran routine.
    //void (*func)(int*, float*, float*);
    void (*func)(int*, float*);
    struct lainode* next;
};
struct lainode* lai_table = NULL;

struct sainode {
    char *name;
    void (*func)(int*, float*);
    struct sainode* next;
};
struct sainode* sai_table = NULL;

struct laiminnode {
    char *name;
    void (*func)(int*, float*);
    struct laiminnode* next;
};
struct laiminnode* laimin_table = NULL;

struct laimaxnode {
    char *name;
    void (*func)(int*, float*);
    struct laimaxnode* next;
};
struct laimaxnode* laimax_table = NULL;

//BOP
// !ROUTINE: registersetlaiattribs
// \label{registersetlaiattribs}
//
//
// !INTERFACE:
void FTN(registersetlaiattribs)(char *j, void (*func)(), int len)
// !DESCRIPTION:
//  Creates an entry in the registry for the routine to
//  read lai data
//
//  The arguments are:
//  \begin{description}
//   \item[j]
//    index of the lai source
//   \end{description}
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
// !ROUTINE: setlaiattribs
// \label{setlaiattribs}
//
// !DESCRIPTION:
// Invokes the routine from the registry for
// reading lai data.
//
//  The arguments are:
//  \begin{description}
//   \item[n]
//    index of the nest
//   \item[j]
//    index of the lai source
//   \item[array]
//    pointer to the lai data
//   \item[marray]
//    pointer to the mask data
//  \end{description}
//
// !INTERFACE:
void FTN(setlaiattribs)(char *j, int len)
//EOP
{

    struct laisetnode* current;

    current = laiset_table;
    while (strcmp(current->name, j) != 0) {
        current = current->next;

        if (current == NULL) {
            printf("****************Error****************************\n");
            printf("setLAIAttribs routine for source %s is not defined\n",j);
            printf("program will seg fault.....\n");
            printf("****************Error****************************\n");
        }
    }
    current->func();
}


//BOP
// !ROUTINE: registerreadlai
// \label{registerreadlai}
//
// !INTERFACE:
//void FTN(registerreadlai)(char *j,void (*func)(int*,float*,float*), int len)
void FTN(registerreadlai)(char *j, void (*func) (int*, float*), int len)
// !DESCRIPTION:
// Makes an entry in the registry for the routine to
// read lai data
//
// The arguments are:
// \begin{description}
// \item[i]
//  index of the domain
//  \item[j]
//  index of the LAI data source
//  \end{description}
//EOP
{
    int len1;
    struct lainode* current;
    struct lainode* pnode;
    // create node

    len1 = len + 1; // ensure that there is space for terminating null
    pnode = (struct lainode*) malloc(sizeof(struct lainode));
    pnode->name = (char*) calloc(len1, sizeof(char));
    strncpy(pnode->name, j, len);
    pnode->func = func;
    pnode->next = NULL;

    if (lai_table == NULL) {
        lai_table = pnode;
    } else {
        current = lai_table;
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
//EMK...Removed unnecessary arguments.
//void FTN(readlai)(char *j, int *n,float *array,float *marray,int len)
void FTN(readlai)(char *j, int *n, float *array, int len)

//
// !DESCRIPTION:
// Invokes the routine from the registry to
// reading lai data
//
// The arguments are:
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//  index of the LAI data source
//  \item[array]
//  pointer to the LAI data
//  \item[marray]
//  pointer to the mask data
//  \end{description}
//EOP
{
    struct lainode* current;

    current = lai_table;
    while (strcmp(current->name, j) != 0 ) {
        current = current->next;

        if (current == NULL) {
            printf("****************Error****************************\n");
            printf("LAI reading routine for source %s is not defined\n", j);
            printf("Please refer to LDT User's Guide for options.\n");
            printf("program will seg fault.....\n");
            printf("****************Error****************************\n");
        }
    }
    //EMK...Third argument not used by any Fortran routine.
    //current->func(n,array,marray);
    current->func(n, array);
}

//BOP
// !ROUTINE: registerreadsai
// \label{registerreadsai}
//
// !INTERFACE:
void FTN(registerreadsai)(char *j, void (*func)(int*, float*), int len)
// !DESCRIPTION:
// Makes an entry in the registry for the routine to
// read sai data
//
// The arguments are:
// \begin{description}
// \item[i]
//  index of the domain
//  \item[j]
//  index of the SAI data source
//  \end{description}
//EOP
{
    int len1;
    struct sainode* current;
    struct sainode* pnode;
    // create node

    len1 = len + 1; // ensure that there is space for terminating null
    pnode = (struct sainode*) malloc(sizeof(struct sainode));
    pnode->name = (char*) calloc(len1, sizeof(char));
    strncpy(pnode->name, j, len);
    pnode->func = func;
    pnode->next = NULL;

    if (sai_table == NULL) {
        sai_table = pnode;
    } else {
        current = sai_table;
        while (current->next != NULL) {
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
void FTN(readsai)(char *j, int *n, float *array, int len)
//
// !DESCRIPTION:
// Invokes the routine from the registry to
// reading sai data
//
// The arguments are:
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//  index of the SAI data source
//  \item[array]
//  pointer to the SAI data
//  \end{description}
//EOP
{

    struct sainode* current;

    current = sai_table;
    while (strcmp(current->name, j) != 0) {
        current = current->next;

        if (current == NULL) {
            printf("****************Error****************************\n");
            printf("SAI reading routine for source %s is not defined\n", j);
            printf("Please refer to LDT User's Guide for options.\n");
            printf("program will seg fault.....\n");
            printf("****************Error****************************\n");
        }
    }
    current->func(n, array);
}


//BOP
// !ROUTINE: registerreadlaimin
// \label{registerreadlaimin}
//
// !INTERFACE:
void FTN(registerreadlaimin)(char *j, void (*func)(int*, float*), int len)
// !DESCRIPTION:
// Makes an entry in the registry for the routine to
// read laimin data
//
// The arguments are:
// \begin{description}
// \item[i]
//  index of the domain
//  \item[j]
//  index of the greenness data source
//  \end{description}
//EOP
{
    int len1;
    struct laiminnode* current;
    struct laiminnode* pnode;
    // create node

    len1 = len + 1; // ensure that there is space for terminating null
    pnode = (struct laiminnode*) malloc(sizeof(struct laiminnode));
    pnode->name = (char*) calloc(len1, sizeof(char));
    strncpy(pnode->name, j, len);
    pnode->func = func;
    pnode->next = NULL;

    if (laimin_table == NULL) {
        laimin_table = pnode;
    } else {
        current = laimin_table;
        while (current->next != NULL) {
            current = current->next;
        }
        current->next = pnode;
    }
}

//BOP
// !ROUTINE: readlaimin
// \label{readlaimin}
//
// !INTERFACE:
void FTN(readlaimin)(char *j, int *n, float *array, int len)
//
// !DESCRIPTION:
// Invokes the routine from the registry to
// reading laimin data
//
// The arguments are:
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//  index of the greenness data source
//  \item[array]
//  pointer to the greenness data
//  \end{description}
//EOP
{
    struct laiminnode* current;

    current = laimin_table;
    while (strcmp(current->name, j) != 0) {
        current = current->next;
        if (current == NULL) {
            printf("****************Error****************************\n");
            printf("Min LAI reading routine for source %s is not defined\n",
                   j);
            printf("Please refer to LDT User's Guide for options.\n");
            printf("program will seg fault.....\n");
            printf("****************Error****************************\n");
        }
    }
    current->func(n, array);
}


//BOP
// !ROUTINE: registerreadlaimax
// \label{registerreadlaimax}
//
// !INTERFACE:
void FTN(registerreadlaimax)(char *j, void (*func)(int*, float*), int len)
// !DESCRIPTION:
// Makes an entry in the registry for the routine to
// read laimax data
//
// The arguments are:
// \begin{description}
// \item[i]
//  index of the domain
//  \item[j]
//  index of the greenness data source
//  \end{description}
//EOP
{
    int len1;
    struct laimaxnode* current;
    struct laimaxnode* pnode;
    // create node

    len1 = len + 1; // ensure that there is space for terminating null
    pnode = (struct laimaxnode*) malloc(sizeof(struct laimaxnode));
    pnode->name = (char*) calloc(len1, sizeof(char));
    strncpy(pnode->name, j, len);
    pnode->func = func;
    pnode->next = NULL;

    if (laimax_table == NULL) {
        laimax_table = pnode;
    } else {
        current = laimax_table;
        while (current->next != NULL) {
            current = current->next;
        }
        current->next = pnode;
    }
}

//BOP
// !ROUTINE: readlaimax
// \label{readlaimax}
//
// !INTERFACE:
void FTN(readlaimax)(char *j, int *n, float *array, int len)
//
// !DESCRIPTION:
// Invokes the routine from the registry to
// reading laimax data
//
// The arguments are:
// \begin{description}
//  \item[n]
//   index of the nest
//  \item[j]
//  index of the greenness data source
//  \item[array]
//  pointer to the greenness data
//  \end{description}
//EOP
{
    struct laimaxnode* current;

    current = laimax_table;
    while (strcmp(current->name, j) != 0) {
        current = current->next;

        if (current == NULL) {
            printf("****************Error****************************\n");
            printf("Max LAI reading routine for source %s is not defined\n",
                   j);
            printf("Please refer to LDT User's Guide for options.\n");
            printf("program will seg fault.....\n");
            printf("****************Error****************************\n");
        }
    }
    current->func(n, array);
}
