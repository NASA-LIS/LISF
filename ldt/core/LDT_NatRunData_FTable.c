//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center Land Data Toolkit (LDT) V7.0 BETA
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
// !MODULE: LDT_natrundata_FTable
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

struct naturerundatasetupnode {
    char *name;
    void (*func)();
    struct naturerundatasetupnode* next;
};
struct naturerundatasetupnode* naturerundatasetup_table = NULL;

struct naturerunsourcenode {
    char *name;
    void (*func)(int*);
    struct naturerunsourcenode* next;
};
struct naturerunsourcenode* naturerunsourcenode_table = NULL;


//BOP
// !ROUTINE: registernaturerunsourcesetup
// \label{registernaturerunsourcesetup}
//
// !INTERFACE:
void FTN(registernaturerunsourcesetup)(char *j, void (*func)(), int len)
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
    struct naturerundatasetupnode* current;
    struct naturerundatasetupnode* pnode;

    // create node
    pnode=(struct naturerundatasetupnode*)
        malloc(sizeof(struct naturerundatasetupnode));
    pnode->name = (char*) malloc(len*sizeof(char));
    strcpy(pnode->name, j);
    pnode->func = func;
    pnode->next = NULL;

    if (naturerundatasetup_table == NULL) {
        naturerundatasetup_table = pnode;
    } else {
        current = naturerundatasetup_table;
        while (current->next != NULL) {
            current = current->next;
        }
        current->next = pnode;
    }
}

//BOP
// !ROUTINE: setupnaturerunsource
// \label{setupnaturerunsource}
//
// !INTERFACE:
void FTN(setupnaturerunsource)(char *j, int len)
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
    struct naturerundatasetupnode* current;

    current = naturerundatasetup_table;
    while (strcmp(current->name, j) != 0 ) {
        current = current->next;
        if (current == NULL) {
            printf("****************Error****************************\n");
            printf("Observation setup routine for %s is not defined\n", j);
            printf("program will seg fault.....\n");
            printf("****************Error****************************\n");
        }
    }
    current->func();
}

//BOP
// !ROUTINE: registerreadnaturerunsource
// \label{registerreadnaturerunsource}
//
// !INTERFACE:
void FTN(registerreadnaturerunsource)(char *j, void (*func)(int*), int len)
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
    struct naturerunsourcenode* current;
    struct naturerunsourcenode* pnode;

    // create node
    pnode = (struct naturerunsourcenode*)
        malloc(sizeof(struct naturerunsourcenode));
    pnode->name = (char*) malloc(len*sizeof(char));
    strcpy(pnode->name, j);
    pnode->func = func;
    pnode->next = NULL;

    if (naturerunsourcenode_table == NULL) {
        naturerunsourcenode_table = pnode;
    } else {
        current = naturerunsourcenode_table;
        while (current->next != NULL) {
            current = current->next;
        }
        current->next = pnode;
    }
}

//BOP
// !ROUTINE: readnaturerunsource
// \label{readnaturerunsource}
//
// !INTERFACE:
void FTN(readnaturerunsource)(char *j, int *n, int len)
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
    struct naturerunsourcenode* current;

    current = naturerunsourcenode_table;
    while (strcmp(current->name, j) != 0) {
        current = current->next;
        if (current == NULL) {
            printf("****************Error****************************\n");
            printf("Read observation source routine for %s is not defined\n",
                   j);
            printf("program will seg fault.....\n");
            printf("****************Error****************************\n");
        }
    }
    current->func(n);
}
