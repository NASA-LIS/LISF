//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ftn_drv.h"
//BOP
//
// !ROUTINE: LIS_memory_managementC
// \label{LIS_memory_managementC}
//
// !DESCRIPTION:
//   This file provides methods for dynamic memory management in LIS 
//   for codes written in the C language.  
//EOP

void lis_log_msgC(char *);

//BOP
// !ROUTINE: lis_calloc
// \label{lis_calloc}
// 
// !DESCRIPTION: 
// 
// Generic call to the C calloc function used in LIS.
//
// !REVISION HISTORY:
// Mar 2004, James Geiger; Initial Specification
// !INTERFACE:
void * lis_calloc(size_t n, size_t size, char * caller)
//EOP
{
   void * ptr;
   char * msg;

   int count;
   int len;

   ptr = calloc(n,size);

   if ( ptr == NULL )
   {
      len = 31 + strlen(caller) + 1;
      msg = (char *) malloc(len);
      count = sprintf(msg,"ERR: %s -- Cannot allocate memory",caller);

      if ( count != len-1 )
      {
         lis_log_msgC("ERR: lis_calloc -- "
                      "Cannot allocate memory for msg string");
      }
      
      lis_log_msgC(msg);
      free(msg);
   }

   return ptr;
}

//BOP
// !ROUTINE: lis_malloc
// \label{lis_malloc}
// 
// !DESCRIPTION: 
// 
// Generic call to the C malloc function used in LIS.
//
// !REVISION HISTORY:
// Mar 2004, James Geiger; Initial Specification
// !INTERFACE:
void * lis_malloc(size_t size, char * caller)
//EOP
{
   void * ptr;
   char * msg;

   int count;
   int len;

   ptr = malloc(size);

   if ( ptr == NULL )
   {
      len = 31 + strlen(caller) + 1;
      msg = (char *) malloc(len);
      count = sprintf(msg,"ERR: %s -- Cannot allocate memory",caller);

      if ( count != len-1 )
      {
         lis_log_msgC("ERR: lis_malloc -- "
                      "Cannot allocate memory for msg string");
      }
      
      lis_log_msgC(msg);
      free(msg);
   }

   return ptr;
}
