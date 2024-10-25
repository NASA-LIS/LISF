//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------

/* Unzip a ".gz", ".z" or ".Z" file, and read a specified number of bytes */
/* into array                                                             */
/* can not use fseek() to jump around for a pipe */ 

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "ftn_drv.h"


int FTN(cm_readzipf)(char *file, char *array, int *recl)
{
  
// CMORPH's .Z files can not be read with zlib. So use pipe/gzip. 
  char command[4096]; 
  FILE *fp; 
  int nrd; 

  snprintf(command, sizeof(command), "/usr/bin/gunzip -c %s", file); 
  fp=popen(command, "r"); 
  if (fp == NULL) { 
     fprintf(stderr, "Zip file open error: %s\n", file); 
     return(-1); 
  }
    
  nrd=fread(array, 1, *recl, fp); 
  if (nrd != *recl) { 
     fprintf(stderr, "Zip file reading error, recl=%d nread=%d\n", *recl, nrd); 
     return(-1); 
  }
  fclose(fp);    
  return(nrd);  
}


  
 
  
