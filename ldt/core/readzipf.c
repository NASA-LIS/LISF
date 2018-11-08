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

/* Unzip a ".gz", ".z" or ".Z" file, and read a specified number of bytes */
/* into array                                                             */
/* can not use fseek() to jump around for a pipe */ 

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <zlib.h>
#include "ftn_drv.h"


int FTN(readzipf)(char *file, char *array, int *recl)
{
  
  int nrd; 
  gzFile fp; 
  
  fp = gzopen(file, "rb"); 
  if (fp == NULL) {
     fprintf(stderr, "Zip file openning error: %s\n", file); 
     return(-1); 
  }
  nrd = gzread(fp, array, *recl); 
  if (nrd != *recl) { 
     fprintf(stderr, "Zip file reading error, recl=%d nread=%d\n", *recl, nrd); 
     gzclose(fp); 
     return(-1); 
  }
  gzclose(fp); 
  return(nrd);  

}



// 3/12/12: Old way of reading. It seems "popen()" is causing problems for parallel 
// runs, as reported by Sujay.  
/* 
  char command[4096]; 
  FILE *fp; 
  int nrd; 

  snprintf(command, sizeof(command), "/usr/bin/gunzip -c %s", file); 
  fp=popen(command, "r"); 

  nrd=fread(array, 1, *recl, fp); 
  if (nrd != *recl) { 
     fprintf(stderr, "Zip file reading error, recl=%d nread=%d", *recl, nrd); 
     return(-1); 
  }
  fclose(fp);    
  return(nrd);  
}

*/ 

  
 
  
