//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------

/*****************************************************************************/
/* NASA/GSFC, Software Systems Support Office, Code 610.3                    */
/*****************************************************************************/
/*                                                                           */
/* FILE: LVT_create_subdirs                                                  */
/*                                                                           */
/* AUTHOR:                                                                   */
/* Eric Kemp, NASA SSSO/SSAI                                                 */
/*                                                                           */
/* DESCRIPTION:                                                              */
/* Function to create hierarchical directories from a Fortran string using   */
/* the POSIX library function 'mkdir'.                                       */
/*                                                                           */
/*****************************************************************************/

/* ANSI C Standard Headers */
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* POSIX header */
#include <sys/stat.h>

/* Name mangling macro header */
#include "ftn.h"

/* Function prototypes */
int FTN(lvt_create_subdirs) (int *length, char *path);

/*****************************************************************************/
/*                                                                           */
/* ROUTINE: LVT_create_subdirs                                               */
/*                                                                           */
/* DESCRIPTION: Receives Fortran string and string length and passes to C    */
/* POSIX library function 'mkdir'. Returns status from this call.            */
/*                                                                           */
/*****************************************************************************/

int FTN(lvt_create_subdirs) (int  *length,
			     char *path) {

  /* Local variables */
  char *work;
  char *ptr;
  int i;
  int length_copy;
  int status;

  /* Sanity check length of string */
  if (*length == 0) {
	puts("ERROR, zero length string sent to LVT_create_subdirs!");
	status = 1;
	return status;
  }

  /* Copy directory path to internal string. Since Fortran doesn't use null 
     terminators, we need to copy this manually and append the terminator. */
  length_copy = *length;
  work = malloc((length_copy+1)*sizeof(char));
  if (work == NULL) {
    perror("ERROR calling malloc from create_subdirs!");
    status = 1;
    return status;
  }
  for (i = 0; i < length_copy; ++i) {
    work[i] = path[i];
  }
  work[length_copy] = '\0';

  /* Using a separate pointer: loop through the character string, temporarily
     change subdirectory tokens to null terminators, and use the mkdir 
     function to create each subdirectory path. If subdirectory already 
     exists, treat this as success. If failure occurs, report. */
  status = 0;
  ptr = work;
  for (;;) {
    if (*ptr == '\0') {     /* At end of the directory string */
      status = mkdir(work, S_IRUSR | S_IWUSR | S_IXUSR);
      if (status != 0) {
        if (errno == EEXIST) {
          status = 0; /* Not really a problem if directory already exists */
        } else {
          perror("ERROR calling mkdir from create_subdirs!");
          status = 1;
        }
      }
      break;

    } else if (*ptr == '/') { /* Subdirectory token */
      *ptr = '\0'; /* Temporarily mark this as end of string */

	   /* If an absolute path is specified, the full string will start with */
	   /* a "/". Thus, the mkdir call below may receive an empty string     */
	   /* since we temporarily replaced the "/" with a null character.      */
	   /* We'll skip that case to prevent fatal problems.                   */

	   if (strlen(work) > 0) { /* NEW LINE */
         status = mkdir(work, S_IRUSR | S_IWUSR | S_IXUSR);
         if (status != 0) {
           if (errno == EEXIST) {
             status = 0; /* Not really a problem if directory already exists */
           } else {
             perror("ERROR calling mkdir!");
             status = 1;
             break;
           }
         }
	   } /* NEW LINE */
      *ptr = '/'; /* Restore subdirectory token */
    }
    ++ptr;
  }

  /* Clean up and exit with status */
  free(work);
  return status;
}
