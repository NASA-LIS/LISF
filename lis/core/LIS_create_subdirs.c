//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center Land Information System (LIS) v7.2
//
// Copyright (c) 2015 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------

/*****************************************************************************/
/* NASA/GSFC, Software Systems Support Office, Code 610.3                    */
/*****************************************************************************/
/*                                                                           */
/* FILE: LIS_create_subdirs                                                  */
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
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* POSIX header */
#include <sys/stat.h>

/* Name mangling macro header */
#include "ftn.h"

/* Function prototypes */
int FTN(lis_create_subdirs) (int *length, char *path);

/*****************************************************************************/
/*                                                                           */
/* ROUTINE: LIS_create_subdirs                                               */
/*                                                                           */
/* DESCRIPTION: Receives Fortran string and string length and passes to C    */
/* POSIX library function 'mkdir'. Returns status from this call.            */
/*                                                                           */
/*****************************************************************************/

int FTN(lis_create_subdirs) (int  *length,
			     char *path) {
    
    /* Local variables */
    char *work;
    char *ptr;
    int i;
    int length_copy;
    int status;
    char *error_string;
    int end_of_line;
    struct stat sb;
    
    /* Sanity check length of string */
    if (*length == 0) {
	puts("ERROR, zero length string sent to LIS_create_subdirs!");
	status = 1;
	return status;
    }
    
    /* Copy directory path to internal string. Since Fortran doesn't use null 
       terminators, we need to copy this manually and append the terminator. */
    length_copy = *length;
    work = NULL;
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
    
    /* Using a separate pointer: loop through the character string and look
       for subdirectory tokens. */
    status = 0;
    end_of_line = 0;
    ptr = work;
    for (;;) {

	/* See if we are at the end of the full directory path. */
	if (*ptr == '\0') {    
	    end_of_line = 1;
	} 

	/* See if we are at the end of a super-directory path */
	if (*ptr == '/') { /* Subdirectory token */
	    *ptr = '\0'; /* Temporarily mark this as end of string */
	}

	if (*ptr == '\0' && strlen(work) > 0) {

	    /* See if the directory already exists */
	    if (stat(work, &sb) == 0 && S_ISDIR(sb.st_mode)) {
		status = 0;
	    } else {
		/* We must create the directory.  Set default permissions to 
		   rwx for user, group, and other.  Permissions can be removed
		   at run time using the shell 'umask' command. */
		status = mkdir(work, 0777);

		/* Handle error */
		if (status != 0) {
		    /* Gracefully handle possible race condition */
		    if (errno == EEXIST) {
			status = 0;
		    } else {
			/* Attempt to report error message */
			error_string = NULL;
			error_string = malloc((length_copy+1+21)*sizeof(char));
			if (error_string != NULL) {
			    sprintf(error_string,
				    "ERROR calling mkdir %s",work);
			    perror(error_string);
			    free(error_string);
			    error_string = NULL;			    
			}
			status = 1;
		    }
		}
	    }
	}

	/* Restore subdirectory token if necessary */
	if (*ptr == '\0' && !end_of_line) {
	    *ptr = '/'; 
	}      

	/* Give up if we had an error */
	if (status != 0) {
	    break;
	}
	
	/* Go to next token, unless we are finished */
	if (end_of_line) {
	    break;
	}
	++ptr;
    }
    
    /* Clean up and return with status */
    free(work);
    fflush(NULL);
    return status;
}

