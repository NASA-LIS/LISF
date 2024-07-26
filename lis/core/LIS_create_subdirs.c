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
/* FILE: LIS_create_subdirs                                                  */
/*                                                                           */
/* AUTHOR:                                                                   */
/* Eric Kemp, NASA SSSO/SSAI                                                 */
/*                                                                           */
/* DESCRIPTION:                                                              */
/* Function to create hierarchical directories from a Fortran string using   */
/* POSIX library calls.                                                      */
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
int FTN(lis_create_subdirs) (int *length, char *path);

/*****************************************************************************/
/*                                                                           */
/* ROUTINE: LIS_create_subdirs                                               */
/*                                                                           */
/* DESCRIPTION: Receives Fortran string and string length and passes to C    */
/* POSIX library calls to create full directory path. Returns error status.  */
/*                                                                           */
/*****************************************************************************/

int FTN(lis_create_subdirs) (int  *length,
			     char *path) 
{
    
    /* Local variables */
    char *work;
    char *ptr;
    int i;
    int length_copy;
    int stat_rc;
    int isdir_rc;
    int mkdir_rc;
    int error;
    char *error_string;
    int end_of_line;
    struct stat sb;
    
    /* Sanity check length of string */
    if (*length == 0) {
        puts("ERROR, zero length string sent to LIS_create_subdirs!");
        error = 1;
        return error;
    }
    
    /* Copy directory path to internal string. Since Fortran doesn't use null 
       terminators, we need to copy this manually and append the terminator. */
    length_copy = *length;
    work = NULL;
    work = malloc((length_copy+1)*sizeof(char));
    if (work == NULL) {
	perror("ERROR calling malloc from create_subdirs!");
	error = 1;
        return error;
    }
    for (i = 0; i < length_copy; ++i) {
        work[i] = path[i];
    }
    work[length_copy] = '\0';
    
    /* Using a separate pointer: loop through the character string and look
       for subdirectory tokens. */
    error = 0;
    end_of_line = 0;
    ptr = work;
    for (;;) {

        /* See if we are at the end of the full directory path. */
        if (*ptr == '\0') {    
            end_of_line = 1;
        } 

        /* See if we are at the end of higher-level directory path */
        if (*ptr == '/') { /* Subdirectory token */
            *ptr = '\0'; /* Temporarily mark this as end of string */
        }

        if (*ptr == '\0' && strlen(work) > 0) {

            stat_rc = 0;
            isdir_rc = 0;
            mkdir_rc = 0;

            /* Get file status of current work path fragment. Uses POSIX
               stat call */
            stat_rc = stat(work, &sb);

            /* Handle errors */
            if (stat_rc == -1) {
                /* ENOENT is returned via errno if component of path doesn't
                   exist.  This isn't an error for our purposes, since we can
                   attempt to create a directory below. */
                if (errno != ENOENT) {
                    /* Attempt to report error message */
                    error_string = NULL;
                    error_string = malloc((length_copy+1+21)*sizeof(char));
                    if (error_string != NULL) {
                        sprintf(error_string,
                                "ERROR calling stat %s", work);
                        perror(error_string);
                        free(error_string);
                        error_string = NULL;			    
                    }
                    error = 1;
                }

            }

            if ( stat_rc == -1 && error == 0) {
                /* At this point, work path doesn't exist but no error
                   occurred.  We will attempt to create the directory. 
                   Set default permissions to rwx for user, group, and
                   other.  Permissions can be removed at run time using
                   the shell 'umask' command. */
                mkdir_rc = mkdir(work, 0777);

                /* Handle error */
                if (mkdir_rc == -1) {
                    /* EEXIST is returned via errno if the directory
                    * already exists. This isn't an error because
                    * we are attempting to create the directory.*/
                    if (errno != EEXIST) {
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
                        error = 1;
                    }
                }            
            }
        }

        /* Restore subdirectory token if necessary */
        if (*ptr == '\0' && !end_of_line) {
            *ptr = '/'; 
        }      

        /* Give up if we had an error */
        if (error) {
            break;
        }
	
        /* Go to next token, unless we are finished */
        if (end_of_line) {
            break;
        }
        ++ptr;
    }
        
    /* Clean up and return with error status */
    free(work);
    fflush(NULL);
    return error;
}

