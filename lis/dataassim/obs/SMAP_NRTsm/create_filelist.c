//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.4
//
// Copyright (c) 2022 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------

/* Standard C headers */
#include <stdio.h>

/* POSIX header */
#include <glob.h>

/* Name mangling macro header */
#include "ftn.h"

/* Function prototypes */
int FTN(create_filelist) (const char *regexp,
			  const char *filelist_name);
//BOP
//
// !ROUTINE: create_filelist
// \label{create_filelist}
//
// !REVISION HISTORY:
//  11 Jan 2019: Eric Kemp; Original implementation.
//
// !INTERFACE:
int FTN(create_filelist) (const char *regexp,
			  const char *filelist_name)

// !DESCRIPTION:
//  Function for searching for files specified by an input regular expression,
//  and then writing the names of those files to a text output file.  The 
//  POSIX C function glob is used under the hood to process the regular
//  expression and create a list of the existing files.
//
//  This is intended for use by the SMAP reader as an alternative to calling
//  'ls' via the 'system' call.  'System' does not work cleanly with all
//  MPI implementations, and 'ls' returns errors to stderr if it cannot find
//  any files; the latter can, in some cases, lead to LIS exhausting memory.
//EOP

{
    /* Local variables */
    glob_t globbuf; 
    FILE * outfile = NULL;
    int ierr;
    int rc;
    int i;

    rc = 0;

    /* Use the POSIX glob function to assemble list of files matching the
       provided regular expression. */
    globbuf.gl_offs = 0;
    ierr = glob(regexp, 0, NULL, &globbuf);

    if (ierr == GLOB_NOMATCH) {
	/* No files were found by glob.  So just create an empty outfile. */
	outfile = fopen(filelist_name,"w");
	fclose(outfile);
    } else if (ierr == 0) {
	/* No error was returned by glob, and at least one file was found.
	   So write the filenames to an output file. */
	outfile = fopen(filelist_name,"w");
	for (i = 0; i < globbuf.gl_pathc; i++) {
	    fprintf(outfile,"%s\n", globbuf.gl_pathv[i]);
	}
	fclose(outfile);
    } else {
	/* An error was reported by glob. We will simply make an empty
	   outfile and return an error to the caller. */
	rc = 1;
	outfile = fopen(filelist_name,"w");
	fclose(outfile);
    }

    /* Clean up globbuf before returning. */
    globfree(&globbuf);    
    return rc;
}
