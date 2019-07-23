#include <stdio.h>
#include <stdlib.h>
#include <vic411_vicNl.h>

static char vcid[] = "$Id: vic411_nrerror.c,v 3.2 1999/08/23 23:59:06 vicadmin Exp $";

void vic411_nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	void _exit();

	fprintf(stderr,"Model run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	_exit(1);
}
