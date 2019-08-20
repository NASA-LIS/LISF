/*
 * Purpose: allocate and free memory for the atmos data struct
 * Usage  : Part of VIC
 * Author : Bart Nijssen
 * E-mail : nijssen@u.washington.edu
 * Created: Fri Aug 27 18:22:42 1999
 * Last Changed: Tue Sep  2 15:11:02 2003 by Keith Cherkauer <cherkaue@u.washington.edu>
 * Notes  :
 */

/****************************************************************************/
/*			  PREPROCESSOR DIRECTIVES                           */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <vic411_vicNl.h>

static char vcid[] = "$Id: vic411_alloc_atmos.c,v 5.5.2.1 2007/01/11 22:31:28 vicadmin Exp $";

/****************************************************************************/
/*			       vic411_alloc_atmos()                                */
/****************************************************************************/
void vic411_alloc_atmos(int nrecs, vic411_atmos_data_struct **atmos)
/*******************************************************************
  vic411_alloc_atmos    

  Modifications:
  01-11-00 Fixed allocation bug                             KAC
  2006-Sep-23 Implemented flexible output configuration; removed
	      LDAS_OUTPUT and OPTIMIZE compile-time vic411_options.  TJB
  2006-Dec-20 All atmos_data arrays are always dynamically allocated now.	TJB

*******************************************************************/
{
  extern vic411_param_set_struct vic411_param_set;

  int i;

  *atmos = (vic411_atmos_data_struct *) calloc(nrecs, sizeof(vic411_atmos_data_struct)); 
  if (*atmos == NULL)
    vic411_vicerror("Memory allocation error in vic411_alloc_atmos().");

  for (i = 0; i < nrecs; i++) {
    (*atmos)[i].prec = (double *) calloc(vic411_NR+1, sizeof(double));
    if ((*atmos)[i].prec == NULL)
      vic411_vicerror("Memory allocation error in vic411_alloc_atmos().");      
    (*atmos)[i].air_temp = (double *) calloc(vic411_NR+1, sizeof(double));
    if ((*atmos)[i].air_temp == NULL)
      vic411_vicerror("Memory allocation error in vic411_alloc_atmos().");
    (*atmos)[i].wind = (double *) calloc(vic411_NR+1, sizeof(double));
    if ((*atmos)[i].wind == NULL)
      vic411_vicerror("Memory allocation error in vic411_alloc_atmos().");
    (*atmos)[i].vpd = (double *) calloc(vic411_NR+1, sizeof(double));
    if ((*atmos)[i].vpd == NULL)
      vic411_vicerror("Memory allocation error in vic411_alloc_atmos().");
    (*atmos)[i].vp = (double *) calloc(vic411_NR+1, sizeof(double));	
    if ((*atmos)[i].vp == NULL)
      vic411_vicerror("Memory allocation error in vic411_alloc_atmos().");
    (*atmos)[i].pressure = (double *) calloc(vic411_NR+1, sizeof(double));
    if ((*atmos)[i].pressure == NULL)
      vic411_vicerror("Memory allocation error in vic411_alloc_atmos().");
    (*atmos)[i].density = (double *) calloc(vic411_NR+1, sizeof(double));	
    if ((*atmos)[i].density == NULL)
      vic411_vicerror("Memory allocation error in vic411_alloc_atmos().");
    (*atmos)[i].shortwave = (double *) calloc(vic411_NR+1, sizeof(double));	
    if ((*atmos)[i].shortwave == NULL)
      vic411_vicerror("Memory allocation error in vic411_alloc_atmos().");
    (*atmos)[i].longwave = (double *) calloc(vic411_NR+1, sizeof(double));	
    if ((*atmos)[i].longwave == NULL)
      vic411_vicerror("Memory allocation error in vic411_alloc_atmos().");
    (*atmos)[i].snowflag = (char *) calloc(vic411_NR+1, sizeof(char));	
    if ((*atmos)[i].snowflag == NULL)
      vic411_vicerror("Memory allocation error in vic411_alloc_atmos().");
  }    			

}

/****************************************************************************/
/*	      		  vic411_free_atmos()                                      */
/****************************************************************************/
void vic411_free_atmos(int nrecs, vic411_atmos_data_struct **atmos)
/***************************************************************************
  Modifications:
  09-02-2003 Added check for LINK_DEBUG global option.  If LINK_DEBUG is
             TRUE atmospheric data is not dynamically allocated, so it
             should not be freed.                                   KAC
  2006-Sep-23 (Port from 4.0.6) Implemented flexible output configuration;
	      removed LDAS_OUTPUT and OPTIMIZE compile-time vic411_options.	TJB
  2006-Dec-20 All atmos_data arrays are always dynamically allocated now.	TJB
***************************************************************************/
{
  int i;

  if (*atmos == NULL)
    return;

  for (i = 0; i < nrecs; i++) {
    free((*atmos)[i].prec);
    free((*atmos)[i].air_temp);
    free((*atmos)[i].wind);
    free((*atmos)[i].vpd);
    free((*atmos)[i].vp);
    free((*atmos)[i].pressure);
    free((*atmos)[i].density);
    free((*atmos)[i].shortwave);
    free((*atmos)[i].longwave);
    free((*atmos)[i].snowflag);
  }

  free(*atmos);
}
