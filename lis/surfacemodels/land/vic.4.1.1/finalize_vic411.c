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
#include "vic411_vicNl.h" 
#include "ftn.h"

//BOP
//
// !ROUTINE: finalize_vic411
// \label{finalize_vic411}
// 
// !REVISION HISTORY: 
//  12 Mar 2012  James Geiger; Initial implementation
// 
// !INTERFACE:
void FTN(finalize_vic411)(int * LIS_NEST, int * LIS_NTILES)
//
// !DESCRIPTION: 
//  This routine frees memory and closes opened files.
//
//  The arguments are:
//  \begin{description}
//  \item[LIS\_NEST] nest index, considered intent(in)
//  \item[LIS\_NTILES] number of sub-grid tiles, considered intent(in)
//  \end{description}
//
//EOP
{
   extern vic411_option_struct         vic411_options;
   extern vic411_out_data_file_struct **vic411_lis_out_data_files;

   int t, filenum;

   vic411_out_data_file_struct    *out_data_files;

/*
 * out_data_files are opened in make_in_outfiles.c.
 * This was commented out in that routine.
 * Thus there are no opened out_data_files.
   for ( t = 0; t < *LIS_NTILES; ++t )
   {
      out_data_files = vic411_lis_out_data_files[t];

      for (filenum=0; filenum<vic411_options.Noutfiles; filenum++)
      {
         fclose(out_data_files[filenum].fh);
         if(vic411_options.COMPRESS) vic411_compress_files(out_data_files[filenum].filename);
      }

   }
*/
}
