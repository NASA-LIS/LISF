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
#include <vicNl.h>

//BOP
//
// !ROUTINE: write_lakecon
// \label{write_lakecon}
// 
// !REVISION HISTORY: 
//  Sept 2011  James Geiger; Initial implementation
// 
// !INTERFACE:
void write_lakecon(lake_con_struct *lake_con)
// !DESCRIPTION: 
//  This routine writes contents of the lake\_con data structure to a log 
//  file, used primarily for debugging.
//
//  The arguments are:
//  \begin{description}
//  \item[lake\_con] given lake\_con data structure
//  \end{description}
//
//EOP
{

   FILE * fp;
   double tmpval;
   int i;

   fp=fopen("lakecon.log","a");
   fprintf(fp,"Lake_con Values:\n");

   //fprintf(fp,"\twetland_veg_class     = %d\n", lake_con->wetland_veg_class);
   fprintf(fp,"\tnumnod     = %d\n", lake_con->numnod);
   for ( i = 0; i < MAX_LAKE_NODES+1; ++i )
   {
      tmpval = lake_con->z[i];
      fprintf(fp,"\tz[%d]     = %f\n", i, tmpval);

      tmpval = lake_con->basin[i];
      fprintf(fp,"\tbasin[%d]     = %f\n", i, tmpval);

      tmpval = lake_con->Cl[i];
      fprintf(fp,"\tCl[%d]     = %f\n", i, tmpval);

      tmpval = lake_con->Cl[i];
      fprintf(fp,"\tCl[%d]     = %f\n", i, tmpval);
   }
   fprintf(fp,"\tb     = %f\n", lake_con->b);
   fprintf(fp,"\tmaxdepth     = %f\n", lake_con->maxdepth);
   fprintf(fp,"\tmindepth     = %f\n", lake_con->mindepth);
   fprintf(fp,"\tmaxvolume     = %f\n", lake_con->maxvolume);
   fprintf(fp,"\tminvolume     = %f\n", lake_con->minvolume);
   fprintf(fp,"\tbpercent     = %f\n", lake_con->bpercent);
   fprintf(fp,"\trpercent     = %f\n", lake_con->rpercent);
   fprintf(fp,"\teta_a     = %f\n", lake_con->eta_a);
   fprintf(fp,"\twfrac     = %f\n", lake_con->wfrac);
   fprintf(fp,"\tdepth_in     = %f\n", lake_con->depth_in);
   fprintf(fp,"\tlake_idx     = %d\n", lake_con->lake_idx);

   fclose(fp);
}
