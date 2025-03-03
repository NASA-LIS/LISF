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
// !ROUTINE: vic412_bin_atmosdata
// \label{vic412_bin_atmosdata}
// 
// !REVISION HISTORY: 
//  Nov 2011  James Geiger; Initial implementation
// 
// !INTERFACE:
void write_bin_atmosdata(atmos_data_struct *atmos, int nrecs, 
                         float lat, float lng)
// !DESCRIPTION: 
//  This routine writes atmospheric data to a file in a binary format.
//
//  The arguments are:
//  \begin{description}
//  \item[atmos] given atmos data structure
//  \item[nrecs] length of atmos data structure (number of time-step records)
//  \item[lat] latitude of corresponding grid-cell
//  \item[lng] longitude of corresponding grid-cell
//  \end{description}
//
//EOP
{
#if LINK_DEBUG
  extern debug_struct debug;
  extern option_struct    options;

  int i;
  int j;

  char latchar[10], lngchar[10], junk[5];
  char fname[120];

  FILE * fp;

  //double dtmp;
  float ftmp;

  //printf("GREP: lat,lon = %f   %f\n",lat,lng);
  sprintf(junk, "%%.%if", options.GRID_DECIMAL);
  sprintf(latchar, junk, lat);
  sprintf(lngchar, junk, lng);

 
  /*  first write all the SNOW_STEP data  - only write if the modelstep !=
      SNOWSTEP */
  if (NR > 0)
  {
     sprintf(fname,"burst_forcing/VIC_snowstep_atmos_%s_%s",latchar,lngchar);
     fp = fopen(fname, "wb");

     for (i = 0; i < nrecs; i++)
     {
        for (j = 0; j < NF; j++)
        {
           //dtmp = (double)i;
           //fwrite(&dtmp,sizeof(double),1,fp);

           //dtmp = (double)j;
           //fwrite(&dtmp,sizeof(double),1,fp);

           ftmp = (float)atmos[i].snowflag[j];
           fwrite(&ftmp,sizeof(float),1,fp);

           ftmp = (float)atmos[i].prec[j];
           fwrite(&ftmp,sizeof(float),1,fp);

           ftmp = (float)atmos[i].air_temp[j];
           fwrite(&ftmp,sizeof(float),1,fp);

           ftmp = (float)atmos[i].wind[j];
           fwrite(&ftmp,sizeof(float),1,fp);

           ftmp = (float)atmos[i].vpd[j];
           fwrite(&ftmp,sizeof(float),1,fp);

           ftmp = (float)atmos[i].vp[j];
           fwrite(&ftmp,sizeof(float),1,fp);

           ftmp = (float)atmos[i].pressure[j];
           fwrite(&ftmp,sizeof(float),1,fp);

           ftmp = (float)atmos[i].density[j];
           fwrite(&ftmp,sizeof(float),1,fp);

           ftmp = (float)atmos[i].shortwave[j];
           fwrite(&ftmp,sizeof(float),1,fp);

           ftmp = (float)atmos[i].longwave[j];
           fwrite(&ftmp,sizeof(float),1,fp);

        }
     }

     fclose(fp);
  }
  

  sprintf(fname, "burst_forcing/VIC_modelstep_atmos_%s_%s",latchar,lngchar);
  fp = fopen(fname, "wb");

  for (i = 0; i < nrecs; i++)
  {
     //dtmp = (double)i;
     //fwrite(&dtmp,sizeof(double),1,fp);

     ftmp = (float)atmos[i].snowflag[NR];
     fwrite(&ftmp,sizeof(float),1,fp);

     ftmp = (float)atmos[i].prec[NR];
     fwrite(&ftmp,sizeof(float),1,fp);

     ftmp = (float)atmos[i].air_temp[NR];
     fwrite(&ftmp,sizeof(float),1,fp);

     ftmp = (float)atmos[i].wind[NR];
     fwrite(&ftmp,sizeof(float),1,fp);

     ftmp = (float)atmos[i].vpd[NR];
     fwrite(&ftmp,sizeof(float),1,fp);

     ftmp = (float)atmos[i].vp[NR];
     fwrite(&ftmp,sizeof(float),1,fp);

     ftmp = (float)atmos[i].pressure[NR];
     fwrite(&ftmp,sizeof(float),1,fp);

     ftmp = (float)atmos[i].density[NR];
     fwrite(&ftmp,sizeof(float),1,fp);

     ftmp = (float)atmos[i].shortwave[NR];
     fwrite(&ftmp,sizeof(float),1,fp);

     ftmp = (float)atmos[i].longwave[NR];
     fwrite(&ftmp,sizeof(float),1,fp);

  }

  fclose(fp);
#endif
}

