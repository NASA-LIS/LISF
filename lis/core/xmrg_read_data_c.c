//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------

#define XMRG_DEBUG 0
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <zlib.h>
#include "xmrg_reverse_byte.h"
#include "ftn.h"
//BOP
// !ROUTINE: xmrg_read_data_c
// \label{xmrg_read_data_c}
// !INTERFACE:
int FTN(xmrg_read_data_c)(char *in_file_name, float *data1d, float *nodata_default)
// !DESCRIPTION:
// This function reads the data records of XMRG/XMRG-like files. It automatically determine
// XMRG and XMRG-like, big endian and little endian, and short integer and real number. 
// Functions invoked by xmrg\_read\_header\_c are
// \begin{description}
//   \item[reverse\_byte\_order] (\ref{reverse_byte_order})
//   \item[reverse\_byte\_order\_short] (\ref{reverse_byte_order_short})
//   \item[gzopen]
//   \item[gzseek]
//   \item[gzread]
//   \item[gzclose]
// \end{description}
//EOP
{
    gzFile  in_file_ptr;
    int   rfchd[4];
    int   numbytes_a[2];
    int   numbytes;
    int   numbytes_cell;

    int   MAXX,MAXY;
    float XORIG,YORIG;
    int   i, j;
    unsigned long k; 
    
    short *onerow_short;
    float *onerow_float;

    float outval;
    short int reversebytes;
    int rec2_len;
    float scale_factor; 
    float nodata_value; 
    void *buf4byte = malloc(4);

    onerow_short = NULL;
    onerow_float = NULL;

    /* set the default no-data value */
    nodata_value = *nodata_default; 

    /* end variable declaration */

    in_file_ptr=gzopen(in_file_name, "rb");
    if (in_file_ptr == NULL)
    {
        printf("xmrg_read_data_c_: can not open file %s for input.\n", in_file_name);
        free(buf4byte);
        return(1) ;
    }

    /* start reading the XMRG file*/

    /* determine if byte reversal is needed */
    gzread(in_file_ptr, &numbytes, sizeof(int)*1);
    if (numbytes != 16)
    {
        reversebytes = 1;
    }
    else
    {
        reversebytes = 0;
    }
    
    /* locate the second record */
    gzseek(in_file_ptr, 24, SEEK_SET);
    gzread(in_file_ptr, &rec2_len, 4);
    
    if(reversebytes==1)
    {
        reverse_byte_order(&rec2_len, 1); 
    }
    

    /*SEEK_SET specifies the position offset from the beginning of the file*/
    gzseek(in_file_ptr, 4, SEEK_SET);
    for(i=0; i<4; i++)
    {
        gzread(in_file_ptr, &rfchd[i], sizeof(int)*1);
    }

    if (reversebytes)
    { 
        reverse_byte_order(rfchd, 4);
    }

     
    /* XMRG-like uses float XOR and YOR*/
    if(rec2_len == 24 || rec2_len == 16)
    {
        /* locate "num of byte"  */
        gzseek(in_file_ptr, 32, SEEK_SET);
        gzread(in_file_ptr, &numbytes_cell, 4);
        if(reversebytes==1)
        {
            reverse_byte_order(&numbytes_cell, 1);
        }
        gzseek(in_file_ptr, 4, SEEK_SET);
        
        gzread(in_file_ptr, buf4byte, 4);
        if(reversebytes==1)
        {
            reverse_byte_order((int*)buf4byte, 1);
        }
        XORIG = *((float *)buf4byte);
        
        gzread(in_file_ptr, buf4byte, 4);
        if(reversebytes==1)
        {
            reverse_byte_order((int*)buf4byte, 1);
        }
        YORIG = *((float *)buf4byte);
    } 
    /* traditional XMRG used integer XOR and YOR */
    else
    {
        numbytes_cell = 2;
        XORIG=rfchd[0];
        YORIG=rfchd[1];
    }
    MAXX=rfchd[2];
    MAXY=rfchd[3];
#if(XMRG_DEBUG)
    /*echo to screen*/
    printf("reverse byte order: %d\n", reversebytes);
    printf("x-coordinate (HRAP) of the lowerleft corner of the RFC rectangle %10.4f\n",XORIG);
    printf("y-coordinate (HRAP) of the lowerleft corner of the RFC rectangle %10.4f\n",YORIG);
    printf("number of HRAP bins along the x-coordinate in the RFC rectangle %d\n",MAXX);
    printf("number of HRAP bins along the y-coordinate in the RFC rectangle %d\n",MAXY);
#endif
    /*each record is preceded and followed by 4 bytes*/
    /*first record is 4+16+4 bytes*/
    gzseek(in_file_ptr, 20, SEEK_SET);

    /*read second FORTRAN record*/
    /*here I am reading an array with two elements instead of 1*/
    /*because I couldn't successfully get the reverse_byte_order*/
    /*routine to work with other approaches*/
    gzread(in_file_ptr, &numbytes_a, sizeof(int)*2);
    if (reversebytes)
        reverse_byte_order(numbytes_a,2);
    numbytes=numbytes_a[1];
  
    /* default scale factor 100 */
    scale_factor = 100.0;

    /***********************************************************
    *   account for all possible header lengths in xmrg files
    ************************************************************/
  
    if ((int) numbytes == 66)
    {
        /* first record (24) plus second record(74=66+8) is 98*/
        gzseek(in_file_ptr, 98, SEEK_SET);
    }
    else if ((int) numbytes==38)
    {
        gzseek(in_file_ptr, 70, SEEK_SET);
    }
    else if ((int) numbytes==37)
    {
        gzseek(in_file_ptr, 69, SEEK_SET);
#if(XMRG_DEBUG)        
        printf("WARNING: SECOND RECORD ONLY HAS 37 BYTES\n");
        printf("SHOULD HAVE 38 BYTES\n");
        printf("Assuming data is still valid. . . \n");
#endif
    }
    else if ((int) numbytes == (MAXX*2))
    {
#if(XMRG_DEBUG)
        printf("Reading pre-1997 format.\n");
#endif
        gzseek(in_file_ptr, 24, SEEK_SET);
    }
    else if ((int) numbytes == 16 )
    {
#if(XMRG_DEBUG)
        printf("Reading XMRG-like (length of 2nd record is 16) \n");
#endif
        /* read scale factor */
        gzseek(in_file_ptr, 28, SEEK_SET);
        gzread(in_file_ptr, buf4byte, 4);
        if(reversebytes==1)
        {
            reverse_byte_order((int*)buf4byte, 1);
        }
        scale_factor = *((int*)buf4byte) * 1.0;

        /* read nodata value */
        gzseek(in_file_ptr, 36, SEEK_SET);
        gzread(in_file_ptr, buf4byte, 4);
        if(reversebytes==1)
        {
            reverse_byte_order((int*)buf4byte, 1);
        }
        nodata_value = *((float*)buf4byte);

        gzseek(in_file_ptr, 48, SEEK_SET);
    }
    else if((int) numbytes == 8 )
    {
#if(XMRG_DEBUG)
        printf("Reading XMRG-like (length of 2nd record is 8)\n");
#endif
        /* read scale factor */
        gzseek(in_file_ptr, 28, SEEK_SET);
        gzread(in_file_ptr, buf4byte, 4);
        if(reversebytes==1)
        {
            reverse_byte_order((int*)buf4byte, 1);
        }
        scale_factor = *((int*)buf4byte) * 1.0; 
        gzseek(in_file_ptr, 40, SEEK_SET);
    }
    else
    {
        /*printf("numbytes %d\n",numbytes);*/
        printf("Error! Header file is in a nonstandard format. Data NOT READ!\n");
        free(buf4byte);
        return(1);
    }

    /* allocate memory for arrays */
    if(numbytes_cell == 2 )
    {
        onerow_short = (short int*) malloc(sizeof(short int*)*MAXX);
    }
    else
    {
        onerow_float = (float*) malloc(sizeof(float*)*MAXX);
    }
    
    /* set the counter for data1d */
    k = 0;  
    for(j=0; j<MAXY; j++)
    {
        gzseek(in_file_ptr, 4, SEEK_CUR);
        /* read one row */
        if(numbytes_cell == 2)
        {
            gzread(in_file_ptr, onerow_short, sizeof(short)*MAXX);
            if(reversebytes==1)
            {
                reverse_byte_order_short(onerow_short, MAXX);
            }
        }
        else
        {
            gzread(in_file_ptr, onerow_float, sizeof(float)*MAXX);
            if(reversebytes==1)
            {
                reverse_byte_order((int*)onerow_float, MAXX);
            }
        }
        
        gzseek(in_file_ptr, 4, SEEK_CUR);

        for(i=0; i<MAXX; i++)
        {
            if(numbytes_cell == 2)
            {
                outval = (float) onerow_short[i];
            }
            else
            {
                outval = onerow_float[i];
            }
            // check nodata value 
            if(fabs(outval-nodata_value)<1e-6)
            {
                data1d[k++] = outval;
            }
            else
            {
                data1d[k++] = outval/scale_factor;
            }
        } 

    } 

    /*free allocated memory*/
    free(onerow_short);
    free(onerow_float);
    free(buf4byte);
    gzclose(in_file_ptr);
    return(0);
}  /**  END OF xmrg_read_data_c  **/

