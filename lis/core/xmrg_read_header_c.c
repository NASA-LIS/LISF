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
#include <math.h>
#include <stdlib.h>
#include <zlib.h>
#include "xmrg_reverse_byte.h"
#include "ftn.h"
//BOP 
// !ROUTINE: xmrg_read_header_c
// \label{xmrg_read_header_c}
// !INTERFACE:
int FTN(xmrg_read_header_c)(char  *in_file_name,
                         float *XOR, 
                         float *YOR,
                         int   *MAXX,
                         int   *MAXY,
                         float *CELL_SIZE, 
                         float *SCALE_FACTOR,
                         float *NODATA_VALUE, 
                         int   *NUM_BYTE,
                         int   *REVERSE_BYTE,
                         int   *REC2_LEN)
// !DESCRIPTION:
// This routine reads the first and the second records of XMRG/XMRG-like file. It retrieves the
// coordinate of the south west corner, numbers of column and rows of XMRG/XMRG-like data. If
// specified, cell size, scale factor and nodata value also will be retrieved. 
// 
// IO functions of LIBZ, such as gzopen, gzread, gzseek and gzclose are used instead one of
// standard IO library. LIBZ functions can deal with compressed (gzipped) and uncompressed
// formats. 
// 
// Functions invoked by xmrg\_read\_header\_c\_ are
// \begin{description}
//   \item[reverse\_byte\_order] (\ref{reverse_byte_order})
//   \item[gzopen]
//   \item[gzseek]
//   \item[gzread]
//   \item[gzclose]
// \end{description}
//EOP
{

    gzFile  in_file_ptr;

    int   rfchd[4];
    int   numbytes;

    int   i,j;
    short int reversebytes;
    void *buf4byte = malloc(4);

    /* end variable declaration */
    
    //printf("file name: %s\n", in_file_name);

    in_file_ptr = gzopen(in_file_name, "rb");
    if (in_file_ptr == NULL)
    {
        printf("xmrg_read_header_c: can not open file %s for input.\n", in_file_name);
        free(buf4byte);
        return(1);
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
    *REVERSE_BYTE = reversebytes;

    /* read the lenght of the second record 4+16+4*/    
    gzseek(in_file_ptr, 24, SEEK_SET);
    gzread(in_file_ptr, REC2_LEN, sizeof(int)*1);
    if(reversebytes==1)
    {
        reverse_byte_order(REC2_LEN, 1);
    }

    /* XMRG-like */
    if(*REC2_LEN == 16 || *REC2_LEN == 8)
    {
        /* locate the first record */
        gzseek(in_file_ptr, 4, SEEK_SET);

        /* read X origin (4 byte float) */
        gzread(in_file_ptr, buf4byte, 4);
        if(reversebytes==1)
        {
            reverse_byte_order((int*)buf4byte, 1);         
        }
        *XOR = *((float *)buf4byte);

        /* read Y origin (4 byte float) */
        gzread(in_file_ptr, buf4byte, 4);
        if(reversebytes==1)
        {
            reverse_byte_order((int*)buf4byte, 1);         
        }
        *YOR = *((float *)buf4byte);
        
        /* read num of columns (4 byte int) */
        gzread(in_file_ptr, buf4byte, 4);
        if(reversebytes==1)
        {
            reverse_byte_order((int*)buf4byte, 1);         
        }
        *MAXX = *((int *)buf4byte);

        /* read num of rows (4 byte int) */
        gzread(in_file_ptr, buf4byte, 4);
        if(reversebytes==1)
        {
            reverse_byte_order((int*)buf4byte, 1);         
        }
        *MAXY = *((int *)buf4byte);

        /* locate the second record */
        gzseek(in_file_ptr, 28, SEEK_SET);
        
        /* read scale factor (4 byte int) */
        gzread(in_file_ptr, buf4byte, 4);
        if(reversebytes==1)
        {
            reverse_byte_order((int*)buf4byte, 1);         
        }
        *SCALE_FACTOR = *((int *)buf4byte); 
        
        /* read num of bytes per value (4 byte int) */
        gzread(in_file_ptr, buf4byte, 4);
        if(reversebytes==1)
        {
            reverse_byte_order((int*)buf4byte, 1);         
        }
        *NUM_BYTE = *((int *)buf4byte);  
        
        if(*REC2_LEN == 16)
        {
            /* read cell size (4 byte float) */
            gzread(in_file_ptr, buf4byte, 4);
            if(reversebytes==1)
            {
                reverse_byte_order((int*)buf4byte, 1);         
            }
            *CELL_SIZE = *((float *)buf4byte);

            /* read nodata value (4 byte float) */
            gzread(in_file_ptr, buf4byte, 4);
            if(reversebytes==1)
            {
                reverse_byte_order((int*)buf4byte, 1);         
            }
            *NODATA_VALUE = *((float *)buf4byte);
        }
        else
        {
            /* default values if not specified */
            *CELL_SIZE    = 1.0;
            *NODATA_VALUE = -9999.0; 
        }
    }
    else /* XMRG */
    {
        /*SEEK_SET specifies the position offset from the beginning of the file*/
        gzseek(in_file_ptr, 4, SEEK_SET);
        for(i=0; i<4; i++)
        {
            gzread(in_file_ptr, &rfchd[i], sizeof(int)*1);
        }

        if (reversebytes == 1)
        {
            reverse_byte_order(rfchd,4);
        }
        *XOR  = 1.0*rfchd[0];
        *YOR  = 1.0*rfchd[1];
        *MAXX = rfchd[2];
        *MAXY = rfchd[3];

        /* default values */
        *NUM_BYTE     = 2; 
        *CELL_SIZE    = 1.0;
        *SCALE_FACTOR = 100;
        *NODATA_VALUE = -9999;
    }

    gzclose(in_file_ptr);
    free(buf4byte);
    return(0);
}  


