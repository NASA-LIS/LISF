/*******************************************************************************
*
*  Name: ztif_read.c
*
*  Read/write data to/from geotiff files
*
*  Updates
*  =======
*  20161215 Initial version....................................Puskar/16WS/WXE
*  20190325 Ported to LDT...Eric Kemp, NASA GSFC/SSAI
*
*******************************************************************************/

#include "LDT_misc.h"
#include "ftn_drv.h"

/* EMK Only compile if LIBGEOTIFF support is enabled */

#ifdef USE_LIBGEOTIFF

#include <stdlib.h> 
#include <stdio.h>

#include "geotiffio.h"
#include "xtiffio.h"
#include "ztif.h"

int
FTN(ztif_read)(int *buffer, char *fname) {
    int i;
    int j;
    int index;
    int status;
    ZTIF ztif;

    if((status = ZTIFOpen(&ztif, fname, "r")) != E_SUCCESS) {
        return status;
    }
    for(j=0; j<ztif.length; j++) {
        if((status = ZTIFReadline(&ztif, j)) != E_SUCCESS) {
            return status;
        }    
        for(i=0; i<ztif.width; i++) {
            buffer[i + (j*ztif.width)] = ((uint8*)ztif.rbuf)[i];
        }
    }
    ZTIFClose(&ztif);
 
    return status;
}

/* EMK Only compile if LIBGEOTIFF support is enabled */
#else

#endif
