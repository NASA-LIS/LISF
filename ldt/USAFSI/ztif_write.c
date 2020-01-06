/*******************************************************************************
*
*  Name: ztif_writ.c
*
*  Write data from buffer to a geotiff file
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
FTN(ztif_write)(int *buffer, char *fname, int *width, int *length) {
    int i;
    int j;
    int index;
    int status;
    int cols = (*width);
    int rows = (*length);
    ZTIF ztif;

    if((status = ZTIFOpen(&ztif, fname, "w")) != E_SUCCESS || 
       (status = ZTIFSetup(&ztif, cols, rows)) != E_SUCCESS) {
        return status;
    }

    for(j=0; j<rows; j++) {
        for(i=0; i<cols; i++) {
            ztif.wbuf[i] = buffer[i + (j*cols)];
        }
        if((status = ZTIFWriteline(&ztif, j) != E_SUCCESS)) {
            break;
        }
    }
    ZTIFClose(&ztif);
 
    return status;
}

/* EMK Only compile if LIBGEOTIFF support is enabled */
#else

#endif
