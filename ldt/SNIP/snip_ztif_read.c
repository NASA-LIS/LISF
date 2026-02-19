/*******************************************************************************
*
*  Name: snip_ztif_read.c
*
*  Read/write data to/from geotiff files
*
*  Updates
*  =======
*  20161215 Initial version....................................Puskar/16WS/WXE
*  20190325 Ported to LDT...Eric Kemp, NASA GSFC/SSAI
*  20250708 Ported to SNIP...Eric Kemp, NASA GSFC/SSAI
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
#include "snip_ztif.h"

/* Function prototype */
int
FTN(snip_ztif_read)(int *buffer, char *fname);

/* Function definition */
int
FTN(snip_ztif_read)(int *buffer, char *fname) {
    int i;
    int j;
    int status;
    SNIP_ZTIF ztif;

    if((status = SNIP_ZTIFOpen(&ztif, fname, "r")) != E_SUCCESS) {
        return status;
    }
    for(j=0; j<ztif.length; j++) {
        if((status = SNIP_ZTIFReadline(&ztif, j)) != E_SUCCESS) {
            return status;
        }
        for(i=0; i<ztif.width; i++) {
            buffer[i + (j*ztif.width)] = ((uint8*)ztif.rbuf)[i];
        }
    }
    SNIP_ZTIFClose(&ztif);

    return status;
}

/* EMK Only compile if LIBGEOTIFF support is enabled */
#else

#endif
