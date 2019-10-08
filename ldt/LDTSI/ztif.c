/*******************************************************************************
*
*  Name: ztif.c
*
*  Utility library for interacting with TIFF/GEOTIFF libraries
*
*  Updates
*  =======
*  20161208 Initial version.....................................Puskar/16WS/WXE
*  20190325 Ported to LDT...Eric Kemp, NASA GSFC/SSAI
*
*******************************************************************************/

#include "LDT_misc.h"
/*EMK Only compile if LIBGEOTIFF support is enabled */

#ifdef USE_LIBGEOTIFF

#include <stdlib.h> 

#include "geotiffio.h"
#include "xtiffio.h"
#include "ztif.h"

int
ZTIFOpen(ZTIF *ztif, char *fname, char *mode) {
    // initialize the ztif object
    ztif->fname = fname;
    ztif->tif = (TIFF*)0;
    ztif->gtif = (GTIF*)0;
    ztif->rbuf = NULL;
    ztif->wbuf = NULL;
    ztif->mode = mode;
    ztif->width = 0;
    ztif->length = 0;
    ztif->size = 0;

    // open the file
    ztif->tif = XTIFFOpen(ztif->fname, ztif->mode);
    if(!ztif->tif) {
            return E_XTIFFOPEN;
    }

    // get the dimensions and allocate a read buffer (if necessary)
    if(ztif->mode[0] == 'r') {
        if(!TIFFGetField(ztif->tif, TIFFTAG_IMAGEWIDTH,  &ztif->width) ||
           !TIFFGetField(ztif->tif, TIFFTAG_IMAGELENGTH, &ztif->length)) {
            return E_TIFFGETFIELD;
        }
        ztif->rbuf = _TIFFmalloc(ztif->width * sizeof(uint8));
        if(ztif->rbuf == NULL) {
            return E_TIFFMALLOC;
        }
        ztif->size = (long)ztif->width * (long)ztif->length;
    }
    return E_SUCCESS;
}

int
ZTIFSetup(ZTIF *ztif, int width, int length) {

    // make sure width and length are greater than 0
    if(width < 1 || length < 1) {
        perror("Width and length must be positive integers!");
        return E_ARGS;
    }

    double scale = 360.0 / width;
    double pixscale[3] = {scale, scale, 0.00};
    double tiepoints[6] = {0, 35999, 0, -179.9975, -89.9975, 0.0};

    // allocate a write buffer
    ztif->wbuf = malloc(width);
    if(ztif->wbuf == NULL) {
        perror("ZTIFSetup");
        return E_MALLOC;
    } 
    ztif->width  = width;
    ztif->length = length;
    ztif->size   = (long)width * (long)length;
    
    TIFFSetField(ztif->tif,TIFFTAG_IMAGEWIDTH,    width);
    TIFFSetField(ztif->tif,TIFFTAG_IMAGELENGTH,   length);
    TIFFSetField(ztif->tif,TIFFTAG_COMPRESSION,   COMPRESSION_PACKBITS);
    TIFFSetField(ztif->tif,TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(ztif->tif,TIFFTAG_ROWSPERSTRIP,  1L);
    TIFFSetField(ztif->tif,TIFFTAG_GEOTIEPOINTS,  6, tiepoints);
    TIFFSetField(ztif->tif,TIFFTAG_GEOPIXELSCALE, 3, pixscale);
    TIFFSetField(ztif->tif,TIFFTAG_PLANARCONFIG,  PLANARCONFIG_CONTIG);
    TIFFSetField(ztif->tif,TIFFTAG_PHOTOMETRIC,   PHOTOMETRIC_MINISBLACK);

    ztif->gtif = GTIFNew(ztif->tif);
    if(!ztif->gtif) {
        return E_GTIFNEW;
    }

    return E_SUCCESS;
}

int
ZTIFSetupRGB(ZTIF *ztif, int width, int length) {
    int status;
    
    // start with default
    if((status = ZTIFSetup(ztif, width, length)) != E_SUCCESS) {
        return status;
    }

    // define colors, default to red
    ushort r_pal[256] = {255};
    ushort g_pal[256] = {0};
    ushort b_pal[256] = {0};

    // black
    r_pal[0] = 0;
    g_pal[0] = 0;
    b_pal[0] = 0;

    // white
    r_pal[1] = 255;
    g_pal[1] = 255;
    b_pal[1] = 255;
    
    // green
    r_pal[2] = 0;
    g_pal[2] = 255;
    b_pal[2] = 0;
    
    // add palette info
    TIFFSetField(ztif->tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_PALETTE);
    TIFFSetField(ztif->tif, TIFFTAG_COLORMAP, r_pal, g_pal, b_pal);

    return E_SUCCESS;
}

int
ZTIFReadline(ZTIF *ztif, int linenum) {
    if(TIFFReadScanline(ztif->tif, ztif->rbuf, linenum, 0) != 1) {
        return E_TIFFREADSCANLINE;
    }
    return E_SUCCESS;
}

int
ZTIFWriteline(ZTIF *ztif, int linenum) {
    if((TIFFWriteScanline(ztif->tif, ztif->wbuf, linenum, 0)) != 1) {
        return E_TIFFWRITESCANLINE;
    }
    return E_SUCCESS;
}

void
ZTIFClose(ZTIF *ztif) {
    if(ztif->mode[0] == 'r') {
        _TIFFfree(ztif->rbuf);
    }
    else {
        free(ztif->wbuf);
    }

    if(ztif->gtif) {
        GTIFWriteKeys(ztif->gtif);
        GTIFFree(ztif->gtif);
    }

    if(ztif->tif) {
        XTIFFClose(ztif->tif);
    }
}

/*EMK Only compile if LIBGEOTIFF support is enabled */
#else

#endif
