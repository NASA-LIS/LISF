/*******************************************************************************
 *
 *   Name: ztif.h
 *
 *   Header file for ZTIF library
 *
 *   Updates
 *   =======
 *   20161208 Initial version ............................... M Puskar 16WS/WXE
 *
 ******************************************************************************/
#ifndef ZTIF_H
#define ZTIF_H

#include "geotiffio.h"
#include "xtiffio.h"

enum _map_values {
    NO_SNOW,
    SNOW,
    NO_DATA = 255
};

enum _io_modes {
    READ,
    WRITE,
};

enum _color_modes {
    GRAYSCALE,
    COLOR,
};

enum _error_codes {
    E_SUCCESS,
    E_USAGE,
    E_XTIFFOPEN,
    E_TIFFMALLOC,
    E_GTIFNEW,
    E_TIFFREADSCANLINE,
    E_TIFFWRITESCANLINE,
    E_MALLOC,
    E_TIFFGETFIELD,
    E_ARGS,
};

typedef struct ZTIF {
    char*          fname;
    TIFF*          tif;
    GTIF*          gtif;
    tdata_t        rbuf;
    unsigned char* wbuf;
    char*          mode;
    int            width;
    int            length;
    long           size;
} ZTIF;


int  ZTIFOpen(ZTIF *ztif, char *fname, char *mode);
int  ZTIFSetup(ZTIF *ztif, int width, int length);
int  ZTIFSetupRGB(ZTIF *ztif, int width, int length);
int  ZTIFReadline(ZTIF *ztif, int linenum);
int  ZTIFWriteline(ZTIF *ztif, int linenum);
void ZTIFClose(ZTIF *ztif);

#endif
