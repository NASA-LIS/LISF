/*******************************************************************************
 *
 *   Name: snip_ztif.h
 *
 *   Header file for ZTIF library
 *
 *   Updates
 *   =======
 *   20161208 Initial version ............................... M Puskar 16WS/WXE
 *   20250708 Ported for SNIP................................ E Kemp GSFC/SSAI
*
 ******************************************************************************/
#ifndef SNIP_ZTIF_H
#define SNIP_ZTIF_H

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

typedef struct SNIP_ZTIF {
    char*          fname;
    TIFF*          tif;
    GTIF*          gtif;
    tdata_t        rbuf;
    unsigned char* wbuf;
    char*          mode;
    int            width;
    int            length;
    long           size;
} SNIP_ZTIF;


int  SNIP_ZTIFOpen(SNIP_ZTIF *ztif, char *fname, const char *mode);
int  SNIP_ZTIFSetup(SNIP_ZTIF *ztif, int width, int length);
int  SNIP_ZTIFSetupRGB(SNIP_ZTIF *ztif, int width, int length);
int  SNIP_ZTIFReadline(SNIP_ZTIF *ztif, int linenum);
int  SNIP_ZTIFWriteline(SNIP_ZTIF *ztif, int linenum);
void SNIP_ZTIFClose(SNIP_ZTIF *ztif);

#endif
