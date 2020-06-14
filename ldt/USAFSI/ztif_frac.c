/*******************************************************************************
*
*  Name: ztif_frac.c
*
*  Get fractional snow from geotiff reference maps
*
*  Updates
*  =======
*  20161208 Initial version .................................. Puskar/16WS/WXE
*  20161219 Rewritten to match ztif library .................. Puskar/16WS/WXE
*  20190325 Ported to LDT...Eric Kemp, NASA GSFC/SSAI
*
*******************************************************************************/

#include "LDT_misc.h"
#include "ftn_drv.h"

/* EMK Only compile of LIBGEOTIFF support is enabled */
#ifdef USE_LIBGEOTIFF

#include <stdlib.h> 
#include <stdio.h>

#include "geotiffio.h"
#include "xtiffio.h"
#include "ztif.h"

// -----------------------------
// private structs
// -----------------------------

struct Parameters {
    char* snomap;
    char* snoage;
    int   max_age;
    int   offset;
    int   grid_width;
    int   grid_length;
    float snow_threshold;
    float bare_threshold;
};

struct Results {
    int *snow;
    int *bare;
    int counter[256];
};

struct ZTIFGroup {
    ZTIF map;
    ZTIF age;
};

struct GridCell {
    int width;
    int length;
    int size;
    int min_snow;
    int min_bare;
};


// -----------------------------
// private functions
// -----------------------------

static int
init_results(struct Results *r, int size) {
    int i;

    printf("... initializing result arrays\n");
    for(i=0; i<256; i++) {
        r->counter[i] = 0;
    }

    r->snow = malloc(size * sizeof(int));
    if(r->snow == NULL) {
        perror("init_results");
        return E_MALLOC;
    }
    r->bare = malloc(size * sizeof(int));
    if(r->bare == NULL) {
        perror("init_results");
        return E_MALLOC;
    }

    for(i=0; i<size; i++) {
        r->snow[i] = 0;
        r->bare[i] = 0;
    }

    return E_SUCCESS;
}

static int
open_files(struct ZTIFGroup *ztif, struct Parameters p) {
    int status;

    printf("... opening reference maps\n");
    if((status = ZTIFOpen(&ztif->map, p.snomap, "r")) != E_SUCCESS ||
       (status = ZTIFOpen(&ztif->age, p.snoage, "r")) != E_SUCCESS ){
        return status;
    }
    printf("... source snomap dimensions: width=%d length=%d\n", ztif->map.width, ztif->map.length);

    // double check map and age are the same size
    if(ztif->map.size != ztif->age.size) {
        fprintf(stderr, "The source map and source age files are not the same size!");
        return E_ARGS;
    }

    return status;
}

static void
init_cell(struct GridCell *cell, ZTIF source, struct Parameters p) {
    printf("... initializing grid cell\n");
    cell->width = source.width / p.grid_width;
    cell->length = source.length / p.grid_length;
    cell->size = cell->width * cell->length;
    cell->min_snow = p.snow_threshold * cell->size;
    cell->min_bare = p.bare_threshold * cell->size;
}

static int
process_grid(struct Results *r, struct ZTIFGroup *ztif, struct Parameters p, struct GridCell cell) {
    int i;
    int j;
    int index;
    int status;
    uint8 map_pixel;
    uint8 age_pixel;

    printf("... processing the grid\n");
    for(j=0; j<ztif->map.length; j++) {

        if((status = ZTIFReadline(&ztif->map, j)) != E_SUCCESS ||
           (status = ZTIFReadline(&ztif->age, j)) != E_SUCCESS ){
            break;
        }

        for(i=0; i<ztif->map.width; i++) {
            map_pixel = NO_DATA;
            age_pixel = ((uint8*)ztif->age.rbuf)[i] + p.offset;
            if(age_pixel <= p.max_age)
                map_pixel = ((uint8*)ztif->map.rbuf)[i];
           
            index = (i/cell.width) + ((j/cell.length)*p.grid_width);
            if(map_pixel == SNOW) {
                r->snow[index]++;
            }
            else if(map_pixel == NO_SNOW) {
                r->bare[index]++;
            }
        }
    }
    return status;
}

static void
fill_buffer(int *buffer, struct Results *r, struct Parameters p, struct GridCell cell) {
    int i;
    int j;
    int index;

    printf("... filling the output buffer\n");
    for(j=0; j<p.grid_length; j++) {
        for(i=0; i<p.grid_width; i++) {
            index = i + (j*p.grid_width);
            buffer[index] = NO_DATA;
            if(r->snow[index] > cell.min_snow) {
                buffer[index] = SNOW;
            }
            else if(r->bare[index] >= cell.min_bare) {
                buffer[index] = NO_SNOW;
            }
            r->counter[buffer[index]]++;
        }
    }
}

static void
print_hrule() {
    printf("------------------------------------------------------------------\n");
}

static void
print_header(char *header) {
    print_hrule();
    printf("%s\n", header);
    print_hrule();
}

static void
print_params(struct Parameters p)  {
    print_header("Running GETFRAC with the following parameters:");
    printf("  reference snomap: %s\n", p.snomap);
    printf("  reference snoage: %s\n", p.snoage);
    printf("  max pixel age   : %d\n", p.max_age);
    printf("  offset          : %d\n", p.offset);
    printf("  target width    : %d\n", p.grid_width);
    printf("  target length   : %d\n", p.grid_length);
    printf("  snow threshold  : %.2f\n", p.snow_threshold);
    printf("  bare threshold  : %.2f\n", p.bare_threshold);
    print_hrule();
}

static void
print_cell(struct GridCell c) {
    print_header("Target grid parameters:");
    printf("      width       : %d\n", c.width);
    printf("      length      : %d\n", c.length);
    printf("      size        : %d\n", c.size);
    printf("      minimum snow: %d\n", c.min_snow);
    printf("      minimum bare: %d\n", c.min_bare);
    print_hrule();
}

static void
print_results(struct Results r) {
    print_header("GETFRAC Results:");
    printf("           Bare   : %d\n", r.counter[0]);
    printf("           Snow   : %d\n", r.counter[1]);
    printf("           No data: %d\n", r.counter[255]);
    print_hrule();
}

// -----------------------------
// public functions
// -----------------------------

int
FTN(ztif_frac)(int *buffer, char *snomap, char *snoage, int *offset, 
	       int *max_age, int *width, int *length, float *snow_threshold, 
	       float *bare_threshold) {

    int status;
    struct Results r;
    struct GridCell cell;
    struct ZTIFGroup ztif;

    // get the parameters
    struct Parameters p = { snomap, snoage, (*max_age), (*offset), (*width), (*length), (*snow_threshold), (*bare_threshold) };
    print_params(p);

    // initialize the result struct
    if((status = init_results(&r, p.grid_width * p.grid_length)) != E_SUCCESS) {
        return status;
    }

    // open the source files
    if((status = open_files(&ztif, p)) != E_SUCCESS) {
        return status;
    }

    // initialize grid cell
    init_cell(&cell, ztif.map, p);

    // print grid info
    print_cell(cell);

    // parse the reference data
    status = process_grid(&r, &ztif, p, cell);

    // free memory and close images
    ZTIFClose(&ztif.map);
    ZTIFClose(&ztif.age);

    // fill the result buffer
    if(status == E_SUCCESS) {
        fill_buffer(buffer, &r, p, cell);
        print_results(r);
    }

    return status;
}

/*EMK Only compile if LIBGEOTIFF support is enabled */
#else

#endif
