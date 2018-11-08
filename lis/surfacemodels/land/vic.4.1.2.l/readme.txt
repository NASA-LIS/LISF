LIS_rc%nt: number of vegetation classes in the landcover dataset
LIS_LMLC(i)%landcover(LIS_rc%lnc(i),LIS_rc%lnr(i), LIS_rc%nt) : data structure for land cover; 
i is the index of nest, 
lnc(i):  Array containing the East-West grid dimension of the running  grid for each processor, for each nest (including the halo regions)
lnr(i):  Array containing the North-South grid dimension of the running grid for each processor, for each nest (including the halo regions)

data structure of tile:
public tiledec
    type tiledec
        integer :: col        !Grid Column of Tile
        integer :: row        !Grid Row of Tile
        integer :: index      !Index of corresponding grid
        integer :: vegt       !Vegetation Type of Tile
        integer :: ensem      !ensemble id for the tile
        integer :: tile_id    !global catchment id for the tile
        integer :: d2g        !local tile count to global tile count
        real    :: fgrd       !Fraction of Grid covered by tile 
        real    :: pens      !ensemble weights
        real    :: com_lon    !center-of-mass longitude of the tile
        real    :: com_lat    !center-of-mass latitude of the tile

col and row are index over the domain for a specific processor.
see create_vegtilespace_latlon.F90
