#define NUM_ATMOS_FORCING ( 10 )
#define ATMOS_SNOWFLAG   ( 0 ) /* TRUE if there is snowfall in any of the snow
                                  bands during the timestep, FALSE otherwise*/
#define ATMOS_PREC       ( 1 ) /* average precipitation in grid cell (mm) */
#define ATMOS_AIR_TEMP   ( 2 ) /* air temperature (C) */
#define ATMOS_WIND       ( 3 ) /* wind speed (m/s) */
#define ATMOS_VPD        ( 4 ) /* atmospheric vapor pressure deficit (kPa) */
#define ATMOS_VP         ( 5 ) /* atmospheric vapor pressure (kPa) */
#define ATMOS_PRESSURE   ( 6 ) /* atmospheric pressure (kPa) */
#define ATMOS_DENSITY    ( 7 ) /* atmospheric density (kg/m^3) */
#define ATMOS_SHORTWAVE  ( 8 ) /* incoming shortwave radiation (W/m^2) */
#define ATMOS_LONGWAVE   ( 9 ) /* incoming longwave radiation (W/m^2) (net 
                                  incoming longwave for water balance model) */

#define VIC_KELVIN ( 273.15 )
#define VIC_EPS    ( 0.62196351 ) /* Ratio of molecular weights: 
                                     M_water_vapor/M_dry_air */
#define VIC_Rd     ( 287 )        /* Gas constant of dry air (J/degC*kg) */
