/*
  Changes for the VIC implementation are preceded by the comment * start
  vic_change * and followed by the comment * end vic_change * */

/* RCS Id String
 * $Id: vic411_mtclim42_vic.h,v 5.4 2004/05/12 01:32:06 tbohn Exp $
 */

/* 
mtclim42.h
constants typedefs, and function prototypes for MTCLIM 4.2

Peter Thornton
NTSG, School of Forestry
University of Montana
5/10/98 

(dim) stands for dimensionless values

Adapted for inclusion in VIC-code:
Bart Nijssen
Sat Aug 21 16:58:43 1999
Last Changed: Fri Apr 25 11:37:15 2003 by Keith Cherkauer <cherkaue@u.washington.edu>
*/



#define TDAYCOEF 0.45     /* daylight air temperature coefficient (dim) */

#define SECPERRAD 13750.9871     /* seconds per radian of hour angle */
#define RADPERDAY 0.017214       /* radians of Earth orbit per julian day */
#define RADPERDEG 0.01745329     /* radians per degree */
#define MINDECL -0.4092797       /* minimum declination (radians) */
#define DAYSOFF 11.25            /* julian day offset of winter solstice */
/* start vic_change */
#define SRADDT 30.0             /* timestep for radiation routine (seconds) */
				/* Note:  Make sure that 3600 % SRADDT == 0 */
/* end vic_change */

#define MA       28.9644e-3      /* (kg mol-1) molecular weight of air */
#define MW       18.0148e-3      /* (kg mol-1) molecular weight of water */
#define R        8.3143          /* (m3 Pa mol-1 K-1) gas law constant */
#define G_STD    9.80665         /* (m s-2) standard gravitational accel. */ 
#define P_STD    101325.0        /* (Pa) standard pressure at 0.0 m elevation */
#define T_STD    288.15          /* (K) standard temp at 0.0 m elevation  */ 
#define CP       1010.0          /* (J kg-1 K-1) specific heat of air */
#define LR_STD   0.0065          /* (-K m-1) standard temperature lapse rate */
/* start vic_change */
#ifndef PI
#define PI       3.14159265
#endif
/* end vic_change */

/****************************
 **                         ** 
 **  STRUCTURE DEFINITIONS  **
 **                         **
 ****************************/
typedef struct
{
  int ndays;             /* number of days of data in input file */
  int indewpt;           /* input dewpoint temperature vic411_flag (0=NO, 1=YES) */
  int outhum;            /* output humidity vic411_flag            (0=VPD, 1=VP) */
  int inyear;            /* input year vic411_flag                 (0=NO, 1=YES) */
} control_struct;

typedef struct
{
  double base_elev;      /* base elevation, meters */
  double base_isoh;      /* base annual precip isohyet, cm */
  double site_lat;       /* site latitude, dec. degrees (- for south) */
  double site_elev;      /* site elevation, meters */
  double site_slp;       /* site slope, degrees */
  double site_asp;       /* site aspect, degrees */
  double site_isoh;      /* site annual precip isohyet, cm */
  double site_ehoriz;    /* site east horizon, degrees */
  double site_whoriz;    /* site west horizon, degrees */
  double tmax_lr;        /* maximum temperature lapse rate, deg C/1000m */
  double tmin_lr;        /* minimum temperature lapse rate, deg C/1000m */
} parameter_struct;

typedef struct
{
  int *year;             /* array of year values */
  int *yday;             /* array of yearday values */
  double *tmax;          /* array of base maximum temperature values */
  double *tmin;          /* array of base minimum temperature values */
  double *prcp;          /* array of base daily precipitation values */
  double *tdew;          /* array of base dewpoint temperature values */
  double *s_tmax;        /* array of site tmax values */
  double *s_tmin;        /* array of site tmin values */
  double *s_tday;        /* array of site daylight temperature values */
  double *s_prcp;        /* array of site prcp values */
  double *s_hum;         /* array of site humidity values (VPD or VP, Pa) */
  double *s_srad;        /* array of site shortwave radiation values */
  double *s_dayl;        /* array of site daylength values */
  /* start vic_change */
  double *s_tskc;	 /* array of cloudiness values */
  /* end vic_change */
} data_struct;

/********************************
 **                             **
 **    FUNCTION PROTOTYPES      **
 **                             **
 ********************************/
int vic411_calc_tair(const control_struct *ctrl, const parameter_struct *p, 
	      data_struct *data);
int vic411_calc_prcp(const control_struct *ctrl, const parameter_struct *p, 
	      data_struct *data);
/* start vic_change */
int vic411_calc_srad_humidity(const control_struct *ctrl, const parameter_struct *p, 
		       data_struct *data, double *tiny_radfract);
/* end vic_change */
/* start vic_change */
int vic411_calc_srad_humidity_iterative(const control_struct *ctrl,
				 const parameter_struct *p, data_struct *data,
				 double *hourly_radfract);
/* end vic_change */
int vic411_data_alloc(const control_struct *ctrl, data_struct *data);
int vic411_data_free(const control_struct *ctrl, data_struct *data);
double vic411_calc_pet(double rad, double ta, double pa, double dayl);
double vic411_atm_pres(double elev);
int vic411_pulled_boxcar(double *input,double *output,int n,int w,int w_flag);

