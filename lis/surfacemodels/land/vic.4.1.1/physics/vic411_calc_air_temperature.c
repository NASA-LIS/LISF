#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vic411_vicNl.h>

static char vcid[] = "$Id: calc_air_temperature.c,v 4.2 2006/10/31 01:18:42 vicadmin Exp $";

/****************************************************************************
  Subroutines developed by Bart Nijssen to estimate the daily temperature
  cycle from maximum and minimum daily temperature measurements.  

  Modifications:
  June 23, 1998 by Keith Cherkauer to be run within the VIC-NL model.
  August 19, 1999 by Bart Nijssen to function in polar regions where
    daylight or darkness may last for 24 hours.
  2006-Oct-26 Shift tminhour and tmaxhour if necessary to remain within
	      the current day. TJB

  ***************************************************************************/

/****************************************************************************/
/*				    vic411_hermite                                 */
/****************************************************************************/
/* calculate the coefficients for the Hermite polynomials */
void vic411_hermite(int n, 
	     double *x, 
	     double *yc1, 
	     double *yc2, 
	     double *yc3, 
	     double *yc4)
{
  int i;
  double dx;
  double divdf1;
  double divdf3;
  
  for (i = 0; i < n-1; i++) {
    dx = x[i+1] - x[i];
    divdf1 = (yc1[i+1] - yc1[i])/dx;
    divdf3 = yc2[i] + yc2[i+1] - 2 * divdf1;
    yc3[i] = (divdf1 - yc2[i] - divdf3)/dx;
    yc4[i] = divdf3/(dx*dx);
  }
}

/**************************************************************************/
/*				    vic411_hermint                               */
/**************************************************************************/
/* use the Hermite polynomials, to find the interpolation function value at 
   xbar */
double vic411_hermint(double xbar, int n, double *x, double *yc1, double *yc2, 
	       double *yc3, double *yc4)
{
  int klo,khi,k;
  double dx;
  double result;

  klo=0;
  khi=n-1;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (x[k] > xbar) khi=k;
    else klo=k;
  }

  dx = xbar - x[klo];
  result = yc1[klo] + dx * (yc2[klo] + dx * (yc3[klo] + dx * yc4[klo]));
  return result;
}

/****************************************************************************/
/*				    vic411_HourlyT                                 */
/****************************************************************************/
void vic411_HourlyT(int Dt, 
	     int ndays, 
	     int *TmaxHour, 
	     double *Tmax, 
	     int *TminHour,
	     double *Tmin, 
	     double *Tair) 
{
  double *x;
  double *Tyc1;
  double *yc2;
  double *yc3;
  double *yc4;
  int i;
  int j;
  int n;
  int hour;
  int nsteps;

  nsteps = HOURSPERDAY/Dt * ndays;

  n     = ndays*2+2;
  x     = (double *) calloc(n, sizeof(double));
  Tyc1  = (double *) calloc(n, sizeof(double));
  yc2   = (double *) calloc(n, sizeof(double));
  yc3   = (double *) calloc(n, sizeof(double));
  yc4   = (double *) calloc(n, sizeof(double));
  if (x == NULL || Tyc1 == NULL || yc2 == NULL || yc3 == NULL || yc4 == NULL)
    vic411_nrerror("Memory allocation failure in vic411_HourlyT()");
  
  /* First fill the x vector with the times for Tmin and Tmax, and fill the 
     Tyc1 with the corresponding temperature and humidity values */
  for (i = 0, j = 1, hour = 0; i < ndays; i++, hour += HOURSPERDAY) {
    if (TminHour[i] < TmaxHour[i]) {
      x[j]       = TminHour[i] + hour;
      Tyc1[j++]  = Tmin[i];
      x[j]       = TmaxHour[i] + hour;
      Tyc1[j++]  = Tmax[i];
    }
    else {
      x[j]       = TmaxHour[i] + hour;
      Tyc1[j++]  = Tmax[i];
      x[j]       = TminHour[i] + hour;
      Tyc1[j++]  = Tmin[i];
    } 
  }
  
  /* To "tie" down the first and last values, repeat those */
  x[0] = x[2] - HOURSPERDAY;
  Tyc1[0] = Tyc1[2];
  x[n-1] = x[n-3] + HOURSPERDAY;
  Tyc1[n-1] = Tyc1[n-3];

  /* we want to preserve maxima and minima, so we require that the first 
     derivative at these points is zero */
  for (i = 0; i < n; i++)
    yc2[i] = 0.;
  
  /* calculate the coefficients for the splines for the temperature */
  vic411_hermite(n, x, Tyc1, yc2, yc3, yc4);
  
  /* interpolate the temperatures */
  for (i = 0, hour = 0; i < nsteps; i++, hour += Dt) {
    Tair[i] = vic411_hermint(hour, n, x, Tyc1, yc2, yc3, yc4);
  }
  
  free(x);   
  free(Tyc1);
  free(yc2);
  free(yc3);
  free(yc4);
  
  return;
}

void vic411_set_max_min_hour(double *hourlyrad, 
		      int ndays, 
		      int *tmaxhour,
		      int *tminhour)
{
  int risehour;
  int sethour;
  int hour;
  int i;
  int j;
 
  /* treat the first day separately */
  risehour = -999;
  sethour = -999;
  if (hourlyrad[0] > 0 && hourlyrad[23] <= 0)
    risehour = 0;
  if (hourlyrad[0] <= 0 && hourlyrad[23] > 0)
    sethour = 0;
  for (hour = 1; hour < 24; hour++) {
    if (hourlyrad[hour] > 0 && hourlyrad[hour-1] <= 0)
      risehour = hour;
    if (hourlyrad[hour] <= 0 && hourlyrad[hour-1] > 0)
      sethour = hour;
  }
  if (risehour == -999 || sethour == -999) {
    /* arbitrarily set the min and max time to 2am and 2pm */
    tminhour[0] = 2;
    tmaxhour[0] = 14;
  }
  else {
    if (risehour > sethour)
      risehour -= 24;
    tmaxhour[0] = 0.67 * (sethour - risehour) + risehour;
    tminhour[0] = risehour - 1;
    // shift tminhour and tmaxhour if necessary to remain within the current day
    if (tminhour[0] < 0) tminhour[0] += 24;
    if (tmaxhour[0] < 0) tmaxhour[0] += 24;
  }

  /* treat remaining days */
  for (i = 1, hour = 24; i < ndays; i++) {
    risehour = -999;
    sethour = -999;
    for (j = 0; j < 24; j++, hour++) {
      if (hourlyrad[hour] > 0 && hourlyrad[hour-1] <= 0)
	risehour = j;
      if (hourlyrad[hour] <= 0 && hourlyrad[hour-1] > 0)
	sethour = j;
    }
    if (risehour == -999 || sethour == -999) {
      /* arbitrarily set the min and max time to 2am and 2pm */
      tminhour[i] = 2;
      tmaxhour[i] = 14;
    }
    else {
      if (risehour > sethour)
	risehour -= 24;
      tmaxhour[i] = 0.67 * (sethour - risehour) + risehour;
      tminhour[i] = risehour - 1;
      // shift tminhour and tmaxhour if necessary to remain within the current day
      if (tminhour[i] < 0) tminhour[i] += 24;
      if (tmaxhour[i] < 0) tmaxhour[i] += 24;
    }
  }
}
