/*  SUBROUTINE COMPUTES A NUMBER OF FROZEN GROUND CONSTANTS  */
/***   V. Koren   7/30/04   ***/

#include "com_header.h"
#include "linux.h"
#include "models.h"

void soilpar1_( float*, float*, float*, float*, float*, float*, float*,
		float*, float*, int*, int*, int*, float*, float*, float* );

/*
 * add one more argument to get smax
 */
void calc_frzcnst(int npix, float *sac_parx, float *frz_parx, 
	    int *nsoil, int *nupl, int *nsac, float *smax)
{

  int i, j, nsl, nup, nsc;
  float stxt, supm, slwm, psisat, brt, swlt, quartz;
  float stype, rtup, rtlw, zsoil[5];
  
  for(i=0; i<npix; i++) {
   stxt = frz_parx[i];
   supm = sac_parx[i] + sac_parx[i+npix];
   slwm = sac_parx[i+8*npix]+sac_parx[i+9*npix]+sac_parx[i+10*npix];

   soilpar1_(&stxt,&supm,&slwm,&smax[ i ],&psisat,&brt,&swlt,&quartz,&stype,
            &nsl,&nup,&nsc,zsoil,&rtup,&rtlw);
            
/*   printf(" calc_frzcnst: par= %f %f %f %f %f %f %f %f %f %f %d %d %d\n",
           rtup,rtlw,smax,psisat,swlt,zsoil[0],zsoil[1],zsoil[2],zsoil[3],
           zsoil[4],nsl,nup,nsc);
*/
           
   nsoil[i] = nsl;
   nupl[i] = nup;
   nsac[i] = nsc;

/*  convert TBOT from Celsius to Kelvin  */
/*
 * Parameter TBOT grid is in Kelvin already. 3/22/09
 */
/*  frz_parx[i+npix] =frz_parx[i+npix] + 273.16;  
 */
#if 0
/*
 *  ZC 4/8/09
 * Try to figure out the unit of TBOT is Celsius or Kelvin
 *   -99 < Celsius  <= 60 
 *   -99 + 273.16 < Kelvin < 273.16 + 50
 *    if -99 < TBOT <= 60, Celsius
 *    if -99 + 273.16 < TBOT <= 273.16 + 60 TBOT
 */
   if ( frz_parx[ i + npix ] > -99.f && frz_parx[ i + npix ] <= 60.f )
   {
      /* TBOT is Celsius, converting to Kelvin */
      frz_parx[i+npix] += 273.16; 
   }
   else if ( frz_parx[ i + npix ] > -99.f + 273.16 
		   && frz_parx[ i + npix ] <= 273.16 + 60.f )
   {
      /* TBOT is Kelvin, keep it as is */
   }
   else
   {
     fprintf( stderr, "TBOT Error: TBOT = %f\n", frz_parx[ i + npix ] );
   }
#endif /*#if 0*/
/*  change ZBOT to negative value as required by HRT1  */ 
   frz_parx[i+4*npix] = -frz_parx[i+4*npix];
/* vk 10/2010 added adjustment to ZBOT if it less than ZSOIL(NSL)  */
    if ( zsoil[nsl-1] <= frz_parx[i+4*npix] ) frz_parx[i+4*npix] = zsoil[nsl-1];
            
   frz_parx[i+5*npix] = rtup;
   frz_parx[i+6*npix] = rtlw;
   frz_parx[i+7*npix] = psisat;
   frz_parx[i+8*npix] = swlt;
   for(j=0; j<nsl; j++) frz_parx[i+(9+j)*npix] = zsoil[j];
   
  }
}               

