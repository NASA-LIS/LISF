/*  SUBROUTINE GENERATE FROZEN GROUND STATES  */
/***   V. Koren   7/30/04   ***/

#include <stdlib.h>
#include "com_header.h"
#include "linux.h"
/* #include "models.h" */

extern void fstfg1_(int *, int *, int *, int *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *);

extern void frzind1_(float *, float *, float *, float *, int *, float *, float *, float *);

int gen_frzst(int npix, int *nsoil, int *nupl, float *sac_parx, 
               float *frz_parx, float *frz_stx, float *sac_stx,
               float *sacst_prv, float *smc, float *sh2o, float dthr,
               float *dtfrz, float *frzd_up, float *frzd_bt, int ivers,
               char* errmsg)
{

  int i, j, k, nx, itxt;
  float brt, smax, xx, sx;
  float zbot, tbot, zx, ts0, fgco[5], fgpm[7], sacst[5], tsoil[5], zsoil[5];
/*  fraction of sand and clay for different soil texture   */
  float sand[12]={0.92,0.82,0.58,0.17,0.09,0.43,0.58,0.10,0.32,0.52,0.06,0.22};
  float clay[12]={0.03,0.06,0.10,0.13,0.05,0.18,0.27,0.34,0.34,0.42,0.47,0.58};
  int error = 0;
  int fill_pixel_states=1, sacstmis;

/*vk2012 Defined Upper (upperbnd) Lower (lowerbnd) bounds of sac/frz storages,
         missing grids (grdmiss), missing values in grids (sacmiss and 
	 tsoilmiss), and tolerance level (tolerance)                         */
  float upperbnd=1.0,lowerbnd=0.0,tolerance=0.0001;
  float sacmiss=-1.0,tsoilmiss=-999.0,grdmiss=-999.0;  
/*  float sacmiss=-1.0,tsoilmiss=-99.0,grdmiss=-999.0;  */
/*vk2012  ----   NOTES ON RDHM ACTIONS ON STATES PROBLEM   ---------------   */
/*vk2012  RDHM takes following actions depending on state problems
    BASIC RELATIVE SAC STATES (range 0-1) UZTWC UZFWC LZTWC LZFSC LZFPC: 
       - If states are missing, it will not perform further simulations 
         at this pixel till states will be defined by user,
       - If states are out of defined here Upper/Lower bounds, it will 
         assign Upper/Lower bound values to respective SAC state and 
	 proceid with simulations; however, if states exceed bounds by 
	 more than defined here 'tolerance' level,it will provide WARNING 
	 message but proceads with simulations,
    SOIL TEMPERATURE STATES ts0 ts1 ts2 ts3 ts4:
       - If temperature at the 1st soil layer (ts0) is missing, it will not
         perform further simulations at this pixel till this state wiil be
         defined by user,
       - If temperature at the 1st soil layer is available but other soil
         layers temperature are missing, it will recover them from tso and
	 the bottom boundary layer temperature and proceid with simulations,
    FROZEN GROUND STATES UZTWH UZFWH LZTWH LZFSH LZFPH:
       - If state is missing (-1.0) at a pixel but grid is available, it will
         not perform further simulations at this pixel till states wiil be
	 defined by user,
       - If state grid is missing (-999.0), it will try to recover states from
         the basic SAC states and soil temperature,
       - If states are above upper bound (1.0), it will take similar action as 
         for the basic SAC state when thet exceed bounds,
       - If states are below lower bound (0.0), it will assign Lower bound 
         value to respective state and proceid with simulations. However, if 
	 state exceed bounds by more than defined tolerance level, it will 
	 provide ERROR message and STOP program.                              */
/*vk2012 ---------------------------------------------------------------------*/


  if( errmsg == NULL )
  {
   fprintf( stderr, "ERORR:gen_frzst.c --- string errmsg is not initialized!" );
    exit(1); 
  }
   
/*  Check SAC states & convert them into mm from relative values  */
  sacstmis=0;
  for(i=0; i<npix; i++) {
   for(j=0; j<6; j++) {
    if(sac_stx[i+j*npix] == sacmiss)
    {
    if(fill_pixel_states != 1) {
     sprintf( errmsg, "ERROR: SAC state #%d is -1.0, sac_stx[ j ] = %f\n !!! NO SIMULATIONS BE PERFORMED FURTHER AT THIS PIXEL\n",
                  j, sac_stx[ i + j * npix ] );
     return 1;
    }
    else {    
     if(j == 0) sac_stx[i+j*npix] = 0.46;
     if(j == 1) sac_stx[i+j*npix] = 0.14;
     if(j == 2) sac_stx[i+j*npix] = 0.56;
     if(j == 3) sac_stx[i+j*npix] = 0.11;
     if(j == 4) sac_stx[i+j*npix] = 0.46;
     if(j == 5) sac_stx[i+j*npix] = 1.00;
     sacstmis=sacstmis+1;
    }                                      
     
    }
     if(sac_stx[i+j*npix] > upperbnd) {
       xx=sac_stx[i+j*npix];
       sx=xx-upperbnd;
       sac_stx[i+j*npix]=upperbnd;
       if(sx > tolerance) {
        sprintf(errmsg, "WARNING: stateN%d > 1.0: value=%g changed to %g\n",
                      j,xx,sac_stx[i+j*npix]);
       }
     }
 
    if(sac_stx[i+j*npix] < lowerbnd) {
       xx=sac_stx[i+j*npix];
       sac_stx[i+j*npix]=lowerbnd;
       if(abs(xx) > tolerance) {
        sprintf(errmsg, "WARNING: stateN%d < 0.0: value=%g changed to %g\n",
                      j,xx,sac_stx[i+j*npix]); 
       }
    }
   }  

   sac_stx[i]=sac_stx[i]*sac_parx[i];
   sac_stx[i+npix]=sac_stx[i+npix]*sac_parx[i+npix];
   sac_stx[i+2*npix]=sac_stx[i+2*npix]*sac_parx[i+8*npix];
   sac_stx[i+3*npix]=sac_stx[i+3*npix]*sac_parx[i+9*npix];
   sac_stx[i+4*npix]=sac_stx[i+4*npix]*sac_parx[i+10*npix];
   sac_stx[i+5*npix]=sac_stx[i+5*npix]*(sac_parx[i]+sac_parx[i+8*npix]);
   if(sac_stx[i+5*npix] < sac_stx[i]) sac_stx[i+5*npix]=sac_stx[i];     

  } /* end i-loop  */
  
  if(ivers == 0) return error;  /*  no frozen states if only SAC  */
  
  for(i=0; i<npix; i++) {
   zbot = -frz_parx[i+4*npix];
   tbot = frz_parx[i+npix];
   ts0 = frz_stx[i];
   nx = 0;

/* VK 2012  Correcting checking missing values */
/*  check missing soil surface temperature     */
   zsoil[0]=frz_parx[i+9*npix];
   if(frz_stx[i] == tsoilmiss) {
 /*    sprintf(errmsg, "ERROR: No soil surface temperature. %d %f\n",
      i,frz_stx[i]);
     printf("!!! NO SIMULATIONS BE PERFORMED FURTHER AT THIS PIXEL\n");
     return 1;            */

    nx=1;
    if(fill_pixel_states == 1) {  
     frz_stx[i] = 15.0;
     ts0=15.0;
    }      
   }

/*  Check/Fill rest soil temperature states if there are no grids  */
   for(j=1; j<nsoil[i]; j++) {
    zsoil[j]=frz_parx[i+(9+j)*npix];    
    if(frz_stx[i+j*npix] == tsoilmiss) nx=nx+1;
   }
     
   if(nx >= 1) {
/*  Interpolate between first and bottom soil layers  */ 
    for(j=0; j<nsoil[i]; j++) {
     if(j == 0) 
      zx = -0.5*zsoil[0];
     else
      zx = -0.5*(zsoil[j-1]+zsoil[j]);
     tsoil[j] = ts0+(tbot-273.16-ts0)*zx/zbot;
     frz_stx[i+j*npix] = tsoil[j];
    }
   }   
   else { 
    for(j=0; j<nsoil[i]; j++) tsoil[j]=frz_stx[i+j*npix];
   }
   if(nx >= 1 && fill_pixel_states == 1) {
/*  Generate all other states w/out checking if fill SAC option*/
   printf("\n"); 
   printf("*** WARNING: Missing SAC states filled by default values:\n");
   printf("Filled SAC:");
   for(k=0; k<5; k++) printf(" %f",sac_stx[k]);
   printf("\n");
   printf("Filled tsoil:");
   for(k=0; k<nsoil[i]; k++) printf(" %f",tsoil[k]);
   printf("\n");
   nx=-9;
   }
   else {
/*  Check/Generate SAC unfrozen water storages if not fill SAC option */
    nx = 0;   
    for(j=0; j<5; j++) {
     k=i+(5+j)*npix;
     if(frz_stx[k] == sacmiss)
     {
      sprintf( errmsg, "ERROR: gen_frzst: frz state # %d is -1.0\n!!! NO SIMULATIONS BE PERFORMED FURTHER AT THIS PIXEL\n", 5 + j ); 
      return 1;
     }
     if(frz_stx[k] > upperbnd) {
       xx=frz_stx[k];
       sx=xx-upperbnd;
       frz_stx[k]=upperbnd;       
       if(sx > tolerance) {
       sprintf(errmsg, "ERROR: gen_frzst: frz_stateN%d > 1.0: value=%g changed to %g sx=%g\n",j,xx,frz_stx[k],sx); 
     error = 2;     
/*vk2012 One possibility not to exclude a pixel for the whole period
         would be replacement of 'error=2' statement by 'nx=nx+1' that
         will lead to regeneration of frozen states using SAC states & 
         tsoil similar to following case 'if(frz_stx[k] == grdmiss)'      */      
      }
     }
     if(frz_stx[k] == grdmiss) {
      nx=nx+1;
     } 
     else {
      if(frz_stx[k] < lowerbnd) {
       xx=frz_stx[k];
       if(abs(xx) > tolerance) {
       sprintf(errmsg, "ERROR: gen_frzst: frz_stateN%d < 0.0: value=%g changed to %g\n",j,xx,frz_stx[k]); 
       sprintf( errmsg, "***************WARNING!!!!!*******WARNING!!!!!*******************\n***************THIS PIXEL ONLY!!!*******************\n******** Consistently missing air temperature at this \n***********pixel can cause this problem!*****\n*********Soil moisture states will have missing values!********\n*****************************************\nFrozen ground sac states is negative!\n" );
       }
      }
     }
    }  

    for(j=0; j<5; j++) sacst[j]=sac_stx[i+j*npix];

/*  recalculate unfrozen water states into mm  */
    if(nx == 0) {
     frz_stx[i+5*npix]=frz_stx[i+5*npix]*sac_parx[i];
     frz_stx[i+6*npix]=frz_stx[i+6*npix]*sac_parx[i+npix];
     frz_stx[i+7*npix]=frz_stx[i+7*npix]*sac_parx[i+8*npix];
     frz_stx[i+8*npix]=frz_stx[i+8*npix]*sac_parx[i+9*npix];
     frz_stx[i+9*npix]=frz_stx[i+9*npix]*sac_parx[i+10*npix];
    }  

   }  /* end else sacstmis  */  

   for(j=0; j<5; j++) sacst[j]=sac_stx[i+j*npix];
   fgco[i]=frz_stx[i+5*npix];
   fgco[i+npix]=frz_stx[i+6*npix];
   fgco[i+2*npix]=frz_stx[i+7*npix];
   fgco[i+3*npix]=frz_stx[i+8*npix];
   fgco[i+4*npix]=frz_stx[i+9*npix];
    
/* fgpm array is a subset of frozen parameter array starting from 
   3rd element, RSMAX                                             */ 
   for(j=0; j<7; j++) fgpm[j]=frz_parx[i+(2+j)*npix];
    
/*  VK  8/2012  corrected itxt to correctly select from sand & clay arrays  */
   itxt=(frz_parx[i]+0.5)-1;
   brt = 15.9*clay[itxt] + 2.91;
   smax = -0.126*sand[itxt] + 0.489; 

   fstfg1_(&nx,nsoil,nupl,&npix,fgco,sacst,tsoil,zsoil,fgpm,&tbot,
          &brt,&smax,sacst_prv,smc,sh2o,&dthr,dtfrz);

/*  estimate Upper-Lower frozen boundaries  */
   frzind1_(smc,sh2o,tsoil,zsoil,nsoil,&fgpm[6],frzd_up,frzd_bt);
   
/*  return frozen ground SAC states  */
   for(j=0; j<5; j++) frz_stx[i+(5+j)*npix]=fgco[j];       

  }  /*  end i-loop  */

  return error;
}  

      
