/*  SUBROUTINE Rescales old states if new soil layer defined  */
/***   V. Koren   09/12   ***/

#include "com_header.h"
#include "linux.h"

extern void fstfg1_(int *, int *, int *, int *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *);

extern void soil_int1_(float *,int *,float *,float *,int *,float *);

extern void scale_sh2o_(float *, float *, float *, int *, int *, float *, float *, float *, float *, float *, float *, float *);

extern void frzind1_(float *, float *, float *, float *, int *, float *, float *, float *);

int scale_frzst(int npix, int *nsoil, int *nupl, float *sac_parx, 
               float *frz_parx, float *frz_stx, float *sac_stx,
               float *sacst_prv, float *smc, float *sh2o, float dthr,
               float *dtfrz, float *frzd_up, float *frzd_bt, int ivers,
               char* errmsg, int *nsoil_old, float *zsoil_old)

/*  VK 08/2012  Scaling states  */
/* Properties nsoil_old and zsoil_old are from old run with different 
   sac parameters, as well as states frz_stx, sac_stx, sacst_prv, smc, sh2o. 
   However, sac_parx, frz_parx, nsoil, and nupl are from new parametric grids.
   At the end of this subroutine, all old states will be rescaled to the new
   parametric data and stored to the same old state arrays.                  */
         	       
{

  int i, j, itxt, nmodold, nx;
  float brt, smax,smcx;
  float zbot, tbot, ts0, fgco[5], fgpm[7], sacst[5], tsoil[5], zsoil[5];
  float dsmodold[6], tsoil_old[6], smc_old[6];
/*  fraction of sand and clay for different soil texture   */
  float sand[12]={0.92,0.82,0.58,0.17,0.09,0.43,0.58,0.10,0.32,0.52,0.06,0.22};
  float clay[12]={0.03,0.06,0.10,0.13,0.05,0.18,0.27,0.34,0.34,0.42,0.47,0.58};
  int error = 0;
  float smup, smlo, suz, slz, rdup, rdlo, deltd, deltx;
  float delt[3], dmax[3], dover[3], dx, splt1, splt2;
  int ndelt[3]; 
  float tolerance=0.0001, sacmiss=-1.0, tsoilmiss=-999.0;

/* !!!   WARNING: below defined hardcoded 'iprint' variable is only for tests.
         if iprint = 0 no output of test information; other value leads to output.
         To print, it should be changed manually and executable regenerated        */
   int iprint=0;
   float sum_sacup, sum_saclo;

/*  Loop through pixels  */
  for(i=0; i<npix; i++) {

 if(iprint != 0) { 
  printf("\n"); 
  printf("UZTWM=%f UZFWM=%f LZTWM=%f LZFSM=%f LZFPM=%f stxt=%f\n",sac_parx[i],sac_parx[i+npix],sac_parx[i+8*npix],sac_parx[i+9*npix],sac_parx[i+10*npix],frz_parx[0]);
  }

/*  Check missing sac states. Skip this pixel if missing  */
   for(j=0; j<6; j++) {
    if(sac_stx[i+j*npix] == sacmiss)
    {
     sprintf( errmsg, "ERROR: SAC state #%d is -1.0, sac_stx[ j ] = %f\n !!! NO SIMULATIONS BE PERFORMED FURTHER AT THIS PIXEL\n",
                  j, sac_stx[ i + j * npix ] );
     return 1;
    }
   } 
    
/*  Relative states will be the same as old one, but real states in mm
    will be changed to match the new parameteric data                 */
/*    printf("sacst_before:\n");
    printf(" %f",sac_stx[i]);
    printf(" %f",sac_stx[i+npix]);
    printf(" %f",sac_stx[i+2*npix]);
    printf(" %f",sac_stx[i+3*npix]);
    printf(" %f",sac_stx[i+4*npix]);
    printf("\n");  */
   sac_stx[i]=sac_stx[i]*sac_parx[i];
   sac_stx[i+npix]=sac_stx[i+npix]*sac_parx[i+npix];
   sac_stx[i+2*npix]=sac_stx[i+2*npix]*sac_parx[i+8*npix];
   sac_stx[i+3*npix]=sac_stx[i+3*npix]*sac_parx[i+9*npix];
   sac_stx[i+4*npix]=sac_stx[i+4*npix]*sac_parx[i+10*npix];
   sac_stx[i+5*npix]=sac_stx[i+5*npix]*(sac_parx[i]+sac_parx[i+8*npix]);
   if(sac_stx[i+5*npix] < sac_stx[i]) sac_stx[i+5*npix]=sac_stx[i];
/*    printf("sacst_after :\n");
    printf(" %f",sac_stx[i]);
    printf(" %f",sac_stx[i+npix]);
    printf(" %f",sac_stx[i+2*npix]);
    printf(" %f",sac_stx[i+3*npix]);
    printf(" %f",sac_stx[i+4*npix]);
    printf("\n");     */
   if(ivers == 0) return error;  /*  no frozen states if only SAC  */

/*  Interpolate old soil temperature states to new soil layers     */
/*  No simulations will be performed if no surface temperature     */  
   zbot = -frz_parx[i+4*npix];
   tbot = frz_parx[i+npix];
   ts0 = frz_stx[i];

   if(frz_stx[i] == tsoilmiss) {
    sprintf(errmsg, "ERROR: No soil surface temperature. %d %f\n",
     i,frz_stx[i]);
    printf("!!! NO SIMULATIONS BE PERFORMED FURTHER AT THIS PIXEL\n");
    return 1;
   }

/*  Interpolate soil temperature states if old one found */
   nmodold=nsoil_old[i]+1;
   for(j=0; j<nsoil[i]; j++) zsoil[j]=-frz_parx[i+(9+j)*npix];
   for(j=0; j<nsoil_old[i]; j++) tsoil_old[j]=frz_stx[i+j*npix];
   tsoil_old[nmodold-1]=tbot-273.16;
   for(j=0; j<nsoil_old[i]; j++) dsmodold[j]=-zsoil_old[j];
   dsmodold[nmodold-1]=zbot;

 if(iprint != 0) {     
   printf("zsoil_old:");
   for(j=0; j<nsoil_old[i]; j++) printf(" %f",zsoil_old[j]);
   printf("\n");
   printf("zsoil_new:");
   for(j=0; j<nsoil[i]; j++) printf(" %f",-zsoil[j]);
   printf("\n");
   printf("tsoil_old:");
   for(j=0; j<nmodold; j++) printf(" %f",tsoil_old[j]);
   printf("\n");                                                              
 }
   soil_int1_(tsoil_old,&nmodold,dsmodold,zsoil,nsoil,tsoil);
   for(j=0; j<nsoil[i]; j++) frz_stx[i+j*npix]=tsoil[j];

  if(iprint != 0) {
   printf("tsoil_new:");
   for(j=0; j<nsoil[i]; j++) printf(" %f",tsoil[j]);
   printf("\n");
   printf("smc___old:");
   for(j=0; j<nsoil_old[i]; j++) printf(" %f",smc[j]);
   printf("\n");
   printf("sh2o__old:");
   for(j=0; j<nsoil[i]; j++) printf(" %f",sh2o[j]);
   printf("\n");
 }                                   

/*  Estimate new  SMC states  */
   for(j=0; j<5; j++) sacst[j]=sac_stx[i+j*npix];
   for(j=0; j<7; j++) fgpm[j]=frz_parx[i+(2+j)*npix];
   itxt=(frz_parx[i]+0.5)-1;
   brt=15.9*clay[itxt]+2.91;
   smax=-0.126*sand[itxt]+0.489;
   nx=0;
   for(j=0; j<nsoil_old[i]; j++) if(smc[j] < 0.0) nx=nx+1;
  /*    printf("soil Properties: stxt=%f smax=%f swlt=%f\n",itxt+1,smax,fgpm[6]);    */
  
   if(nx == 0) {
/* Old SMC states available; Interpolate & adjust them to get new SMC   */
    smc_old[0]=smc[1];
    smcx=smc[0]; 
    for(j=1; j<nsoil_old[i]; j++) smc_old[j]=smc[j];
    smc_old[nmodold-1]=smc[nsoil_old[i]-1];   
    soil_int1_(smc_old,&nmodold,dsmodold,zsoil,nsoil,smc);
    smc[0]=smcx;

 if(iprint != 0) {
  printf("smcint   :");
  for(j=0; j<nsoil[i]; j++) printf(" %f",smc[j]);
  printf("\n");
 }
                           
    for(j=0; j<nsoil[i]; j++) zsoil[j]=-zsoil[j];
    suz=sacst[0]+sacst[1];
    slz=sacst[2]+sacst[3]+sacst[4];
    smup=0.0;
    smlo=0.0;
    for(j=1; j<nsoil[i]; j++) {
     if(j < nupl[i]) {
      smup=smup+1000.0*(smc[j]-fgpm[6])*(zsoil[j-1]-zsoil[j]);
     }
     else {
      smlo=smlo+1000.0*(smc[j]-fgpm[6])*(zsoil[j-1]-zsoil[j]);
     }    
    }  
    rdup=smup-suz;
    rdlo=smlo-slz;
    deltx=0.0;
    deltd=sacst[0]*rdup/suz;
    sacst[0]=sacst[0]+deltd; 
    sacst[1]=sacst[1]+(rdup-deltd);
    if(sacst[1] < 0.0) {
     deltx=-sacst[1];
     sacst[1]=0.0;
    }
    if(sacst[1] > sac_parx[i+npix]) {
     deltx=sacst[1]-sac_parx[i+npix];
     sacst[1]=sac_parx[i+npix];
    } 
    sacst[0]=sacst[0]+deltx;
    if(sacst[0] > sac_parx[i]) {
     deltx=sacst[0] - sac_parx[i];
     sacst[0]=sac_parx[i];
     sacst[1]=sacst[1]+deltx;
     deltd=sacst[1]-sac_parx[i+npix];
     if(deltd > 0.0) {
      if(deltd < tolerance) {
       sacst[1] = sac_parx[i+npix];
      }
      else {
      printf("ERROR: Upper zone overfilled uztwc=%f uzfwc=%f deltx=%f\n",sacst[0],sacst[1],deltx);
      printf("      sacpar: uztwm=%f uzfwm=%f\n",sac_parx[i],sac_parx[i+npix]);
      printf("start: suz=%f smup=%f rdup=%f swlt=%f smax=%f zbot=%f\n",suz,smup,rdup,fgpm[6],smax,zbot);
/*       exit (1);  */
      }
     }
     if(sacst[1] < 0.0 || sacst[0] < 0.0) {
     printf("WARNING: SM below wilting point uztwc=%f uzfwc=%f\n",sacst[1],sacst[0]);
      sacst[1]=0.0;
      sacst[0]=0.0;
     }
    }

/*  Lower zone SAC states adjustment  */
    for(j=0; j<3; j++) {
     ndelt[j]=0;
     delt[j]=sacst[j+2]*rdlo/slz;
     dmax[j]=sac_parx[i+(j+8)*npix]-sacst[j+2];
     dover[j]=delt[j]-dmax[j];
     if(dover[j] > tolerance) {
      ndelt[j]=-1;
      delt[j]=dmax[j];
     }
     else {
      dover[j]=0.0;
     }
    }

/* 1st state over case  */
/*  printf("before lower: nsoil=%d nupl=%d\n",nsoil[i],nupl[i]);
  printf("states: lztwc=%f lzfsc=%f lzfpc=%f\n",sacst[2],sacst[3],sacst[4]);
  printf("slz=%f smlo=%f rdlo=%f rdlo/slz=%f\n",slz,smlo,rdlo,rdlo/slz);
  for(j=0; j<3; j++) printf("ndelt=%d delt=%f dover=%f dmax=%f",ndelt[j],delt[j],dover[j],dmax[j]);
  printf("\n");       */
 
  
    if(ndelt[0] == -1) {
     if(ndelt[1] == -1) {
      delt[2]=delt[2]+dover[0]+dover[1];
     }
     else {
      if(ndelt[2] == -1) {
       delt[1]=delt[1]+dover[0]+dover[2];
      }
      else {
       splt1=dover[0]*delt[1]/(delt[1]+delt[2]);
       splt2=dover[0]-splt1;
       delt[1]=delt[1]+splt1;
       delt[2]=delt[2]+splt2;
   /*       printf("splt: splt1=%f splt2=%f dover[0]=%f delt1=%f delt2=%f\n",splt1,splt2,dover[0],delt[1],delt[2]); */
       if(delt[1] > dmax[1]) {
        dx=delt[1]-dmax[1];
        delt[1]=dmax[1];
        delt[2]=delt[2]+dx;
   /*        printf("inif>: dx=%f delt1=%f delt2=%f\n",dx,delt[1],delt[2]);   */
       }
       else {
        if(delt[2] >dmax[2]) {
         dx=delt[2]-dmax[2];
         delt[2]=dmax[2];
         delt[1]=delt[1]+dx;
     /*      printf("inif>: dx=%f delt1=%f delt2=%f\n",dx,delt[1],delt[2]);  */
        }   
       }
      }
     }
    }
    else {    

/*  1st state not over case  */
     if(ndelt[1] == 0) {
      if(dover[2] > tolerance) {
       splt1=dover[2]*delt[0]/(delt[0]+delt[1]);
       splt2=dover[2]-splt1;
       delt[0]=delt[0]+splt1;
       delt[1]=delt[1]+splt2;
       if(delt[0] > dmax[0]) {
        dx=delt[0]-dmax[0];
        delt[0]=dmax[0];
        delt[1]=delt[1]+dx;
       }
       if(delt[1] > dmax[1]) {
        dx=delt[1]-dmax[1];
        delt[1]=dmax[1];
        delt[0]=delt[0]+dx;
       }
      }
     }
     else {
      if(ndelt[2] == 0) {
       splt1=dover[1]*delt[0]/(delt[0]+delt[2]);
       splt2=dover[1]-splt1;
       delt[0]=delt[0]+splt1;
       delt[2]=delt[2]+splt2;
       if(delt[0] > dmax[0]) {
        dx=delt[0]-dmax[0];
        delt[0]=dmax[0];
        delt[2]=delt[2]+dx;
       }
       if(delt[2] > dmax[2]) {
        dx=delt[2]-dmax[2];
        delt[2]=dmax[2];
        delt[0]=delt[0]+dx;
       }
       else {
        delt[0]=delt[0]+dover[1]+dover[2];
       }
      }
     }
    } 
    for(j=0; j<3; j++) sacst[j+2]=sacst[j+2]+delt[j];

/*  generate sh2o from smc_new and tsoil_new   */
    scale_sh2o_(smc,tsoil,fgpm,&nsoil[i],&nupl[i],zsoil,
                       sacst,&smax,&brt,&tbot,fgco,sh2o);

    for(j=0; j<5; j++) sacst_prv[j]=sacst[j];
   } 

/*  No old SMC states: Generate new SMC states using old SAC states  */
   else {
    nx = 1;
    smc[1] = -1.0;
    sh2o[1] = -1.0;
    fstfg1_(&nx,nsoil,nupl,&npix,fgco,sacst,tsoil,zsoil,fgpm,&tbot,
          &brt,&smax,sacst_prv,smc,sh2o,&dthr,dtfrz);
   }

/*  rescaling done; fill output states  */
   sac_stx[i]=sacst[0];
   sac_stx[i+npix]=sacst[1];
   sac_stx[i+2*npix]=sacst[2];
   sac_stx[i+3*npix]=sacst[3];
   sac_stx[i+4*npix]=sacst[4];
   sac_stx[i+5*npix]=sac_stx[i+5*npix]*(sac_parx[i]+sac_parx[i+8*npix]);
   if(sac_stx[i+5*npix] < sac_stx[i]) sac_stx[i+5*npix]=sac_stx[i];   

/*  estimate Upper-Lower frozen boundaries  */
   frzind1_(smc,sh2o,tsoil,zsoil,nsoil,&fgpm[6],frzd_up,frzd_bt);

/*  return SAC states into output states  */
   for(j=0; j<5; j++) frz_stx[i+(5+j)*npix]=fgco[j];
  }  
 
 if(iprint != 0) {
  printf("smc___new:");
  for(j=0; j<nsoil[i]; j++) printf(" %f",smc[j]);
  printf("\n");
  printf("sh2o__new:");
  for(j=0; j<nsoil[i]; j++) printf(" %f",sh2o[j]);
  printf("\n");
  printf("fgco__new:");
  for(j=0; j<nsoil[i]; j++) printf(" %f",fgco[j]);
  printf("\n");
  printf("sacst_new:");
  for(j=0; j<5; j++) printf(" %f",sacst[j]);
  printf("\n");
  printf("sacst_prv:");
  for(j=0; j<5; j++) printf(" %f",sacst_prv[j]);
  printf("\n"); 
  sum_sacup=frz_stx[5]+frz_stx[6];
  sum_saclo=frz_stx[7]+frz_stx[8]+frz_stx[9];
  printf("frzsmsum: sacup=%f smup=%f saclo=%f smlo=%f\n",sum_sacup,smup,sum_saclo,smlo);
  sum_sacup=sacst_prv[0]+sacst_prv[1];
  sum_saclo=sacst_prv[2]+sacst_prv[3]+sacst_prv[4];
  printf("prvsmsum: sacup=%f smup=%f saclo=%f smlo=%f\n",sum_sacup,smup,sum_saclo,smlo);
 }
/*  sum_sacup=sac_stx[0]+sac_stx[1];
  sum_saclo=sac_stx[2]+sac_stx[3]+sac_stx[4];
  printf("sacsmsum: sacup=%f smup=%f saclo=%f smlo=%f\n",sum_sacup,smup,sum_saclo,smlo); 
  sum_sacup=abs(sum_sacup-smup);
  sum_saclo=abs(sum_saclo-smlo);
  if(sum_sacup >tolerance || sum_saclo > tolerance) printf("ERROR_SCALE: sum_sacup=%f sum_saclo=%f\n",sum_sacup,sum_saclo);
  sum_sacup=sacst[0]-sac_parx[0];
  smup=sacst[1]-sac_parx[1];
  sum_saclo=sacst[2]-sac_parx[8];
  smlo=sacst[3]-sac_parx[9];
  dx=sacst[4]-sac_parx[10];
  if((sum_sacup > tolerance || sacst[0]<0.0) ||( smup>tolerance || sacst[1]<0.0) || (sum_saclo>tolerance || sacst[2]<0.0) ||(smlo>tolerance || sacst[3]<0.0) || (dx>tolerance || sacst[4]<0.0)) {
   printf("STERROR: duztc=%f uztc=%f duzfc=%f uzfc=%f dlztc=%f lztc=%f dlzfc=%f lzfc=%f dlzfpc=%f lzfpc=%f fpmax=%f\n",sum_sacup,sacst[0],smup,sacst[1],sum_saclo,sacst[2],smlo,sacst[3],dx,sacst[4],sac_parx[10]);
  }            */    
  
  return error;
}

      
