!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: read_geosit
! \label{read_geosit}
!
! !REVISION HISTORY:
! 18 Mar 2015: James Geiger, initial code (based on merra-land)
! 20 Apr 2023: David Mocko,  initial code (based on merra2)
!
! !INTERFACE:
      subroutine read_geosit(n,order,month,findex,                     &
                           slvname,flxname,lfoname,radname,            &
                           geositforc,ferror)
! !USES:
      use LIS_coreMod,       only : LIS_rc, LIS_domain, LIS_masterproc
      use LIS_logMod
      use LIS_FORC_AttributesMod
      use LIS_metforcingMod, only : LIS_forc
      use geosit_forcingMod, only : geosit_struc
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
      use netcdf
#endif

      implicit none
! !ARGUMENTS:
      integer, intent(in)          :: n
      integer, intent(in)          :: order
      integer, intent(in)          :: month
      integer, intent(in)          :: findex
      character(len=*), intent(in) :: slvname
      character(len=*), intent(in) :: flxname
      character(len=*), intent(in) :: lfoname
      character(len=*), intent(in) :: radname
      real, intent(inout)          :: geositforc(geosit_struc(n)%nvars,&
                                            LIS_rc%lnc(n)*LIS_rc%lnr(n))
      integer, intent(out)         :: ferror

! !DESCRIPTION:
!  For the given time, reads parameters from GEOS-IT data, transforms
!  into 9 LIS forcing parameters and interpolates to the LIS domain. \newline
!
! GEOS-IT FORCING VARIABLES: \newline
!  1. tair    Temperature of air at 2-m, 10-m, or LML [$K$] \newline
!  2. qair    Specific humidity of air at 2-m, 10-m, or LML [$kg/kg$] \newline
!  3. swgdn   Downward shortwave radiation at the ground [$W/m^2$] \newline
!  4. lwgab   Downward longwave radiation at the ground [$W/m^2$] \newline
!  5. uwind   Zonal wind at 2-m, 10-m, or LML [$m/s$] \newline
!  6. vwind   Meridional wind at 2-m, 10-m, or LML [$m/s$] \newline
!  7. ps      Instantaneous Surface Pressure [$Pa$] \newline
!  8. pretot  Total precipitation [$mm/s$] \newline
!  9. precon  Convective precipitation [$mm/s$] \newline
! 10. presno  Precipitation falling as snow [$mm/s$] \newline
! 11. swland  Net shortwave radiation at the ground [$W/m^2$] \newline
! 12. pardr   Surface downward PAR beam flux [$W/m^2$] \newline
! 13. pardf   Surface downward PAR diffuse flux [$W/m^2$] \newline
! 14. hlml    Height of center of lowest model layer (LML) [$m$] \newline
!
!  The arguments are:
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous
!    1 hourly instance, order=2, read the next 1 hourly instance)
!  \item[n]
!    index of the nest
!  \item[name]
!    name of the 1 hour GEOS-IT analysis file
!  \item[tscount]
!    time step count
!  \item[ferror]
!    return error code (0 indicates success)
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
!  \end{description}
!EOP

      integer   :: ftn_slv,ftn_flx,ftn_lfo,ftn_rad
      integer   :: tmpId,qId,uId,vId,psId
      integer   :: prectotId,precconId,swgdnId,lwgabId
      integer   :: precsnoId,hlmlID
      integer   :: swlandId,pardrId,pardfId
      integer   :: nr_index,nc_index
      logical   :: file_exists,file_exists1
      integer   :: c,r,t,k,iret
      integer   :: mo
      logical   :: read_lnd

      real      :: tair(geosit_struc(n)%ncold,geosit_struc(n)%nrold)
      real      :: qair(geosit_struc(n)%ncold,geosit_struc(n)%nrold)
      real      :: uwind(geosit_struc(n)%ncold,geosit_struc(n)%nrold)
      real      :: vwind(geosit_struc(n)%ncold,geosit_struc(n)%nrold)
      real      :: ps(geosit_struc(n)%ncold,geosit_struc(n)%nrold)
      real      :: prectot(geosit_struc(n)%ncold,geosit_struc(n)%nrold)
      real      :: precsno(geosit_struc(n)%ncold,geosit_struc(n)%nrold)
      real      :: preccon(geosit_struc(n)%ncold,geosit_struc(n)%nrold)
      real      :: swgdn(geosit_struc(n)%ncold,geosit_struc(n)%nrold)
      real      :: lwgab(geosit_struc(n)%ncold,geosit_struc(n)%nrold)
      real      :: swland(geosit_struc(n)%ncold,geosit_struc(n)%nrold)
      real      :: pardr(geosit_struc(n)%ncold,geosit_struc(n)%nrold)
      real      :: pardf(geosit_struc(n)%ncold,geosit_struc(n)%nrold)
      real      :: hlml(geosit_struc(n)%ncold,geosit_struc(n)%nrold)
! __________________________________________________________________________

#if (defined USE_NETCDF3)
      write(LIS_logunit,*) "[ERR] GEOS-IT reader requires NetCDF4"
      call LIS_endrun()
#endif

#if (defined USE_NETCDF4)
      ferror = 0
      nr_index = geosit_struc(n)%nrold
      nc_index = geosit_struc(n)%ncold
      mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)

! Read single layer file (*slv) fields:
      inquire(file=slvname,exist=file_exists)
      if (file_exists) then
         write(LIS_logunit,*)                                          &
            '[INFO] Reading GEOS-IT file: ',trim(slvname)
         call LIS_verify(nf90_open(path=trim(slvname),                 &
                         mode=NF90_NOWRITE,ncid=ftn_slv),              &
                         'nf90_open failed for slvfile in read_geosit')

         call LIS_verify(nf90_inq_varid(ftn_slv,'PS',psId),            &
                         'nf90_inq_varid failed for ps in read_geosit')
         call LIS_verify(nf90_get_var(ftn_slv,psId,ps),                &
                         'nf90_get_var failed for ps in read_geosit')

! If using the GEOS-IT lowest model level forcing (*flx):
         if (geosit_struc(n)%uselml.eq.1) then
            inquire(file=flxname,exist=file_exists1)
            if (.not.file_exists1) then
               write(LIS_logunit,*)                                    &
                               '[ERR] ',trim(flxname)//' does not exist'
               call LIS_endrun()
            endif

            write(LIS_logunit,*)                                       &
               '[INFO] Reading GEOS-IT file: ',trim(flxname)
            write(LIS_logunit,*) ' (for lowest model level fields)'
            call LIS_verify(nf90_open(path=trim(flxname),              &
                            mode=NF90_NOWRITE,ncid=ftn_flx),           &
                          'nf90_open failed for flxfile in read_geosit')

            call LIS_verify(nf90_inq_varid(ftn_flx,'TLML',tmpId),      &
                    'nf90_inq_varid failed for tlml in read_geosit')
            call LIS_verify(nf90_get_var(ftn_flx,tmpId,tair),          &
                    'nf90_get_var failed for tlml in read_geosit')

            call LIS_verify(nf90_inq_varid(ftn_flx,'QLML',qId),        &
                    'nf90_inq_varid failed for qlml in read_geosit')
            call LIS_verify(nf90_get_var(ftn_flx,qId,qair),            &
                    'nf90_get_var failed for qlml in read_geosit')

            call LIS_verify(nf90_inq_varid(ftn_flx,'ULML',uId),        &
                    'nf90_inq_varid failed for ulml in read_geosit')
            call LIS_verify(nf90_get_var(ftn_flx,uId,uwind),           &
                    'nf90_get_var failed for ulml in read_geosit')

            call LIS_verify(nf90_inq_varid(ftn_flx,'VLML',vId),        &
                    'nf90_inq_varid failed for vlml in read_geosit')
            call LIS_verify(nf90_get_var(ftn_flx,vId,vwind), &
                    'nf90_get_var failed for vlml in read_geosit')

            call LIS_verify(nf90_close(ftn_flx),                       &
                    '[WARN] Failed to close flxfile in read_geosit.')

! Else use the single-layer fields (e.g., 2-m ref. height):
         else
            call LIS_verify(nf90_inq_varid(ftn_slv,'T2M',tmpId),       &
                    'nf90_inq_varid failed for t2m in read_geosit')
            call LIS_verify(nf90_get_var(ftn_slv,tmpId,tair),          &
                    'nf90_get_var failed for t2m in read_geosit')

            call LIS_verify(nf90_inq_varid(ftn_slv,'QV2M',qId),        &
                    'nf90_inq_varid failed for qv2m in read_geosit')
            call LIS_verify(nf90_get_var(ftn_slv,qId,qair),            &
                    'nf90_get_var failed for qv2m in read_geosit')

! Use either 2-m winds or 10-m winds:
            if (geosit_struc(n)%use2mwind.eq.1) then
               call LIS_verify(nf90_inq_varid(ftn_slv,'U2M',uId),      &
                       'nf90_inq_varid failed for u2m in read_geosit')
               call LIS_verify(nf90_get_var(ftn_slv,uId,uwind),        &
                       'nf90_get_var failed for u2m in read_geosit')

               call LIS_verify(nf90_inq_varid(ftn_slv,'V2M',vId),      &
                       'nf90_inq_varid failed for v2m in read_geosit')
               call LIS_verify(nf90_get_var(ftn_slv,vId,vwind),        &
                       'nf90_get_var failed for v2m in read_geosit')
            else
               call LIS_verify(nf90_inq_varid(ftn_slv,'U10M',uId),     &
                       'nf90_inq_varid failed for u10m in read_geosit')
               call LIS_verify(nf90_get_var(ftn_slv,uId,uwind),        &
                       'nf90_get_var failed for u10m in read_geosit')

               call LIS_verify(nf90_inq_varid(ftn_slv,'V10M',vId),     &
                       'nf90_inq_varid failed for v10m in read_geosit')
               call LIS_verify(nf90_get_var(ftn_slv,vId,vwind),        &
                       'nf90_get_var failed for v10m in read_geosit')
            endif
         endif

         call LIS_verify(nf90_close(ftn_slv),                          &
                 '[WARN] Failed to close slvfile in read_geosit.')

         call interp_geosit_var(n,findex,month,tair,  1,.false.,geositforc)
         call interp_geosit_var(n,findex,month,qair,  2,.false.,geositforc)
         call interp_geosit_var(n,findex,month,uwind, 5,.false.,geositforc)
         call interp_geosit_var(n,findex,month,vwind, 6,.false.,geositforc)
         call interp_geosit_var(n,findex,month,ps,    7,.false.,geositforc)
      else
         write(LIS_logunit,*) '[ERR] ',trim(slvname)//' does not exist'
         call LIS_endrun()
      endif

! Read in the flux file fields (*flx):
      inquire(file=flxname,exist=file_exists)
      if (file_exists) then
         write(LIS_logunit,*)                                          &
            '[INFO] Reading GEOS-IT file: ',trim(flxname)
         call LIS_verify(nf90_open(path=trim(flxname),                 &
                         mode=NF90_NOWRITE,ncid=ftn_flx),              &
                         'nf90_open failed for flxfile in read_geosit')

         call LIS_verify(nf90_inq_varid(ftn_flx,'PRECTOT',prectotId),  &
                 'nf90_inq_varid failed for prectot in read_geosit')
         call LIS_verify(nf90_get_var(ftn_flx,prectotId,prectot),      &
                 'nf90_get_var failed for prectot in read_geosit')

         call interp_geosit_var(n,findex,month,prectot,8,.true.,geositforc)

! Read in the convective precipitation, if selected in the forcing table:
         if (LIS_FORC_CRainf%selectOpt.eq.1) then
            call LIS_verify(nf90_inq_varid(ftn_flx,'PRECCON',precconId),&
                    'nf90_inq_varid failed for preccon in read_geosit')
            call LIS_verify(nf90_get_var(ftn_flx,precconId,preccon),   &
                    'nf90_get_var failed for preccon in read_geosit')

            call interp_geosit_var(n,findex,month,preccon,9,.true.,geositforc)
         endif

! Read in the snowfall field, if selected in the forcing input table:
         if (LIS_FORC_Snowf%selectOpt.eq.1) then
            call LIS_verify(nf90_inq_varid(ftn_flx,'PRECSNO',precsnoId),&
                    'nf90_inq_varid failed for precsno in read_geosit')
            call LIS_verify(nf90_get_var(ftn_flx,precsnoId,precsno),   &
                    'nf90_get_var failed for precsno in read_geosit')

            call interp_geosit_var(n,findex,month,precsno,10,.true.,geositforc)
         endif

! Read in Forcing Height, if selected in the forcing input table:
         if (LIS_FORC_Forc_Hgt%selectOpt.eq.1) then
            call LIS_verify(nf90_inq_varid(ftn_flx,'HLML',hlmlId),     &
                    'nf90_inq_varid failed for hlml in read_geosit')
            call LIS_verify(nf90_get_var(ftn_flx,hlmlId,hlml),         &
                    'nf90_get_var failed for hlml in read_geosit')

            call interp_geosit_var(n,findex,month,hlml,14,.false.,geositforc)
         endif

         call LIS_verify(nf90_close(ftn_flx),                          &
                 '[WARN] Failed to close flxfile in read_geosit.')

      else
         write(LIS_logunit,*) '[ERR] ',trim(flxname)//' does not exist'
         call LIS_endrun()
      endif

! Read in the radiation file fields (*rad):
      inquire(file=radname,exist=file_exists)
      if (file_exists) then
         write(LIS_logunit,*)                                          &
            '[INFO] Reading GEOS-IT file: ',trim(radname)
         call LIS_verify(nf90_open(path=trim(radname),                 &
                         mode=NF90_NOWRITE,ncid=ftn_rad),              &
                         'nf90_open failed for radfile in read_geosit')

         call LIS_verify(nf90_inq_varid(ftn_rad,'SWGDN',swgdnId),      &
                 'nf90_inq_varid failed for swgdn in read_geosit')
         call LIS_verify(nf90_get_var(ftn_rad,swgdnId,swgdn),          &
                 'nf90_get_var failed for swgdn in read_geosit')

         call LIS_verify(nf90_inq_varid(ftn_rad,'LWGAB',lwgabId),      &
                 'nf90_inq_varid failed for lwgab in read_geosit')
         call LIS_verify(nf90_get_var(ftn_rad,lwgabId,lwgab),          &
                 'nf90_get_var failed for lwgab in read_geosit')

         call interp_geosit_var(n,findex,month,swgdn,3,.false.,geositforc)
         call interp_geosit_var(n,findex,month,lwgab,4,.false.,geositforc)

         call LIS_verify(nf90_close(ftn_rad),                          &
                 '[WARN] Failed to close radfile in read_geosit.')
      else
         write(LIS_logunit,*) '[ERR] ',trim(radname)//' does not exist'
         call LIS_endrun()
      endif

! Checks: For reading in the surface layer file fields (*lfo):
      read_lnd = .false.
      if ((LIS_FORC_Pardr%selectOpt.eq.1).or.                          &
          (LIS_FORC_Pardf%selectOpt.eq.1).or.                          &
          (LIS_FORC_SWnet%selectOpt.eq.1)) then
         read_lnd = .true.
      endif

! Read in the surface layer fields (*lfo):
      if (read_lnd) then
         inquire(file=lfoname,exist=file_exists)
         if (file_exists) then
            write(LIS_logunit,*)                                       &
               '[INFO] Reading GEOS-IT file: ',trim(radname)
            call LIS_verify(nf90_open(path=trim(radname),              &
                            mode=NF90_NOWRITE,ncid=ftn_rad),           &
                          'nf90_open failed for radfile in read_geosit')

            if (LIS_FORC_SWnet%selectOpt.eq.1) then
               call LIS_verify(nf90_inq_varid(ftn_lfo,'SWLAND',swlandId),&
                      'nf90_inq_varid failed for swland in read_geosit')
               call LIS_verify(nf90_get_var(ftn_lfo,swlandId,swland),  &
                      'nf90_get_var failed for swland in read_geosit')

               call interp_geosit_var(n,findex,month,swland,11,.false.,geositforc)
            endif

            if (LIS_FORC_Pardr%selectOpt.eq.1) then
               call LIS_verify(nf90_inq_varid(ftn_lfo,'PARDR',pardrId),&
                       'nf90_inq_varid failed for pardr in read_geosit')
               call LIS_verify(nf90_get_var(ftn_lfo,pardrId,pardr),    &
                       'nf90_get_var failed for pardr in read_geosit')

               call interp_geosit_var(n,findex,month,pardr,12,.false.,geositforc)
            endif

            if (LIS_FORC_Pardf%selectOpt.eq.1) then
               call LIS_verify(nf90_inq_varid(ftn_lfo,'PARDF',pardfId),&
                       'nf90_inq_varid failed for pardf in read_geosit')
               call LIS_verify(nf90_get_var(ftn_lfo,pardfId,pardf),    &
                       'nf90_get_var failed for pardf in read_geosit')

               call interp_geosit_var(n,findex,month,pardf,13,.false.,geositforc)
            endif

            call LIS_verify(nf90_close(ftn_lfo),                       &
                    '[WARN] Failed to close lfofile in read_geosit.')

         else
            write(LIS_logunit,*) '[ERR] ',trim(lfoname)//              &
                                 ' does not exist'
            call LIS_endrun()
         endif
      endif
#endif

      end subroutine read_geosit

!BOP
!
! !ROUTINE: interp_geosit_var
! \label{interp_geosit_var}
!
! !INTERFACE:
      subroutine interp_geosit_var(n,findex,month,input_var,var_index, &
                                 pcp_flag,geositforc)

! !USES:
      use LIS_coreMod
      use LIS_logMod
      use LIS_spatialDownscalingMod
      use geosit_forcingMod, only : geosit_struc
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
      use netcdf
#endif

      implicit none
! !ARGUMENTS:
      integer, intent(in)    :: n
      integer, intent(in)    :: findex
      integer, intent(in)    :: month
      real,    intent(in)    :: input_var(geosit_struc(n)%ncold,       &
                                          geosit_struc(n)%nrold)
      integer, intent(in)    :: var_index
      logical, intent(in)    :: pcp_flag
      real,    intent(inout) :: geositforc(geosit_struc(n)%nvars,      &
                                          LIS_rc%lnc(n)*LIS_rc%lnr(n))

! !DESCRIPTION:
!  This subroutine spatially interpolates a GEOS-IT field
!  to the LIS running domain
!
!EOP
      integer   :: t,c,r,k,iret
      integer   :: doy
      integer   :: ftn
      integer   :: pcp1Id,pcp2Id,pcp3Id,pcp4Id,pcp5Id,pcp6Id
      real      :: f (geosit_struc(n)%ncold*geosit_struc(n)%nrold)
      logical*1 :: lb(geosit_struc(n)%ncold*geosit_struc(n)%nrold)
      logical*1 :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
      integer   :: input_size
      logical   :: scal_read_flag
! _____________________________________________________________

      input_size = geosit_struc(n)%ncold*geosit_struc(n)%nrold

!-----------------------------------------------------------------------
! Apply downscaling
!-----------------------------------------------------------------------
      lb = .true.
      do r = 1,geosit_struc(n)%nrold
         do c = 1,geosit_struc(n)%ncold
            k = c+(r-1)*geosit_struc(n)%ncold
            f(k) = input_var(c,r)
            if (f(k).eq.1.e+15) then 
               f(k)  = LIS_rc%udef
               lb(k) = .false. 
            endif
         enddo
      enddo

      if (pcp_flag.and.(LIS_rc%pcp_downscale(findex).ne.0)) then
! input_data becomes the ratio field.
         call LIS_generatePcpClimoRatioField(n,findex,"GEOS-IT",       &
                                             month,input_size,f,lb)
      endif

      if (pcp_flag.and.                                                &
          (trim(LIS_rc%met_interp(findex)).eq."budget-bilinear")) then

         call conserv_interp(LIS_rc%gridDesc(n,:),lb,f,lo,             &
                     geositforc(var_index,:),                          &
                     geosit_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),   &
                     LIS_domain(n)%lat,LIS_domain(n)%lon,              &
                     geosit_struc(n)%w112,geosit_struc(n)%w122,        &
                     geosit_struc(n)%w212,geosit_struc(n)%w222,        &
                     geosit_struc(n)%n112,geosit_struc(n)%n122,        &
                     geosit_struc(n)%n212,geosit_struc(n)%n222,        &
                     LIS_rc%udef,iret)

      elseif ((trim(LIS_rc%met_interp(findex)).eq."bilinear").or.      &
              (trim(LIS_rc%met_interp(findex)).eq."budget-bilinear")) then
         call bilinear_interp(LIS_rc%gridDesc(n,:),lb,f,lo,            &
                      geositforc(var_index,:),                         &
                      geosit_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),  &
                      LIS_domain(n)%lat,LIS_domain(n)%lon,             &
                      geosit_struc(n)%w111,geosit_struc(n)%w121,       &
                      geosit_struc(n)%w211,geosit_struc(n)%w221,       &
                      geosit_struc(n)%n111,geosit_struc(n)%n121,       &
                      geosit_struc(n)%n211,geosit_struc(n)%n221,       &
                      LIS_rc%udef,iret)

      elseif (trim(LIS_rc%met_interp(findex)).eq."neighbor") then
         call neighbor_interp(LIS_rc%gridDesc(n,:),lb,f,lo,            &
                      geositforc(var_index,:),                         &
                      geosit_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),  &
                      LIS_domain(n)%lat,LIS_domain(n)%lon,             &
                      geosit_struc(n)%n113,LIS_rc%udef,iret)

      else
         write(LIS_logunit,*) '[ERR] Spatial interpolation option '//  &
                               trim(LIS_rc%met_interp(findex))//       &
                              ' not supported for GEOS-IT.'
         call LIS_endrun()
      endif

      if (pcp_flag.and.(LIS_rc%pcp_downscale(findex).ne.0)) then
         call LIS_pcpClimoDownscaling(n,findex,month,                  &
                 LIS_rc%lnc(n)*LIS_rc%lnr(n),geositforc(var_index,:),lo)
      endif

      end subroutine interp_geosit_var

