!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: read_merra2
! \label{read_merra2}
! 
! !REVISION HISTORY:
! 18 Mar 2015: James Geiger, initial code (based on merra-land)
!
! !INTERFACE:
subroutine read_merra2(n, order, month, findex,          &
                       slvname, flxname, lfoname, radname, &
                       merraforc, ferror)
! !USES:
  use LIS_coreMod,       only : LIS_rc, LIS_domain, LIS_masterproc
  use LIS_logMod
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod, only : LIS_forc
  use merra2_forcingMod, only : merra2_struc
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
  real, intent(inout)          :: merraforc(merra2_struc(n)%nvars, 24, &
       LIS_rc%lnc(n)*LIS_rc%lnr(n))
  integer, intent(out)         :: ferror          

!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  MERRA2 data, transforms into 9 LIS forcing 
!  parameters and interpolates to the LIS domain. \newline
!
! merra2 FORCING VARIABLES (unless noted, fields are 1-hr upstream averaged): \newline
!  1. T 2m    Temperature interpolated to 2 metres [$K$] \newline
!  2. q 2m    Instantaneous specific humidity interpolated to 2 metres[$kg/kg$] \newline
!  3. radswg  Downward shortwave flux at the ground [$W/m^2$] \newline
!  4. lwgdwn  Downward longwave radiation at the ground [$W/m^2$] \newline
!  5. u 10m   Instantaneous zonal wind interpolated to 10 metres [$m/s$] \newline
!  6. v 10m   Instantaneous meridional wind interpolated to 10 metres[$m/s$] \newline
!  7. ps      Instantaneous Surface Pressure [$Pa$] \newline
!  8. preacc  Total precipitation [$mm/s$] \newline
!  9. precon  Convective precipitation [$mm/s$] \newline
! 10. albedo  Surface albedo (0-1) \newline
!
!  The arguments are: 
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    1 hourly instance, order=2, read the next 1 hourly instance)
!  \item[n]
!    index of the nest
!  \item[name]
!    name of the 1 hour MERRA2 analysis file
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
  
  integer   :: ftn_slv, ftn_flx, ftn_lfo,ftn_rad
  integer   :: tmpId, qId, uId, vId, psId
  integer   :: prectotId, precconId, swgdnId, lwgabId, emisId
  integer   :: precsnoId, hlmlID
  integer   :: swlandId, pardrId, pardfId
  integer   :: nr_index, nc_index
  logical   :: file_exists, file_exists1
  integer   :: c,r,t,k,iret
  integer   :: mo
  logical   :: read_lnd

  real      :: tair(merra2_struc(n)%ncold, merra2_struc(n)%nrold,24)
  real      :: qair(merra2_struc(n)%ncold, merra2_struc(n)%nrold,24)
  real      :: uwind(merra2_struc(n)%ncold, merra2_struc(n)%nrold,24)
  real      :: vwind(merra2_struc(n)%ncold, merra2_struc(n)%nrold,24)
  real      :: ps(merra2_struc(n)%ncold, merra2_struc(n)%nrold,24)
  real      :: prectot(merra2_struc(n)%ncold, merra2_struc(n)%nrold,24)
  real      :: precsno(merra2_struc(n)%ncold, merra2_struc(n)%nrold,24)
  real      :: preccon(merra2_struc(n)%ncold, merra2_struc(n)%nrold,24)
  real      :: prectot_flx(merra2_struc(n)%ncold, merra2_struc(n)%nrold,24)
  real      :: preccon_flx(merra2_struc(n)%ncold, merra2_struc(n)%nrold,24)
  real      :: swgdn(merra2_struc(n)%ncold, merra2_struc(n)%nrold,24)
  real      :: lwgab(merra2_struc(n)%ncold, merra2_struc(n)%nrold,24)
  real      :: swland(merra2_struc(n)%ncold, merra2_struc(n)%nrold,24)
  real      :: pardr(merra2_struc(n)%ncold, merra2_struc(n)%nrold,24)
  real      :: pardf(merra2_struc(n)%ncold, merra2_struc(n)%nrold,24)
  real      :: hlml(merra2_struc(n)%ncold, merra2_struc(n)%nrold,24)
!  real      :: emis(merra2_struc(n)%ncold, merra2_struc(n)%nrold,24)
! __________________________________________________________________________

#if (defined USE_NETCDF3) 
  write(LIS_logunit,*) "[ERR] MERRA2 reader requires NetCDF4"
  call LIS_endrun()
#endif

#if (defined USE_NETCDF4) 
  ferror = 0 
  nr_index = merra2_struc(n)%nrold
  nc_index = merra2_struc(n)%ncold
  mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)

! Read single layer file (*slv) fields:
  inquire(file=slvname,exist=file_exists) 
  if(file_exists) then 
     write(LIS_logunit,*)'[INFO] Reading MERRA-2 file (bookend,', order,' ... ',trim(slvname)
     call LIS_verify(nf90_open(path=trim(slvname), mode=NF90_NOWRITE, &
          ncid=ftn_slv), 'nf90_open failed for slvfile in read_merra2')
     call LIS_verify(nf90_inq_varid(ftn_slv,'PS',psId), &
          'nf90_inq_varid failed for ps in read_merra2')

     call LIS_verify(nf90_get_var(ftn_slv,psId, ps), &
          'nf90_get_var failed for ps in read_merra2') 
     
   ! If using the MERRA2 lowest model level forcing (*flx): 
     if(merra2_struc(n)%uselml.eq.1) then 
        inquire(file=flxname,exist=file_exists1)
        if(.not.file_exists1) then 
           write(LIS_logunit,*) '[ERR] ',trim(flxname)//' does not exist'
           call LIS_endrun()
        endif
        write(LIS_logunit,*) '[INFO] Reading MERRA-2 file (bookend,',order,' ... '
        write(LIS_logunit,*)  trim(flxname), &
            ' (for lowest model level fields)'
        call LIS_verify(nf90_open(path=trim(flxname), mode=NF90_NOWRITE, &
             ncid=ftn_flx), 'nf90_open failed for flxfile in read_merra2')

        call LIS_verify(nf90_inq_varid(ftn_flx,'TLML',tmpId), &
             'nf90_inq_varid failed for tlml in read_merra2')
        call LIS_verify(nf90_inq_varid(ftn_flx,'QLML',qId), &
             'nf90_inq_varid failed for qlml in read_merra2')
        call LIS_verify(nf90_inq_varid(ftn_flx,'ULML',uId), &
             'nf90_inq_varid failed for ulml in read_merra2')
        call LIS_verify(nf90_inq_varid(ftn_flx,'VLML',vId), &
             'nf90_inq_varid failed for vlml in read_merra2')

        call LIS_verify(nf90_get_var(ftn_flx,tmpId, tair), &
             'nf90_get_var failed for tlml in read_merra2')
        call LIS_verify(nf90_get_var(ftn_flx,qId, qair), &
             'nf90_get_var failed for qlml in read_merra2')
        call LIS_verify(nf90_get_var(ftn_flx,uId, uwind), &
             'nf90_get_var failed for ulml in read_merra2')
        call LIS_verify(nf90_get_var(ftn_flx,vId, vwind), &
             'nf90_get_var failed for vlml in read_merra2')
        call LIS_verify(nf90_close(ftn_flx), &
             'failed to close flxfile in read_merra2')

   ! Else use the single-layer fields (e.g., 2m ref. height):
     else
        call LIS_verify(nf90_inq_varid(ftn_slv,'T2M',tmpId), &
             'nf90_inq_varid failed for t2m in read_merra2')
        call LIS_verify(nf90_inq_varid(ftn_slv,'QV2M',qId), &
             'nf90_inq_varid failed for qv2m in read_merra2')

        call LIS_verify(nf90_get_var(ftn_slv,tmpId, tair), &
             'nf90_get_var failed for t2ml in read_merra2')
        call LIS_verify(nf90_get_var(ftn_slv,qId, qair), &
             'nf90_get_var failed for qv2m in read_merra2')

        if(merra2_struc(n)%use2mwind.eq.1) then
            call LIS_verify(nf90_inq_varid(ftn_slv,'U2M',uId), &
                 'nf90_inq_varid failed for u2m in read_merra2')
            call LIS_verify(nf90_inq_varid(ftn_slv,'V2M',vId), &
                 'nf90_inq_varid failed for v2m in read_merra2')
            call LIS_verify(nf90_get_var(ftn_slv,uId, uwind), &
                 'nf90_get_var failed for u2m in read_merra2')
            call LIS_verify(nf90_get_var(ftn_slv,vId, vwind), &
                 'nf90_get_var failed for v2m in read_merra2')
        else
            call LIS_verify(nf90_inq_varid(ftn_slv,'U10M',uId), & 
                 'nf90_inq_varid failed for u10m in read_merra2')
            call LIS_verify(nf90_inq_varid(ftn_slv,'V10M',vId), &
                 'nf90_inq_varid failed for v10m in read_merra2')
            call LIS_verify(nf90_get_var(ftn_slv,uId, uwind), &
                 'nf90_get_var failed for u10m in read_merra2')
            call LIS_verify(nf90_get_var(ftn_slv,vId, vwind), &
                 'nf90_get_var failed for v10m in read_merra2')
        endif
           
     endif

     call LIS_verify(nf90_close(ftn_slv), &
          'failed to close slvfile in read_merra2')

     call interp_merra2_var(n,findex,month,tair,  1, .false., merraforc)
     call interp_merra2_var(n,findex,month,qair,  2, .false., merraforc)
     call interp_merra2_var(n,findex,month,uwind, 5, .false., merraforc)
     call interp_merra2_var(n,findex,month,vwind, 6, .false., merraforc)
     call interp_merra2_var(n,findex,month,ps,    7, .false., merraforc)
  else
     write(LIS_logunit,*) '[ERR] ',trim(slvname)//' does not exist'
     call LIS_endrun()
  endif

! Read in the flux file fields (*flx):
  inquire(file=flxname,exist=file_exists) 
  if(file_exists) then 

     write(LIS_logunit,*) '[INFO] Reading MERRA-2 file (bookend,',order,' ... '
     write(LIS_logunit,*) trim(flxname),' (for corrected precipitation fields)'
     call LIS_verify(nf90_open(path=trim(flxname), mode=NF90_NOWRITE, &
          ncid=ftn_flx),'nf90_open failed for flxfile in read_merra2')
     
   ! Read in the *corrected* total precipitation field:
     if(merra2_struc(n)%usecorr.eq.1) then 
        call LIS_verify(nf90_inq_varid(ftn_flx,'PRECTOTCORR',prectotId), &
             'nf90_inq_varid failed for prectotcorr (flx) in read_merra2')

        call LIS_verify(nf90_get_var(ftn_flx,prectotId, prectot), &
             'nf90_get_var failed for prectotcorr (flx) in read_merra2')

     else  ! Original total precip:
        call LIS_verify(nf90_inq_varid(ftn_flx,'PRECTOT',prectotId), &
             'nf90_inq_varid failed for prectot (flx) in read_merra2')

        call LIS_verify(nf90_get_var(ftn_flx,prectotId, prectot), &
             'nf90_get_var failed for prectot (flx) in read_merra2')
     endif

     call interp_merra2_var(n,findex,month,prectot,  8, .true.,merraforc)

   ! Read in the convective precipitation, if selected in the forcing table:
     if( LIS_FORC_CRainf%selectOpt.eq.1 .and. &
          merra2_struc(n)%usecorr.eq.0 ) then
        call LIS_verify(nf90_inq_varid(ftn_flx,'PRECCON',precconId), &
             'nf90_inq_varid failed for preccon (flx) in read_merra2')

        call LIS_verify(nf90_get_var(ftn_flx,precconId, preccon), &
             'nf90_get_var failed for preccon (flx) in read_merra2')

        call interp_merra2_var(n,findex,month,preccon,  9, .true.,merraforc)
     endif

   ! Read in the snowfall field, if selected in the forcing input table:
     if( LIS_FORC_Snowf%selectOpt.eq.1 .and. &
          merra2_struc(n)%usecorr.eq.0 ) then
        call LIS_verify(nf90_inq_varid(ftn_flx,'PRECSNO',precsnoId), &
             'nf90_inq_varid failed for precsno (flx) in read_merra2')  

        call LIS_verify(nf90_get_var(ftn_flx,precsnoId, precsno), &
            'nf90_get_var failed for precsno (flx) in read_merra2')

        call interp_merra2_var(n,findex,month,precsno,  10, .true.,merraforc)
     endif

   ! Read in Forcing Height, if selected in the forcing input table:
     if (LIS_FORC_Forc_Hgt%selectOpt.eq.1) then
        call LIS_verify(nf90_inq_varid(ftn_flx,'HLML',hlmlId), &
             'nf90_inq_varid failed for hlml (flx) in read_merra2')

        call LIS_verify(nf90_get_var(ftn_flx,hlmlId, hlml), &
             'nf90_get_var failed for hlml (flx) in read_merra2')

        call interp_merra2_var(n,findex,month,hlml,  14, .false.,merraforc)
     endif

     call LIS_verify(nf90_close(ftn_flx), &
          'failed to close flxfile in read_merra2')
     
  else
     write(LIS_logunit,*) '[ERR] ',trim(flxname)//' does not exist'
     call LIS_endrun()
  endif
  
! Read in the radiation file fields (*rad):
  inquire(file=radname,exist=file_exists) 
  if(file_exists) then 
     write(LIS_logunit,*) '[INFO] Reading MERRA-2 file (bookend,',order,' ... ',trim(radname)
     call LIS_verify(nf90_open(path=trim(radname), mode=NF90_NOWRITE, &
          ncid=ftn_rad), 'nf90_open failed in read_merra2')
     call LIS_verify(nf90_inq_varid(ftn_rad,'SWGDN',swgdnId), &
          'nf90_inq_varid failed for swgdn in read_merra2')
     call LIS_verify(nf90_inq_varid(ftn_rad,'LWGAB',lwgabId), &
          'nf90_inq_varid failed for lwgab in read_merra2')
!     call LIS_verify(nf90_inq_varid(ftn_rad,'emis',emisId), &
!          'nf90_inq_varid failed for emis in read_merra2')
    
     call LIS_verify(nf90_get_var(ftn_rad,swgdnId,swgdn), &
          'nf90_get_var failed for swgdn in read_merra2')
     call LIS_verify(nf90_get_var(ftn_rad,lwgabId,lwgab), &
          'nf90_get_var failed for lwgab in read_merra2')
!     call LIS_verify(nf90_get_var(ftn_rad,emisId,emis), &
!          'nf90_get_var failed for emis in read_merra2')

     call interp_merra2_var(n,findex,month,swgdn, 3, .false.,merraforc)
     call interp_merra2_var(n,findex,month,lwgab, 4, .false.,merraforc)

     call LIS_verify(nf90_close(ftn_rad), &
          'failed to close rad file in read_merra2')
  else
     write(LIS_logunit,*) '[ERR] ',trim(radname)//' does not exist'
     call LIS_endrun()
  endif

! Checks: For reading in the surface layer file fields (*lfo):
  read_lnd = .false. 
  if ( LIS_FORC_Pardr%selectOpt.eq.1.or.&
       LIS_FORC_Pardf%selectOpt.eq.1.or.&
       LIS_FORC_SWnet%selectOpt.eq.1 ) then
     read_lnd = .true.
  endif
  if( merra2_struc(n)%usecorr.eq.1 ) then
    if( LIS_FORC_CRainf%selectOpt.eq.1 .or. &
       LIS_FORC_Snowf%selectOpt.eq.1 ) then
      read_lnd = .true.
    endif
  endif

! Read in the surface layer fields (*lfo):
  if(read_lnd) then 
     inquire(file=lfoname,exist=file_exists) 
     if(file_exists) then 
       write(LIS_logunit,*) '[INFO] Reading MERRA-2 file (bookend,',order,&
             ' ... ',trim(lfoname)
       call LIS_verify(nf90_open(path=trim(lfoname), mode=NF90_NOWRITE, &
            ncid=ftn_lfo), 'nf90_open failed in read_merra2')

       ! Read *corrected* convective precipitation, if selected: 
       if( LIS_FORC_CRainf%selectOpt.eq.1 .and. &
          merra2_struc(n)%usecorr.eq.1 ) then
         call LIS_verify(nf90_inq_varid(ftn_lfo,'PRECCUCORR',precconId), &
             'nf90_inq_varid failed for preccucorr (lfo) in read_merra2')

         call LIS_verify(nf90_get_var(ftn_lfo,precconId, preccon), &
             'nf90_get_var failed for preccucorr (lfo) in read_merra2')

         call interp_merra2_var(n,findex,month,preccon, 9,.true.,merraforc)
       endif

       ! Read *corrected* snowfall field, if selected: 
       if( LIS_FORC_Snowf%selectOpt.eq.1 .and. &
          merra2_struc(n)%usecorr.eq.1 ) then
         call LIS_verify(nf90_inq_varid(ftn_flx,'PRECSNOCORR',precsnoId), &
             'nf90_inq_varid failed for precsnocorr (lfo) in read_merra2')

         call LIS_verify(nf90_get_var(ftn_flx,precsnoId, precsno), &
             'nf90_get_var failed for precsnocorr (lfo) in read_merra2')

         call interp_merra2_var(n,findex,month,precsno, 10,.true.,merraforc)
       endif

       if(LIS_FORC_SWnet%selectOpt.eq.1) then 
          call LIS_verify(nf90_inq_varid(ftn_lfo,'SWLAND',swlandId), &
               'nf90_inq_varid failed for swland in read_merra2')

          call LIS_verify(nf90_get_var(ftn_lfo,swlandId,swland), &
               'nf90_get_var failed for swland in read_merra2')
          
          call interp_merra2_var(n,findex,month,swland,11,.false.,merraforc)
       endif
        
        if(LIS_FORC_Pardr%selectOpt.eq.1) then 
           call LIS_verify(nf90_inq_varid(ftn_lfo,'PARDR',pardrId), &
                'nf90_inq_varid failed for pardr in read_merra2')

           call LIS_verify(nf90_get_var(ftn_lfo,pardrId,pardr), &
                'nf90_get_var failed for pardr in read_merra2')
           
           call interp_merra2_var(n,findex,month,pardr,12,.false.,merraforc)
        endif

        if(LIS_FORC_Pardf%selectOpt.eq.1) then 
           call LIS_verify(nf90_inq_varid(ftn_lfo,'PARDF',pardfId), &
                'nf90_inq_varid failed for pardf in read_merra2')

           call LIS_verify(nf90_get_var(ftn_lfo,pardfId,pardf), &
                'nf90_get_var failed for pardf in read_merra2')
           call interp_merra2_var(n,findex,month,pardf,13,.false.,merraforc)
           
        endif

        call LIS_verify(nf90_close(ftn_lfo),&
             'nf90_close failed for lfofile in read_merra2')
     else
        write(LIS_logunit,*) '[ERR] ',trim(lfoname)//' does not exist'
        call LIS_endrun()
     endif

  endif
  
#endif
end subroutine read_merra2

!BOP
! 
! !ROUTINE: interp_merra2_var
! \label{interp_merra2_var}
! 
! !INTERFACE: 
subroutine interp_merra2_var(n,findex,month, input_var,  var_index, &
     pcp_flag, merraforc)

! !USES: 
  use LIS_coreMod
  use LIS_logMod
  use LIS_spatialDownscalingMod
  use merra2_forcingMod, only : merra2_struc
#if(defined USE_NETCDF3 || defined USE_NETCDF4)      
  use netcdf
#endif
  implicit none

! !ARGUMENTS:
  integer, intent(in)    :: n
  integer, intent(in)    :: findex
  integer, intent(in)    :: month
  real,    intent(in)    :: input_var(merra2_struc(n)%ncold, &
       merra2_struc(n)%nrold, 24)
  integer, intent(in)    :: var_index
  logical, intent(in)    :: pcp_flag
  real,    intent(inout) :: merraforc(merra2_struc(n)%nvars, &
       24, LIS_rc%lnc(n)*LIS_rc%lnr(n))
  !
! !DESCRIPTION: 
!  This subroutine spatially interpolates a MERRA2 field
!  to the LIS running domain
! 
!EOP

  integer   :: t,c,r,k,iret
  integer   :: doy
  integer   :: ftn
  integer   :: pcp1Id, pcp2Id, pcp3Id, pcp4Id,pcp5Id, pcp6Id
  real      :: f (merra2_struc(n)%ncold*merra2_struc(n)%nrold)
  logical*1 :: lb(merra2_struc(n)%ncold*merra2_struc(n)%nrold)
  logical*1 :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  integer   :: input_size
  logical   :: scal_read_flag
! _____________________________________________________________

  input_size = merra2_struc(n)%ncold*merra2_struc(n)%nrold

!-----------------------------------------------------------------------    
! Apply corrections
!-----------------------------------------------------------------------  
  
  scal_read_flag = .false. 

  if(merra2_struc(n)%pcpscal_cmo.ne.LIS_rc%mo) then 
     scal_read_flag = .true. 
  endif
     
  if ( pcp_flag .and. merra2_struc(n)%usescalef==1.and.&
       scal_read_flag ) then 
     
!     call finddoy(doy, LIS_rc%yr, LIS_rc%mo, LIS_rc%da)
     
#if(defined USE_NETCDF3 || defined USE_NETCDF4)      
     call LIS_verify(nf90_open(path=merra2_struc(n)%scaleffile,&
          mode=nf90_nowrite,ncid=ftn),&
          'failed to open MERRA2 precip scaling factor input file')
     
     call LIS_verify(nf90_inq_varid(ftn,"REF_XRANGE",pcp1id),&
          'nf90_inq_varid failed for REF_XRANGE')
     call LIS_verify(nf90_get_var(ftn,pcp1id,merra2_struc(n)%refxrange,&
          start=(/1,1,LIS_rc%mo,1/),&
          count=(/merra2_struc(n)%ncold,merra2_struc(n)%nrold,1,&
          merra2_struc(n)%nbins/)),&
          'nf90_get_var failed for REF_XRANGE')
     
     call LIS_verify(nf90_inq_varid(ftn,"REF_CDF",pcp2id),&
          'nf90_inq_varid failed for REF_CDF')
     call LIS_verify(nf90_get_var(ftn,pcp2id,merra2_struc(n)%refcdf,&
          start=(/1,1,LIS_rc%mo,1/),&
          count=(/merra2_struc(n)%ncold,merra2_struc(n)%nrold,1,&
          merra2_struc(n)%nbins/)),&
          'nf90_get_var failed for REF_CDF')
     
     call LIS_verify(nf90_inq_varid(ftn,"MERRA2_XRANGE",pcp3id),&
          'nf90_inq_varid failed for MERRA2_XRANGE')
     call LIS_verify(nf90_get_var(ftn,pcp3id,merra2_struc(n)%merraxrange,&
          start=(/1,1,LIS_rc%mo,1/),&
          count=(/merra2_struc(n)%ncold,merra2_struc(n)%nrold,1,&
          merra2_struc(n)%nbins/)),&
          'nf90_get_var failed for MERRA2_XRANGE')
     
     call LIS_verify(nf90_inq_varid(ftn,"MERRA2_CDF",pcp4id),&
          'nf90_inq_varid failed for MERRA2_CDF')
     call LIS_verify(nf90_get_var(ftn,pcp4id,merra2_struc(n)%merracdf,&
          start=(/1,1,LIS_rc%mo,1/),&
          count=(/merra2_struc(n)%ncold,merra2_struc(n)%nrold,1,&
          merra2_struc(n)%nbins/)),&
          'nf90_get_var failed for MERRA2_CDF')

     if(merra2_struc(n)%usepcpsampling.gt.0) then 
        !REF MEAN, STDEV
        call LIS_verify(nf90_inq_varid(ftn,"REF_MEAN",pcp5id),&
             'nf90_inq_varid failed for REF_MEAN')
        call LIS_verify(nf90_get_var(ftn,pcp5id,merra2_struc(n)%refmean,&
             start=(/1,1,LIS_rc%mo/),&
             count=(/merra2_struc(n)%ncold,merra2_struc(n)%nrold,1/)),&
             'nf90_get_var failed for REF_MEAN')

        call LIS_verify(nf90_inq_varid(ftn,"REF_STDEV",pcp6id),&
             'nf90_inq_varid failed for REF_STDEV')
        call LIS_verify(nf90_get_var(ftn,pcp6id,merra2_struc(n)%refstdev,&
             start=(/1,1,LIS_rc%mo/),&
             count=(/merra2_struc(n)%ncold,merra2_struc(n)%nrold,1/)),&
             'nf90_get_var failed for REF_STDEV')
     endif
#endif
  endif

  do t=1,24
     lb = .true.
     do r=1,merra2_struc(n)%nrold
        do c=1,merra2_struc(n)%ncold
           k= c+(r-1)*merra2_struc(n)%ncold
           f(k) = input_var(c,r,t)
           if ( f(k) == 1.e+15 ) then 
              f(k)  = LIS_rc%udef
              lb(k) = .false. 
           endif
        enddo
     enddo
     if ( pcp_flag .and. merra2_struc(n)%usescalef==1 ) then 

        call rescaleWithCDFmatching(&
             merra2_struc(n)%ncold,&
             merra2_struc(n)%nrold,&
             merra2_struc(n)%nbins,&
             merra2_struc(n)%refxrange,&
             merra2_struc(n)%merraxrange,&
             merra2_struc(n)%refcdf,&
             merra2_struc(n)%merracdf,&
             f)

     endif
!-----------------------------------------------------------------------    
! Apply downscaling
!-----------------------------------------------------------------------    
     
     if(pcp_flag.and.LIS_rc%pcp_downscale(findex).ne.0) then 
        !input_data becomes the ratio field. 
        call LIS_generatePcpClimoRatioField(n,findex,"MERRA2",&
             month, & 
             input_size, &
             f, &
             lb)     
     endif
          
     if(pcp_flag.and.&
         trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 

        call conserv_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
             merraforc(var_index,t,:), &
             merra2_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),& 
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             merra2_struc(n)%w112,merra2_struc(n)%w122,&
             merra2_struc(n)%w212,merra2_struc(n)%w222,&
             merra2_struc(n)%n112,merra2_struc(n)%n122,&
             merra2_struc(n)%n212,merra2_struc(n)%n222,&
             LIS_rc%udef, iret)

     elseif(trim(LIS_rc%met_interp(findex)).eq."bilinear".or.&
             trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 
        call bilinear_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
             merraforc(var_index,t,:), &
             merra2_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n), & 
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             merra2_struc(n)%w111,merra2_struc(n)%w121,&
             merra2_struc(n)%w211,merra2_struc(n)%w221,&
             merra2_struc(n)%n111,merra2_struc(n)%n121,&
             merra2_struc(n)%n211,merra2_struc(n)%n221,&
             LIS_rc%udef, iret)

     elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then 
        call neighbor_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
             merraforc(var_index,t,:),merra2_struc(n)%mi,&
             LIS_rc%lnc(n)*LIS_rc%lnr(n),&
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             merra2_struc(n)%n113,LIS_rc%udef,iret)
     else
        write(LIS_logunit,*) '[ERR] Spatial interpolation option '//&
                             trim(LIS_rc%met_interp(findex))//&
                             ' not supported for MERRA2'
        call LIS_endrun()
     endif

!Interpolate the refmean fields
     if( pcp_flag .and. scal_read_flag.and. &
          merra2_struc(n)%usescalef==1.and. &
          merra2_struc(n)%usepcpsampling.eq.1) then 

        lb = .false. 
        f = -9999.0
        do r=1,merra2_struc(n)%nrold
           do c=1,merra2_struc(n)%ncold
              k= c+(r-1)*merra2_struc(n)%ncold
              f(k) = merra2_struc(n)%refmean(c,r,1)
              if(f(k).ne.-9999.0) then 
                 lb(k) = .true.
              endif
           enddo
        enddo
        
        if(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 

           call conserv_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
                merra2_struc(n)%refmean_ip, &
                merra2_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),& 
                LIS_domain(n)%lat, LIS_domain(n)%lon,&
                merra2_struc(n)%w112,merra2_struc(n)%w122,&
                merra2_struc(n)%w212,merra2_struc(n)%w222,&
                merra2_struc(n)%n112,merra2_struc(n)%n122,&
                merra2_struc(n)%n212,merra2_struc(n)%n222,&
                LIS_rc%udef, iret)
           
        elseif(trim(LIS_rc%met_interp(findex)).eq."bilinear".or.&
             trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then

           call bilinear_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
                merra2_struc(n)%refmean_ip, &
                merra2_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n), & 
                LIS_domain(n)%lat, LIS_domain(n)%lon,&
                merra2_struc(n)%w111,merra2_struc(n)%w121,&
                merra2_struc(n)%w211,merra2_struc(n)%w221,&
                merra2_struc(n)%n111,merra2_struc(n)%n121,&
                merra2_struc(n)%n211,merra2_struc(n)%n221,&
                LIS_rc%udef, iret)

        elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then 
           call neighbor_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
                merra2_struc(n)%refmean_ip,merra2_struc(n)%mi,&
                LIS_rc%lnc(n)*LIS_rc%lnr(n),&
                LIS_domain(n)%lat, LIS_domain(n)%lon,&
                merra2_struc(n)%n113,LIS_rc%udef,iret)
        endif
     endif

!Interpolate the refstdev fields
     if( pcp_flag .and. scal_read_flag.and. &
          merra2_struc(n)%usescalef==1.and. &
          merra2_struc(n)%usepcpsampling.eq.1) then 

        lb = .false. 
        f = -9999.0
        do r=1,merra2_struc(n)%nrold
           do c=1,merra2_struc(n)%ncold
              k= c+(r-1)*merra2_struc(n)%ncold
              f(k) = merra2_struc(n)%refstdev(c,r,1)
              if(f(k).ne.-9999.0) then 
                 lb(k) = .true.
              endif
           enddo
        enddo
        
        if(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 

           call conserv_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
                merra2_struc(n)%refstdev_ip, &
                merra2_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),& 
                LIS_domain(n)%lat, LIS_domain(n)%lon,&
                merra2_struc(n)%w112,merra2_struc(n)%w122,&
                merra2_struc(n)%w212,merra2_struc(n)%w222,&
                merra2_struc(n)%n112,merra2_struc(n)%n122,&
                merra2_struc(n)%n212,merra2_struc(n)%n222,&
                LIS_rc%udef, iret)
           
        elseif(trim(LIS_rc%met_interp(findex)).eq."bilinear".or.&
             trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then

           call bilinear_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
                merra2_struc(n)%refstdev_ip, &
                merra2_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n), & 
                LIS_domain(n)%lat, LIS_domain(n)%lon,&
                merra2_struc(n)%w111,merra2_struc(n)%w121,&
                merra2_struc(n)%w211,merra2_struc(n)%w221,&
                merra2_struc(n)%n111,merra2_struc(n)%n121,&
                merra2_struc(n)%n211,merra2_struc(n)%n221,&
                LIS_rc%udef, iret)

        elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then 
           call neighbor_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
                merra2_struc(n)%refstdev_ip,merra2_struc(n)%mi,&
                LIS_rc%lnc(n)*LIS_rc%lnr(n),&
                LIS_domain(n)%lat, LIS_domain(n)%lon,&
                merra2_struc(n)%n113,LIS_rc%udef,iret)
        endif
     endif

     if( pcp_flag.and.LIS_rc%pcp_downscale(findex).ne.0 ) then 

        call LIS_pcpClimoDownscaling(n, findex, month,&
             LIS_rc%lnc(n)*LIS_rc%lnr(n), merraforc(var_index,t,:), lo)
        
     endif
     
  enddo
  
end subroutine interp_merra2_var

#if 0 
!BOP
!
! !ROUTINE: finddoy
! \label{finddoy} 
! 
! !INTERFACE:
  subroutine finddoy(doy,yr,mo,da)

    implicit none
! !ARGUMENTS:
    integer,intent(in)   :: yr,mo,da
    integer,intent(out)  :: doy
!
! !DESCRIPTION:
! 
!  Determines the time, time in GMT, and the day of the year
!  based on the value of year, month, day of month, hour of 
!  the day, minute and second. This method is the inverse of 
!  time2date.
!
!   NOTE: This routine has been known to give round off error
!   problems when attempting to retrieving minutes and seconds
!   from the given number. Use at your own risk!
! 
!  The arguments are: 
!  \begin{description}
!  \item[yr]
!    year
!  \item[mo]
!   month
!  \item[da]
!   day of the month
!  \item[hr]
!   hour of day
!  \item[mn]
!   minute
!  \item[ss]
!   second
!  \item[time]
!   lis time
!  \item[gmt]
!   time in GMT
!  \item[doy]
!   day of the year
!  \end{description}
!EOP
    integer :: yrdays,days(13),k
!    data days /31,28,31,30,31,30,31,31,30,31,30,31,30/
    data days /31,29,31,30,31,30,31,31,30,31,30,31,30/

    if((mod(yr,4).eq.0.and.mod(yr,100).ne.0) &     !correct for leap year
         .or.(mod(yr,400).eq.0))then             !correct for y2k
       yrdays=366                  
    else
       yrdays=365
    endif
    
    doy=0
    do k=1,(mo-1)
       doy=doy+days(k)
    enddo
    doy=doy+da
    
!    if(yrdays.eq.366.and.mo.gt.2)doy=doy+1
    
    return
  end subroutine finddoy
#endif

!BOP
! 
! !ROUTINE: rescaleWithCDFmatching
! \label{rescaleWithCDFmatching}
!
! !INTERFACE:
  subroutine rescaleWithCDFmatching(&
       nc,            & 
       nr,            & 
       nbins,         & 
       ref_xrange,    & 
       merra_xrange,  &       
       ref_cdf,       &
       merra_cdf,     &
       out_value)

    use LIS_logMod
       
    implicit none
! 
! !ARGUMENTS: 
    integer             :: nc
    integer             :: nr
    integer             :: nbins
    real                :: ref_xrange(nc,nr,1,nbins)
    real                :: merra_xrange(nc,nr,1,nbins)
    real                :: ref_cdf(nc,nr,1,nbins)
    real                :: merra_cdf(nc,nr,1,nbins)
    real                :: out_value(nc*nr)
!
! !DESCRIPTION: 
! 
!   This routine rescales the input merraervation data to the ref's
!   climatology so that the cumulative distribution functions (CDFs)
!   of the merraervations and the ref match (for each grid point). 
! 
!   Ref: Reichle and Koster, 2004, Bias reduction in short records of 
!   satellite soil moisture, Geophys. Res. Lett. 31, L19501, 
!   doi:10.1029/2004GL020938. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]               index of the nest
!  \item[nbins]           number of bins used to compute the ref and merra CDFs
!  \item[max\_merra\_value] maximum allowable value of merraervation
!  \item[min\_merra\_value] minimum allowable value of merraervation
!  \item[ref\_xrange]   x-axis values corresponding to the ref CDF
!  \item[merra\_xrange]     x-axis values corresponding to the merra CDF
!  \item[ref\_cdf]      y-axis (CDF) values corresponding to the ref CDF
!  \item[merra\_cdf]        y-axis (CDF) values corresponding to the merra CDF
!  \item[merra\_value]      merraervation value to be rescaled. 
! \end{description}
!EOP

    real                 :: max_merra_value
    real                 :: min_merra_value

    real                 :: ref_delta(nc,nr)
    real                 :: merra_delta(nc,nr)

    integer              :: c,r,i
    integer              :: binval
    integer              :: col,row
    real                 :: cdf_merraval
    real                 :: merra_tmp, merra_in
    integer,dimension(1) :: index_25 , index_75
    real                 :: Lb_xrange, Ub_xrange, iqr_merra, iqr_ref

    ref_delta   = 0 
    merra_delta = 0
   
    do r=1,nr
       do c=1,nc
          if(ref_xrange(c,r,1,2).ne.-9999.0.and.&
               ref_xrange(c,r,1,1).ne.-9999.0.and.&
               merra_xrange(c,r,1,2).ne.-9999.0.and.&
               merra_xrange(c,r,1,1).ne.-9999.0) then 
             ref_delta(c,r) = ref_xrange(c,r,1,2)-ref_xrange(c,r,1,1)
             merra_delta(c,r)   = merra_xrange(c,r,1,2)-merra_xrange(c,r,1,1)
          endif
       enddo
    enddo
    do r=1,nr
       do c=1,nc
          if(ref_delta(c,r).gt.1E-14.and.merra_delta(c,r).gt.1E-14) then 
             min_merra_value = ref_xrange(c,r,1,1)
             max_merra_value = ref_xrange(c,r,1,nbins)
             index_25 = minloc(abs(merra_cdf(c,r,1,:) - 0.25))
             index_75 = minloc(abs(merra_cdf(c,r,1,:) - 0.75))
             Lb_xrange = merra_xrange(c,r,1,index_25(1))
             Ub_xrange = merra_xrange(c,r,1,index_75(1))
             iqr_merra = Ub_xrange - Lb_xrange
             
             index_25 = minloc(abs(ref_cdf(c,r,1,:) - 0.25))
             index_75 = minloc(abs(ref_cdf(c,r,1,:) - 0.75))
             Lb_xrange = ref_xrange(c,r,1,index_25(1))
             Ub_xrange = ref_xrange(c,r,1,index_75(1))
             iqr_ref = Ub_xrange - Lb_xrange
             
! In a normal distribution 50% of data will fall between +/- 0.67448 sigma (1.134896 sigma)
! and 99.7% of data are between the +/-3 sigma (6 sigma) By dividing these two values we
! get 1.134896 / 6 = 0.189. That means IQR is about the 0.19 of the 
! x-range (~ dynamic range of the variable).
! We can say if IQR is less than 0.05 of the x-range then CDF is too 
! steep and it is better to ignore that for CDF matching.  
! NOTE: more tests are needed to determine the best threshold

!             if( iqr_merra .lt. 0.05 * (merra_xrange(c,r,1,nbins)-merra_xrange(c,r,1,1)) .or. &
!                  iqr_ref .lt. 0.05 * (ref_xrange(c,r,1,nbins)-ref_xrange(c,r,1,1)) ) then 
!                merra_delta(c,r) = 0
!             endif
             if(out_value(c+(r-1)*nc).lt.merra_xrange(c,r,1,nbins)) then 
             
                if(out_value(c+(r-1)*nc).ne.-9999.0) then 
                   merra_in = out_value(c+(r-1)*nc)
!                print*,c,r,merra_in,merra_delta(c,r)
                   if(merra_in.gt.1E-14) then 
!                      print*, c,r,out_value(c+(r-1)*nc),merra_xrange(c,r,1,1),&
!                           merra_xrange(c,r,1,nbins),&
!                           merra_delta(c,r)                   
                      binval = nint((out_value(c+(r-1)*nc)-merra_xrange(c,r,1,1))/&
                           merra_delta(c,r))+1
                      
                      if(binval.gt.nbins) binval = nbins
                      if(binval.le.0) binval = 1
                      cdf_merraval = merra_cdf(c,r,1,binval)
                      if(cdf_merraval.gt.1.0) cdf_merraval = 1.0
                      i=1
                      do while((ref_cdf(c,r,1,i).lt.cdf_merraval).and.&
                           (i.le.nbins))
                         i = i+1
                         if(i.gt.nbins) exit
                      enddo
                      if(i.gt.nbins) i = i-1
                      merra_tmp = ref_xrange(c,r,1,i)
                      
                      if(merra_tmp.gt.max_merra_value) then 
                         !                      print*, 'max const',merra_tmp
                         merra_tmp = -9999.0
                      endif
                      
                      if(merra_tmp.lt.min_merra_value) then 
                         merra_tmp = -9999.0
                      endif
                      out_value(c+(r-1)*nc) = merra_tmp
                   endif
                   if(out_value(c+(r-1)*nc).lt.min_merra_value.and.&
                        out_value(c+(r-1)*nc).ne.-9999.0) then 
                      write(LIS_logunit,*) '[ERR] Problem in CDF scaling '
                      call LIS_endrun()
                   endif
                else
                   out_value(c+(r-1)*nc) = -9999.0
                endif
                if(merra_in.ne.-9999.0.and.out_value(c+(r-1)*nc).eq.-9999.0) then 
                   print*,'mismatch', c,r,merra_in, out_value(c+(r-1)*nc)
                   stop
                endif
             endif
          endif
    enddo
 enddo
end subroutine rescaleWithCDFmatching

