!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: read_merra2
! \label{read_merra2}
! 
! !REVISION HISTORY:
! 18 Mar 2015: James Geiger, initial code (based on merra-land)
! 12 Nov 2015: KR Arsenault, added to LDT
!
! !INTERFACE:
subroutine read_merra2(n, order, findex,          &
                       slvname, flxname, lfoname, radname, &
                       merraforc, ferror)
! !USES:
  use LDT_coreMod,       only : LDT_rc, LDT_domain, LDT_masterproc
  use LDT_logMod
  use LDT_FORC_AttributesMod
  use LDT_metforcingMod, only : LDT_forc
  use merra2_forcingMod, only : merra2_struc
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in)          :: n
  integer, intent(in)          :: order
  integer, intent(in)          :: findex
  character(len=*), intent(in) :: slvname
  character(len=*), intent(in) :: flxname
  character(len=*), intent(in) :: lfoname
  character(len=*), intent(in) :: radname
  real, intent(inout)          :: merraforc(LDT_rc%met_nf(findex), 24, &
                                            LDT_rc%lnc(n)*LDT_rc%lnr(n))
  integer, intent(out)         :: ferror          

!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  MERRA2 data, transforms into 9 LDT forcing 
!  parameters and interpolates to the LDT domain. \newline
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
!  9. precon  Convective precipatation [$mm/s$] \newline
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

  real      :: tair(merra2_struc(n)%nc, merra2_struc(n)%nr,24)
  real      :: qair(merra2_struc(n)%nc, merra2_struc(n)%nr,24)
  real      :: uwind(merra2_struc(n)%nc, merra2_struc(n)%nr,24)
  real      :: vwind(merra2_struc(n)%nc, merra2_struc(n)%nr,24)
  real      :: ps(merra2_struc(n)%nc, merra2_struc(n)%nr,24)
  real      :: prectot(merra2_struc(n)%nc, merra2_struc(n)%nr,24)
  real      :: precsno(merra2_struc(n)%nc, merra2_struc(n)%nr,24)
  real      :: preccon(merra2_struc(n)%nc, merra2_struc(n)%nr,24)
  real      :: prectot_flx(merra2_struc(n)%nc, merra2_struc(n)%nr,24)
  real      :: preccon_flx(merra2_struc(n)%nc, merra2_struc(n)%nr,24)
  real      :: swgdn(merra2_struc(n)%nc, merra2_struc(n)%nr,24)
  real      :: lwgab(merra2_struc(n)%nc, merra2_struc(n)%nr,24)
  real      :: swland(merra2_struc(n)%nc, merra2_struc(n)%nr,24)
  real      :: pardr(merra2_struc(n)%nc, merra2_struc(n)%nr,24)
  real      :: pardf(merra2_struc(n)%nc, merra2_struc(n)%nr,24)
  real      :: hlml(merra2_struc(n)%nc, merra2_struc(n)%nr,24)
!  real      :: emis(merra2_struc(n)%nc, merra2_struc(n)%nr,24)
! __________________________________________________________________________

#if (defined USE_NETCDF3) 
  write(LDT_logunit,*) "[ERR] MERRA2 reader requires NetCDF4"
  call LDT_endrun()
#endif

#if (defined USE_NETCDF4) 
  ferror = 0 
  nr_index = merra2_struc(n)%nr
  nc_index = merra2_struc(n)%nc
  mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)

! Read single layer file (*slv) fields:
  inquire(file=slvname,exist=file_exists) 
  if(file_exists) then 
     write(LDT_logunit,*) '[INFO] Reading MERRA-2 file (bookend,',order,' ... ',trim(slvname)
     call LDT_verify(nf90_open(path=trim(slvname), mode=NF90_NOWRITE, &
          ncid=ftn_slv), 'nf90_open failed for slvfile in read_merra2')
     call LDT_verify(nf90_inq_varid(ftn_slv,'PS',psId), &
          'nf90_inq_varid failed for ps in read_merra2')

     call LDT_verify(nf90_get_var(ftn_slv,psId, ps), &
          'nf90_get_var failed for ps in read_merra2') 
     
   ! If using the MERRA2 lowest model level forcing (*flx): 
     if(merra2_struc(n)%uselml.eq.1) then 
        inquire(file=flxname,exist=file_exists1)
        if(.not.file_exists1) then 
           write(LDT_logunit,*)'[ERR] ',trim(flxname)//' does not exist'
           call LDT_endrun()
        endif
        write(LDT_logunit,*) '[INFO] Reading MERRA-2 file (bookend,',order,' ... ' 
        write(LDT_logunit,*)  trim(flxname), &
            ' (for lowest model level fields)'
        call LDT_verify(nf90_open(path=trim(flxname), mode=NF90_NOWRITE, &
             ncid=ftn_flx), 'nf90_open failed for flxfile in read_merra2')

        call LDT_verify(nf90_inq_varid(ftn_flx,'TLML',tmpId), &
             'nf90_inq_varid failed for tlml in read_merra2')
        call LDT_verify(nf90_inq_varid(ftn_flx,'QLML',qId), &
             'nf90_inq_varid failed for qlml in read_merra2')
        call LDT_verify(nf90_inq_varid(ftn_flx,'ULML',uId), &
             'nf90_inq_varid failed for ulml in read_merra2')
        call LDT_verify(nf90_inq_varid(ftn_flx,'VLML',vId), &
             'nf90_inq_varid failed for vlml in read_merra2')

        call LDT_verify(nf90_get_var(ftn_flx,tmpId, tair), &
             'nf90_get_var failed for tlml in read_merra2')
        call LDT_verify(nf90_get_var(ftn_flx,qId, qair), &
             'nf90_get_var failed for qlml in read_merra2')
        call LDT_verify(nf90_get_var(ftn_flx,uId, uwind), &
             'nf90_get_var failed for ulml in read_merra2')
        call LDT_verify(nf90_get_var(ftn_flx,vId, vwind), &
             'nf90_get_var failed for vlml in read_merra2')
        call LDT_verify(nf90_close(ftn_flx), &
             'failed to close flxfile in read_merra2')

   ! Else use the single-layer fields (e.g., 2m ref. height):
     else
        call LDT_verify(nf90_inq_varid(ftn_slv,'T2M',tmpId), &
             'nf90_inq_varid failed for t2m in read_merra2')
        call LDT_verify(nf90_inq_varid(ftn_slv,'QV2M',qId), &
             'nf90_inq_varid failed for qv2m in read_merra2')
        call LDT_verify(nf90_inq_varid(ftn_slv,'U10M',uId), &
             'nf90_inq_varid failed for u10m in read_merra2')
        call LDT_verify(nf90_inq_varid(ftn_slv,'V10M',vId), &
             'nf90_inq_varid failed for v10m in read_merra2')

        call LDT_verify(nf90_get_var(ftn_slv,tmpId, tair), &
             'nf90_get_var failed for t2ml in read_merra2')
        call LDT_verify(nf90_get_var(ftn_slv,qId, qair), &
             'nf90_get_var failed for qv2m in read_merra2')
        call LDT_verify(nf90_get_var(ftn_slv,uId, uwind), &
             'nf90_get_var failed for u10m in read_merra2')
        call LDT_verify(nf90_get_var(ftn_slv,vId, vwind), &
             'nf90_get_var failed for v10m in read_merra2')
     endif

     call LDT_verify(nf90_close(ftn_slv), &
          'failed to close slvfile in read_merra2')

     call interp_merra2_var(n,findex,tair,  1, .false., merraforc)
     call interp_merra2_var(n,findex,qair,  2, .false., merraforc)
     call interp_merra2_var(n,findex,uwind, 5, .false., merraforc)
     call interp_merra2_var(n,findex,vwind, 6, .false., merraforc)
     call interp_merra2_var(n,findex,ps,    7, .false., merraforc)
  else
     write(LDT_logunit,*) '[ERR] ',trim(slvname)//' does not exist'
     call LDT_endrun()
  endif

! Read in the flux file fields (*flx):
  inquire(file=flxname,exist=file_exists) 
  if(file_exists) then 

     write(LDT_logunit,*) '[INFO] Reading MERRA-2 file (bookend,',order,' ... '
     write(LDT_logunit,*) trim(flxname),' (for corrected precipitation fields)'
     call LDT_verify(nf90_open(path=trim(flxname), mode=NF90_NOWRITE, &
          ncid=ftn_flx),'nf90_open failed for flxfile in read_merra2')
     
   ! Read in the corrected total precipitation fields:
     if(merra2_struc(n)%usecorr.eq.1) then 
        call LDT_verify(nf90_inq_varid(ftn_flx,'PRECTOTCORR',prectotId), &
             'nf90_inq_varid failed for prectotcorr (flx) in read_merra2')

        call LDT_verify(nf90_get_var(ftn_flx,prectotId, prectot), &
             'nf90_get_var failed for prectotcorr (flx) in read_merra2')

     else  ! Original total precip:
        call LDT_verify(nf90_inq_varid(ftn_flx,'PRECTOT',prectotId), &
             'nf90_inq_varid failed for prectot (flx) in read_merra2')

        call LDT_verify(nf90_get_var(ftn_flx,prectotId, prectot), &
             'nf90_get_var failed for prectot (flx) in read_merra2')
     endif

     call interp_merra2_var(n,findex,prectot,  8, .true.,merraforc)

   ! Read in the convective precipitation, if selected in the forcing table:
     if( LDT_FORC_CRainf%selectOpt.eq.1 .and. &
          merra2_struc(n)%usecorr.eq.0 ) then
        call LDT_verify(nf90_inq_varid(ftn_flx,'PRECCON',precconId), &
             'nf90_inq_varid failed for preccon (flx) in read_merra2')

        call LDT_verify(nf90_get_var(ftn_flx,precconId, preccon), &
             'nf90_get_var failed for preccon (flx) in read_merra2')

        call interp_merra2_var(n,findex,preccon,  9, .true.,merraforc)
     endif

   ! Read in the snowfall field, if selected in the forcing input table:
     if( LDT_FORC_Snowf%selectOpt.eq.1 .and. &
          merra2_struc(n)%usecorr.eq.0 ) then
        call LDT_verify(nf90_inq_varid(ftn_flx,'PRECSNO',precsnoId), &
             'nf90_inq_varid failed for precsno (flx) in read_merra2')  

        call LDT_verify(nf90_get_var(ftn_flx,precsnoId, precsno), &
          'nf90_get_var failed for precsno (flx) in read_merra2')

        call interp_merra2_var(n,findex,precsno,  10, .true.,merraforc)
     endif

   ! Read in Forcing Height, if selected in the forcing input table:
     if (LDT_FORC_Forc_Hgt%selectOpt.eq.1) then
        call LDT_verify(nf90_inq_varid(ftn_flx,'HLML',hlmlId), &
             'nf90_inq_varid failed for hlml (flx) in read_merra2')

        call LDT_verify(nf90_get_var(ftn_flx,hlmlId, hlml), &
             'nf90_get_var failed for hlml (flx) in read_merra2')

        call interp_merra2_var(n,findex,hlml,  14, .false.,merraforc)
     endif

     call LDT_verify(nf90_close(ftn_flx), &
          'failed to close flxfile in read_merra2')
     
  else
     write(LDT_logunit,*) '[ERR] ',trim(flxname)//' does not exist'
     call LDT_endrun()
  endif
  
! Read in the radiation file fields (*rad):
  inquire(file=radname,exist=file_exists) 
  if(file_exists) then 
     write(LDT_logunit,*) '[INFO] Reading MERRA-2 file (bookend,',order,' ... ',trim(radname)
     call LDT_verify(nf90_open(path=trim(radname), mode=NF90_NOWRITE, &
          ncid=ftn_rad), 'nf90_open failed in read_merra2')
     call LDT_verify(nf90_inq_varid(ftn_rad,'SWGDN',swgdnId), &
          'nf90_inq_varid failed for swgdn in read_merra2')
     call LDT_verify(nf90_inq_varid(ftn_rad,'LWGAB',lwgabId), &
          'nf90_inq_varid failed for lwgab in read_merra2')
!     call LDT_verify(nf90_inq_varid(ftn_rad,'emis',emisId), &
!          'nf90_inq_varid failed for emis in read_merra2')
    
     call LDT_verify(nf90_get_var(ftn_rad,swgdnId,swgdn), &
          'nf90_get_var failed for swgdn in read_merra2')
     call LDT_verify(nf90_get_var(ftn_rad,lwgabId,lwgab), &
          'nf90_get_var failed for lwgab in read_merra2')
!     call LDT_verify(nf90_get_var(ftn_rad,emisId,emis), &
!          'nf90_get_var failed for emis in read_merra2')

     call interp_merra2_var(n,findex,swgdn, 3, .false.,merraforc)
     call interp_merra2_var(n,findex,lwgab, 4, .false.,merraforc)

     call LDT_verify(nf90_close(ftn_rad), &
          'failed to close lfofile in read_merra2')
  else
     write(LDT_logunit,*) '[ERR] ',trim(radname)//' does not exist'
     call LDT_endrun()
  endif

  ! Checks: For reading in the surface layer file fields (*lfo):
  read_lnd = .false. 
  if ( LDT_FORC_Pardr%selectOpt.eq.1.or.&
       LDT_FORC_Pardf%selectOpt.eq.1.or.&
       LDT_FORC_SWnet%selectOpt.eq.1 ) then
     read_lnd = .true.
  endif
  if( merra2_struc(n)%usecorr.eq.1 ) then
    if( LDT_FORC_CRainf%selectOpt.eq.1 .or. &
       LDT_FORC_Snowf%selectOpt.eq.1 ) then
      read_lnd = .true.
    endif
  endif

! Read in the surface layer fields (*lfo):
  if(read_lnd) then 
     inquire(file=lfoname,exist=file_exists) 
     if(file_exists) then 
        write(LDT_logunit,*) '[INFO] Reading MERRA-2 file (bookend,',order,&
        ' ... ',trim(lfoname)
        call LDT_verify(nf90_open(path=trim(lfoname), mode=NF90_NOWRITE, &
             ncid=ftn_lfo), 'nf90_open failed in read_merra2')

        ! Read *corrected* convective precipitation, if selected: 
        if( LDT_FORC_CRainf%selectOpt.eq.1 .and. &
          merra2_struc(n)%usecorr.eq.1 ) then
          call LDT_verify(nf90_inq_varid(ftn_lfo,'PRECCUCORR',precconId), &
              'nf90_inq_varid failed for preccucorr (lfo) in read_merra2')

          call LDT_verify(nf90_get_var(ftn_lfo,precconId, preccon), &
              'nf90_get_var failed for preccucorr (lfo) in read_merra2')

          call interp_merra2_var(n,findex,preccon, 9,.true.,merraforc)
        endif

        ! Read *corrected* snowfall field, if selected: 
        if( LDT_FORC_Snowf%selectOpt.eq.1 .and. &
          merra2_struc(n)%usecorr.eq.1 ) then
          call LDT_verify(nf90_inq_varid(ftn_flx,'PRECSNOCORR',precsnoId), &
              'nf90_inq_varid failed for precsnocorr (lfo) in read_merra2')

          call LDT_verify(nf90_get_var(ftn_flx,precsnoId, precsno), &
              'nf90_get_var failed for precsnocorr (lfo) in read_merra2')

          call interp_merra2_var(n,findex,precsno, 10,.true.,merraforc)
        endif

        if(LDT_FORC_SWnet%selectOpt.eq.1) then 
           call LDT_verify(nf90_inq_varid(ftn_lfo,'SWLAND',swlandId), &
                'nf90_inq_varid failed for swland in read_merra2')

           call LDT_verify(nf90_get_var(ftn_lfo,swlandId,swland), &
                'nf90_get_var failed for swland in read_merra2')
           
           call interp_merra2_var(n,findex,swland,11,.false.,merraforc)
        endif
        
        if(LDT_FORC_Pardr%selectOpt.eq.1) then 
           call LDT_verify(nf90_inq_varid(ftn_lfo,'PARDR',pardrId), &
                'nf90_inq_varid failed for pardr in read_merra2')

           call LDT_verify(nf90_get_var(ftn_lfo,pardrId,pardr), &
                'nf90_get_var failed for pardr in read_merra2')
           
           call interp_merra2_var(n,findex,pardr,12,.false.,merraforc)
        endif

        if(LDT_FORC_Pardf%selectOpt.eq.1) then 
           call LDT_verify(nf90_inq_varid(ftn_lfo,'PARDF',pardfId), &
                'nf90_inq_varid failed for pardf in read_merra2')

           call LDT_verify(nf90_get_var(ftn_lfo,pardfId,pardf), &
                'nf90_get_var failed for pardf in read_merra2')
           call interp_merra2_var(n,findex,pardf,13,.false.,merraforc)
           
        endif

        call LDT_verify(nf90_close(ftn_lfo),&
             'nf90_close failed for lfofile in read_merra2')
     else
        write(LDT_logunit,*) '[ERR] ',trim(lfoname)//' does not exist'
        call LDT_endrun()
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
subroutine interp_merra2_var(n,findex, input_var,  var_index, &
     pcp_flag, merraforc)

! !USES: 
  use LDT_coreMod
  use LDT_logMod
  use merra2_forcingMod, only : merra2_struc

  implicit none

! !ARGUMENTS:
  integer, intent(in)    :: n
  integer, intent(in)    :: findex
  real,    intent(in)    :: input_var(merra2_struc(n)%nc, &
                                      merra2_struc(n)%nr, 24)
  integer, intent(in)    :: var_index
  logical, intent(in)    :: pcp_flag
  real,    intent(inout) :: merraforc(LDT_rc%met_nf(findex), &
                                      24, LDT_rc%lnc(n)*LDT_rc%lnr(n))
!
! !DESCRIPTION: 
!  This subroutine spatially interpolates a MERRA2 field
!  to the LDT running domain
! 
!EOP

  integer   :: t,c,r,k,iret
  real      :: f (merra2_struc(n)%nc*merra2_struc(n)%nr)
  logical*1 :: lb(merra2_struc(n)%nc*merra2_struc(n)%nr)
  logical*1 :: lo(LDT_rc%lnc(n)*LDT_rc%lnr(n))
! _____________________________________________________________

  do t=1,24
     lb = .true.
     do r=1,merra2_struc(n)%nr
        do c=1,merra2_struc(n)%nc
           k= c+(r-1)*merra2_struc(n)%nc
           f(k) = input_var(c,r,t)
           if ( f(k) == 1.e+15 ) then 
              f(k)  = LDT_rc%udef
              lb(k) = .false. 
           endif
        enddo
     enddo
          
     if(pcp_flag.and.&
         trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then 

        call conserv_interp(LDT_rc%gridDesc(n,:),lb,f,lo,&
             merraforc(var_index,t,:), &
             merra2_struc(n)%mi,LDT_rc%lnc(n)*LDT_rc%lnr(n),& 
             LDT_domain(n)%lat, LDT_domain(n)%lon,&
             merra2_struc(n)%w112,merra2_struc(n)%w122,&
             merra2_struc(n)%w212,merra2_struc(n)%w222,&
             merra2_struc(n)%n112,merra2_struc(n)%n122,&
             merra2_struc(n)%n212,merra2_struc(n)%n222,&
             LDT_rc%udef, iret)

     elseif(trim(LDT_rc%met_gridtransform(findex)).eq."bilinear".or.&
             trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then 
        call bilinear_interp(LDT_rc%gridDesc(n,:),lb,f,lo,&
             merraforc(var_index,t,:), &
             merra2_struc(n)%mi,LDT_rc%lnc(n)*LDT_rc%lnr(n), & 
             LDT_domain(n)%lat, LDT_domain(n)%lon,&
             merra2_struc(n)%w111,merra2_struc(n)%w121,&
             merra2_struc(n)%w211,merra2_struc(n)%w221,&
             merra2_struc(n)%n111,merra2_struc(n)%n121,&
             merra2_struc(n)%n211,merra2_struc(n)%n221,&
             LDT_rc%udef, iret)

     elseif(trim(LDT_rc%met_gridtransform(findex)).eq."neighbor") then 
        call neighbor_interp(LDT_rc%gridDesc(n,:),lb,f,lo,&
             merraforc(var_index,t,:),merra2_struc(n)%mi,&
             LDT_rc%lnc(n)*LDT_rc%lnr(n),&
             LDT_domain(n)%lat, LDT_domain(n)%lon,&
             merra2_struc(n)%n113,LDT_rc%udef,iret)

     ! If want to match the incoming domain / grid:
     elseif( LDT_rc%met_gridtransform(findex) == "none" ) then
        merraforc(var_index,t,:) = f 

     else
        write(LDT_logunit,*) '[ERR] Spatial interpolation option '//&
                             trim(LDT_rc%met_gridtransform(findex))//&
                             ' not supported for MERRA2'
        call LDT_endrun()
     endif
  enddo
  
end subroutine interp_merra2_var
