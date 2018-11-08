!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: read_merraland
! \label{read_merraland}
! 
! !REVISION HISTORY:
! 12 Oct 2009: Eric Kemp, Initial code
! 22 Jul 2010: David Mocko, changed to hourly forcing
!  5 Apr 2013: Sujay Kumar, updated for the direct use of files from GES-DISC
!
! !INTERFACE:      
subroutine read_merraland(n,order,findex,&
     slvname, flxname, radname, mldname, lndname, &
     merraforc,ferror)
! !USES:
  use LIS_coreMod, only : LIS_rc,LIS_domain, LIS_masterproc
  use LIS_logMod
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod, only: LIS_forc
  use merraland_forcingMod, only : merraland_struc
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
  character(len=*), intent(in) :: radname
  character(len=*), intent(in) :: mldname
  character(len=*), intent(in) :: lndname
  real                         :: merraforc(&
       merraland_struc(n)%nvars,24,&
       LIS_rc%lnc(n)*LIS_rc%lnr(n))
  integer, intent(out)         :: ferror          

!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  MERRA-Land data, transforms into 9 LIS forcing 
!  parameters and interpolates to the LIS domain. \newline
!
! MERRALAND FORCING VARIABLES (unless noted, fields are 1-hr upstream averaged): \newline
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
!    name of the 1 hour MERRA-Land analysis file
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
  
  integer   :: ftn_slv, ftn_flx, ftn_rad, ftn_mld, ftn_lnd
  integer   :: tmpId, qId, uId, vId, psId
  integer   :: prectotId, precconId, swgdnId, lwgabId, emisId
  integer   :: prectot2Id, preccon2Id, precsnoId, hlmlID
  integer   :: swlandId, pardrId, pardfId
  integer   :: nr_index, nc_index
  logical   :: file_exists, file_exists1
  integer   :: c,r,t,k,iret
  integer   :: mo
  logical   :: read_lnd

  real      :: tair(merraland_struc(n)%ncold, merraland_struc(n)%nrold,24)
  real      :: qair(merraland_struc(n)%ncold, merraland_struc(n)%nrold,24)
  real      :: uwind(merraland_struc(n)%ncold, merraland_struc(n)%nrold,24)
  real      :: vwind(merraland_struc(n)%ncold, merraland_struc(n)%nrold,24)
  real      :: ps(merraland_struc(n)%ncold, merraland_struc(n)%nrold,24)
  real      :: prectot(merraland_struc(n)%ncold, merraland_struc(n)%nrold,24)
  real      :: precsno(merraland_struc(n)%ncold, merraland_struc(n)%nrold,24)
  real      :: preccon(merraland_struc(n)%ncold, merraland_struc(n)%nrold,24)
  real      :: prectot_flx(merraland_struc(n)%ncold, merraland_struc(n)%nrold,24)
  real      :: preccon_flx(merraland_struc(n)%ncold, merraland_struc(n)%nrold,24)
  real      :: swgdn(merraland_struc(n)%ncold, merraland_struc(n)%nrold,24)
  real      :: lwgab(merraland_struc(n)%ncold, merraland_struc(n)%nrold,24)
  real      :: swland(merraland_struc(n)%ncold, merraland_struc(n)%nrold,24)
  real      :: pardr(merraland_struc(n)%ncold, merraland_struc(n)%nrold,24)
  real      :: pardf(merraland_struc(n)%ncold, merraland_struc(n)%nrold,24)
  real      :: hlml(merraland_struc(n)%ncold, merraland_struc(n)%nrold,24)
!  real      :: emis(merraland_struc(n)%ncold, merraland_struc(n)%nrold,24)

#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  ferror = 0 
  nr_index = merraland_struc(n)%nrold
  nc_index = merraland_struc(n)%ncold
  mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)

  inquire(file=slvname,exist=file_exists) 
  if(file_exists) then 
     write(LIS_logunit,*) 'Reading .. ', trim(slvname)
     call LIS_verify(nf90_open(path=trim(slvname), mode=NF90_NOWRITE, &
          ncid=ftn_slv), 'nf90_open failed for slvfile in read_merraland')
     call LIS_verify(nf90_inq_varid(ftn_slv,'ps',psId), &
          'nf90_inq_varid failed for ps in read_merraland')

     call LIS_verify(nf90_get_var(ftn_slv,psId, ps), &
          'nf90_get_var failed for ps in read_merraland') 
     
     if(merraland_struc(n)%uselml.eq.1) then 
        inquire(file=flxname,exist=file_exists1)
        if(.not.file_exists1) then 
           write(LIS_logunit,*) trim(flxname)//' does not exist'
           call LIS_endrun()
        endif
        write(LIS_logunit,*) 'Reading .. ', trim(flxname)
        call LIS_verify(nf90_open(path=trim(flxname), mode=NF90_NOWRITE, &
             ncid=ftn_flx), 'nf90_open failed for flxfile in read_merraland')

        call LIS_verify(nf90_inq_varid(ftn_flx,'tlml',tmpId), &
             'nf90_inq_varid failed for tlml in read_merraland')
        call LIS_verify(nf90_inq_varid(ftn_flx,'qlml',qId), &
             'nf90_inq_varid failed for qlml in read_merraland')
        call LIS_verify(nf90_inq_varid(ftn_flx,'ulml',uId), &
             'nf90_inq_varid failed for ulml in read_merraland')
        call LIS_verify(nf90_inq_varid(ftn_flx,'vlml',vId), &
             'nf90_inq_varid failed for vlml in read_merraland')

        call LIS_verify(nf90_get_var(ftn_flx,tmpId, tair), &
             'nf90_get_var failed for tlml in read_merraland')
        call LIS_verify(nf90_get_var(ftn_flx,qId, qair), &
             'nf90_get_var failed for qlml in read_merraland')
        call LIS_verify(nf90_get_var(ftn_flx,uId, uwind), &
             'nf90_get_var failed for ulml in read_merraland')
        call LIS_verify(nf90_get_var(ftn_flx,vId, vwind), &
             'nf90_get_var failed for vlml in read_merraland')

     else
        call LIS_verify(nf90_inq_varid(ftn_slv,'t2m',tmpId), &
             'nf90_inq_varid failed for t2m in read_merraland')
        call LIS_verify(nf90_inq_varid(ftn_slv,'qv2m',qId), &
             'nf90_inq_varid failed for qv2m in read_merraland')
        call LIS_verify(nf90_inq_varid(ftn_slv,'u10m',uId), &
             'nf90_inq_varid failed for u10m in read_merraland')
        call LIS_verify(nf90_inq_varid(ftn_slv,'v10m',vId), &
             'nf90_inq_varid failed for v10m in read_merraland')

        call LIS_verify(nf90_get_var(ftn_slv,tmpId, tair), &
             'nf90_get_var failed for t2ml in read_merraland')
        call LIS_verify(nf90_get_var(ftn_slv,qId, qair), &
             'nf90_get_var failed for qv2m in read_merraland')
        call LIS_verify(nf90_get_var(ftn_slv,uId, uwind), &
             'nf90_get_var failed for u10m in read_merraland')
        call LIS_verify(nf90_get_var(ftn_slv,vId, vwind), &
             'nf90_get_var failed for v10m in read_merraland')
     endif

     call LIS_verify(nf90_close(ftn_slv), &
          'failed to close slvfile in read_merraland')

     call interp_merraland_var(n,findex,tair,  1, .false., merraforc)
     call interp_merraland_var(n,findex,qair,  2, .false., merraforc)
     call interp_merraland_var(n,findex,uwind, 5, .false., merraforc)
     call interp_merraland_var(n,findex,vwind, 6, .false., merraforc)
     call interp_merraland_var(n,findex,ps,    7, .false., merraforc)
  else
     write(LIS_logunit,*) trim(slvname)//' does not exist'
     call LIS_endrun()
  endif

  inquire(file=mldname,exist=file_exists) 
  if(file_exists) then 

     write(LIS_logunit,*) 'Reading .. ', trim(mldname)
     call LIS_verify(nf90_open(path=trim(mldname), mode=NF90_NOWRITE, &
          ncid=ftn_mld),'nf90_open failed for mldfile in read_merraland')

     call LIS_verify(nf90_inq_varid(ftn_mld,'prectot',prectotId), &
          'nf90_inq_varid failed for prectot (mld) in read_merraland')

     call LIS_verify(nf90_get_var(ftn_mld,prectotId, prectot), &
          'nf90_get_var failed for prectot (mld) in read_merraland')

     if (LIS_FORC_Forc_Hgt%selectOpt.eq.1) then
        call LIS_verify(nf90_inq_varid(ftn_flx,'hlml',hlmlId), &
             'nf90_inq_varid failed for hlml (flx) in read_merraland')

        call LIS_verify(nf90_get_var(ftn_flx,hlmlId, hlml), &
             'nf90_get_var failed for hlml (flx) in read_merraland')
     endif
     
     if (LIS_FORC_CRainf%selectOpt.eq.1) then
        call LIS_verify(nf90_inq_varid(ftn_flx,'prectot',prectot2Id), &
             'nf90_inq_varid failed for prectot (flx) in read_merraland')

        call LIS_verify(nf90_inq_varid(ftn_flx,'preccon',preccon2Id), &
             'nf90_inq_varid failed for preccon (flx) in read_merraland')

        call LIS_verify(nf90_get_var(ftn_flx,prectot2Id, prectot_flx), &
             'nf90_get_var failed for prectot (flx) in read_merraland')

        call LIS_verify(nf90_get_var(ftn_flx,preccon2Id, preccon_flx), &
             'nf90_get_var failed for preccon (flx) in read_merraland')
     endif

     if (LIS_FORC_Snowf%selectOpt.eq.1) then
        call LIS_verify(nf90_inq_varid(ftn_mld,'precsno',precsnoId), &
             'nf90_inq_varid failed for precsno (mld) in read_merraland')  

        call LIS_verify(nf90_get_var(ftn_mld,precsnoId, precsno), &
          'nf90_get_var failed for precsno (mld) in read_merraland')
     endif

     if (LIS_FORC_CRainf%selectOpt.eq.1) then
        do t=1,24
           do r=1,nr_index
              do c=1,nc_index
                 if ( prectot_flx(c,r,t) /= 1.e+15 .and. &
                      preccon_flx(c,r,t) /= 1.e+15 .and. &
                      prectot(c,r,t)     /= 1.e+15 .and. &
                      prectot_flx(c,r,t) /= 1.e+15 ) then
                    if(prectot_flx(c,r,t).gt.0.00000001) then 
                       preccon(c,r,t) = preccon_flx(c,r,t)* &
                            (prectot(c,r,t)/prectot_flx(c,r,t))
                       preccon(c,r,t) = min(preccon(c,r,t),prectot(c,r,t))
                       preccon(c,r,t) = max(preccon(c,r,t),0.0)
                    else
                       preccon(c,r,t) = 0.0
                    endif
                 else
                    preccon(c,r,t) = 0.0
                 endif
              enddo
           enddo
        enddo
     endif

     if (LIS_FORC_Snowf%selectOpt.eq.1) then
        do t=1,24
           do r=1,nr_index
              do c=1,nc_index
                 if ( prectot(c,r,t) /= 1.e+15 .and. &
                      precsno(c,r,t) /= 1.e+15 ) then
                    prectot(c,r,t) = prectot(c,r,t) - precsno(c,r,t)
                    if(prectot(c,r,t).lt.0) then 
                       !prectot(c,r,t) = precsno(c,r,t)
                       prectot(c,r,t) = 0.0
                    endif
                 endif
              enddo
           enddo
        enddo

        call interp_merraland_var(n,findex,precsno,  10, .true.,merraforc)
     endif

     call interp_merraland_var(n,findex,prectot,  8, .true.,merraforc)
     
     if (LIS_FORC_CRainf%selectOpt.eq.1) then
        call interp_merraland_var(n,findex,preccon,  9, .true.,merraforc)
     endif

     if (LIS_FORC_Forc_Hgt%selectOpt.eq.1) then
        call interp_merraland_var(n,findex,hlml,  14, .false.,merraforc)
     endif
  else
     write(LIS_logunit,*) trim(mldname)//' does not exist'
     call LIS_endrun()
  endif
  
  inquire(file=radname,exist=file_exists) 
  if(file_exists) then 
     write(LIS_logunit,*) 'Reading .. ', trim(radname)
     call LIS_verify(nf90_open(path=trim(radname), mode=NF90_NOWRITE, &
          ncid=ftn_rad), 'nf90_open failed in read_merraland')
     call LIS_verify(nf90_inq_varid(ftn_rad,'swgdn',swgdnId), &
          'nf90_inq_varid failed for swgdn in read_merraland')
     call LIS_verify(nf90_inq_varid(ftn_rad,'lwgab',lwgabId), &
          'nf90_inq_varid failed for lwgab in read_merraland')
!     call LIS_verify(nf90_inq_varid(ftn_rad,'emis',emisId), &
!          'nf90_inq_varid failed for emis in read_merraland')
    
     call LIS_verify(nf90_get_var(ftn_rad,swgdnId,swgdn), &
          'nf90_get_var failed for swgdn in read_merraland')
     call LIS_verify(nf90_get_var(ftn_rad,lwgabId,lwgab), &
          'nf90_get_var failed for lwgab in read_merraland')
!     call LIS_verify(nf90_get_var(ftn_rad,emisId,emis), &
!          'nf90_get_var failed for emis in read_merraland')

     call interp_merraland_var(n,findex,swgdn, 3, .false.,merraforc)
     call interp_merraland_var(n,findex,lwgab, 4, .false.,merraforc)
  else
     write(LIS_logunit,*) trim(radname)//' does not exist'
     call LIS_endrun()
  endif

  read_lnd = .false. 
  if (LIS_FORC_Pardr%selectOpt.eq.1.or.&
       LIS_FORC_Pardf%selectOpt.eq.1.or.&
       LIS_FORC_SWnet%selectOpt.eq.1) then
     read_lnd = .true.
  endif

  if(read_lnd) then 
     inquire(file=lndname,exist=file_exists) 
     if(file_exists) then 
        write(LIS_logunit,*) 'Reading .. ', trim(lndname)
        call LIS_verify(nf90_open(path=trim(lndname), mode=NF90_NOWRITE, &
             ncid=ftn_lnd), 'nf90_open failed in read_merraland')
        if(LIS_FORC_SWnet%selectOpt.eq.1) then 
           call LIS_verify(nf90_inq_varid(ftn_lnd,'swland',swlandId), &
                'nf90_inq_varid failed for swland in read_merraland')

           call LIS_verify(nf90_get_var(ftn_lnd,swlandId,swland), &
                'nf90_get_var failed for swland in read_merraland')
           
           call interp_merraland_var(n,findex,swland,11,.false.,merraforc)
        endif
        
        if(LIS_FORC_Pardr%selectOpt.eq.1) then 
           call LIS_verify(nf90_inq_varid(ftn_lnd,'pardr',pardrId), &
                'nf90_inq_varid failed for pardr in read_merraland')

           call LIS_verify(nf90_get_var(ftn_lnd,pardrId,pardr), &
                'nf90_get_var failed for pardr in read_merraland')
           
           call interp_merraland_var(n,findex,pardr,12,.false.,merraforc)
        endif

        if(LIS_FORC_Pardf%selectOpt.eq.1) then 
           call LIS_verify(nf90_inq_varid(ftn_lnd,'pardf',pardfId), &
                'nf90_inq_varid failed for pardf in read_merraland')

           call LIS_verify(nf90_get_var(ftn_lnd,pardfId,pardf), &
                'nf90_get_var failed for pardf in read_merraland')
           call interp_merraland_var(n,findex,pardf,13,.false.,merraforc)
           
        endif

        call LIS_verify(nf90_close(ftn_lnd),&
             'nf90_close failed for lndfile in read_merraland')
     else
        write(LIS_logunit,*) trim(lndname)//' does not exist'
        call LIS_endrun()
     endif

  endif
  
  call LIS_verify(nf90_close(ftn_rad),&
       'nf90_close failed for radfile in read_merraland')
  call LIS_verify(nf90_close(ftn_mld),&
       'nf90_close failed for mldfile in read_merraland')
  call LIS_verify(nf90_close(ftn_flx),&
       'nf90_close failed for flxfile in read_merraland')
#endif
end subroutine read_merraland

!BOP
! 
! !ROUTINE: interp_merraland_var
! \label{interp_merraland_var}
! 
! !INTERFACE: 
subroutine interp_merraland_var(n,findex, input_var,  var_index, &
     pcp_flag, merraforc)
! !USES: 
  use LIS_coreMod
  use LIS_logMod
  use merraland_forcingMod, only : merraland_struc

  implicit none

! !ARGUMENTS:
  integer, intent(in)    :: n
  integer, intent(in)    :: findex
  real,    intent(in)    :: input_var(merraland_struc(n)%ncold, &
       merraland_struc(n)%nrold, 24)
  integer, intent(in)    :: var_index
  logical, intent(in)    :: pcp_flag
  real,    intent(inout) :: merraforc(merraland_struc(n)%nvars, &
       24, LIS_rc%lnc(n)*LIS_rc%lnr(n))
!
! !DESCRIPTION: 
!  This subroutine spatially interpolates a MERRA-Land field
!  to the LIS running domain
! 
!EOP

  integer   :: t,c,r,k,iret
  real      :: f (merraland_struc(n)%ncold*merraland_struc(n)%nrold)
  logical*1 :: lb(merraland_struc(n)%ncold*merraland_struc(n)%nrold)
  logical*1 :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))

  do t=1,24
     lb = .true.
     do r=1,merraland_struc(n)%nrold
        do c=1,merraland_struc(n)%ncold
           k= c+(r-1)*merraland_struc(n)%ncold
           f(k) = input_var(c,r,t)
           if ( f(k) == 1.e+15 ) then 
              f(k)  = LIS_rc%udef
              lb(k) = .false. 
           endif
        enddo
     enddo
          
     if(pcp_flag.and.&
          trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 
        call conserv_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
             merraforc(var_index,t,:), &
             merraland_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),& 
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             merraland_struc(n)%w112,merraland_struc(n)%w122,&
             merraland_struc(n)%w212,merraland_struc(n)%w222,&
             merraland_struc(n)%n112,merraland_struc(n)%n122,&
             merraland_struc(n)%n212,merraland_struc(n)%n222,&
             LIS_rc%udef, iret)
     elseif(trim(LIS_rc%met_interp(findex)).eq."bilinear".or.&
          trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 
        call bilinear_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
             merraforc(var_index,t,:), &
             merraland_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n), & 
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             merraland_struc(n)%w111,merraland_struc(n)%w121,&
             merraland_struc(n)%w211,merraland_struc(n)%w221,&
             merraland_struc(n)%n111,merraland_struc(n)%n121,&
             merraland_struc(n)%n211,merraland_struc(n)%n221,&
             LIS_rc%udef, iret)
     elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then 
        call neighbor_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
             merraforc(var_index,t,:),merraland_struc(n)%mi,&
             LIS_rc%lnc(n)*LIS_rc%lnr(n),&
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             merraland_struc(n)%n113,LIS_rc%udef,iret)
     else
        write(LIS_logunit,*) 'Spatial interpolation option '//&
                             trim(LIS_rc%met_interp(findex))//&
                             ' not supported for MERRA-Land'
        call LIS_endrun()
     endif
  enddo
  
end subroutine interp_merraland_var
