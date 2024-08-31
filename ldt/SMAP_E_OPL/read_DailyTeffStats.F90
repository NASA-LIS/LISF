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
! !ROUTINE: read_DailyTeffStats
! \label{read_DailyTeffStats}
! 
! !REVISION HISTORY: 
!  12 JAN 2022: Yonghwan Kwon, Initial Specification
! 
! !INTERFACE:
subroutine read_DailyTeffStats(doy)
! !USES:
  use LDT_logMod
  use LDT_smap_e_oplMod

  implicit none
! !ARGUMENTS: 
  integer            :: doy
!EOP
  logical            :: file_exists_ref, file_exists_lis

  ! read reference Teff daily mean and std dev
  inquire(file=trim(SMAPeOPL%dailystats_ref),exist=file_exists_ref)
  if(file_exists_ref) then
     write(LDT_logunit,*) '[INFO] Reading reference Teff daily mean and std dev: ',&
                          trim(SMAPeOPL%dailystats_ref)
     call read_meanstddev(SMAPeOPL%dailystats_ref,doy,&
                          SMAPeOPL%ngrid,&
                          SMAPeOPL%mu_6am_ref,SMAPeOPL%mu_6pm_ref,&
                          SMAPeOPL%sigma_6am_ref,SMAPeOPL%sigma_6pm_ref,&
                          SMAPeOPL%grid_col,SMAPeOPL%grid_row)
     write(LDT_logunit,*) '[INFO] Finished reading ',trim(SMAPeOPL%dailystats_ref)
  endif
    
  ! read lis Teff daily mean and std dev
  inquire(file=trim(SMAPeOPL%dailystats_lis),exist=file_exists_lis)
  if(file_exists_lis) then
     write(LDT_logunit,*) '[INFO] Reading LIS Teff daily mean and std dev: ',&
                          trim(SMAPeOPL%dailystats_lis)
     call read_meanstddev(SMAPeOPL%dailystats_lis,doy,&
                          SMAPeOPL%ngrid,&
                          SMAPeOPL%mu_6am_lis,SMAPeOPL%mu_6pm_lis,&
                          SMAPeOPL%sigma_6am_lis,SMAPeOPL%sigma_6pm_lis,&
                          SMAPeOPL%grid_col,SMAPeOPL%grid_row)
     write(LDT_logunit,*) '[INFO] Finished reading ',trim(SMAPeOPL%dailystats_lis)
  endif

end subroutine read_DailyTeffStats

!BOP
! 
! !ROUTINE: getattributes
! \label{getattributes}
!
! !INTERFACE:
subroutine getattributes(fname,ntimes,ngrid)
! 
! !USES:
  use LDT_logMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
!
  character (len=*)  :: fname
  integer            :: ntimes,ngrid
!EOP
  integer            :: ios,nid,gid

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
  call LDT_verify(ios,'Error opening file '//trim(fname))

  !get ngrid
  ios = nf90_inq_dimid(nid, 'ngrid',gid)
  call LDT_verify(ios, 'Error nf90_inq_dimid: ngrid')

  ios = nf90_inquire_dimension(nid, gid, len=ngrid)
  call LDT_verify(ios, 'Error nf90_inquire_dimension: ngrid')

  !get ntimes
  ios = nf90_get_att(nid, NF90_GLOBAL, &
                 'temporal_resolution_CDF', &
                 ntimes)
  call LDT_verify(ios, 'Error nf90_get_att: temporal_resolution_CDF')

  ios = nf90_close(ncid=nid)
  call LDT_verify(ios,'Error closing file '//trim(fname))
#endif

end subroutine getattributes


!BOP
! 
! !ROUTINE: read_meanstddev
! \label{read_meanstddev}
!
! !INTERFACE:
subroutine read_meanstddev(fname,doy,&
                           ngrid,&
                           mu_6am,mu_6pm,&
                           sigma_6am,sigma_6pm,&
                           grid_col,grid_row)
! 
! !USES:
  use LDT_logMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
!
  character (len=*)  :: fname
  integer            :: doy,ngrid
  real               :: mu_6am(ngrid,1,1), mu_6pm(ngrid,1,1)
  real               :: sigma_6am(ngrid,1,1), sigma_6pm(ngrid,1,1)
  integer            :: grid_col(ngrid), grid_row(ngrid)
!EOP
  integer            :: ios, nid
  integer            :: mu6amid, mu6pmid
  integer            :: sigma6amid, sigma6pmid
  integer            :: gridcolid, gridrowid

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
  call LDT_verify(ios,'Error opening file '//trim(fname))

  !get mean and std dev
  ios = nf90_inq_varid(nid,'SoilTeff_mu_6am',mu6amid)
  call LDT_verify(ios, 'Error nf90_inq_varid: SoilTeff_mu_6am')

  ios = nf90_inq_varid(nid,'SoilTeff_mu_6pm',mu6pmid)
  call LDT_verify(ios, 'Error nf90_inq_varid: SoilTeff_mu_6pm')

  ios = nf90_inq_varid(nid,'SoilTeff_sigma_6am',sigma6amid)
  call LDT_verify(ios, 'Error nf90_inq_varid: SoilTeff_sigma_6am')

  ios = nf90_inq_varid(nid,'SoilTeff_sigma_6pm',sigma6pmid)
  call LDT_verify(ios, 'Error nf90_inq_varid: SoilTeff_sigma_6pm')

  ios = nf90_inq_varid(nid,'grid_col',gridcolid)
  call LDT_verify(ios, 'Error nf90_inq_varid: grid_col')

  ios = nf90_inq_varid(nid,'grid_row',gridrowid)
  call LDT_verify(ios, 'Error nf90_inq_varid: grid_row')

  ios = nf90_get_var(nid, mu6amid, mu_6am, &
        start=(/1,doy,1/), &
        count=(/ngrid,1,1/))
  call LDT_verify(ios, 'Error nf90_get_var: SoilTeff_mu_6am')

  ios = nf90_get_var(nid, mu6pmid, mu_6pm, &
        start=(/1,doy,1/), &
        count=(/ngrid,1,1/))
  call LDT_verify(ios, 'Error nf90_get_var: SoilTeff_mu_6pm')

  ios = nf90_get_var(nid, sigma6amid, sigma_6am, &
        start=(/1,doy,1/), &
        count=(/ngrid,1,1/))
  call LDT_verify(ios, 'Error nf90_get_var: SoilTeff_sigma_6am')

  ios = nf90_get_var(nid, sigma6pmid, sigma_6pm, &
        start=(/1,doy,1/), &
        count=(/ngrid,1,1/))
  call LDT_verify(ios, 'Error nf90_get_var: SoilTeff_sigma_6pm')

  ios = nf90_get_var(nid, gridcolid, grid_col, &
        start=(/1/), &
        count=(/ngrid/))
  call LDT_verify(ios, 'Error nf90_get_var: grid_col')

  ios = nf90_get_var(nid, gridrowid, grid_row, &
        start=(/1/), &
        count=(/ngrid/))
  call LDT_verify(ios, 'Error nf90_get_var: grid_row')

  ios = nf90_close(ncid=nid)
  call LDT_verify(ios,'Error closing file '//trim(fname))
#endif

end subroutine read_meanstddev

!BOP
! 
! !ROUTINE: get_doy
! \label{get_doy}
!
! !INTERFACE:
subroutine get_doy(mo,da,doy)
! 
! !USES:

  implicit none

! !ARGUMENTS:
  integer    :: mo, da, doy
  integer    :: imo

!EOP

doy = 0
do imo = 1,mo
   if(imo.lt.mo) then
      if(imo.eq.1.or.&
       imo.eq.3.or.&
       imo.eq.5.or.&
       imo.eq.7.or.&
       imo.eq.8.or.&
       imo.eq.10.or.&
       imo.eq.12) then
         doy = doy + 31
      elseif(imo.eq.2) then
         doy = doy + 28
      else
         doy = doy + 30
      endif
   elseif(imo.eq.mo) then
      if(imo.eq.2) then
         if(da.eq.29) then
            da = 28
         endif
      endif
      doy = doy + da
   endif
enddo

end subroutine get_doy

!BOP
! 
! !ROUTINE: get_UTC
! \label{get_UTC}
!
! !INTERFACE:
subroutine get_UTC(n,TIMEsec,UTChr)
! 
! !USES:
  use LDT_coreMod
  use LDT_logMod, only: LDT_logunit

  implicit none

! !ARGUMENTS:
  integer, intent(in) :: n
  real*8              :: TIMEsec(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real                :: UTChr(LDT_rc%lnc(n),LDT_rc%lnr(n))

!EOP
  integer             :: ilat, ilon, imo, ida
  real*8              :: TIMEday
  real                :: TIMEhr
  real                :: UTCyr, UTCmo, UTCda
  integer             :: count_yr, dayremove

  do ilat=1,LDT_rc%lnr(n)
     do ilon=1,LDT_rc%lnc(n)

        if (TIMEsec(ilon,ilat).gt.0) then
           !write(LDT_logunit,*) 'EMK: ilon,ilat, TIMEsec= ', &
           !     ilon, ilat, TIMEsec(ilon,ilat)

           TIMEday = TIMEsec(ilon,ilat)/DBLE(60*60*24)
           !write(*,*) 'TIMEday= ', TIMEday

           count_yr = 0
           do while (TIMEday.ge.366)
              if(mod(count_yr,4).eq.0) then
                 TIMEday = TIMEday - 366
              else
                 TIMEday = TIMEday - 365
              endif
              count_yr = count_yr + 1
           enddo
           UTCyr = 2000 + count_yr

           !write(*,*) 'UTCyr= ', UTCyr
           !write(*,*) 'TIMEday= ', TIMEday          
 
           imo = 1
           do while (TIMEday.gt.0)
              if(imo.eq.1.or.&
               imo.eq.3.or.&
               imo.eq.5.or.&
               imo.eq.7.or.&
               imo.eq.8.or.&
               imo.eq.10.or.&
               imo.eq.12) then
                 TIMEday = TIMEday - 31
                 dayremove = 31
              elseif(imo.eq.2) then
                 if(mod(UTCyr,4.).eq.0) then
                    TIMEday = TIMEday - 29
                    dayremove = 29
                 else
                    TIMEday = TIMEday - 28
                    dayremove = 28
                 endif
              else
                 TIMEday = TIMEday - 30
                 dayremove = 30
              endif

              imo = imo + 1
           enddo  
          
           if(TIMEday.lt.0) then
              TIMEday = TIMEday + dayremove
              imo = imo - 1
           endif
           UTCmo = imo
           UTCda = ceiling(TIMEday)
           TIMEhr = (TIMEday - UTCda + 1)*24
           UTChr(ilon,ilat) = TIMEhr + 12  !reference time: Jan.1, 2000 at 12pm.

           if(UTChr(ilon,ilat).gt.24) then
              UTCda = UTCda + 1
              UTChr(ilon,ilat) = UTChr(ilon,ilat) - 24
           endif

           !write(LDT_logunit,*) 'EMK: ilon,ilat,TIMEsec,UTChr= ', &
           !     ilon, ilat, TIMEsec(ilon,ilat), UTChr(ilon,ilat)

        else
           UTChr(ilon,ilat) = LDT_rc%udef 
        endif

     enddo
  enddo

end subroutine get_UTC
