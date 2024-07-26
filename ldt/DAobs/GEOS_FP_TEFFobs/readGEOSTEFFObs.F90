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
! !ROUTINE: readGEOSTEFFObs
! \label{readGEOSTEFFObs}
! 
! !REVISION HISTORY: 
!  29 NOV 2021: Yonghwan Kwon, Initial Specification
! 
! !INTERFACE: 
subroutine readGEOSTEFFObs(n)
! !USES:
  use ESMF
  use LDT_coreMod
  use LDT_logMod
  use LDT_DAobsDataMod
  use LDT_timeMgrMod
  use GEOSTEFF_obsMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for 
! GEOS FP soil temperature data
!
!EOP

  logical           :: file_exists
  integer           :: c,r
  character*100     :: fname
  real              :: tsoil01obs(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real              :: tsoil02obs(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real              :: kk, cc_6am, cc_6pm    !parameters for calculating effective soil temperature
                                                   !from soil layer temperature at 6 AM and 6 PM local time
  real              :: lon
  real              :: gmt
  real              :: lhour
  integer           :: zone

  GEOSTeffobs(n)%teffobs = LDT_rc%udef
  tsoil01obs = LDT_rc%udef
  tsoil02obs = LDT_rc%udef

  call create_GEOSsoilT_filename(GEOSTeffobs(n)%odir, &
       LDT_rc%yr, LDT_rc%mo, LDT_rc%da, LDT_rc%hr, fname)

  inquire(file=trim(fname),exist=file_exists)
  if(file_exists) then
     write(LDT_logunit,*) '[INFO] Reading ',trim(fname)
     call read_GEOSTeff_data(n, fname, tsoil01obs, tsoil02obs)
     write(LDT_logunit,*) '[INFO] Finished reading ',trim(fname)

     kk = 1.007
     cc_6am = 0.246;    !Descending
     cc_6pm = 1.000;    !Ascending

     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
           
           !calculate effective soil temperature
           lon = LDT_domain(n)%lon(c+(r-1)*LDT_rc%lnc(n))
           gmt = LDT_rc%hr
           call LDT_gmt2localtime(gmt, lon, lhour, zone)

           if(lhour.eq.6) then
              if(tsoil01obs(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0.and.&
                 tsoil02obs(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then
                 GEOSTeffobs(n)%teffobs(c,r) = &
                        kk * (cc_6am * tsoil01obs(c+(r-1)*LDT_rc%lnc(n)) + &
                        (1 - cc_6am) * tsoil02obs(c+(r-1)*LDT_rc%lnc(n)))
              endif
           elseif(lhour.eq.18) then
              if(tsoil01obs(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0.and.&
                 tsoil02obs(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then
                 GEOSTeffobs(n)%teffobs(c,r) = &
                        kk * (cc_6pm * tsoil01obs(c+(r-1)*LDT_rc%lnc(n)) + &
                        (1 - cc_6pm) * tsoil02obs(c+(r-1)*LDT_rc%lnc(n)))
              endif 
           endif

        enddo
     enddo
  endif

  call LDT_logSingleDAobs(n,LDT_DAobsData(n)%teff_obs,&
       GEOSTeffobs(n)%teffobs,vlevel=1)

end subroutine readGEOSTEFFObs

!BOP
! 
! !ROUTINE: read_GEOSTeff_data
! \label{read_GEOSTeff_data}
!
! !INTERFACE:
subroutine read_GEOSTeff_data(n, fname, tsoil01obs_ip, tsoil02obs_ip)
! 
! !USES:

  use LDT_coreMod
  use LDT_logMod
  use LDT_timeMgrMod
  use GEOSTEFF_obsMod

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n
  character (len=*)             :: fname
  real                          :: tsoil01obs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real                          :: tsoil02obs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
!
! !DESCRIPTION:
!  This subroutine reads the GEOS FP soil temperature file 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the GEOS FP file
!  \item[tsoil01obs_ip]   First layer soil temperature data processed to the LIS domain
!  \item[tsoil02obs_ip]   Second layer soil temperature data processed to the LIS domain
! \end{description}
!
!
!EOP

  integer,  parameter     :: nc=1152, nr=721
  real                    :: tsoil01(GEOSTeffobs(n)%nc,GEOSTeffobs(n)%nr,1)
  real                    :: tsoil02(GEOSTeffobs(n)%nc,GEOSTeffobs(n)%nr,1)
  real                    :: tsoil01_in(GEOSTeffobs(n)%nc*GEOSTeffobs(n)%nr)
  real                    :: tsoil02_in(GEOSTeffobs(n)%nc*GEOSTeffobs(n)%nr)
  logical*1               :: tsoil01_data_b(GEOSTeffobs(n)%nc*GEOSTeffobs(n)%nr)
  logical*1               :: tsoil02_data_b(GEOSTeffobs(n)%nc*GEOSTeffobs(n)%nr)
  logical*1               :: tsoil01obs_b_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  logical*1               :: tsoil02obs_b_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  integer                 :: tsoil01id, tsoil02id
  integer                 :: ios, nid
  integer                 :: c,r

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
  call LDT_verify(ios,'Error opening file '//trim(fname))

  ios = nf90_inq_varid(nid,'TSOIL1',tsoil01id)
  call LDT_verify(ios, 'Error nf90_inq_varid: TSOIL1')

  ios = nf90_inq_varid(nid,'TSOIL2',tsoil02id)
  call LDT_verify(ios, 'Error nf90_inq_varid: TSOIL2')

  !values 

  tsoil01_data_b = .false.
  tsoil02_data_b = .false.

  ios = nf90_get_var(nid, tsoil01id, tsoil01, &
        start=(/1,1,1/), &
        count=(/GEOSTeffobs(n)%nc,GEOSTeffobs(n)%nr,1/))
  call LDT_verify(ios, 'Error nf90_get_var: tsoil01')

  ios = nf90_get_var(nid, tsoil02id, tsoil02, &
        start=(/1,1,1/), &
        count=(/GEOSTeffobs(n)%nc,GEOSTeffobs(n)%nr,1/))
  call LDT_verify(ios, 'Error nf90_get_var: tsoil02')

  ios = nf90_close(ncid=nid)
  call LDT_verify(ios,'Error closing file '//trim(fname))

#endif

  do r=1,GEOSTeffobs(n)%nr
     do c=1,GEOSTeffobs(n)%nc

        if (tsoil01(c,r,1).gt.273.15.and.tsoil01(c,r,1).lt.1000) then
           tsoil01_in(c+(r-1)*GEOSTeffobs(n)%nc) = tsoil01(c,r,1)
           tsoil01_data_b(c+(r-1)*GEOSTeffobs(n)%nc) = .true.
        else
           tsoil01_in(c+(r-1)*GEOSTeffobs(n)%nc) = LDT_rc%udef
           tsoil01_data_b(c+(r-1)*GEOSTeffobs(n)%nc) = .false.
        endif

        if (tsoil02(c,r,1).gt.273.15.and.tsoil02(c,r,1).lt.1000) then
           tsoil02_in(c+(r-1)*GEOSTeffobs(n)%nc) = tsoil02(c,r,1)
           tsoil02_data_b(c+(r-1)*GEOSTeffobs(n)%nc) = .true.
        else
           tsoil02_in(c+(r-1)*GEOSTeffobs(n)%nc) = LDT_rc%udef
           tsoil02_data_b(c+(r-1)*GEOSTeffobs(n)%nc) = .false.
        endif

     enddo
  enddo

!--------------------------------------------------------------------------
! Interpolate to the LDT running domain
!--------------------------------------------------------------------------
  if(LDT_isLDTatAfinerResolution(n,0.3125)) then
     call bilinear_interp(LDT_rc%gridDesc(n,:),&
          tsoil01_data_b, tsoil01_in, tsoil01obs_b_ip, tsoil01obs_ip, &
          GEOSTeffobs(n)%nc*GEOSTeffobs(n)%nr, &
          LDT_rc%lnc(n)*LDT_rc%lnr(n), &
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          GEOSTeffobs(n)%w11,GEOSTeffobs(n)%w12,&
          GEOSTeffobs(n)%w21,GEOSTeffobs(n)%w22,&
          GEOSTeffobs(n)%n11,GEOSTeffobs(n)%n12,&
          GEOSTeffobs(n)%n21,GEOSTeffobs(n)%n22,LDT_rc%udef,ios)

     call bilinear_interp(LDT_rc%gridDesc(n,:),&
          tsoil02_data_b, tsoil02_in, tsoil02obs_b_ip, tsoil02obs_ip, &
          GEOSTeffobs(n)%nc*GEOSTeffobs(n)%nr, &
          LDT_rc%lnc(n)*LDT_rc%lnr(n), &
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          GEOSTeffobs(n)%w11,GEOSTeffobs(n)%w12,&
          GEOSTeffobs(n)%w21,GEOSTeffobs(n)%w22,&
          GEOSTeffobs(n)%n11,GEOSTeffobs(n)%n12,&
          GEOSTeffobs(n)%n21,GEOSTeffobs(n)%n22,LDT_rc%udef,ios)
  else
     call upscaleByAveraging(GEOSTeffobs(n)%nc*GEOSTeffobs(n)%nr,&
          LDT_rc%lnc(n)*LDT_rc%lnr(n), &
          LDT_rc%udef, GEOSTeffobs(n)%n11,&
          tsoil01_data_b,tsoil01_in, tsoil01obs_b_ip, tsoil01obs_ip)

     call upscaleByAveraging(GEOSTeffobs(n)%nc*GEOSTeffobs(n)%nr,&
          LDT_rc%lnc(n)*LDT_rc%lnr(n), &
          LDT_rc%udef, GEOSTeffobs(n)%n11,&
          tsoil02_data_b,tsoil02_in, tsoil02obs_b_ip, tsoil02obs_ip)
  endif

end subroutine read_GEOSTeff_data


!BOP
! !ROUTINE: create_GEOSsoilT_filename
! \label{create_GEOSsoilT_filename}
! 
! !INTERFACE:
subroutine create_GEOSsoilT_filename(ndir,yr,mo,da,hr,filename)
! !USES:

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr, mo, da, hr
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the GEOS soilT filename
!  based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the GEOS soilT data directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[hr]  current hour
!  \item[filename] Generated GEOS soilT filename
! \end{description}
!EOP

  character*4             :: fyr
  character*2             :: fmo
  character*2             :: fda
  character*2             :: fhr

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr

  filename = trim(ndir)//'/Y'//trim(fyr)//'/M'//trim(fmo)//'/D'//trim(fda)//&
             '/GEOS.fp.asm.inst1_2d_smp_Nx.'//&
             trim(fyr)//trim(fmo)//trim(fda)//'_'//trim(fhr)//'00.V01.nc4'

end subroutine create_GEOSsoilT_filename
