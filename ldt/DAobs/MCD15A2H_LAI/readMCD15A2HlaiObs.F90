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
! !ROUTINE: readMCD15A2HlaiObs
! \label{readMCD15A2HlaiObs}
! 
! !REVISION HISTORY: 
!  12 Nov 2020: Wanshu Nie, Initial Specification
! 
! !INTERFACE: 
subroutine readMCD15A2HlaiObs(n)
! !USES:   
  use ESMF
  use LDT_coreMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_DAobsDataMod
  use MCD15A2Hlai_obsMod
  use map_utils

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for the MODIS MCD15A2H
! LAI data product.
!
!EOP

  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r,i,j
  character(len=LDT_CONST_PATH_LEN)     :: fname
  character(len=LDT_CONST_PATH_LEN)     :: climofile
  real              :: laiobs(LDT_rc%lnc(n)*LDT_rc%lnr(n))

!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations. 
!-----------------------------------------------------------------------
  MCD15A2Hlaiobs(n)%laiobs = LDT_rc%udef
  laiobs= LDT_rc%udef

  call create_MCD15A2Hlai_filename(MCD15A2Hlaiobs(n)%odir, &
       MCD15A2Hlaiobs(n)%version,&
       LDT_rc%yr, LDT_rc%doy, fname, climofile)

  inquire(file=trim(fname),exist=file_exists)
  if(file_exists) then

     write(LDT_logunit,*) '[INFO] Reading ',trim(fname)
     call read_MCD15A2Hlai_data(n, fname, climofile, laiobs)
     write(LDT_logunit,*) '[INFO] Finished reading ',trim(fname)

     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
           if(laiobs(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then
              MCD15A2Hlaiobs(n)%laiobs(c,r) = laiobs(c+(r-1)*LDT_rc%lnc(n))
           endif
        enddo
     enddo
  endif

  call LDT_logSingleDAobs(n,LDT_DAobsData(n)%lai_obs,&
       MCD15A2Hlaiobs(n)%laiobs,vlevel=1)

end subroutine readMCD15A2HlaiObs

!BOP
! 
! !ROUTINE: read_MCD15A2Hlai_data
! \label{read_MCD15A2Hlai_data}
!
! !INTERFACE:
subroutine read_MCD15A2Hlai_data(n, fname, climofile,laiobs_ip)
! 
! !USES:   

  use LDT_coreMod
  use LDT_logMod
  use LDT_timeMgrMod
  use MCD15A2Hlai_obsMod

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n
  character (len=*)             :: fname
  character (len=*)             :: climofile
  real                          :: laiobs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real*8                        :: cornerlat(2), cornerlon(2)
!
! !DESCRIPTION: 
!  This subroutine reads the MODIS MCD15A2H LAI file and applies the data
!  quality flags to filter the data. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the MODIS MCD15A2H LAI file
!  \item[laiobs\_ip]   LAI data processed to the LDT domain
! \end{description}
!
!
!EOP

! !USES:   
  integer,  parameter     :: nc=86400, nr=43200
  integer                 :: lat_off, lon_off
  integer    :: lai(MCD15A2Hlaiobs(n)%nc,MCD15A2Hlaiobs(n)%nr)
  integer    :: flag(MCD15A2Hlaiobs(n)%nc,MCD15A2Hlaiobs(n)%nr)
  real       :: lai_flagged(MCD15A2Hlaiobs(n)%nc,MCD15A2Hlaiobs(n)%nr)
  real                    :: lai_in(MCD15A2Hlaiobs(n)%nc*MCD15A2Hlaiobs(n)%nr)
  logical*1               :: lai_data_b(MCD15A2Hlaiobs(n)%nc*MCD15A2Hlaiobs(n)%nr)
  logical*1               :: laiobs_b_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real                    :: laiobs_climo_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  integer                 :: c,r,t
  integer                 :: nid
  integer                 :: laiid, flagid
  integer                 :: ios
  integer, dimension(nf90_max_var_dims) :: dimIDs
  integer                                :: numLons, numLats


  ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
  call LDT_verify(ios,'Error opening file '//trim(fname))

  ios = nf90_inq_varid(nid, 'Lai_500m',laiid)
  call LDT_verify(ios, 'Error nf90_inq_varid: Lai_500m')

  ios = nf90_inq_varid(nid, 'FparLai_QC',flagid)
  call LDT_verify(ios, 'Error nf90_inq_varid: flag')

  !values

  cornerlat(1)=MCD15A2Hlaiobs(n)%gridDesci(4)
  cornerlon(1)=MCD15A2Hlaiobs(n)%gridDesci(5)
  cornerlat(2)=MCD15A2Hlaiobs(n)%gridDesci(7)
  cornerlon(2)=MCD15A2Hlaiobs(n)%gridDesci(8)

  lai_data_b = .false.

  lat_off = nint((cornerlat(1)+89.9979167)/0.00416667)+1
  lon_off = nint((cornerlon(1)+179.9979167)/0.00416667)+1

  ios = nf90_get_var(nid, laiid, lai, &
       start=(/lon_off,lat_off/), &
       count=(/MCD15A2Hlaiobs(n)%nc,MCD15A2Hlaiobs(n)%nr/))

  call LDT_verify(ios, 'Error nf90_get_var: Lai_500m')

  ios = nf90_get_var(nid, flagid, flag, &
       start=(/lon_off,lat_off/), &
       count=(/MCD15A2Hlaiobs(n)%nc,MCD15A2Hlaiobs(n)%nr/))

  call LDT_verify(ios, 'Error nf90_get_var: flag')

  ios = nf90_close(ncid=nid)
  call LDT_verify(ios,'Error closing file '//trim(fname))

  do r=1, MCD15A2Hlaiobs(n)%nr
     do c=1, MCD15A2Hlaiobs(n)%nc

        if(MCD15A2Hlaiobs(n)%qcflag.eq.1) then !apply QC flag

          if(lai(c,r).gt.0.and.lai(c,r).le.100) then
             if (MOD(flag(c,r),2) ==0.and.flag(c,r).le.62) then
                lai_flagged(c,r) =&
                   lai(c,r)*0.1
             else
               lai_flagged(c,r) = LDT_rc%udef
             endif
          else
            lai_flagged(c,r) = LDT_rc%udef
          endif

        else  ! no QC flag applied                

           if(lai(c,r).gt.0.and.lai(c,r).le.100) then
              lai_flagged(c,r) =&
                   lai(c,r)*0.1
           else
              lai_flagged(c,r) = LDT_rc%udef
           endif
        endif
     end do
  end do


  do r=1, MCD15A2Hlaiobs(n)%nr
     do c=1, MCD15A2Hlaiobs(n)%nc
        lai_in(c+(r-1)*MCD15A2Hlaiobs(n)%nc) = lai_flagged(c,r)
        if(lai_flagged(c,r).ne.LDT_rc%udef) then
           lai_data_b(c+(r-1)*MCD15A2Hlaiobs(n)%nc) = .true.
        else
           lai_data_b(c+(r-1)*MCD15A2Hlaiobs(n)%nc) = .false.
        endif
     enddo
  enddo


  if(LDT_isLDTatAfinerResolution(n,0.00416667)) then

!--------------------------------------------------------------------------
! Interpolate to the LDT running domain
!-------------------------------------------------------------------------- 

     call bilinear_interp(LDT_rc%gridDesc(n,:),&
          lai_data_b, lai_in, laiobs_b_ip, laiobs_ip, &
          MCD15A2Hlaiobs(n)%nc*MCD15A2Hlaiobs(n)%nr, &
          LDT_rc%lnc(n)*LDT_rc%lnr(n), &
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          MCD15A2Hlaiobs(n)%w11,MCD15A2Hlaiobs(n)%w12,&
          MCD15A2Hlaiobs(n)%w21,MCD15A2Hlaiobs(n)%w22,&
          MCD15A2Hlaiobs(n)%n11,MCD15A2Hlaiobs(n)%n12,&
          MCD15A2Hlaiobs(n)%n21,MCD15A2Hlaiobs(n)%n22,LDT_rc%udef,ios)
  else
     call upscaleByAveraging(MCD15A2Hlaiobs(n)%nc*MCD15A2Hlaiobs(n)%nr,&
          LDT_rc%lnc(n)*LDT_rc%lnr(n), &
          LDT_rc%udef, MCD15A2Hlaiobs(n)%n11,&
          lai_data_b,lai_in, laiobs_b_ip, laiobs_ip)

  endif

  if(MCD15A2Hlaiobs(n)%climofill.eq.1) then

     write(LDT_logunit,*) '[INFO] Opening climo file ',trim(climofile)
     ios = nf90_open(path=trim(climofile),mode=NF90_NOWRITE,ncid=nid)


     call LDT_verify(ios,'Error opening file '//trim(climofile))

     ios = nf90_inq_varid(nid, 'Lai_500m',laiid)
     call LDT_verify(ios, 'Error nf90_inq_varid: Lai_500m')


     cornerlat(1)=MCD15A2Hlaiobs(n)%gridDesci(4)
     cornerlon(1)=MCD15A2Hlaiobs(n)%gridDesci(5)
     cornerlat(2)=MCD15A2Hlaiobs(n)%gridDesci(7)
     cornerlon(2)=MCD15A2Hlaiobs(n)%gridDesci(8)

     lai_data_b = .false.

     lat_off = nint((cornerlat(1)+89.9979167)/0.00416667)+1
     lon_off = nint((cornerlon(1)+179.9979167)/0.00416667)+1

     ios = nf90_get_var(nid, laiid, lai, &
          start=(/lon_off,lat_off/), &
          count=(/MCD15A2Hlaiobs(n)%nc,MCD15A2Hlaiobs(n)%nr/))


     call LDT_verify(ios, 'Error nf90_get_var: Lai_500m')

     ios = nf90_close(ncid=nid)
     call LDT_verify(ios,'Error closing file '//trim(fname))

     do r=1, MCD15A2Hlaiobs(n)%nr
        do c=1, MCD15A2Hlaiobs(n)%nc

           if(lai(c,r).gt.0.and.lai(c,r).le.100) then
              lai_flagged(c,r) = lai(c,r)*0.1
           else
              lai_flagged(c,r) = LDT_rc%udef
           endif

        end do
     end do


     do r=1, MCD15A2Hlaiobs(n)%nr
        do c=1, MCD15A2Hlaiobs(n)%nc
           lai_in(c+(r-1)*MCD15A2Hlaiobs(n)%nc) = lai_flagged(c,r)
           if(lai_flagged(c,r).ne.LDT_rc%udef) then
              lai_data_b(c+(r-1)*MCD15A2Hlaiobs(n)%nc) = .true.
           else
              lai_data_b(c+(r-1)*MCD15A2Hlaiobs(n)%nc) = .false.
           endif
        enddo
     enddo


 if(LDT_isLDTatAfinerResolution(n,0.00416667)) then

!--------------------------------------------------------------------------
! Interpolate to the LDT running domain
!-------------------------------------------------------------------------- 

     call bilinear_interp(LDT_rc%gridDesc(n,:),&
          lai_data_b, lai_in, laiobs_b_ip, laiobs_climo_ip, &
          MCD15A2Hlaiobs(n)%nc*MCD15A2Hlaiobs(n)%nr, &
          LDT_rc%lnc(n)*LDT_rc%lnr(n), &
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          MCD15A2Hlaiobs(n)%w11,MCD15A2Hlaiobs(n)%w12,&
          MCD15A2Hlaiobs(n)%w21,MCD15A2Hlaiobs(n)%w22,&
          MCD15A2Hlaiobs(n)%n11,MCD15A2Hlaiobs(n)%n12,&
          MCD15A2Hlaiobs(n)%n21,MCD15A2Hlaiobs(n)%n22,LDT_rc%udef,ios)
  else


     call upscaleByAveraging(MCD15A2Hlaiobs(n)%nc*MCD15A2Hlaiobs(n)%nr,&
          LDT_rc%lnc(n)*LDT_rc%lnr(n), &
          LDT_rc%udef, MCD15A2Hlaiobs(n)%n11,&
          lai_data_b,lai_in, laiobs_b_ip, laiobs_climo_ip)


  endif


     do t=1,LDT_rc%lnc(n)*LDT_rc%lnr(n)

        if(laiobs_ip(t).eq.-9999.0.and.laiobs_climo_ip(t).ne.-9999.0) then

           laiobs_ip(t) = laiobs_climo_ip(t)
        endif
     enddo
  endif

end subroutine read_MCD15A2Hlai_data

!BOP
! !ROUTINE: create_MCD15A2Hlai_filename
! \label{create_MCD15A2Hlai_filename}
! 
! !INTERFACE: 
subroutine create_MCD15A2Hlai_filename(ndir, version, yr, doy, filename, climofile)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  character(len=*)  :: climofile
  character(len=*)  :: version
  integer           :: yr, doy
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the MCD15A2H LAI filename
!  based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the MCD15A2H LAI data directory
!  \item[yr]  current year
!  \item[mo]  current doy
!  \item[filename] Generated MCD15A2H LAI filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=3) :: fdoy

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fdoy, fmt='(i3.3)') doy

  if(version.eq."006") then
     filename = trim(ndir)//'/'//trim(fyr)//'/MCD15A2H.006_LAI_'//&
          trim(fyr)//trim(fdoy)//'.nc4'
  endif
  climofile = trim(ndir)//'/MCD15A2H.006_LAI_YYYY'//&
       trim(fdoy)//'.nc4'

end subroutine create_MCD15A2Hlai_filename

