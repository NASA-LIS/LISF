!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readMCD15A2Hobs
! \label{readMCD15A2Hobs}
!
! !INTERFACE:
subroutine readMCD15A2Hobs(source)
! 
! !USES:   
  use LVT_coreMod,    only : LVT_rc
  use LVT_logMod
  use LVT_histDataMod
  use MCD15A2H_obsMod, only : mcd15a2hobs
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
  integer,  intent(in)      :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  
!  This subroutine provides the data reader for the 500m MCD15A2H LAI
!  data . The routine also spatially
!  interpolates the MCD15A2H data to the analysis grid 
!  and resolution. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  20 Oct 2020: Sujay Kumar, Initial Specification
! 
!EOP

  integer, parameter           :: nc = 36000, nr=15000
  integer                      :: lat_off, lon_off
  character*100                :: filename
  real*8                       :: cornerlat(2), cornerlon(2)
  logical                      :: file_exists
  logical*1, allocatable       :: lb(:)
  real, allocatable            :: lai1(:,:), lai_flagged(:,:),lai_in(:)
  integer, allocatable         :: flag(:,:)
  real                         :: lai(LVT_rc%lnc, LVT_rc%lnr)
  integer                      :: ftn, c, r
  integer                      :: laiid, flagid
  logical                      :: alarmcheck
  real                         :: timenow

  lai = -9999.0
  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) + &
     LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)
  if(mcd15a2hobs(source)%startflag.or.alarmCheck) then 
     
     mcd15a2hobs(source)%startflag = .false. 
     
     call create_mcd15a2hobs_filename(mcd15a2hobs(source)%odir,&
          LVT_rc%dyr(source), LVT_rc%ddoy(source),&
          filename)
     
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
     inquire(file=trim(filename),exist=file_exists)
     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading MCD15A2H file ',trim(filename)
        

        allocate(lai1(mcd15a2hobs(source)%nc,mcd15a2hobs(source)%nr))
        allocate(flag(mcd15a2hobs(source)%nc,mcd15a2hobs(source)%nr))
        allocate(lai_flagged(mcd15a2hobs(source)%nc,mcd15a2hobs(source)%nr))

        call LVT_verify(nf90_open(path=trim(filename),mode=NF90_NOWRITE, &
           ncid=ftn), 'Error in nf90_open')
        
        call LVT_verify(nf90_inq_varid(ftn,'Lai_500m',laiid),&
           'nf90_inq_varid failed for Lai_500m')

        call LVT_verify(nf90_inq_varid(ftn,'FparLai_QC',flagid),&
           'nf90_inq_varid failed for FparLai_Qc')

        cornerlat(1)=MCD15A2Hobs(source)%gridDesci(4)
        cornerlon(1)=MCD15A2Hobs(source)%gridDesci(5)
        cornerlat(2)=MCD15A2Hobs(source)%gridDesci(7)
        cornerlon(2)=MCD15A2Hobs(source)%gridDesci(8)
  
        lat_off = nint((cornerlat(1)+89.99979167)/0.00416667)+1
        lon_off = nint((cornerlon(1)+179.9979167)/0.00416667)+1
           
        call LVT_verify(nf90_get_var(ftn,laiid,lai1, &
           start=(/lon_off, lat_off/),&
           count=(/mcd15a2hobs(source)%nc, mcd15a2hobs(source)%nr/)),&
           'Error in nf90_get_var: Lai_500m')
           
        call LVT_verify(nf90_get_var(ftn, flagid, flag, &
           start=(/lon_off,lat_off/), &
           count=(/mcd15a2hobs(source)%nc, mcd15a2hobs(source)%nr/)),&
           'Error in nf90_get_var: flag')
                
        call LVT_verify(nf90_close(ftn), 'Error in nf90_close')
#endif        
        do r=1,mcd15a2hobs(source)%nr
           do c=1,mcd15a2hobs(source)%nc
              if(lai1(c,r).gt.0.or.lai1(c,r).le.100) then 
                 if (MOD(flag(c,r),2) ==0.and.flag(c,r).le.62) then
                    lai_flagged(c,r) =&
                         lai1(c,r)*0.1
                 else
                    lai_flagged(c,r) = LVT_rc%udef
                 endif
              else
                 lai_flagged(c,r) = LVT_rc%udef
              endif
           enddo
        enddo
        deallocate(lai1)
        deallocate(flag)

        allocate(lb(mcd15a2hobs(source)%nc*mcd15a2hobs(source)%nr))
        allocate(lai_in(mcd15a2hobs(source)%nc*mcd15a2hobs(source)%nr))
        
        lb = .true. 

        do r=1, mcd15a2hobs(source)%nr
           do c=1, mcd15a2hobs(source)%nc
              lai_in(c+(r-1)*mcd15a2hobs(source)%nc) = lai_flagged(c,r)
              if(lai_flagged(c,r).ne.LVT_rc%udef) then
                 lb(c+(r-1)*mcd15a2hobs(source)%nc) = .true.
              else
                 lb(c+(r-1)*mcd15a2hobs(source)%nc) = .false.
              endif
           enddo
        enddo

        call interp_mcd15a2h(source, &
             mcd15a2hobs(source)%nc, mcd15a2hobs(source)%nr, &
             lai_in, lb, lai)

        deallocate(lb)
        deallocate(lai_in)
     else
        lai = -9999
     endif
  else
     lai = -9999.0
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_LAI,source,lai,&
       vlevel=1,units="-")


end subroutine readMCD15A2Hobs

!BOP
! 
! !ROUTINE: interp_mcd15a2h
!  \label{interp_mcd15a2h}
!
! !INTERFACE: 
subroutine interp_mcd15a2h(source, nc, nr, var_input, lb, var_output)
! 
! !USES:   
  use LVT_coreMod
  use MCD15A2H_obsMod, only : mcd15a2hobs

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!   This subroutine spatially interpolates the MCD15A2H variable to the
!   analysis grid and resolution.
! 
!   The arguments are: 
!   \begin{description}
!    \item[nc]      number of columns in the input (MCD15A2H) grid
!    \item[nr]      number of rows in the input (MCD15A2H) grid
!    \item[var_input] input variable to be interpolated
!    \item[lb]        input bitmap (true//false)
!    \item[var_output] resulting interplated field
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS:   
  integer            :: source
  integer            :: nc
  integer            :: nr
  real               :: var_input(nc*nr)
  logical*1          :: lb(nc*nr)
  real               :: var_output(LVT_rc%lnc, LVT_rc%lnr)
!EOP
  integer            :: iret
  integer            :: c,r
  logical*1          :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real               :: go(LVT_rc%lnc*LVT_rc%lnr)

  var_output = LVT_rc%udef
  if(LVT_isAtAfinerResolution(0.00416667)) then
     call bilinear_interp(LVT_rc%gridDesc,lb, var_input,&
          lo,go,nc*nr,LVT_rc%lnc*LVT_rc%lnr,&
          mcd15a2hobs(source)%rlat,mcd15a2hobs(source)%rlon,&
          mcd15a2hobs(source)%w11,mcd15a2hobs(source)%w12,&
          mcd15a2hobs(source)%w21,mcd15a2hobs(source)%w22,&
          mcd15a2hobs(source)%n11,mcd15a2hobs(source)%n12,&
          mcd15a2hobs(source)%n21,mcd15a2hobs(source)%n22,LVT_rc%udef,iret)
     
  else
     call upscaleByAveraging(&
          mcd15a2hobs(source)%nc*mcd15a2hobs(source)%nr, &
          LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
          mcd15a2hobs(source)%n11, lb, &
          var_input, lo, go)
  endif

  do r=1, LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(go(c+(r-1)*LVT_rc%lnc).ne.-9999.0) then 
           var_output(c,r) = go(c+(r-1)*LVT_rc%lnc)/100.0
        endif
     enddo
  enddo
  
end subroutine interp_mcd15a2h

!BOP
! 
! !ROUTINE: create_mcd15a2hobs_filename
! \label{create_mcd15a2hobs_filename}
!
! !INTERFACE: 
subroutine create_mcd15a2hobs_filename(odir, yr, doy, filename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine creates the MCD15A2H filename based on the given 
!  date (year, month, day and hour)
! 
!  The arguments are: 
!  \begin{description}
!   \item[odir]      MCD15A2H base directory
!   \item[yr]        year of data
!   \item[mo]        month of data
!   \item[da]        day of data
!   \item[filename]  Name of the MCD15A2Hv-c6 file
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS:  
  character(len=*)  :: odir
  integer           :: yr
  integer           :: doy
  character(len=*)  :: filename
! 
!EOP

  character*4       :: fyr
  character*3       :: fdoy

  write(unit=fyr,  fmt='(i4.4)') yr
  write(unit=fdoy, fmt='(i3.3)') doy

  filename = trim(odir)//'/MCD15A2H.006_LAI_'//trim(fyr)//trim(fdoy)//'.nc4'

end subroutine create_mcd15a2hobs_filename
  
