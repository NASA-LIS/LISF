!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! 
! !ROUTINE: readFMISWEobs
! \label{readFMISWEobs}
!
! !INTERFACE: 
subroutine readFMISWEobs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,    only : LVT_rc
  use LVT_histDataMod
  use LVT_timeMgrMod, only : LVT_tick
  use LVT_logMod,     only : LVT_verify
  use FMISWE_obsMod,    only : fmisweobs
  use LVT_logMod,  only : LVT_verify, LVT_endrun, LVT_logunit,&
       LVT_getNextUnitNumber, LVT_releaseUnitNumber
  
  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in)  :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for FMI SWE data. 
! LVT expects the data to be organized as daily timestamped files,
! Each reported observation is assumed to be time averaged. 
!
! The data is then interpolated using the neighbor approach
! to the LIS output grid. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  16 Feb 2008: Sujay Kumar, Initial Specification
! 
!EOP
  integer                :: ftn
  character*100          :: fname
  integer                :: doy
  real                   :: gmt
  integer                :: hr, mn, ss,ts
  real*8                 :: time
  logical                :: file_exists
  integer                :: c,r
  integer                :: iret
  logical*1              :: lb(fmisweobs(source)%nc*fmisweobs(source)%nr)
  logical*1              :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                   :: swe(LVT_rc%lnc, LVT_rc%lnr)
  real                   :: swe1(fmisweobs(source)%nc,fmisweobs(source)%nr)
  
  swe = LVT_rc%udef
  if(fmisweobs(source)%startFlag) then 
    fmisweobs(source)%yr = LVT_rc%dyr(source)
    fmisweobs(source)%mo = LVT_rc%dmo(source)
    fmisweobs(source)%da = LVT_rc%dda(source)

!set time back by one day.
    hr = 0
    mn = 0 
    ss = 0 
    ts = -86400
    call LVT_tick(time, doy, gmt, fmisweobs(source)%yr, fmisweobs(source)%mo, fmisweobs(source)%da, &
         hr,mn,ss, ts) 
  endif


  if(fmisweobs(source)%startFlag.or.(fmisweobs(source)%da.ne.&
       LVT_rc%dda(source).and.LVT_rc%dhr(source).eq.0)) then !new day

     fmisweobs(source)%startFlag = .false. 

     call create_FMISWE_filename(fmisweobs(source)%odir, fmisweobs(source)%yr, fmisweobs(source)%mo, &
          fmisweobs(source)%da, fname)
     inquire(file=trim(fname), exist=file_exists) 
     if(file_exists) then 
        ftn = LVT_getNextUnitNumber()
        write(LVT_logunit,*) '[INFO] Reading FMISWE ',trim(fname)
        open(ftn,file=trim(fname),form='unformatted')
        read(ftn) swe1

        lb = .false. 
        do r=1,fmisweobs(source)%nr
           do c=1,fmisweobs(source)%nc
              if(swe1(c,r).ge.0) then 
                 lb(c+(r-1)*fmisweobs(source)%nc) = .true.                  
              endif
           enddo
        enddo

        call neighbor_interp(LVT_rc%gridDesc,lb,swe1,lo,swe,&
             fmisweobs(source)%nc*fmisweobs(source)%nr,LVT_rc%lnc*LVT_rc%lnr,&
             fmisweobs(source)%rlat,fmisweobs(source)%rlon,&
             fmisweobs(source)%n11,LVT_rc%udef,iret)

       call LVT_releaseUnitNumber(ftn)
     else
        write(LVT_logunit,*) '[ERR] File ',trim(fname), ' does not exist'
        swe = LVT_rc%udef
     endif
! advance by one day  
     hr = 0
     mn = 0 
     ss = 0 
     ts = 86400
     call LVT_tick(time, doy, gmt, fmisweobs(source)%yr, &
          fmisweobs(source)%mo, fmisweobs(source)%da, &
          hr,mn,ss,ts)
  else
     swe = LVT_rc%udef
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_SWE,source,swe,vlevel=1,units="m")

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(swe(c,r).ne.-9999.0) then 
           swe(c,r) = swe(c,r)*1000.0
        endif
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_SWE,source,swe,vlevel=1,units="kg/m2")

end subroutine readFMISWEobs



!BOP
! 
! !ROUTINE: create_FMISWE_filename
! \label{create_FMISWE_filename}
!
! !INTERFACE: 
subroutine create_FMISWE_filename(odir, yr, mo, da, fmiswename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
  character(len=*), intent(in)  :: odir
  integer,          intent(in)  :: yr
  integer,          intent(in)  :: mo
  integer,          intent(in)  :: da
! 
! !OUTPUT PARAMETERS:
  character(len=*), intent(out) :: fmiswename
!
! !DESCRIPTION: 
! This routine creates a filename for the FMISWE station
! 
!  The arguments are: 
!  \begin{description}
!   \item[odir] FMI snow data base directory
!   \item[yr]   year of data 
!   \item[mo]   month of data 
!   \item[da]   day of data 
!   \item[fmiswename]  Name of the FMI snow file  
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  character*4             :: fyr
  character*2             :: fmo
  character*2             :: fda
  
  write(fyr, '(i4.4)' ) yr
  write(fmo, '(i2.2)' ) mo
  write(fda, '(i2.2)' ) da

  fmiswename = trim(odir)//'/FMI_SWE_'//trim(fyr)//trim(fmo)//trim(fda)//'.bin'
  
end subroutine create_FMISWE_filename
