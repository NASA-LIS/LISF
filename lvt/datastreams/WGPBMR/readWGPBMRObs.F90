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
! !ROUTINE: readWGPBMRObs
! \label{readWGPBMRObs}
!
! !INTERFACE: 
subroutine readWGPBMRObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,      only : LVT_rc
  use LVT_histDataMod
  use LVT_logMod,       only : LVT_logunit, LVT_getNextUnitNumber, &
       LVT_releaseUnitNumber, LVT_endrun, LVT_verify
  use WGPBMRobsMod,      only : wgpbmrobs

  implicit none
!
! !INPUT PARAMETERS: 
  integer,   intent(in)  :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! This subroutine provides the data reader for Walnut Gulch PBMR data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  16 Feb 2008: Sujay Kumar, Initial Specification
! 
!EOP

  integer             :: ftn 
  integer             :: c,r
  logical             :: file_exists
  character*100       :: wgpbmrfilename
  real                :: dummy(LVT_rc%lnc, LVT_rc%lnr)
  real                :: smc(LVT_rc%lnc, LVT_rc%lnr)
  real                :: gridDesc(6)
  integer             :: istat

  smc     = LVT_rc%udef

  call create_WGPBMR_filename(wgpbmrobs(source)%odir, wgpbmrobs(source)%stnid, wgpbmrfilename, &
       LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source), LVT_rc%dhr(source), LVT_rc%dmn(source), LVT_rc%dss(source))
  
  inquire(file=trim(wgpbmrfilename), exist=file_exists) 

  if(file_exists.and.LVT_rc%dhr(source).eq.12.and.LVT_rc%dmn(source).eq.0) then 
     write(LVT_logunit,*) '[INFO] Reading WGPBMR file ',trim(wgpbmrfilename)
     ftn=LVT_getNextUnitNumber()
     open(ftn,file=trim(wgpbmrfilename),form='unformatted',status='old',&
          access='direct',recl=4, iostat = istat)
     
     gridDesc(1) = 12
     gridDesc(2) = 3507393.0
     gridDesc(3) = 586018.0
     gridDesc(4) = 660
     gridDesc(5) = 333
     gridDesc(6) = 40.0

     write(LVT_logunit,*) '[ERR] readWgpbmrobs(Source) needs to be updated '
     call LVT_endrun()
!     call LVT_readData(ftn,gridDesc,smc)

     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(smc(c,r).ne.-9999.0) then 
              smc(c,r) = smc(c,r)/100.0
           endif
        enddo
     enddo
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist, source, smc,vlevel=1,units="m3/m3")

  do c=2,LVT_rc%nsmlayers
     call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist, source, dummy,vlevel=c,units="m3/m3")
  enddo

end subroutine readWGPBMRObs

!BOP
! 
! !ROUTINE: create_WGPBMR_filename
! \label(create_WGPBMR_filename)
!
! !INTERFACE:
subroutine create_WGPBMR_filename(odir, stnid, wgpbmrname, &
     yr, mo, da, hr, mn, ss)
! 
! !USES:   
  use LVT_timeMgrMod, only : LVT_date2time  

  implicit none
!
! !INPUT PARAMETERS: 
  character(len=*), intent(in)  :: odir
  integer,          intent(in)  :: stnid
! 
! !OUTPUT PARAMETERS:
  character(len=*), intent(out) :: wgpbmrname
!
! !DESCRIPTION: 
! 
! This routine creates a filename for the WGPBMR station
! 
!  The arguments are: 
!  \begin{description}
!   \item[stnid] Station ID 
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  integer                       :: yr, mo, da, hr, mn, ss
!EOP
  
  character*1       :: fsite
  character*3       :: fdoy
  real*8            :: time
  integer           :: doy
  real              :: gmt
 
  call LVT_date2time(time, doy, gmt, yr, mo, da, hr, mn, ss)
  write(unit=fdoy,fmt='(i3.3)') doy
  write(unit=fsite,fmt='(i1.1)') stnid
  wgpbmrname = trim(odir)//'/pbmrsm'//trim(fdoy)//'_site'//trim(fsite)//'.1gd4r'
  
end subroutine create_WGPBMR_filename
