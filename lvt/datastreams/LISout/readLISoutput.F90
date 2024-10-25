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
! !ROUTINE: readLISoutput
! \label(readLISoutput)
!
! !INTERFACE:
subroutine readLISoutput(source)
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_logMod
  use LVT_fileIOMod
  use LVT_timeMgrMod
  use LVT_statsDataMod
  use LVT_LISoutputHandlerMod

  implicit none
  
  integer, intent(in)    :: source
!
! !INPUT PARAMETERS: 

! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  
  character*500                    :: fname1 
  character*500, allocatable       :: fname(:)

  logical                          :: file_exists
  real                             :: obsData(LVT_rc%lnc, LVT_rc%lnr)
  
  integer                          :: t
  integer                          :: ftn
  type(LVT_metadataEntry), pointer :: obs
  integer                          :: c,r,k
  integer                          :: nfiles
  integer                          :: status
  integer                          :: yr,mo,da,hr,mn,ss
  type(ESMF_Time)                  :: cTime, stTime, enTime

  if(LVT_rc%chkTS) then 
     if(LVT_LIS_rc(source)%anlys_data_class.eq."LSM") then 
        if(LVT_rc%smoothObs.eq.1) then 

           call ESMF_TimeSet(cTime, yy = LVT_rc%dyr(source), &
                mm = LVT_rc%dmo(source), &
                dd = LVT_rc%dda(source), &
                h  = LVT_rc%dhr(source), &
                m  = LVT_rc%dmn(source), & 
                s  = LVT_rc%dss(source), &
                calendar = LVT_calendar, & 
                rc = status)
           call LVT_verify(status)
           
           stTime = cTime - LVT_obsSmTwL
           enTime = cTime + LVT_obsSmTwL
           
           nfiles = (enTime - stTime)/LVT_obsSmTwI + 1
           allocate(fname(nfiles))

           cTime = stTime
           k = 1
           do while (cTime.le.enTime) 
              call ESMF_TimeGet(cTime, yy = yr, mm=mo, dd=da,&
                   h = hr, m = mn, s =ss, calendar=LVT_calendar, &
                   rc=status)
              call LVT_verify(status)
              call LVT_create_output_filename(LVT_LIS_rc(source)%nest, &
                   source,fname(k), &
                   yr,mo,da,hr,mn,ss, &
                   'SURFACEMODEL',&
                   writeint = LVT_LIS_rc(source)%ts, &
                   wout = LVT_LIS_rc(source)%format,&
                   style=LVT_LIS_rc(source)%style, &
                   odir=LVT_LIS_rc(source)%odir)

              k = k + 1
              cTime = cTime + LVT_obsSmTwI
           enddo
        else
           call LVT_create_output_filename(LVT_LIS_rc(source)%nest, &
                source,fname1, &
                'SURFACEMODEL',&
! EMK Bug fix to LIS file name
!             writeint = LVT_rc%ts, wout = LVT_LIS_rc(source)%format,&
                writeint = LVT_LIS_rc(source)%ts, &
                wout = LVT_LIS_rc(source)%format,&
                style=LVT_LIS_rc(source)%style, &
                odir=LVT_LIS_rc(source)%odir)
        endif
     elseif(LVT_LIS_rc(source)%anlys_data_class.eq."Routing") then 
        call LVT_create_output_filename(LVT_LIS_rc(source)%nest, &
             source, fname1, &
             'ROUTING',&
             writeint = LVT_rc%ts, wout = LVT_LIS_rc(source)%format,&
             style=LVT_LIS_rc(source)%style,&
             odir=LVT_LIS_rc(source)%odir)
     elseif(LVT_LIS_rc(source)%anlys_data_class.eq."RTM") then 
        call LVT_create_output_filename(LVT_LIS_rc(source)%nest,&
             source, fname1, &
             'RTM',&
             writeint = LVT_rc%ts, wout = LVT_LIS_rc(source)%format,&
             style=LVT_LIS_rc(source)%style, &
             odir=LVT_LIS_rc(source)%odir)
     elseif(LVT_LIS_rc(source)%anlys_data_class.eq."Irrigation") then 
        call LVT_create_output_filename(LVT_LIS_rc(source)%nest, &
             source, fname1, &
             'IRRIGATION',&
             writeint = LVT_rc%ts, wout = LVT_LIS_rc(source)%format,&
             style=LVT_LIS_rc(source)%style,&
             odir=LVT_LIS_rc(source)%odir)
     endif
     
     if(LVT_rc%smoothObs.eq.1) then 

        do k=1,nfiles
           write(LVT_logunit,*) '[INFO] Reading LIS output ',trim(fname(k))
        enddo
        
        call LVT_readLISModelOutput(nfiles,&
             fname,source, &
             LVT_LIS_rc(source)%format,&
             wopt=LVT_LIS_rc(source)%wopt)

     else
        inquire(file=trim(fname1),exist=file_exists)
        
        if(file_exists) then 
        ! The following line is modified by Shugong Wang: LVT should be LIS
        !write(LVT_logunit,*) 'Reading LVT output ',trim(fname1)
           write(LVT_logunit,*) '[INFO] Reading LIS output ',trim(fname1)
           call LVT_readLISModelOutput(trim(fname1),source, &
                LVT_LIS_rc(source)%format,&
                wopt=LVT_LIS_rc(source)%wopt)
        else
           write(LVT_logunit,*) '[WARN] LIS file ',&
                trim(fname1),' does not exist'
        endif
        
     endif
  endif

end subroutine readLISoutput

