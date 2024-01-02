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
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
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
!BOP
! 
! !ROUTINE: readCMC_SNWDObs
! \label{readCMC_SNWDObs}
! 
! !REVISION HISTORY: 
!  23 APR 2010: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readCMC_SNWDObs(source)
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_logMod
  use LVT_timeMgrMod
  use CMCSNWD_obsMod,    only : cmcsnwdobs
  use map_utils

  implicit none
! !ARGUMENTS: 

! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for CMC snow depth data. 
! LVT expects the data to be organized per calendar year, with 
! each file containing a daily data. Each reported observation is
! assumed to be time averaged. 
!
! The yearly file is read until the data for the current time is 
! found. The data is then interpolated using the neighbor approach
! to the LIS output grid. 
! 
!EOP
  integer                :: source
  integer                :: yr, mo, da, hr
  integer                :: i,j,t,c,r
  integer                :: stn_col, stn_row
  real                   :: col,row
  character*100          :: cmcsnwdfilename
  logical                :: file_exists
  logical                :: readflag
  integer                :: ftn, ios
  integer                :: status
  type(ESMF_Time)        :: cmcsnwdtime, cmcsnwdtime1
  integer                :: stnindex,tind
  real                   :: offset
  real                   :: snwd_data
  integer                :: iret
  logical*1              :: lb(cmcsnwdobs(source)%nc*cmcsnwdobs(source)%nr)
  logical*1              :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                   :: snwd(cmcsnwdobs(source)%nc, cmcsnwdobs(source)%nr)
  real                   :: snwd1d(cmcsnwdobs(source)%nc*cmcsnwdobs(source)%nr)
  real                   :: snwd_ip(LVT_rc%lnc,LVT_rc%lnr)
  
  snwd_ip  = LVT_rc%udef

  call ESMF_TimeSet(cmcsnwdtime1,yy=LVT_rc%dyr(source), &
       mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), &
       h=LVT_rc%dhr(source), m=LVT_rc%dmn(source), &
       s = LVT_rc%dss(source), calendar=LVT_calendar, rc=status)
  call LVT_verify(status, 'cmcsnwdtime1 set failed')

  offset = (cmcsnwdtime1-cmcsnwdobs(source)%starttime)/&
       cmcsnwdobs(source)%timestep

  if((nint(offset)-offset).eq.0) then !only when LIS time matches the observation
  
     call create_CMCSNWD_filename(cmcsnwdobs(source)%odir, LVT_rc%dyr(source), &
          LVT_rc%dmo(source), LVT_rc%dda(source), cmcsnwdfilename)
     
     inquire(file=trim(cmcsnwdfilename),exist=file_exists)
     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading CMCSNWD data ',&
             trim(cmcsnwdfilename)
        
        ftn=LVT_getNextUnitNumber()
        open(ftn,file=trim(cmcsnwdfilename),form='formatted')
        
        readflag = .true. 
        do while(readflag) 
           read(ftn,fmt='(1X,I4.4,3I2.2)', iostat=ios) yr, mo, da, hr
           do r=1,cmcsnwdobs(source)%nr
              read(ftn,fmt='(706F6.2)', iostat=ios) snwd(:,r)
           enddo
           
           if(yr.eq.LVT_rc%dyr(source).and.mo.eq.&
                LVT_rc%dmo(source).and.da.eq.LVT_rc%dda(source)&
                .and.hr.eq.LVT_rc%dhr(source)) then 
              lb = .true. 
              do r=1,cmcsnwdobs(source)%nr
                 do c=1,cmcsnwdobs(source)%nc
                    snwd1d(c+(r-1)*cmcsnwdobs(source)%nc) = snwd(c,r)/100 !cm to m
                    if(snwd(c,r).lt.0.or.snwd(c,r).gt.990) &
                         lb(c+(r-1)*cmcsnwdobs(source)%nc) = .false. 
                 enddo
              end do
              call neighbor_interp(LVT_rc%gridDesc,lb,snwd1d,lo,snwd_ip,&
                   cmcsnwdobs(source)%nc*cmcsnwdobs(source)%nr,&
                   LVT_rc%lnc*LVT_rc%lnr,&
                   cmcsnwdobs(source)%rlat,cmcsnwdobs(source)%rlon,&
                   cmcsnwdobs(source)%n11,LVT_rc%udef,iret)
              readflag = .false. 
           endif
        enddo
        call LVT_releaseUnitNumber(ftn)
        write(LVT_logunit,*) '[INFO] Successfully processed ',&
        trim(cmcsnwdfilename)
     endif
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_snowdepth,source,&
       snwd_ip,vlevel=1,units="m")

end subroutine readCMC_SNWDObs

!BOP
! 
! !ROUTINE: create_CMCSNWD_filename
!  \label(create_CMCSNWD_filename)
!
! !INTERFACE: 
subroutine create_CMCSNWD_filename(odir, yr,mo,da, cmcsnwdname)
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
  character(len=*), intent(out) :: cmcsnwdname
!
! !DESCRIPTION: 
! 
! This routine creates a filename for the CMC snow analysis file
! 
!  The arguments are: 
!  \begin{description}
!   \item[odir] CMC snow data base directory
!   \item[yr]   year of data 
!   \item[cmcsnwdname]  Name of the CMC snow analysis file  
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

  cmcsnwdname = trim(odir)//'/'//trim(fyr)//'/cmc_analysis_'&
       //trim(fyr)//trim(fmo)//trim(fda)//'.txt'
  
end subroutine create_CMCSNWD_filename
