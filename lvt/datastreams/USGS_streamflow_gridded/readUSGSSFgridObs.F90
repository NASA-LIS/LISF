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
! !ROUTINE: readUSGSSFgridObs
! \label{readUSGSSFgridObs}
!
! !INTERFACE: 
  subroutine readUSGSSFgridObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_logMod
  use USGSSFgrid_obsMod, only : USGSSFgridobs
  use map_utils

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in) :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  11 May 2011: Sujay Kumar, Initial Specification
! 
!EOP
!--------------------------------------------------------------------------



  character*100       :: filename
  logical             :: file_exists
  integer             :: data_index,ftn,k,c,r
  real                :: q(LVT_rc%lnc,LVT_rc%lnr)
  integer             :: days(12)
  data days /31,28,31,30,31,30,31,31,30,31,30,31/

  q = LVT_rc%udef
!every new year read the data, for each station and store it in memory

  if(USGSSFgridobs(source)%startFlag.eq.-1) then 
     USGSSFgridobs(source)%startFlag = 1
     USGSSFgridobs(source)%q = LVT_rc%udef

     filename = trim(USGSSFgridobs(source)%odir)//"/USGS_Q.bin"
     write(LVT_logunit,*) &
          '[INFO] reading USGS streamflow file ',trim(filename)
     inquire(file=trim(filename), exist=file_exists)    
     if(file_exists) then            
        ftn = LVT_getNextUnitNumber()
        open(ftn,file=trim(filename),form='unformatted')
        do k=1,USGSSFgridobs(source)%nts
           read(ftn) USGSSFgridobs(source)%q(:,:,k)
        enddo
        call LVT_releaseUnitNumber(ftn)
     endif
  endif

  if((USGSSFgridobs(source)%mo.ne.LVT_rc%d_nmo(source)).or.&
       LVT_rc%resetFlag(source)) then

     LVT_rc%resetFlag(source) = .false. 
     USGSSFgridobs(source)%mo = LVT_rc%d_nmo(source)
     call findUSGSSFgridDataIndex(LVT_rc%dyr(source),&
          LVT_rc%dmo(source),data_index)

     q(:,:) = USGSSFgridobs(source)%q(:,:,data_index)
     do r=1,USGSSFgridobs(source)%nr
        do c=1,USGSSFgridobs(source)%nc
           if(q(c,r).lt.0.0) then 
              q(c,r) = LVT_rc%udef
           else
              q(c,r) = q(c,r)/(86400.0* days(LVT_rc%dmo(source)))  !mm/s
           endif
        enddo
     enddo
     
  else
     q = LVT_rc%udef
          
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_RUNOFF,source, q,vlevel=1,units="kg/m2s")

end subroutine readUSGSSFgridObs

subroutine findUSGSSFgridDataIndex(eyr,emo,data_index)
  
  use ESMF
  use LVT_timeMgrMod

  integer         :: eyr
  integer         :: emo
  integer         :: data_index

  integer         :: syr
  integer         :: smo
  integer         :: yr
  integer         :: mo
  integer         :: status
  logical         :: togo
  type(ESMF_Time) :: startTime, endTime, currTime

  data_index = 1
  
  syr = 1979
  smo = 1

  togo = .true. 
  
  call ESMF_TimeSet(startTime, yy=syr,&
       mm= smo,&
       calendar = LVT_calendar, &
       rc=status)

  call ESMF_TimeSet(endTime, yy=eyr,&
       mm= emo,&
       calendar = LVT_calendar, &
       rc=status)
  
  yr = syr
  mo = smo

  do while(togo) 
     call ESMF_TimeSet(currTime, yy=yr,&
          mm= mo,&
          calendar = LVT_calendar, &
          rc=status)

     if(currTime.ge.endTime) then 
        togo = .false. 
     else
        mo = mo +1
        do while(mo.gt.12)
           mo = mo-12
           yr = yr +1
        enddo
        data_index = data_index + 1
     endif
  end do
end subroutine findUSGSSFgridDataIndex
