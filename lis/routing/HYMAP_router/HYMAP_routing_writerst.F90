!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine HYMAP_routing_writerst(n)


  ! Augusto Getirana - 11/15/2011

  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_fileIOMod
  use LIS_timeMgrMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use HYMAP_routingMod

  implicit none
  
  integer, intent(in)   :: n 
  
  character(len=LIS_CONST_PATH_LEN) :: filename
  integer               :: ftn
  integer               :: status
  logical               :: alarmCheck

  if(LIS_masterproc) then
     alarmCheck = LIS_isAlarmRinging(LIS_rc,&
          "HYMAP router restart alarm")
     if(alarmCheck.or.(LIS_rc%endtime ==1)) then 
        ftn = LIS_getNextUnitNumber()
        call LIS_create_output_directory('ROUTING')
        call LIS_create_restart_filename(n,filename,&
             'ROUTING','HYMAP_router',&
             wformat="binary")
        write(LIS_logunit,*) '[INFO] Writing routing restart ',trim(filename)

        open(ftn,file=trim(filename), form='unformatted')

        write(ftn) HYMAP_routing_struc(n)%rivsto
        write(ftn) HYMAP_routing_struc(n)%fldsto
        write(ftn) HYMAP_routing_struc(n)%rnfsto
        write(ftn) HYMAP_routing_struc(n)%bsfsto

        call LIS_releaseUnitNumber(ftn)
        
     endif
  endif
end subroutine HYMAP_routing_writerst
