!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine NLDAS_routing_writerst(n)

  use ESMF
  use LIS_coreMod,      only : LIS_rc, LIS_masterproc
  use LIS_timeMgrMod,   only : LIS_isAlarmRinging
  use LIS_logMod,       only : LIS_logunit, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber
  use LIS_fileIOMod,   only : LIS_create_output_directory, &
       LIS_create_restart_filename
  use NLDAS_routingMod, only : NLDAS_routing_struc

  implicit none
  
  integer, intent(in)   :: n 
  
  character*100         :: filename
  integer               :: ftn
  integer               :: status
  logical               :: alarmCheck

  if(LIS_masterproc) then
     alarmCheck = LIS_isAlarmRinging(LIS_rc,&
          "NLDAS router restart alarm")
     if(alarmCheck.or.(LIS_rc%endtime ==1)) then 
        ftn = LIS_getNextUnitNumber()
        
        call LIS_create_output_directory('ROUTING')
        call LIS_create_restart_filename(n,filename,&
             'ROUTING','NLDAS_router',&
             wformat="binary")
        write(LIS_logunit,*) 'Writing routing restart ',trim(filename)

        open(ftn,file=trim(filename), form='unformatted')
     
        write(ftn) NLDAS_routing_struc(n)%runoff_intern
        write(ftn) NLDAS_routing_struc(n)%runoff_trans
        call LIS_releaseUnitNumber(ftn)
        
     endif
  endif
end subroutine NLDAS_routing_writerst
