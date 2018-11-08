!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.1
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine HYMAP2_routing_writerst(n)


  ! Augusto Getirana - 11/15/2011
  ! 19 Jan 2016: Augusto Getirana, Inclusion of four Local Inertia variables

  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_fileIOMod
  use LIS_timeMgrMod
  use HYMAP2_routingMod

  implicit none
  
  integer, intent(in)   :: n 
  
  character*100         :: filename
  integer               :: ftn
  integer               :: status
  logical               :: alarmCheck

  if(LIS_masterproc) then
     alarmCheck = LIS_isAlarmRinging(LIS_rc,&
          "HYMAP2 router restart alarm")
     if(alarmCheck.or.(LIS_rc%endtime ==1)) then 
        ftn = LIS_getNextUnitNumber()
        call LIS_create_output_directory('ROUTING')
        call LIS_create_restart_filename(n,filename,&
             'ROUTING','HYMAP2_router',&
             wformat="binary")
        write(LIS_logunit,*) 'Writing routing restart ',trim(filename)

        open(ftn,file=trim(filename), form='unformatted')
		
        write(ftn) HYMAP2_routing_struc(n)%rivsto
        write(ftn) HYMAP2_routing_struc(n)%fldsto
        write(ftn) HYMAP2_routing_struc(n)%rnfsto
        write(ftn) HYMAP2_routing_struc(n)%bsfsto
        write(ftn) HYMAP2_routing_struc(n)%rivout_pre
        write(ftn) HYMAP2_routing_struc(n)%rivdph_pre
        write(ftn) HYMAP2_routing_struc(n)%fldout_pre
        write(ftn) HYMAP2_routing_struc(n)%flddph_pre

        call LIS_releaseUnitNumber(ftn)
        
     endif
  endif
end subroutine HYMAP2_routing_writerst
