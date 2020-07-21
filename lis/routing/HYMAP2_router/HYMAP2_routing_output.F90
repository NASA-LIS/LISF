!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.1
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine HYMAP2_routing_output(n)
  
  ! Augusto Getirana - 11/15/2011
  
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_histDataMod
  use LIS_historyMod
  use LIS_fileIOMod
  use HYMAP2_routingMod

  use LIS_mpiMod
  implicit none
  
  integer, intent(in)   :: n 
  
  character(len=12)     :: cdate1
  integer               :: iret
  character*100         :: filename
  character*100         :: name
  integer               :: ftn
  integer               :: mo, da
  logical               :: open_stats
  logical               :: alarmCheck
  integer               :: status

  alarmCheck = .false. 
  if ( LIS_rc%time >= LIS_histData(n)%time ) then
!------------------------------------------------------------------
! Number of variables to be outputted from Noah LSM 
!------------------------------------------------------------------
  if ( LIS_rc%wsingle .ne. 1 ) then 
!------------------------------------------------------------------
! Writes bundled output
!------------------------------------------------------------------
     if(trim(LIS_rc%wopt).ne."none") then
        if(LIS_rc%output_at_specifictime.eq.1) then 
           if(LIS_histData(n)%month.eq.-1) then 
              mo = LIS_rc%mo
           else
              mo = LIS_histData(n)%month
           endif
           if(LIS_histData(n)%day.eq.-1) then 
              da = LIS_rc%da
           else
              da = LIS_histData(n)%day
           endif

           if(LIS_rc%mo.eq.mo.and.&
                LIS_rc%da.eq.da.and.&
                LIS_rc%hr.eq.LIS_histData(n)%hour.and.&
                LIS_rc%mn.eq.LIS_histData(n)%min.and.&
                LIS_rc%ss.eq.LIS_histData(n)%sec) then 
              alarmCheck = .true. 
           endif
        else
           alarmCheck = LIS_isAlarmRinging(LIS_rc,&
                "HYMAP2 router output alarm")
        endif
        if(alarmCheck) then 

           open_stats = .false.
           if(LIS_masterproc) then 
              HYMAP2_routing_struc(n)%numout=HYMAP2_routing_struc(n)%numout+1    
              call LIS_create_output_directory('ROUTING')
              call LIS_create_output_filename(n, filename, &
                   model_name='ROUTING', &
                   writeint=HYMAP2_routing_struc(n)%outInterval)

!-----------------------------------------------------------------------
! Open statistical output file
!-----------------------------------------------------------------------
              if(HYMAP2_routing_struc(n)%fileopen.eq.0)then
                 call LIS_create_stats_filename(n, name,'HYMAP2_routerstats')
                 HYMAP2_routing_struc(n)%fileopen=1
                 open_stats = .true.
              endif
           endif
     
!-----------------------------------------------------------------------
! Write Output 
!-----------------------------------------------------------------------
           ! Grib expects soil layers to be in cm.
           ! lyrthk = (/100.0,300.0,600.0,1000.0/) mm.
           call LIS_writeModelOutput(n,filename, name, open_stats,  &
                outInterval=HYMAP2_routing_struc(n)%outInterval,     &
                nsoillayers = 1,lyrthk = (/1.0/), nsoillayers2 = 1, &
                group=2)
        endif
     endif
  endif
endif

end subroutine HYMAP2_routing_output
