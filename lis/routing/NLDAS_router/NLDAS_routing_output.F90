!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine NLDAS_routing_output(n)

  use ESMF
  use LIS_coreMod,      only : LIS_rc, LIS_masterproc
  use LIS_timeMgrMod,   only : LIS_isAlarmRinging
  use LIS_logMod,       only : LIS_logunit, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber
  use LIS_histDataMod,  only : LIS_histData
  use LIS_historyMod,   only : LIS_writeModelOutput
  use LIS_fileIOMod,   only : LIS_create_output_directory, & 
       LIS_create_stats_filename, LIS_create_output_filename
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use NLDAS_routingMod, only : NLDAS_routing_struc

  implicit none
  
  integer, intent(in)   :: n 
  
  character(len=12)     :: cdate1
  integer               :: iret
  character(len=LIS_CONST_PATH_LEN) :: filename
  character*100         :: name
  integer               :: ftn
  integer               :: mo, da
  logical               :: open_stats
  logical               :: alarmCheck

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
                "NLDAS router output alarm")
        endif
        if(alarmCheck) then 

           open_stats = .false.
           if(LIS_masterproc) then 
              NLDAS_routing_struc(n)%numout=NLDAS_routing_struc(n)%numout+1    
              call LIS_create_output_directory('ROUTING')

!-----------------------------------------------------------------------
! Open statistical output file
!-----------------------------------------------------------------------
              if(NLDAS_routing_struc(n)%fileopen.eq.0)then
                 call LIS_create_stats_filename(n, name,'NLDAS_routerstats')
                 NLDAS_routing_struc(n)%fileopen=1
                 open_stats = .true.
              endif
           endif

           call LIS_create_output_filename(n, filename, &
                model_name='ROUTING', &
                writeint=NLDAS_routing_struc(n)%outInterval)

!-----------------------------------------------------------------------
! Write Output 
!-----------------------------------------------------------------------
           ! Grib expects soil layers to be in cm.
           ! lyrthk = (/100.0,300.0,600.0,1000.0/) mm.

           call LIS_writeModelOutput(n,filename, name, open_stats,  &
                outInterval=NLDAS_routing_struc(n)%outInterval,     &
                nsoillayers = 1,lyrthk = (/1.0/), nsoillayers2 = 1, &
                group=2)
                         
        endif
     endif
  endif
endif

end subroutine NLDAS_routing_output
