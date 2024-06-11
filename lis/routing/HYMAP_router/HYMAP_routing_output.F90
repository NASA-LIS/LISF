!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine HYMAP_routing_output(n)
  
  ! Augusto Getirana - 11/15/2011
  
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_histDataMod
  use LIS_historyMod
  use LIS_fileIOMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use HYMAP_routingMod
#if ( defined USE_PFIO )
      use LIS_PFIO_historyMod
#endif 

!
! !HISTORY:
! 01 Sep 2023 Jules Kouatchou; Introduce preprocessing directives for calls
!             of HISTORY related subroutines with and without PFIO components.

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
  integer               :: vcol_id

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
                "HYMAP router output alarm")
        endif

        if(alarmCheck) then 

#if ( defined USE_PFIO )
           HYMAP_routing_struc(n)%numout=HYMAP_routing_struc(n)%numout+1    
           IF (PFIO_bundle%first_time(n, 1, PFIO_ROUTING_idx)) THEN
              ! Create the file metadata ONCE.
              call LIS_PFIO_create_file_metadata(n, PFIO_ROUTING_idx, HYMAP_routing_struc(n)%outInterval, &
                                   1, (/1.0/), group=2)
              PFIO_bundle%first_time(n, :, PFIO_ROUTING_idx) = .FALSE.
           ENDIF
           call LIS_create_output_directory('ROUTING')
           call LIS_create_output_filename(n, outfile,&
                            model_name = 'ROUTING',&
                            writeint=HYMAP_routing_struc(n)%outInterval)
           vcol_id = MOD(PFIO_bundle%counter(n, PFIO_ROUTING_idx)-1, LIS_rc%n_vcollections) + 1
           CALL LIS_PFIO_write_data(n, PFIO_ROUTING_idx, vcol_id, &
                                outfile, HYMAP_routing_struc(n)%outInterval)
    
           PFIO_bundle%counter(n, PFIO_ROUTING_idx) = PFIO_bundle%counter(n, PFIO_ROUTING_idx) + 1
#else
           open_stats = .false.
           if(LIS_masterproc) then 
              HYMAP_routing_struc(n)%numout=HYMAP_routing_struc(n)%numout+1    
              call LIS_create_output_directory('ROUTING')

!-----------------------------------------------------------------------
! Open statistical output file
!-----------------------------------------------------------------------
              if(HYMAP_routing_struc(n)%fileopen.eq.0)then
                 call LIS_create_stats_filename(n, name,'HYMAP_routerstats')
                 HYMAP_routing_struc(n)%fileopen=1
                 open_stats = .true.
              endif
           endif

           call LIS_create_output_filename(n, filename, &
                model_name='ROUTING', &
                writeint=HYMAP_routing_struc(n)%outInterval)

!-----------------------------------------------------------------------
! Write Output 
!-----------------------------------------------------------------------
           ! Grib expects soil layers to be in cm.
           ! lyrthk = (/100.0,300.0,600.0,1000.0/) mm.

           call LIS_writeModelOutput(n,filename, name, open_stats,  &
                outInterval=HYMAP_routing_struc(n)%outInterval,     &
                nsoillayers = 1,lyrthk = (/1.0/), nsoillayers2 = 1, &
                group=2)
                         
#endif
        endif
     endif
  endif
endif

end subroutine HYMAP_routing_output
