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
! !MODULE: GRDC_obsMod
! \label(GRDC_obsMod)
!
! !INTERFACE:
module GRDC_obsMod
! 
! !USES:   
  use ESMF

  implicit none

  PRIVATE 

  PUBLIC :: GRDC_obsinit
  PUBLIC :: GRDCobs

  type, public :: GRDCobsdec
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!
!  The Global Runoff Data Centre (GRDC), a repository for !  the world's river discharge data. The Global Runoff Database 
!  is a unique collection of river discharge data collected at 
!  daily or monthly intervals from more than 9,300 stations in 
!  160 countries. This adds up to around 400,000 station-years 
!  with an average record length of 43 years. The GRDC provides 
!  discharge data and data products for non-commercial applications.
!
!  The GRDC operates under the auspices of the World Meteorological 
!  Organisation (WMO), and the German Federal Institute of Hydrology 
!  (Bundesanstalt für Gewässerkunde or BfG) hosts the GRDC in Koblenz. 
!
!  http://www.bafg.de/GRDC/EN/Home/homepage_node.html
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  13 May 2011   Sujay Kumar  Initial Specification
! 
!EOP

     integer                     :: version
     character*100               :: odir
     integer                     :: use_daily
     integer                     :: n_stns
     integer                     :: nts
     character*100, allocatable  :: stn_id(:)
     real,          allocatable  :: stnlat(:)
     real,          allocatable  :: stnlon(:)
     integer                     :: yr
     integer                     :: mo
     real,          allocatable  :: q(:,:)
     type(ESMF_Time)             :: startTime
     type(ESMF_TimeInterval)     :: timestep
  end type GRDCobsdec

  type(GRDCobsdec), allocatable :: GRDCobs(:)

contains
  
!BOP
! 
! !ROUTINE: GRDC_obsInit
! \label{GRDC_obsInit}
!
! !INTERFACE: 
 subroutine GRDC_obsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_logMod
    use LVT_timeMgrMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
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

    character*100 :: stnlist_file
    character*200 :: currentLine
    character*20  :: data_freq
    integer       :: iloc
    integer       :: ftn, k, status

    if(.not.allocated(GRDCobs)) then 
       allocate(GRDCobs(LVT_rc%nDataStreams))
    endif

!------------------------------------------------------------------------------
! Read any runtime specifications from the lvt.config file. 
!------------------------------------------------------------------------------

    write(LVT_logunit,*) '[INFO] Initializing GRDC streamflow data reader ...'

    call ESMF_ConfigGetAttribute(LVT_config, GRDCobs(i)%odir, &
         label='GRDC data directory:',rc=status)
    call LVT_verify(status, 'GRDC data directory: not defined')

    write(LVT_logunit,*) '[INFO] Reading from GRDC data directory ... '
    write(LVT_logunit,*) trim(GRDCobs(i)%odir)
  
    call ESMF_ConfigGetAttribute(LVT_config, stnlist_file, &
         label='GRDC station list file:',rc=status)
    call LVT_verify(status, 'GRDC station list file: not defined')

    call ESMF_ConfigGetAttribute(LVT_config, data_freq, &
         label='GRDC frequency of data:',rc=status)
    call LVT_verify(status, 'GRDC frequency of data: not defined')

    if(data_freq.eq."daily") then 
       GRDCobs(i)%use_daily = 1
    elseif(data_freq.eq."monthly") then 
       GRDCobs(i)%use_daily = 0
    endif

    call ESMF_ConfigGetAttribute(LVT_config, GRDCobs(i)%version, &
         label='GRDC file version:',rc=status)
    call LVT_verify(status, 'GRDC file version: not defined')
    write(LVT_logunit,*) '[INFO] Reading in GRDC file format version :: ',&
          GRDCobs(i)%version 
       

    ! Open station list file, and read number of station header:
    ftn = LVT_getNextUnitNumber()
    open(ftn, file=trim(stnlist_file), form='formatted')
    read(ftn,*)
    read(ftn,*) GRDCobs(i)%n_stns
    read(ftn,*) 

    allocate(GRDCobs(i)%stn_id(GRDCobs(i)%n_stns))
    allocate(GRDCobs(i)%stnlat(GRDCobs(i)%n_stns))
    allocate(GRDCobs(i)%stnlon(GRDCobs(i)%n_stns))

!------------------------------------------------------------------------------
! For each station, this reads the site name, station name, and station
! position in lat, lon coordinates.
!------------------------------------------------------------------------------
    do k=1,GRDCobs(i)%n_stns
       read(ftn,'(a)') currentLine
       iloc=index(currentLine,",")
       read(currentLine(1:iloc-1),*) GRDCobs(i)%stn_id(k)
       currentLine = currentLine(iloc+1:len(currentLine))
       
       iloc = Index(currentLine, ",")
       READ(currentLine(1: iloc - 1), *) GRDCobs(i)%stnlat(k)
       currentLine = currentLine(iloc + 1: Len(currentLine))

       iloc = Index(currentLine, ",")
       READ(currentLine(1: iloc - 1), *) GRDCobs(i)%stnlon(k)
       currentLine = currentLine(iloc + 1: Len(currentLine))
       
       write(LVT_logunit,*) '[INFO] ', trim(GRDCobs(i)%stn_id(k)), &
            GRDCobs(i)%stnlat(k), &
            GRDCobs(i)%stnlon(k)
    end do
    call LVT_releaseUnitNumber(ftn)

    ! Assign number of yearly timesteps and time step interval:
    if(GRDCobs(i)%use_daily.eq.1) then
       GRDCobs(i)%nts = 366
       call LVT_update_timestep(LVT_rc, 86400)
       call ESMF_TimeIntervalSet(GRDCobs(i)%timestep,s=86400,rc=status)
       call LVT_verify(status,"ESMF_TimeIntervalSet failed in GRDC_obsInit")

    else
       GRDCobs(i)%nts = 12
       call LVT_update_timestep(LVT_rc, 2592000)
    endif

    ! Allocate and initialize streamflow, "q", variable:

    allocate(GRDCobs(i)%q(GRDCobs(i)%n_stns,GRDCobs(i)%nts))
    GRDCobs(i)%q = -9999.0

    GRDCobs(i)%yr = -1
    GRDCobs(i)%mo = LVT_rc%mo

  end subroutine GRDC_obsinit


end module GRDC_obsMod
