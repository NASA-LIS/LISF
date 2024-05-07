!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
!BOP
! 
! !MODULE: Ameriflux_obsMod
!  \label(Ameriflux_obsMod)
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
!  This is for defining new observational plugins into
!  LVT. The rountine checks to see if the current year is a leap year, since
!  this would affect the amount of space needed in arrays. It then extracts
!  station names, id's, and locations from a file. Finally, this initializes
!  the obs entry log so that it can be used by the next rountine.
!
!  This file uses the lvt.config.ctrl, located in the same folder as the LVT
!  executable. The config controls such data as start year/month/day/...second,
!  the location of output files, and the location of files the program uses
!  (such as the stations file, and the ts_locations.txt file)
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  18 Apr 2009   Sujay Kumar  Initial Specification
!  30 Jun 2010   Teodor Georgiev Adapted for Ameriflux data
! 
!EOP
! 
! 
!
module Ameriflux_obsMod

  use ESMF

  implicit none

  PRIVATE 

  PUBLIC :: Ameriflux_obsinit
  PUBLIC :: Amerifluxobs

  type, public :: Amerifluxobsdec
     character*100 :: odir
     character*50            :: version
     integer                 :: n_stns
     character*100, allocatable  :: site_name(:)
     character*100, allocatable  :: stn_name(:)
     real,          allocatable  :: stnlat(:)
     real,          allocatable  :: stnlon(:)
     real,          allocatable  :: sfsm_wt(:,:)
     real,          allocatable  :: rzsm_wt(:,:)
     real,          allocatable  :: sfst_wt(:,:)
     real,          allocatable  :: rzst_wt(:,:)
     real                    :: change
     logical                 :: startflag
     integer                 :: yr
     integer                 :: nsmlayers
     integer                 :: nstlayers

     type(ESMF_Time)         :: starttime
     type(ESMF_TimeInterval) :: timestep

     real,          allocatable  :: Qle(:,:)        ! Latent heat flux
     real,          allocatable  :: Qh(:,:)         ! Sensible heat flux
     real,          allocatable  :: Qg(:,:)         ! Ground Heat Flux
     real,          allocatable  :: Ta(:,:)
     real,          allocatable  :: sfsm(:,:)
     real,          allocatable  :: rzsm(:,:)
     real,          allocatable  :: sfst(:,:)
     real,          allocatable  :: rzst(:,:)
     real,          allocatable  :: precip(:,:)
  end type Amerifluxobsdec

  type(Amerifluxobsdec), allocatable :: Amerifluxobs(:)

contains
  
!BOP
! 
! !ROUTINE: Ameriflux_obsInit
! \label{Ameriflux_obsInit}
!
! !INTERFACE: 
 subroutine Ameriflux_obsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_logMod
    use LVT_timeMgrMod

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine should contain code to initialize and set up the data 
!  structures required for reading the specific data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! 
! !ARGUMENTS: 
    integer,   intent(IN) :: i 
! 
!
!EOP
    character*100 :: stnlist_file
    integer       :: ftn 
    integer       :: k,iloc
    integer       :: arrayLen
    character*200 :: currentLine
    real          :: swcd(2), tsd(2)
    integer       :: status

    if(.not.allocated(Amerifluxobs)) then 
       allocate(Amerifluxobs(LVT_rc%nDataStreams))
    endif
!------------------------------------------------------------------------------
! Read any runtime specifications from the lvt.config file. 
!------------------------------------------------------------------------------

    call ESMF_ConfigGetAttribute(LVT_config, Amerifluxobs(i)%odir, &
         label='Ameriflux observation directory:',rc=status)
    call LVT_verify(status, 'Ameriflux observation directory: not defined')
  
    call ESMF_ConfigGetAttribute(LVT_config, Amerifluxobs(i)%version, &
         label='Ameriflux data level:',rc=status)
    call LVT_verify(status, 'Ameriflux data level: not defined')

    call ESMF_ConfigGetAttribute(LVT_config, stnlist_file, &
         label='Ameriflux station list file:',rc=status)
    call LVT_verify(status, 'Ameriflux station list file: not defined')

!------------------------------------------------------------------------------
! Initialize any other variables.
! This is done using the Ameriflux_stns.txt, created specifically for this 
! purpose. The file is found in AMERIFLUX_DATA. It follows the straighforward
! format of:
! nstations, followed by a new line containing the number of stations.
! next is a line describing that stations are listed by name, location, and
! coordinates. On the fourth line, the station list begins. For example, the
! header of the file when it was written looked like this:
! #nstns
! 76
! #stnname; location name; lat; lon
! ARM_SGP_Burn; USARb; 35.5497; -98.0402
! Since there were 76 stations, 75 follwed the one listed
!------------------------------------------------------------------------------
    ftn = LVT_getNextUnitNumber()
    open(ftn, file=trim(stnlist_file), form='formatted')
    read(ftn,*)
    read(ftn,*) Amerifluxobs(i)%n_stns
    read(ftn,*) 

    allocate(Amerifluxobs(i)%site_name(Amerifluxobs(i)%n_stns))
    allocate(Amerifluxobs(i)%stn_name(Amerifluxobs(i)%n_stns))
    allocate(Amerifluxobs(i)%stnlat(Amerifluxobs(i)%n_stns))
    allocate(Amerifluxobs(i)%stnlon(Amerifluxobs(i)%n_stns))
    allocate(Amerifluxobs(i)%sfsm_wt(Amerifluxobs(i)%n_stns,2))
    allocate(Amerifluxobs(i)%rzsm_wt(Amerifluxobs(i)%n_stns,2))
    allocate(Amerifluxobs(i)%sfst_wt(Amerifluxobs(i)%n_stns,2))
    allocate(Amerifluxobs(i)%rzst_wt(Amerifluxobs(i)%n_stns,2))

!---------------------------------------------------------------------------
! Leap year test, following rules by Microsoft. The depth of this test may not
! be entirely necessary due to how recent the data is, but it is correct.
!---------------------------------------------------------------------------   
    if (Mod(LVT_rc%yr, 4) == 0) Then            !Step1
        if (Mod(LVT_rc%yr, 100) == 0) Then      !Step2
           if (Mod(LVT_rc%yr, 400) == 0) Then   !Step3
              arrayLen = 17568                  !Step4, leap LVT_rc%yr
           else                           
              arrayLen = 17520                  !Step5, not leap LVT_rc%yr
           end if
        else
           arrayLen = 17568                     !Step4, leap LVT_rc%yr
        end if
     else
        arrayLen = 17520                        !Step5, not leap LVT_rc%yr
     End if
!svk : to be safe:
     arrayLen = 17570

    allocate(Amerifluxobs(i)%Qle(Amerifluxobs(i)%n_stns, arrayLen))
    allocate(Amerifluxobs(i)%Qh(Amerifluxobs(i)%n_stns, arrayLen))
    allocate(Amerifluxobs(i)%Qg(Amerifluxobs(i)%n_stns, arrayLen))
    allocate(Amerifluxobs(i)%Ta(Amerifluxobs(i)%n_stns, arrayLen))
    allocate(Amerifluxobs(i)%sfsm(Amerifluxobs(i)%n_stns, arrayLen))
    allocate(Amerifluxobs(i)%sfst(Amerifluxobs(i)%n_stns, arrayLen))
    allocate(Amerifluxobs(i)%precip(Amerifluxobs(i)%n_stns, arrayLen))

    Amerifluxobs(i)%Qle = LVT_rc%udef
    Amerifluxobs(i)%Qh = LVT_rc%udef
    Amerifluxobs(i)%Qg = LVT_rc%udef
    Amerifluxobs(i)%Ta = LVT_rc%udef
    Amerifluxobs(i)%sfsm = LVT_rc%udef
    Amerifluxobs(i)%sfst = LVT_rc%udef
    Amerifluxobs(i)%precip = LVT_rc%udef

    write(LVT_logunit,*) '[INFO] Processing Ameriflux stations'
!------------------------------------------------------------------------------
! For each station, this reads the site name, station name, and station
! position in lat, lon coordinates.
!------------------------------------------------------------------------------
    do k=1,Amerifluxobs(i)%n_stns
       read(ftn,'(a)') currentLine
       iloc = Index(currentLine, ";")
       READ(currentLine(1: iloc - 1), *) Amerifluxobs(i)%site_name(k)
       currentLine = currentLine(iloc + 1: Len(currentLine))

       iloc = Index(currentLine, ";")
       READ(currentLine(1: iloc - 1), *) Amerifluxobs(i)%stn_name(k)
       currentLine = currentLine(iloc + 1: Len(currentLine))

       iloc = Index(currentLine, ";")
       READ(currentLine(1: iloc - 1), *) Amerifluxobs(i)%stnlat(k)
       currentLine = currentLine(iloc + 1: Len(currentLine))

       iloc = Index(currentLine, ";")
       READ(currentLine(1: iloc - 1), *) Amerifluxobs(i)%stnlon(k)
       currentLine = currentLine(iloc + 1: Len(currentLine))

       iloc = Index(currentLine, ";")
       READ(currentLine(1: iloc - 1), *) SWCD(1)
       currentLine = currentLine(iloc + 1: Len(currentLine))

       iloc = Index(currentLine, ";")
       READ(currentLine(1: iloc - 1), *) SWCD(2)
       currentLine = currentLine(iloc + 1: Len(currentLine))

       iloc = Index(currentLine, ";") 
       READ(currentLine(1: iloc - 1), *) tsd(1)
       currentLine = currentLine(iloc + 1: Len(currentLine))

       READ(currentLine, *) tsd(2)
       
       if(swcd(1).gt.0.and.swcd(2).gt.0) then 
          Amerifluxobs(i)%nsmlayers = 2
       elseif(swcd(1).gt.0.and.swcd(2).le.0) then 
          Amerifluxobs(i)%nsmlayers = 1
       elseif(swcd(1).le.0.and.swcd(2).gt.0) then 
          Amerifluxobs(i)%nsmlayers = 1
       else
          Amerifluxobs(i)%nsmlayers = 0
       endif

       if(tsd(1).gt.0.and.tsd(2).gt.0) then 
          Amerifluxobs(i)%nstlayers = 2
       elseif(swcd(1).gt.0.and.swcd(2).le.0) then 
          Amerifluxobs(i)%nstlayers = 1
       elseif(swcd(1).le.0.and.swcd(2).gt.0) then 
          Amerifluxobs(i)%nstlayers = 1
       else
          Amerifluxobs(i)%nstlayers = 0
       endif

       !print*, Amerifluxobs(i)%site_name(k), ' ', Amerifluxobs(i)%stn_name(k), ' ', Amerifluxobs(i)%stnlat(k), ' ', &
!& Amerifluxobs(i)%stnlon(k), ' ', Amerifluxobs(i)%SWC1D(k), ' ', Amerifluxobs(i)%SWC2D(k)
       
       
       if(Amerifluxobs(i)%nsmlayers.ne.0) then 
          call compute_vinterp_weights(&
               Amerifluxobs(i)%nsmlayers,&
               LVT_rc%lis_sf_d, &
               LVT_rc%lis_rz_d, &
               swcd(1:Amerifluxobs(i)%nsmlayers)/100.0, &
               Amerifluxobs(i)%sfsm_wt(k,:), &
               Amerifluxobs(i)%rzsm_wt(k,:))
       else
             Amerifluxobs(i)%sfsm_wt(k,:) = 0.0
             Amerifluxobs(i)%rzsm_wt(k,:) = 0.0
       endif

       if(Amerifluxobs(i)%nstlayers.ne.0) then 
          call compute_vinterp_weights(&
               Amerifluxobs(i)%nstlayers,&
               LVT_rc%lis_sf_d, &
               LVT_rc%lis_rz_d, &
               tsd(1:Amerifluxobs(i)%nstlayers)/100.0, &
               Amerifluxobs(i)%sfst_wt(k,:), &
               Amerifluxobs(i)%rzst_wt(k,:))
       else
             Amerifluxobs(i)%sfst_wt(k,:) = 0.0
             Amerifluxobs(i)%rzst_wt(k,:) = 0.0
       endif
       write(LVT_logunit,*) '[INFO] ',k, Amerifluxobs(i)%site_name(k), &
            Amerifluxobs(i)%stn_name(k), &
            Amerifluxobs(i)%stnlat(k), Amerifluxobs(i)%stnlon(k)
    end do
    
    call LVT_releaseUnitNumber(ftn)
!------------------------------------------------------------------------------
! Sets the start flag to true since this is the start of data. The flag is then
! set to false after the initial file is read so that a switch in year will not
! be mistaken for a new data set.
!------------------------------------------------------------------------------

    Amerifluxobs(i)%startflag = .true. 
    Amerifluxobs(i)%yr = -1

    call ESMF_TimeIntervalSet(Amerifluxobs(i)%timestep, s=1800, rc=status)
    call LVT_verify(status, 'error in setting timestep (Amerifluxobs)')

    call LVT_update_timestep(LVT_rc, 1800)

  end subroutine Ameriflux_obsinit

end module Ameriflux_obsMod
