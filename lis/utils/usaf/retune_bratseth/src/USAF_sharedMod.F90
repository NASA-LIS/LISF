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
!
! MODULE: USAF_sharedMod
!
! This contains miscellaneous routines shared by procOBA_NWP and procOBA_Sat.
!
! REVISION HISTORY:
! 26 Oct 2020:  Eric Kemp.  Initial specification.
!
module USAF_sharedMod

  implicit none
  private

  public :: is_ssmi
  public :: is_cmorph
  public :: is_gage
  public :: is_geoprecip
  public :: is_imerg
  public :: endrun
  public :: pick_proc
  public :: fetch_blacklist

contains

  logical function is_ssmi(net)
    implicit none
    character(len=10) :: net
    logical :: answer
    answer = .false.
    if (trim(net) .eq. "SSMI") answer = .true.
    is_ssmi = answer
  end function is_ssmi

  logical function is_cmorph(net)
    implicit none
    character(len=10) :: net
    logical :: answer
    answer = .false.
    if (trim(net) .eq. "CMORPH") answer = .true.
    is_cmorph = answer
  end function is_cmorph

  logical function is_geoprecip(net)
    implicit none
    character(len=10) :: net
    logical :: answer
    answer = .false.
    if (trim(net) .eq. "GEOPRECIP") answer = .true.
    is_geoprecip = answer
  end function is_geoprecip

  logical function is_imerg(net)
    implicit none
    character(len=10) :: net
    logical :: answer
    answer = .false.
    if (trim(net) .eq. "IMERG") answer = .true.
    is_imerg = answer
  end function is_imerg

  logical function is_gage(net)
    implicit none
    character(len=10) :: net
    logical :: answer
    answer = .false.
    if (trim(net) .eq. "AMIL") answer = .true.
    if (trim(net) .eq. "CANA") answer = .true.
    if (trim(net) .eq. "FAA") answer = .true.
    if (trim(net) .eq. "ICAO") answer = .true.
    if (trim(net) .eq. "WMO") answer = .true.
    !Skip MOBL and SUPERGAGE. MOBL may move around with time, while
    !SUPERGAGE may include MOBL data.
    !if (trim(net) .eq. "MOBL") answer = .true.
    !if (trim(net) .eq. "SUPERGAGE") answer = .true.
    ! Handle reformatted CDMS data that are missing the network type.
    if (trim(net) .eq. "CDMS") answer = .true.
    is_gage = answer
  end function is_gage

  subroutine endrun(status)
    use esmf
    implicit none
    integer,intent(in) :: status
    call ESMF_LogFlush()
    if (status .ne. 0) then
       call esmf_finalize(endflag=ESMF_END_ABORT)
       stop 1
    else
       call esmf_finalize()
       stop
    endif
  end subroutine endrun

  subroutine pick_proc(id, id_incr, numprocs)
    implicit none
    integer,intent(inout) :: id
    integer,intent(inout) :: id_incr
    integer,intent(in) :: numprocs
    id = id + id_incr
    if (id .ge. numprocs) then
       id = numprocs - 1
       id_incr = -1
    else if (id .lt. 0) then
       id = 0
       id_incr = 1
    end if
  end subroutine pick_proc

  ! Construct list of stations to blacklist
  subroutine fetch_blacklist(blacklist_file, blacklist_stns, nstns)

    use esmf
    
    ! Defaults
    implicit none

    ! Arguments
    character(len=*), intent(in) :: blacklist_file
    character(len=9), allocatable, intent(inout) :: blacklist_stns(:)
    integer, intent(inout) :: nstns

    ! Local variables
    logical :: ex
    character(len=255) :: line
    integer :: istat
    integer :: i, j

    ! Check if file exists
    inquire(file=trim(blacklist_file), exist=ex)
    if (.not. ex) then
       call ESMF_LogWrite( &
            'Cannot find blacklist file '//trim(blacklist_file), &
            ESMF_LOGMSG_ERROR)
       call endrun(1)
    end if

    ! Open the file
    open(unit=15, &
         file=trim(blacklist_file), &
         status="OLD", &
         iostat=istat)

    ! Count the number of stations
    nstns = 0
    do
       read(unit=15, fmt='(A)', end=100) line
       if (line(1:1) .eq. "#") cycle
       nstns = nstns + 1
    end do
100 continue

    ! If no stations were found, we are done.
    if (nstns .eq. 0) then
       close(unit=15)
       return
    end if

    ! Rewind the file
    rewind(unit=15)

    ! Allocate space for the blacklisted stations
    allocate(blacklist_stns(nstns))
    blacklist_stns(:) = ''

    ! Loop through the file, adding the stations to the list
    i = 0
    do
       read(unit=15, fmt='(A)', end=200) line
       if (line(1:1) .eq. "#") cycle
       i = i + 1
       do j = 1, 9
          if (line(j:j) .eq. ' ') exit
       end do
       blacklist_stns(i) = trim(line(1:j))
    end do
200 continue

    ! Finish up
    close(unit=10)

  end subroutine fetch_blacklist

end module USAF_sharedMod
