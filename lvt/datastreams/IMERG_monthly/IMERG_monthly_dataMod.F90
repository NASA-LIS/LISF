!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! NOTE:  Currently only V07A IMERG Final Run data are supported
module IMERG_monthly_dataMod

  ! Defaults
  implicit none
  private

  ! Public routines
  public :: IMERG_monthly_datainit

  ! Public types
  public :: imergmonthlydata

  type, public :: imergmonthlydatadec
     character*100 :: odir
     character*10  :: imergver, imergprd
     real          :: datares
     real, allocatable           :: rlat(:)
     real, allocatable           :: rlon(:)

     ! This is only used with upscale averaging
     integer, allocatable        :: n11(:)

     ! These are only used with budget interpolation
     integer, allocatable        :: n112(:,:)
     integer, allocatable        :: n122(:,:)
     integer, allocatable        :: n212(:,:)
     integer, allocatable        :: n222(:,:)
     real,    allocatable        :: w112(:,:)
     real,    allocatable        :: w122(:,:)
     real,    allocatable        :: w212(:,:)
     real,    allocatable        :: w222(:,:)

     integer                     :: nc
     integer                     :: nr
  end type imergmonthlydatadec

  type(imergmonthlydatadec), allocatable :: imergmonthlydata(:)

contains

  subroutine IMERG_monthly_datainit(i)

    ! Imports
    use ESMF
    use LVT_coreMod, only: LVT_rc, LVT_config, LVT_isAtAFinerResolution
    use LVT_logMod, only: LVT_logunit, LVT_verify, LVT_endrun

    ! Defaults
    implicit none

    ! Arguments
    integer, intent(in) :: i

    ! Local variables
    integer              :: status
    real                 :: gridDesci(50)

    ! External routines
    external :: conserv_interp_input
    external :: upscaleByAveraging_input

    if (.not. allocated(imergmonthlydata)) then
       allocate(imergmonthlydata(LVT_rc%nDataStreams))
    endif

    ! Get top level IMERG data directory
    call ESMF_ConfigGetAttribute(LVT_Config, imergmonthlydata(i)%odir, &
         label='IMERG monthly data directory:', rc=status)
    call LVT_verify(status, 'IMERG monthly data directory: not defined')

    ! Get IMERG product
    call ESMF_ConfigFindLabel(LVT_config,"IMERG monthly product:", &
         rc=status)
    call LVT_verify(status, 'IMERG monthly product: not defined')
    call ESMF_ConfigGetAttribute(LVT_config, &
         imergmonthlydata(i)%imergprd, rc=status)
    if (imergmonthlydata(i)%imergprd .ne. 'final') then
       write(LVT_logunit,*)'[ERR] Invalid IMERG monthly product specified'
       write(LVT_logunit,*)"[ERR] Currently only 'final' is supported"
       call LVT_endrun()
    end if

    ! Get IMERG version
    call ESMF_ConfigFindLabel(LVT_config, "IMERG monthly version:", &
         rc=status)
    call LVT_verify(status, 'IMERG monthly version: not defined')
    call ESMF_ConfigGetAttribute(LVT_config, &
         imergmonthlydata(i)%imergver, rc=status)
    if (imergmonthlydata(i)%imergver .ne. 'V07A') then
       write(LVT_logunit,*)'[ERR] Invalid IMERG monthly version specified'
       write(LVT_logunit,*)"[ERR] Currently only 'V07A' is supported"
       call LVT_endrun()
    end if

    ! Allocate arrays on LVT grid
    allocate(imergmonthlydata(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(imergmonthlydata(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))

    ! Set IMERG grid and map projection information.
    gridDesci(:) =    0
    gridDesci(1) =    0     ! Lat/lon
    gridDesci(2) = 3600     ! Points along latitude circle
    gridDesci(3) = 1800     ! Points along longitude circle
    gridDesci(4) =  -89.95  ! Latitude of first grid point
    gridDesci(5) = -179.95  ! Longitude of first grid point
    gridDesci(6) =  128
    gridDesci(7) =   89.95  ! Latitude of last grid point
    gridDesci(8) =  179.95  ! Longitude of last grid point
    gridDesci(9) =    0.1   ! Longitudinal direction increment
    gridDesci(10) =   0.1   ! Latitude direction increment
    gridDesci(20) =  64

    ! Set up interpolation data
    imergmonthlydata(i)%datares = 0.1
    imergmonthlydata(i)%nc = 3600
    imergmonthlydata(i)%nr = 1800

    ! Use budget-bilinear interpolation if IMERG resolution (0.1 deg)
    ! is coarser than the analysis grid. Use upscale averaging if
    ! IMERG is finer than the analysis grid.
    if (LVT_isAtAFinerResolution(imergmonthlydata(i)%datares)) then

       ! Used only with budget interpolation
       allocate(imergmonthlydata(i)%n112(LVT_rc%lnc*LVT_rc%lnr,25))
       allocate(imergmonthlydata(i)%n122(LVT_rc%lnc*LVT_rc%lnr,25))
       allocate(imergmonthlydata(i)%n212(LVT_rc%lnc*LVT_rc%lnr,25))
       allocate(imergmonthlydata(i)%n222(LVT_rc%lnc*LVT_rc%lnr,25))
       allocate(imergmonthlydata(i)%w112(LVT_rc%lnc*LVT_rc%lnr,25))
       allocate(imergmonthlydata(i)%w122(LVT_rc%lnc*LVT_rc%lnr,25))
       allocate(imergmonthlydata(i)%w212(LVT_rc%lnc*LVT_rc%lnr,25))
       allocate(imergmonthlydata(i)%w222(LVT_rc%lnc*LVT_rc%lnr,25))
       imergmonthlydata(i)%n112 = 0
       imergmonthlydata(i)%n122 = 0
       imergmonthlydata(i)%n212 = 0
       imergmonthlydata(i)%n222 = 0
       imergmonthlydata(i)%w112 = 0
       imergmonthlydata(i)%w122 = 0
       imergmonthlydata(i)%w212 = 0
       imergmonthlydata(i)%w222 = 0
       call conserv_interp_input(gridDesci, LVT_rc%gridDesc, &
            LVT_rc%lnc*LVT_rc%lnr, &
            imergmonthlydata(i)%rlat, imergmonthlydata(i)%rlon, &
            imergmonthlydata(i)%n112, imergmonthlydata(i)%n122, &
            imergmonthlydata(i)%n212, imergmonthlydata(i)%n222, &
            imergmonthlydata(i)%w112, imergmonthlydata(i)%w122, &
            imergmonthlydata(i)%w212, imergmonthlydata(i)%w222)
    else

       ! Used only with upscale averaging
       allocate(imergmonthlydata(i)%n11(imergmonthlydata(i)%nc * &
            imergmonthlydata(i)%nr))
       imergmonthlydata(i)%n11 = 0
       call upscaleByAveraging_input(gridDesci, LVT_rc%gridDesc,&
            imergmonthlydata(i)%nc*imergmonthlydata(i)%nr, &
            LVT_rc%lnc*LVT_rc%lnr, &
            imergmonthlydata(i)%n11)
    end if
  end subroutine IMERG_monthly_datainit
end module IMERG_monthly_dataMod
