!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !MODULE: GLASSlaiobsMod
! \label(GLASSlaiobsMod)
!
! !INTERFACE:
module GLASSlaiobsMod
! 
! !USES:   
  use ESMF

  implicit none

  PRIVATE 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  7 Mar 2015   Sujay Kumar  Initial Specification
! 
!EOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GLASSlaiobsinit !Initializes structures for reading MOD16A2 data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GLASSlaiobs !Object to hold GLASSlai observation attributes
!EOP

  type, public :: glasslaidec
     character*100           :: odir
     character*100           :: source
     integer                 :: nc, nr
     real                    :: gridDesc(50)
     logical                 :: startFlag
     real                    :: datares

     real,    allocatable :: rlat(:)
     real,    allocatable :: rlon(:)

     integer, allocatable :: n11(:)
     integer, allocatable :: n12(:)
     integer, allocatable :: n21(:)
     integer, allocatable :: n22(:)
     real,    allocatable :: w11(:)
     real,    allocatable :: w12(:)
     real,    allocatable :: w21(:)
     real,    allocatable :: w22(:)

  end type glasslaidec
     
  type(glasslaidec), save :: GLASSlaiObs(2)

contains
  
!BOP
! 
! !ROUTINE: GLASSlaiobsInit
! \label{GLASSlaiobsInit}
!
! !INTERFACE: 
  subroutine GLASSlaiobsinit(i)
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
!   This subroutine initializes and sets up the data structures required
!   for reading the GIMMSAVHRR NDVI data, including the computation of spatial 
!   interpolation weights. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: status

    call ESMF_ConfigGetAttribute(LVT_Config, GLASSlaiobs(i)%odir, &
         label='GLASS LAI data directory:',rc=status)
    call LVT_verify(status, 'GLASS LAI data directory: not defined')

! source = "AVHRR" or "MODIS"
    call ESMF_ConfigGetAttribute(LVT_Config, GLASSlaiobs(i)%source, &
         label='GLASS LAI data source:',rc=status)
    call LVT_verify(status, 'GLASS LAI data source: not defined')

    call LVT_update_timestep(LVT_rc, 86400)

    glasslaiobs(i)%gridDesc = 0
        
    glasslaiobs(i)%nc = 7200
    glasslaiobs(i)%nr = 3600

    !filling the items needed by the interpolation library
    glasslaiobs(i)%gridDesc(1) = 0  
    glasslaiobs(i)%gridDesc(2) = glasslaiobs(i)%nc
    glasslaiobs(i)%gridDesc(3) = glasslaiobs(i)%nr
    glasslaiobs(i)%gridDesc(4) = -89.975
    glasslaiobs(i)%gridDesc(5) = -179.975
    glasslaiobs(i)%gridDesc(7) = 89.975
    glasslaiobs(i)%gridDesc(8) = 179.975
    glasslaiobs(i)%gridDesc(6) = 128
    glasslaiobs(i)%gridDesc(9) = 0.05
    glasslaiobs(i)%gridDesc(10) = 0.05
    glasslaiobs(i)%gridDesc(20) = 64

    glasslaiobs(i)%datares  = 0.05

    if(LVT_isAtAfinerResolution(glasslaiobs(i)%datares)) then
       
       allocate(glasslaiobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
       allocate(glasslaiobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
       allocate(glasslaiobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
       allocate(glasslaiobs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
       allocate(glasslaiobs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
       allocate(glasslaiobs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))
       allocate(glasslaiobs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
       allocate(glasslaiobs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
       allocate(glasslaiobs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
       allocate(glasslaiobs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))
       
       call bilinear_interp_input(glasslaiobs(i)%gridDesc, &
            LVT_rc%gridDesc, &
            LVT_rc%lnc*LVT_rc%lnr, &
            glasslaiobs(i)%rlat, &
            glasslaiobs(i)%rlon, &
            glasslaiobs(i)%n11,&
            glasslaiobs(i)%n12,&
            glasslaiobs(i)%n21,&
            glasslaiobs(i)%n22,&
            glasslaiobs(i)%w11,&
            glasslaiobs(i)%w12,&
            glasslaiobs(i)%w21,&
            glasslaiobs(i)%w22)
    else
       allocate(glasslaiobs(i)%n11(glasslaiobs(i)%nc*glasslaiobs(i)%nr))
       call upscaleByAveraging_input(glasslaiobs(i)%gridDesc,&
            LVT_rc%gridDesc,glasslaiobs(i)%nc*glasslaiobs(i)%nr,&
            LVT_rc%lnc*LVT_rc%lnr,glasslaiobs(i)%n11)
    endif

    glasslaiobs(i)%startFlag = .false. 

  end subroutine GLASSlaiobsinit


end module GLASSlaiobsMod
