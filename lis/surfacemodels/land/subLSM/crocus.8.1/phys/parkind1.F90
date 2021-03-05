MODULE PARKIND1
!
!     *** Define usual kinds for strong typing ***
!
IMPLICIT NONE
SAVE
!
!     Integer Kinds
!     -------------
!
INTEGER, PARAMETER :: JPIT = SELECTED_INT_KIND(2)
INTEGER, PARAMETER :: JPIS = SELECTED_INT_KIND(4)
INTEGER, PARAMETER :: JPIM = SELECTED_INT_KIND(9)
INTEGER, PARAMETER :: JPIB = SELECTED_INT_KIND(12)

!Special integer type to be used for sensative adress calculations
!should be *8 for a machine with 8byte adressing for optimum performance
#ifdef ADDRESS64
INTEGER, PARAMETER :: JPIA = JPIB
#else
INTEGER, PARAMETER :: JPIA = JPIM
#endif

!
!     Real Kinds
!     ----------
!
INTEGER, PARAMETER :: JPRT = SELECTED_REAL_KIND(2,1)
INTEGER, PARAMETER :: JPRS = SELECTED_REAL_KIND(4,2)
INTEGER, PARAMETER :: JPRM = SELECTED_REAL_KIND(6,37)
#ifdef PARKIND1_SINGLE
INTEGER, PARAMETER :: JPRB = SELECTED_REAL_KIND(6,37)
#else
INTEGER, PARAMETER :: JPRB = SELECTED_REAL_KIND(13,300)
#endif

! Double real for C code and special places requiring 
!    higher precision. 
INTEGER, PARAMETER :: JPRD = SELECTED_REAL_KIND(13,300)


! Logical Kinds for RTTOV....

INTEGER, PARAMETER :: JPLM = JPIM   !Standard logical type

END MODULE PARKIND1
