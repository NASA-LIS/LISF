! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

MODULE flake_configure

!------------------------------------------------------------------------------
!
! Description:
!
!  Switches and reference values of parameters 
!  that configure the lake model FLake are set.
!
!
! Current Code Owner: DWD, Dmitrii Mironov
!  Phone:  +49-69-8062 2705
!  Fax:    +49-69-8062 3721
!  E-mail: dmitrii.mironov@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.00       2005/11/17 Dmitrii Mironov 
!  Initial release 
! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

USE data_parameters , ONLY : &
  ireals                   , & ! KIND-type parameter for real variables 
  iintegers                    ! KIND-type parameter for "normal" integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations

LOGICAL, PARAMETER :: &
  lflk_botsed_use   = .TRUE.         ! .TRUE. indicates that the bottom-sediment scheme is used
                                     ! to compute the depth penetrated by the thermal wave, 
                                     ! the temperature at this depth and the bottom heat flux.
                                     ! Otherwise, the heat flux at the water-bottom sediment interface
                                     ! is set to zero, the depth penetrated by the thermal wave 
                                     ! is set to a reference value defined below,
                                     ! and the temperature at this depth is set to 
                                     ! the temperature of maximum density of the fresh water.

REAL (KIND = ireals), PARAMETER :: &
  rflk_depth_bs_ref = 10.0_ireals    ! Reference value of the depth of the thermally active
                                     ! layer of bottom sediments [m].
                                     ! This value is used to (formally) define
                                     ! the depth penetrated by the thermal wave
                                     ! in case the bottom-sediment scheme is not used.

!==============================================================================

END MODULE flake_configure

