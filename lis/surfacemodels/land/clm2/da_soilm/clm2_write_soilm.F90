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
! !ROUTINE: clm2_write_soilm
! \label{clm2_write_soilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine clm2_write_soilm(ftn, n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod,    only : LIS_rc
  use LIS_historyMod, only : LIS_writevar_restart
  use clm2_lsmMod,     only : clm2_struc
  use clm2_varcon,     only : denh2o, denice
  use clm2_varpar,     only : nlevsoi

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: ftn
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the soilmoisture related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  real,allocatable   :: soilm(:)
  real,allocatable   :: eff_porosity(:)
  real,allocatable   :: vol_ice(:), vol_liq(:)
  integer :: m, t
 
  allocate(soilm(LIS_rc%ntiles(n)))
  allocate(eff_porosity(LIS_rc%ntiles(n)))
  allocate(vol_ice(LIS_rc%ntiles(n)))
  allocate(vol_liq(LIS_rc%ntiles(n)))

  do m=1,nlevsoi 
     do t=1,LIS_rc%ntiles(n)
        vol_ice(m) = min(clm2_struc(n)%clm(t)%watsat(m), &
             clm2_struc(n)%clm(t)%h2osoi_ice(m)/&
             (clm2_struc(n)%clm(t)%dz(m)*denice))
        eff_porosity(m) = clm2_struc(n)%clm(t)%watsat(m)-vol_ice(m)
        vol_liq(m) = min(eff_porosity(m), &
             clm2_struc(n)%clm(t)%h2osoi_liq(m)/&
             (clm2_struc(n)%clm(t)%dz(m)*denh2o))
        soilm(t)= (vol_liq(m)+vol_ice(m))
     enddo
     call LIS_writevar_restart(ftn,n,soilm)
  enddo

  deallocate(soilm)
  deallocate(eff_porosity)
  deallocate(vol_ice)
  deallocate(vol_liq)

end subroutine clm2_write_soilm

