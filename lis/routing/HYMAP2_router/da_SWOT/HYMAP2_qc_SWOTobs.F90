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
! !ROUTINE: HYMAP2_qc_SWOTobs
! \label{HYMAP2_qc_SWOTobs}
!
! !REVISION HISTORY:
! 15 Apr 24: Yeosang Yoon; Initial specification;
!                          copied from HYMAP2_qc_WLobs
!
! !INTERFACE:
subroutine HYMAP2_qc_SWOTobs(n,k,OBS_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod,  only : LIS_verify
  use LIS_constantsMod, only : LIS_CONST_TKFRZ
  use LIS_DAobservationsMod
  use HYMAP2_routingMod

  implicit none

! !ARGUMENTS:
  integer, intent(in)      :: n
  integer, intent(in)      :: k
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION:
!
!
!  The arguments are:
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF state container for observations \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)          :: obs_wl_field
  real,       pointer       :: wlobs(:)
  integer                   :: i,t
  real                      :: wl(HYMAP2_routing_struc(n)%nseqall)
  real                      :: wl_obsspace(LIS_rc%obs_ngrid(k))
  integer                   :: status

  external :: HYMAP2_convertRoutingSpaceToObsSpace

  call ESMF_StateGet(OBS_State,"Observation01",obs_wl_field,rc=status)
  call LIS_verify(status,"ESMF_StateGet failed in HYMAP2_qc_SWOTobs")
  call ESMF_FieldGet(obs_wl_field,localDE=0,farrayPtr=wlobs,rc=status)
  call LIS_verify(status,"ESMF_FieldGet failed in HYMAP2_qc_SWOTobs")

  do i=1,HYMAP2_routing_struc(n)%nseqall
     wl(i) = HYMAP2_routing_struc(n)%rivelv(i)
  enddo

  call HYMAP2_convertRoutingSpaceToObsSpace(n, k, &
       wl, wl_obsspace)

  do t = 1,LIS_rc%obs_ngrid(k)
     if(wlobs(t).le.wl_obsspace(t)) then
        wlobs(t) = LIS_rc%udef
     endif
  enddo

end subroutine HYMAP2_qc_SWOTobs
