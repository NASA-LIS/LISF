!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: HYMAP2_qcWL
! \label{HYMAP2_qcWL}
!
! !REVISION HISTORY:
!  07 Nov 2019: Sujay Kumar, Initial specification
!
! !INTERFACE:
subroutine HYMAP2_qcWL(n, Routing_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use HYMAP2_routingMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: Routing_State
!
! !DESCRIPTION:
!
!  QCs the water level related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[Routing\_State] ESMF State container for Routing state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)       :: sfcelvField
  integer                :: i,t,m
  integer                :: status
  real, pointer          :: sfcelv(:)
  real                   :: sfcelvmax
  real                   :: sfcelvmin


  call ESMF_StateGet(Routing_State,"Surface elevation",sfcelvField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet for Surface elevation failed in HYMAP2_qcWL")
 
  call ESMF_FieldGet(sfcelvField,localDE=0,farrayPtr=sfcelv,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet for Surface elevation failed in HYMAP2_qcWL")

  do i=1,HYMAP2_routing_struc(n)%nseqall
     do m=1,LIS_rc%nensem(n)
        t = (i-1)*LIS_rc%nensem(n)+m
! make sure that the surface elevations are never below the river bed
! elevatiion
        if(sfcelv(t).lt.0) then 
!           print*, 'qc issue ',i,m, sfcelv(t)
        endif
     enddo
  enddo


end subroutine HYMAP2_qcWL

