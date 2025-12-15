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
! !ROUTINE: HYMAP3_getWL
! \label{HYMAP3_getWL}
!
! !REVISION HISTORY:
!  07 Nov 2019: Sujay Kumar, Initial specification
!
! !INTERFACE:
subroutine HYMAP3_getWL(n, Routing_State)

! !USES:
  use ESMF
  use HYMAP3_modelMod
  use HYMAP3_routingMod
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify

  implicit none

! !ARGUMENTS:
  integer, intent(in)          :: n
  type(ESMF_State), intent(in) :: Routing_State
!
! !DESCRIPTION:
!
!  Returns the water level DA related state prognostic variables for
!  data assimilation
!
!  The arguments are:
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[Routing\_State] ESMF State container for Routing state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)       :: sfcelevField
  integer                :: t,i,m
  integer                :: status
  real, pointer          :: sfcelev(:)
  real*8                :: elevtn
  real*8                :: fldhgt(HYMAP3_routing_struc(n)%nz)
  real*8                   :: fldstomax(HYMAP3_routing_struc(n)%nz)
  real*8                   :: grarea
  real*8                   :: rivstomax
  real*8                   :: rivlen
  real*8                   :: rivwth
  real*8                   :: elv
  real*8                 :: rivelv
  real*8                 :: vol

  call ESMF_StateGet(Routing_State,"Surface elevation",sfcelevField, &
       rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm1 in HYMAP3_getWL')

  call ESMF_FieldGet(sfcelevField,localDE=0,farrayPtr=sfcelev,rc=status)
  call LIS_verify(status, &
       'ESMF_FieldGet failed for sfcelev in HYMAP3_getWL')

  do i=1,HYMAP3_routing_struc(n)%nseqall
     do m=1,LIS_rc%nensem(n)
        t = (i-1)*LIS_rc%nensem(n)+m
        elevtn=dble(HYMAP3_routing_struc(n)%elevtn(i))
        fldhgt=dble(HYMAP3_routing_struc(n)%fldhgt(i,:))
        fldstomax = dble(HYMAP3_routing_struc(n)%fldstomax(i,:,m))
        rivelv=dble(HYMAP3_routing_struc(n)%rivelv(i))
        grarea = dble(HYMAP3_routing_struc(n)%grarea(i))
        rivstomax = dble(HYMAP3_routing_struc(n)%rivstomax(i,m))
        rivelv=dble(HYMAP3_routing_struc(n)%rivelv(i))
        rivlen = dble(HYMAP3_routing_struc(n)%rivlen(i))
        rivwth = dble(HYMAP3_routing_struc(n)%rivwth(i,m))
        vol = dble(HYMAP3_routing_struc(n)%rivsto(i,m)+&
             HYMAP3_routing_struc(n)%fldsto(i,m))

        call HYMAP3_get_elevation_profile(&
             HYMAP3_routing_struc(n)%nz,&
             elevtn,&
             fldhgt,&
             fldstomax,&
             grarea,&
             rivstomax,&
             rivelv,&
             rivlen,&
             rivwth,&
             elv,&
             vol)
        sfcelev(t) = real(elv) - HYMAP3_routing_struc(n)%rivelv(i)
     enddo

  enddo

end subroutine HYMAP3_getWL

