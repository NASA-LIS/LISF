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
! !ROUTINE: HYMAP2_getWL
! \label{HYMAP2_getWL}
!
! !REVISION HISTORY:
!  07 Nov 2019: Sujay Kumar, Initial specification
!
! !INTERFACE:
subroutine HYMAP2_getWL(n, Routing_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use HYMAP2_routingMod
  use HYMAP2_modelMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: Routing_State
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
  character*100          :: lsm_state_objs(4)

#if 0 
  integer                :: ielevtn
  integer                :: ifldhgt(HYMAP2_routing_struc(n)%nz)
  real                   :: fldstomax(HYMAP2_routing_struc(n)%nz)
  real                   :: grarea
  real                   :: rivstomax
  real                   :: rivlen
  real                   :: rivwth
  integer                :: ielv
  integer                :: irivelv
  real*8                 :: vol

  real*8                 :: grarea1
  real*8                 :: fldstomax1(HYMAP2_routing_struc(n)%nz)
  real*8                 :: rivstomax1
  real*8                 :: rivlen1
  real*8                 :: rivwth1
#endif
  real*8                :: elevtn
  real*8                :: fldhgt(HYMAP2_routing_struc(n)%nz)
  real*8                   :: fldstomax(HYMAP2_routing_struc(n)%nz)
  real*8                   :: grarea
  real*8                   :: rivstomax
  real*8                   :: rivlen
  real*8                   :: rivwth
  real*8                   :: elv
  real*8                 :: rivelv
  real*8                 :: vol

  real*8                 :: grarea1
  real*8                 :: fldstomax1(HYMAP2_routing_struc(n)%nz)
  real*8                 :: rivstomax1
  real*8                 :: rivlen1
  real*8                 :: rivwth1


  call ESMF_StateGet(Routing_State,"Surface elevation",sfcelevField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sfcelev in HYMAP2_getWL')

  call ESMF_FieldGet(sfcelevField,localDE=0,farrayPtr=sfcelev,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sfcelev in HYMAP2_getWL')

  do i=1,HYMAP2_routing_struc(n)%nseqall
     do m=1,LIS_rc%nensem(n)
        t = (i-1)*LIS_rc%nensem(n)+m
        elevtn=dble(HYMAP2_routing_struc(n)%elevtn(i))
        fldhgt=dble(HYMAP2_routing_struc(n)%fldhgt(i,:))
        fldstomax = dble(HYMAP2_routing_struc(n)%fldstomax(i,:,m))
        rivelv=dble(HYMAP2_routing_struc(n)%rivelv(i))
        grarea = dble(HYMAP2_routing_struc(n)%grarea(i))
        rivstomax = dble(HYMAP2_routing_struc(n)%rivstomax(i,m))
        rivelv=dble(HYMAP2_routing_struc(n)%rivelv(i))
        rivlen = dble(HYMAP2_routing_struc(n)%rivlen(i))
        rivwth = dble(HYMAP2_routing_struc(n)%rivwth(i,m))
        vol = dble(HYMAP2_routing_struc(n)%rivsto(i,m)+&
             HYMAP2_routing_struc(n)%fldsto(i,m))

        call HYMAP2_get_elevation_profile(&
             HYMAP2_routing_struc(n)%nz,&
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
        sfcelev(t) = real(elv) - HYMAP2_routing_struc(n)%rivelv(i)
     enddo
     
  enddo

end subroutine HYMAP2_getWL
