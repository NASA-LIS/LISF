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
! !ROUTINE: HYMAP2_setSWOT
!  \label{HYMAP2_setSWOT}
!
! !REVISION HISTORY:
! 15 Apr 24: Yeosang Yoon; Initial specification;
!                          copied from HYMAP2_setWL
!
! !INTERFACE:
subroutine HYMAP2_setSWOT(n, Routing_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use HYMAP2_routingMod
  use HYMAP2_modelMod

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n
  type(ESMF_State)       :: Routing_State
!
! !DESCRIPTION:
!
!  This routine assigns the water level prognostic variables to the HyMAP2
!  model space.
!
!  The arguments are:
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[Routing\_State] ESMF State container for Routing state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)       :: sfcelvField
  integer                :: t,i,m
  integer                :: status
  real, pointer          :: sfcelv(:)

  real*8                 :: elevtn
  real*8                 :: fldhgt(HYMAP2_routing_struc(n)%nz)
  real*8                 :: rivelv
  real*8                 :: grarea
  real*8                 :: fldstomax(HYMAP2_routing_struc(n)%nz)
  real*8                 :: rivstomax
  real*8                 :: rivlen
  real*8                 :: rivwth
  real*8                 :: elv
  real*8                 :: vol

  logical                :: diffCheck(HYMAP2_routing_struc(n)%nseqall)
  logical                :: ensCheck(HYMAP2_routing_struc(n)%nseqall)

  external :: reorderEnsForOutliers

  call ESMF_StateGet(Routing_State,"Surface elevation",sfcelvField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sfcelvField in HYMAP2_getSWOT')

  call ESMF_FieldGet(sfcelvField,localDE=0,farrayPtr=sfcelv,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sfcelv in HYMAP2_getSWOT')

  ensCheck = .true.
  diffCheck = .false.

  do i=1,HYMAP2_routing_struc(n)%nseqall
     do m=1,LIS_rc%nensem(n)
        t = (i-1)*LIS_rc%nensem(n)+m
        sfcelv(t)  = sfcelv(t) + HYMAP2_routing_struc(n)%rivelv(i)
     enddo
  enddo

  do i=1,HYMAP2_routing_struc(n)%nseqall
     do m=1,LIS_rc%nensem(n)
        t = (i-1)*LIS_rc%nensem(n)+m
        if(sfcelv(t).le.HYMAP2_routing_struc(n)%rivelv(i)) then
           ensCheck(i) = .false.
        endif

        if(sfcelv(t).ne.sfcelv(i*LIS_rc%nensem(n))) then
           diffCheck(i) = .true.
        endif
     enddo
     if(.not.ensCheck(i).and.diffCheck(i)) then
        call reorderEnsForOutliers(&
             LIS_rc%nensem(n),&
             sfcelv((i-1)*LIS_rc%nensem(n)+1:i*LIS_rc%nensem(n)),&
             HYMAP2_routing_struc(n)%rivelv(i))
     endif
  enddo

  do i=1,HYMAP2_routing_struc(n)%nseqall
     do m=1,LIS_rc%nensem(n)
        t = (i-1)*LIS_rc%nensem(n)+m
        HYMAP2_routing_struc(n)%sfcelv(i,m) = sfcelv(t)
     enddo
  enddo

  !update surface level and then update the storage.
  do i=1,HYMAP2_routing_struc(n)%nseqall
     do m=1,LIS_rc%nensem(n)
        t = (i-1)*LIS_rc%nensem(n)+m

        elevtn=dble(HYMAP2_routing_struc(n)%elevtn(i))
        fldhgt=dble(HYMAP2_routing_struc(n)%fldhgt(i,:))
        fldstomax=dble(HYMAP2_routing_struc(n)%fldstomax(i,:,m))
        rivelv=dble(HYMAP2_routing_struc(n)%rivelv(i))
        grarea=dble(HYMAP2_routing_struc(n)%grarea(i))
        rivstomax=dble(HYMAP2_routing_struc(n)%rivstomax(i,m))
        rivlen=dble(HYMAP2_routing_struc(n)%rivlen(i))
        rivwth=dble(HYMAP2_routing_struc(n)%rivwth(i,m))
        elv=dble(HYMAP2_routing_struc(n)%sfcelv(i,m))

        call HYMAP2_get_volume_profile(&
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
        HYMAP2_routing_struc(n)%rivsto(i,m) = real(vol)
        HYMAP2_routing_struc(n)%fldsto(i,m) = 0.0
     enddo
  enddo

end subroutine HYMAP2_setSWOT
