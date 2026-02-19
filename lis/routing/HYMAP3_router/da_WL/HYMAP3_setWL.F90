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
! !ROUTINE: HYMAP3_setWL
!  \label{HYMAP3_setWL}
!
! !REVISION HISTORY:
! 07 Nov 2019: Sujay Kumar, Initial specification
!
! !INTERFACE:
subroutine HYMAP3_setWL(n, Routing_State)

! !USES:
  use ESMF
  use HYMAP3_modelMod
  use HYMAP3_routingMod
  use LIS_coreMod
  use LIS_logMod

  implicit none

! !ARGUMENTS:
  integer, intent(in)                :: n
  type(ESMF_State), intent(in)       :: Routing_State
!
! !DESCRIPTION:
!
!  This routine assigns the water level prognostic variables to the HyMAP3
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
  real*8                 :: fldhgt(HYMAP3_routing_struc(n)%nz)
  real*8                 :: rivelv
  real*8                 :: grarea
  real*8                 :: fldstomax(HYMAP3_routing_struc(n)%nz)
  real*8                 :: rivstomax
  real*8                 :: rivlen
  real*8                 :: rivwth
  real*8                 :: elv
  real*8                 :: vol

  logical                :: diffCheck(HYMAP3_routing_struc(n)%nseqall)
  logical                :: ensCheck(HYMAP3_routing_struc(n)%nseqall)

  external :: HYMAP3_reorderEnsForOutliers

  call ESMF_StateGet(Routing_State,"Surface elevation",sfcelvField, &
       rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm1 in HYMAP3_getWL')

  call ESMF_FieldGet(sfcelvField,localDE=0,farrayPtr=sfcelv,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sfcelv in HYMAP3_getWL')

  ensCheck = .true.
  diffCheck = .false.

  do i=1,HYMAP3_routing_struc(n)%nseqall
     do m=1,LIS_rc%nensem(n)
        t = (i-1)*LIS_rc%nensem(n)+m
        sfcelv(t)  = sfcelv(t) + HYMAP3_routing_struc(n)%rivelv(i)
     enddo
  enddo

  do i=1,HYMAP3_routing_struc(n)%nseqall
     do m=1,LIS_rc%nensem(n)
        t = (i-1)*LIS_rc%nensem(n)+m
        if(sfcelv(t).le.HYMAP3_routing_struc(n)%rivelv(i)) then
           ensCheck(i) = .false.
        endif

        if(sfcelv(t).ne.sfcelv(i*LIS_rc%nensem(n))) then
           diffCheck(i) = .true.
        endif
     enddo
     if(.not.ensCheck(i).and.diffCheck(i)) then
        call HYMAP3_reorderEnsForOutliers(&
             LIS_rc%nensem(n),&
             sfcelv((i-1)*LIS_rc%nensem(n)+1:i*LIS_rc%nensem(n)),&
             HYMAP3_routing_struc(n)%rivelv(i))
     endif
  enddo

  do i=1,HYMAP3_routing_struc(n)%nseqall
     do m=1,LIS_rc%nensem(n)
        t = (i-1)*LIS_rc%nensem(n)+m
        HYMAP3_routing_struc(n)%sfcelv(i,m) = &
             sfcelv(t)
     enddo
  enddo

  !update surface level and then update the storage.

  do i=1,HYMAP3_routing_struc(n)%nseqall
     do m=1,LIS_rc%nensem(n)
        t = (i-1)*LIS_rc%nensem(n)+m

        elevtn=dble(HYMAP3_routing_struc(n)%elevtn(i))
        fldhgt=dble(HYMAP3_routing_struc(n)%fldhgt(i,:))
        fldstomax = dble(HYMAP3_routing_struc(n)%fldstomax(i,:,m))
        rivelv=dble(HYMAP3_routing_struc(n)%rivelv(i))
        grarea = dble(HYMAP3_routing_struc(n)%grarea(i))
        rivstomax = dble(HYMAP3_routing_struc(n)%rivstomax(i,m))
        rivlen = dble(HYMAP3_routing_struc(n)%rivlen(i))
        rivwth = dble(HYMAP3_routing_struc(n)%rivwth(i,m))
        elv = dble(HYMAP3_routing_struc(n)%sfcelv(i,m))

        call HYMAP3_get_volume_profile(&
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

        HYMAP3_routing_struc(n)%rivsto(i,m) = real(vol)
        HYMAP3_routing_struc(n)%fldsto(i,m) = 0.0
     enddo
  enddo

end subroutine HYMAP3_setWL

subroutine HYMAP3_reorderEnsForOutliers(nensem, statevec, minvalue)

  implicit none

  integer, intent(in)              :: nensem
  real, intent(out)                :: statevec(nensem)
  real, intent(in)                 :: minvalue

  real                 :: minvT, maxvT, minvG, maxvG
  integer              :: k
  real                 :: spread_total, spread_good, spread_ratio

  !Ensemble spread (total and with 'good' ensemble members
  minvT = 1E10
  maxvT = -1E10
  minvG = 1E10
  maxvG = -1E10

  do k=1,nensem

     if(statevec(k).lt.minvT) then
        minvT = statevec(k)
     endif
     if(statevec(k).gt.maxvT) then
        maxvT = statevec(k)
     endif

     if(statevec(k).gt.minvalue) then
        if(statevec(k).lt.minvG) then
           minvG = statevec(k)
        endif
        if(statevec(k).gt.maxvG) then
           maxvG = statevec(k)
        endif
     endif
  enddo

  if(minvG.eq.1E10.and.maxvG.eq.-1E10) then
     statevec = minvalue
  else
     spread_total = (maxvT - minvT)
     spread_good  = (maxvG - minvG)

     spread_ratio = spread_good/spread_total

     !rescale the ensemble
     do k=1,nensem-1
        statevec(k) = statevec(nensem) + &
             (statevec(k) - statevec(nensem))*spread_ratio
     enddo
  endif

end subroutine HYMAP3_reorderEnsForOutliers
