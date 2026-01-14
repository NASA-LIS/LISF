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
! !ROUTINE: HYMAP2_updateWL
!  \label{HYMAP2_updateWL}
!
! !REVISION HISTORY:
!  07 Nov 2019: Sujay Kumar, Initial Specification
!
! !INTERFACE:
subroutine HYMAP2_updateWL(n, Routing_State, Routing_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_routingMod
  use HYMAP2_routingMod
  use HYMAP2_daWL_Mod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: Routing_State
  type(ESMF_State)       :: Routing_Incr_State
!
! !DESCRIPTION:
!  
!  This routine updates the water level prognostic variables 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[Routing\_State] ESMF State container for Routing state variables \newline
!  \item[Routing\_Incr\_State] ESMF State container for Routing increment state variables \newline
!  \end{description}
!
!EOP

  type(ESMF_Field)       :: sfcelevField
  type(ESMF_Field)       :: sfcelevIncrField
  real, pointer          :: sfcelev(:)
  real, pointer          :: sfcelevIncr(:)
  real, allocatable      :: sfcelevIncr_tmp(:)
  integer, allocatable   :: nsfcelevIncr_tmp(:)
  integer                :: t,i,i1,m
  integer                :: ix,iy
  integer                :: c,r,c1,c2,r1,r2,t1
  integer                :: siteid
  integer                :: status
  real                   :: maxdistance, weight
  real                   :: localWeight(LIS_rc%lnc(n),LIS_rc%lnr(n))

  call ESMF_StateGet(Routing_State,"Surface elevation",sfcelevField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Surface elevation failed in HYMAP2_updateWL")

  call ESMF_FieldGet(sfcelevField,localDE=0,farrayPtr=sfcelev,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Surface elevation failed in HYMAP2_updateWL")

  call ESMF_StateGet(Routing_Incr_State,"Surface elevation",sfcelevIncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Surface elevation failed in HYMAP2_updateWL")

  call ESMF_FieldGet(sfcelevIncrField,localDE=0,farrayPtr=sfcelevIncr,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Surface elevation failed in HYMAP2_updateWL")

  allocate(sfcelevIncr_tmp(HYMAP2_routing_struc(n)%nseqall*LIS_rc%nensem(n)))
  allocate(nsfcelevIncr_tmp(HYMAP2_routing_struc(n)%nseqall*LIS_rc%nensem(n)))
  sfcelevIncr_tmp = 0.0
  nsfcelevIncr_tmp = 0

  do i=1,HYMAP2_routing_struc(n)%nseqall
     do m=1,LIS_rc%nensem(n)
        t = (i-1)*LIS_rc%nensem(n)+m
        if (HYMAP2_daWL_struc(n)%useLocalUpd.eq.1) then 
           if(abs(sfcelevIncr(t)).gt.0) then 
              localweight = -9999.0

              ix = HYMAP2_routing_struc(n)%seqx(i)
              iy = HYMAP2_routing_struc(n)%seqy(i)
              
!              call HYMAP2_map_l2g_index(n,i,siteid)
              siteid = HYMAP2_dawl_struc(n)%sites(ix,iy)
              localweight(:,:) = &
                   HYMAP2_daWL_struc(n)%localWeight(&
                   LIS_ews_halo_ind(n,LIS_localPet+1):&
                   LIS_ewe_halo_ind(n,LIS_localPet+1), &
                   LIS_nss_halo_ind(n,LIS_localPet+1):&
                   LIS_nse_halo_ind(n,LIS_localPet+1), &
                   siteid)
              
              c1=max(1,ix-HYMAP2_daWL_struc(n)%localupdDX)
              c2=min(LIS_rc%lnc(n),ix+HYMAP2_daWL_struc(n)%localupdDX)
              r1=max(1,iy-HYMAP2_daWL_struc(n)%localupdDX)
              r2=min(LIS_rc%lnr(n),iy+HYMAP2_daWL_struc(n)%localupdDX) 

              maxdistance = 0.0
              do r=r1,r2
                 do c=c1,c2
                    if(localweight(c,r).gt.maxdistance) then
                       maxdistance = localweight(c,r)
                    endif
                 enddo
              enddo
              

              do r=r1,r2
                 do c=c1,c2
                    i1 = LIS_routing(n)%gindex(c,r)
                    if(i1.gt.0) then 
                       t1 = (i1-1)*LIS_rc%nensem(n)+m

                       if(sfcelev(t1).ne.-9999.0.and.&
                            localweight(c,r).ne.-9999.0.and.&
                            localweight(ix,iy).ne.-9999.0) then 
                          weight = exp(-localweight(c,r)**2/&
                               (2*maxdistance**2))
                          sfcelevIncr_tmp(t1) = sfcelevIncr(t)*weight
                          nsfcelevIncr_tmp(t1) = nsfcelevIncr_tmp(t1) + 1
!                          sfcelev(t1) = sfcelev(t1) + &
!                               sfcelevIncr(t)*weight
                       endif
                    endif
                 enddo
              enddo

           endif
        else
           sfcelevIncr_tmp(t) = sfcelevIncr(t)
           nsfcelevIncr_tmp(t) = nsfcelevIncr_tmp(t) + 1
!           sfcelev(t) = sfcelev(t) + sfcelevIncr(t)
        endif
     enddo
  enddo

  do i=1,HYMAP2_routing_struc(n)%nseqall
     do m=1,LIS_rc%nensem(n)
        t = (i-1)*LIS_rc%nensem(n)+m
        if(nsfcelevIncr_tmp(t).gt.0) then 
          sfcelev(t) = sfcelev(t)+sfcelevIncr_tmp(t)/&
               nsfcelevIncr_tmp(t)

       endif
    enddo
 enddo
 
 deallocate(sfcelevIncr_tmp)
 deallocate(nsfcelevIncr_tmp)

end subroutine HYMAP2_updateWL
