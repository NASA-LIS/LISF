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
! !ROUTINE: HYMAP2_getWLpred
! \label{HYMAP2_getWLpred}
!
! !REVISION HISTORY:
! 07 Nov 2019: Sujay Kumar, Initial Specification
!
! !INTERFACE:
subroutine HYMAP2_getWLpred(n, k,obs_pred)
! !USES:
  use ESMF
  use LIS_constantsMod
  use LIS_coreMod
  use LIS_dataAssimMod
  use LIS_DAobservationsMod
  use HYMAP2_routingMod
!EOP

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: k
  real                   :: obs_pred(LIS_rc%obs_ngrid(k),LIS_rc%nensem(n))
!
! !DESCRIPTION:
!
!  Returns the water level obs pred (model's estimate of 
!  observations) for data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[obs\_pred] model's estimate of observations \newline
!  \end{description}
!EOP
  integer                :: i,t,m
  real                   :: wl(HYMAP2_routing_struc(n)%nseqall*LIS_rc%nensem(n))


  do i=1,HYMAP2_routing_struc(n)%nseqall
     do m=1,LIS_rc%nensem(n)
        t=(i-1)*LIS_rc%nensem(n)+m
        wl(t) = HYMAP2_routing_struc(n)%sfcelv(i,m)
     enddo
  enddo

  call HYMAP2_convertPatchSpaceToObsEnsSpace(n,k,&
       wl,&
       obs_pred)
  
end subroutine HYMAP2_getWLpred

!BOP
! !ROUTINE: HYMAP2_convertPatchSpaceToObsEnsSpace
! \label{HYMAP2_convertPatchSpaceToObsEnsSpace}
! 
! !INTERFACE:
  subroutine HYMAP2_convertPatchSpaceToObsEnsSpace(&
       n,&
       k,&
       pvar, &
       ovar)

!USES:
    use HYMAP2_routingMod
    use LIS_DAobservationsMod
    use LIS_coreMod

! !ARGUMENTS: 
    integer,          intent(in) :: n 
    integer,          intent(in) :: k
    real                         :: pvar(HYMAP2_routing_struc(n)%nseqall*LIS_rc%nensem(n))
    real                         :: ovar(LIS_rc%obs_ngrid(k),LIS_rc%nensem(n))
!
! !DESCRIPTION: 
!
!  This routine converts a variable in the patch space to the observation
!  ensemble grid space. If the observation space is at a coarser resolution than 
!  the patch space, then the variable is upscaled. On the other hand, the
!  variable is spatially interpolated to the observation space. These 
!  transformations are done separately for each ensemble member. 
! 
! The arguments are:
!  \begin{description}
!   \item [n]
!     index of the current nest
!   \item [k]
!     index of the DA instance
!   \item [patch\_index]
!     index of the patch to which the variable belong to
!   \item [pvar]
!     variable in the patch space
!   \item [ovar]
!     variable in the observation ensemble space
!  \end{description}
!  
!EOP

    integer                      :: c,r,t,i,m,g,gid
    real                         :: lis_gvar(LIS_rc%lnc(n)*LIS_rc%lnr(n))
    integer                      :: nlis_gvar(LIS_rc%lnc(n)*LIS_rc%lnr(n))
    logical*1                    :: li(LIS_rc%lnc(n)*LIS_rc%lnr(n))
    logical*1                    :: lo(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
    real                         :: obs_gvar(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
    integer                      :: iret

    ovar = LIS_rc%udef
    do m=1,LIS_rc%nensem(n)
       lis_gvar  = 0.0
       nlis_gvar = 0
       obs_gvar = LIS_rc%udef

       do i=1,HYMAP2_routing_struc(n)%nseqall
          t = (i-1)*LIS_rc%nensem(n)+m
          c = HYMAP2_routing_struc(n)%seqx(i)
          r = HYMAP2_routing_struc(n)%seqy(i)
          gid = c+(r-1)*LIS_rc%lnc(n)
          lis_gvar(gid)  = lis_gvar(gid) + pvar(t)
          nlis_gvar(gid) = nlis_gvar(gid) + 1
       enddo       

       li = .false. 
       do g=1,LIS_rc%lnc(n)*LIS_rc%lnr(n)
          if(nlis_gvar(g).ne.0) then 
             lis_gvar(g)  = lis_gvar(g)/ &
                  nlis_gvar(g) 
             li(g) = .true. 
          else
             lis_gvar(g) = LIS_rc%udef
          endif
       enddo

       if(LIS_isatAfinerResolution(n,LIS_obs_domain(n,k)%datares)) then     
          call upscaleByAveraging(&
               LIS_rc%lnc(n)*LIS_rc%lnr(n), &
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
               LIS_rc%udef, &
               LIS_obs_domain(n,k)%nbr_index, &
               li, &
               lis_gvar, &
               lo, &
               obs_gvar)
       else          
          call neighbor_interp(LIS_rc%obs_gridDesc(k,:), &
               li, lis_gvar, lo, obs_gvar,&
               LIS_rc%lnc(n)*LIS_rc%lnr(n), &
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
               LIS_obs_domain(n,k)%rlat, &
               LIS_obs_domain(n,k)%rlon, &
               LIS_obs_domain(n,k)%nbr_index, &
               LIS_rc%udef,iret)
       endif
       
       do r=1,LIS_rc%obs_lnr(k)
          do c=1,LIS_rc%obs_lnc(k)          
             if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                ovar(LIS_obs_domain(n,k)%gindex(c,r),m) = & 
                     obs_gvar(c+(r-1)*LIS_rc%obs_lnc(k))                
             endif
          enddo
       enddo
    enddo
 
  end subroutine HYMAP2_convertPatchSpaceToObsEnsSpace
