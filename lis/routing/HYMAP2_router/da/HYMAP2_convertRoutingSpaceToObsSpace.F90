!BOP
! !ROUTINE: HYMAP2_convertRoutingSpaceToObsSpace
! \label{HYMAP2_convertRoutingSpaceToObsSpace}
!
! !INTERFACE:
  subroutine HYMAP2_convertRoutingSpaceToObsSpace( &
       n, &
       k, &
       mvar, &
       ovar)

!USES:
    use HYMAP2_routingMod
    use LIS_DAobservationsMod
    use LIS_coreMod

! !ARGUMENTS:
    integer,          intent(in) :: n
    integer,          intent(in) :: k
    real                         :: mvar(HYMAP2_routing_struc(n)%nseqall)
    real                         :: ovar(LIS_rc%obs_ngrid(k))
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
!   \item [mvar]
!     variable in the patch space
!   \item [ovar]
!     variable in the observation ensemble space
!  \end{description}
!
!EOP

    integer                      :: c,r,i,g,gid
    real                         :: lis_gvar(LIS_rc%lnc(n)*LIS_rc%lnr(n))
    integer                      :: nlis_gvar(LIS_rc%lnc(n)*LIS_rc%lnr(n))
    logical*1                    :: li(LIS_rc%lnc(n)*LIS_rc%lnr(n))
    logical*1                    :: lo(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
    real                         :: obs_gvar(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
    integer                      :: iret

    external :: upscaleByAveraging
    external :: neighbor_interp

    ovar = LIS_rc%udef
    lis_gvar  = 0.0
    nlis_gvar = 0
    obs_gvar = LIS_rc%udef

    do i=1,HYMAP2_routing_struc(n)%nseqall
       c = HYMAP2_routing_struc(n)%seqx(i)
       r = HYMAP2_routing_struc(n)%seqy(i)
       gid = c+(r-1)*LIS_rc%lnc(n)
       lis_gvar(gid)  = lis_gvar(gid) + mvar(i)
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
             ovar(LIS_obs_domain(n,k)%gindex(c,r)) = &
                  obs_gvar(c+(r-1)*LIS_rc%obs_lnc(k))
          endif
       enddo
    enddo

  end subroutine HYMAP2_convertRoutingSpaceToObsSpace
