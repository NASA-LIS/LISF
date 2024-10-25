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
! 
! !ROUTINE: polfixs
!  \label{polfixs}
!
! !REVISION HISTORY:
!   04-10-96  Mark Iredell; Initial Specification
!
! !INTERFACE:
subroutine polfixs(nm,nx,km,rlat,rlon,ib,lo,go)
  implicit none
! !ARGUMENTS: 
  integer         :: nm
  integer         :: nx
  integer         :: km
  real            :: rlat(nm)
  real            :: rlon(nm)
  integer         :: ib(km)
  logical*1       :: lo(nx,km)
  real            :: go(nx,km)

!
! !DESCRIPTION: 
! This subroutine averages multiple pole scalar values
! on a latitude/longitude grid.  bitmaps may be averaged too.
!
!  The arguments are:
!  \begin{description}
!  \item[nm]
!    number of grid points
!  \item[nx]
!    leading dimensition of fields
!  \item[rlat]
!    latitudes in degrees
!  \item[rlon]
!    longitudes in degrees
!  \item[ib]
!    integer bitmap flags
!  \item[lo]
!    logical bitmaps 
!  \item[go]
!    returned scalar value
!  \end{description}
!        
!EOP  
  integer         :: n, k
  real            :: tsp, gnp, gsp, wsp, tnp, wnp
  real, PARAMETER :: rlatnp=89.9995, rlatsp=-89.9995

  do k=1,km
     wnp=0.0
     gnp=0.0
     tnp=0.0
     wsp=0.0
     gsp=0.0
     tsp=0.0
     !  average multiple pole values
     do n=1,nm
        if(rlat(n).ge.rlatnp) then
           wnp=wnp+1
           if(ib(k).eq.0.or.lo(n,k)) then
              gnp=gnp+go(n,k)
              tnp=tnp+1
           endif
        elseif(rlat(n).le.rlatsp) then
           wsp=wsp+1
           if(ib(k).eq.0.or.lo(n,k)) then
              gsp=gsp+go(n,k)
              tsp=tsp+1
           endif
        endif
     enddo
     !  distribute average values back to multiple poles
     if(wnp.gt.1) then
        if(tnp.ge.wnp/2) then
           gnp=gnp/tnp
        else
           gnp=0.
        endif
        do n=1,nm
           if(rlat(n).ge.rlatnp) then
              if(ib(k).ne.0) lo(n,k)=tnp.ge.wnp/2
              go(n,k)=gnp
           endif
        enddo
     endif
     if(wsp.gt.1) then
        if(tsp.ge.wsp/2) then
           gsp=gsp/tsp
        else
           gsp=0.
        endif
        do n=1,nm
           if(rlat(n).le.rlatsp) then
              if(ib(k).ne.0) lo(n,k)=tsp.ge.wsp/2
              go(n,k)=gsp
           endif
        enddo
     endif
  enddo
end subroutine polfixs
