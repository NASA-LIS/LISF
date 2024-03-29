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
! !INTERFACE:
subroutine polfixs(nm,nx,km,rlat,rlon,ib,lo,go)
  implicit none
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
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
! !FILES USED:
!
! !REVISION HISTORY:
!   04-10-96  Mark Iredell; Initial Specification
! 
!EOP
!BOP
! 
! !ARGUMENTS: 
  integer         :: nm
  integer         :: nx
  integer         :: km
  real            :: rlat(nm)
  real            :: rlon(nm)
  integer         :: ib
  logical*1       :: lo(nx)
  real            :: go(nx)

!
!        
!EOP  
  integer         :: n, k
  real            :: tsp, gnp, gsp, wsp, tnp, wnp
  real, PARAMETER :: rlatnp=89.9995, rlatsp=-89.9995

!  do k=1,km
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
           if(ib.eq.0.or.lo(n)) then
              gnp=gnp+go(n)
              tnp=tnp+1
           endif
        elseif(rlat(n).le.rlatsp) then
           wsp=wsp+1
           if(ib.eq.0.or.lo(n)) then
              gsp=gsp+go(n)
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
              if(ib.ne.0) lo(n)=tnp.ge.wnp/2
              go(n)=gnp
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
              if(ib.ne.0) lo(n)=tsp.ge.wsp/2
              go(n)=gsp
           endif
        enddo
     endif
!  enddo
end subroutine polfixs
