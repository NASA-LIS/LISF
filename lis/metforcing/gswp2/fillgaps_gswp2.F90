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
! !ROUTINE: fillgaps_gswp2
! \label{fillgaps_gswp2}
! 
! !INTERFACE:
subroutine fillgaps_gswp2(n,ip,varfield)
! !USES:
  use LIS_coreMod,       only : LIS_rc,LIS_domain
  use LIS_logMod,        only : LIS_logunit, LIS_endrun
  use gswp2_forcingMod,   only : gswp2_struc

  implicit none
! !USES: 
  integer, intent(in)    :: n
  integer, intent(in)    :: ip
  real, intent(inout)    :: varfield(LIS_rc%lnc(n),LIS_rc%lnr(n)) 
! !DESCRIPTION:
!   This subroutine fills in invalid grid points introduced due to 
!   mismatches in landmasks in LIS and GSWP2. 
!   This routine assumes that the undef
!   or invalid value is the LIS undefined value. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[ip]
!    interpolation option
!  \item[varfield]
!    updated output field
!  \end{description}
!
!EOP
  integer                :: c,r
  logical                :: foundPt
  integer                :: i,j,str,enr,stc,enc,kk
  integer                :: try

  try = 0 
  if(ip.eq.1) then 
     if(gswp2_struc(n)%fillflag1) then !This will be done once 
        gswp2_struc(n)%smask1 = 0
        do r=1,LIS_rc%lnr(n)
           do c=1,LIS_rc%lnc(n)
              if((LIS_domain(n)%gindex(c,r).ne.-1).and.&
                   varfield(c,r).eq.LIS_rc%udef) then !mismatch
                 gswp2_struc(n)%smask1(c,r) = 1
              endif
           enddo
        enddo
        gswp2_struc(n)%fillflag1 = .false. 
     endif
     
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(gswp2_struc(n)%smask1(c,r).eq.1) then 
              foundPt = .false.
              kk = 1
              try = 0 
              do while(.not.foundPt) 
                 try = try +1
                 str = max(r-kk,1)
                 enr = min(r+kk,LIS_rc%lnr(n))
                 stc = max(c-kk,1)
                 enc = min(c+kk,LIS_rc%lnc(n))
                 do j=str,enr
                    do i=stc,enc
                       if(LIS_domain(n)%gindex(i,j).ne.-1&
                            .and.gswp2_struc(n)%smask1(i,j).ne.1) then 
                          varfield(c,r) = varfield(i,j)
                          foundPt = .true.
                          exit
                       endif
                    enddo
                 enddo
                 kk = kk+1
                 if(try.gt.100) then 
                    write(LIS_logunit,*) 'GSWP2 fillgaps failed, stopping..',try,kk,c,r
                    call LIS_endrun()
                 endif
              enddo
           endif
        enddo
     enddo
  else
     write(LIS_logunit,*) 'GSWP2 fillgaps is implemented only for '
     write(LIS_logunit,*) 'bilinear interpolation method. Program stopping ....'
     call LIS_endrun()
  endif

end subroutine fillgaps_gswp2
