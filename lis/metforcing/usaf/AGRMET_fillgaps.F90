!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: AGRMET_fillgaps
!  \label{AGRMET_fillgaps}
! 
! !INTERFACE:
subroutine AGRMET_fillgaps(n,ip,varfield)
! !USES:
  use LIS_coreMod,       only : LIS_rc, LIS_domain
  use LIS_logMod,        only : LIS_logunit, LIS_endrun, LIS_abort
  use AGRMET_forcingMod, only : agrmet_struc

  implicit none
! !USES: 
  integer, intent(in)    :: n
  integer, intent(in)    :: ip
  real, intent(inout)    :: varfield(LIS_rc%lnc(n),LIS_rc%lnr(n)) 
! !DESCRIPTION:
!   This subroutine fills in invalid grid points introduced due to 
!   reprojection from PS to lat/lon. This routine assumes that the undef
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
  character*255                 :: message     (20)


  try = 0 
  if(ip.eq.1) then 
     if(agrmet_struc(n)%fillflag1) then !This will be done once 
        agrmet_struc(n)%smask1 = 0
        do r=1,LIS_rc%lnr(n)
           do c=1,LIS_rc%lnc(n)
              if((LIS_domain(n)%gindex(c,r).ne.-1).and.&
                   varfield(c,r).eq.LIS_rc%udef) then !mismatch
                 agrmet_struc(n)%smask1(c,r) = 1
              endif
           enddo
        enddo
        agrmet_struc(n)%fillflag1 = .false. 
     endif
     
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(agrmet_struc(n)%smask1(c,r).eq.1) then 
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
                            .and.agrmet_struc(n)%smask1(i,j).ne.1) then 
                          varfield(c,r) = varfield(i,j)
                          foundPt = .true.
                          exit
                       endif
                    enddo
                 enddo
                 kk = kk+1
                 if(try.gt.100) then 
!                    write(LIS_logunit,*) 'AGRMET fillgaps failed, stopping..',try,kk,c,r
                    ! EMK...Force abort with logging!
                    write(LIS_logunit,*) &
                         '[ERR] AGRMET fillgaps failed, stopping..',kk,c,r
                    flush(LIS_logunit)
                    message(1) = 'program: LIS'
                    message(2) = ' routine: AGRMET_fillgaps'
                    message(3) = ' Cannot fill gap in forcing data!'
                    call LIS_abort(message)
                    call LIS_endrun()
                 endif
              enddo
           endif
        enddo
     enddo
  elseif(ip.eq.2) then 
     if(agrmet_struc(n)%fillflag2) then !This will be done once 
        agrmet_struc(n)%smask2 = 0
        do r=1,LIS_rc%lnr(n)
           do c=1,LIS_rc%lnc(n)
              if((LIS_domain(n)%gindex(c,r).ne.-1).and.&
                   varfield(c,r).eq.LIS_rc%udef) then !mismatch
                 agrmet_struc(n)%smask2(c,r) = 1
              endif
           enddo
        enddo
        agrmet_struc(n)%fillflag2 = .false. 
     endif
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(agrmet_struc(n)%smask2(c,r).eq.1) then 
              foundPt = .false.
              kk = 1
              do while(.not.foundPt) 
                 str = max(r-kk,1)
                 enr = min(r+kk,LIS_rc%lnr(n))
                 stc = max(c-kk,1)
                 enc = min(c+kk,LIS_rc%lnc(n))
                 do j=str,enr
                    do i=stc,enc
                       if(LIS_domain(n)%gindex(i,j).ne.-1&
                            .and.agrmet_struc(n)%smask2(i,j).ne.1) then 
                          varfield(c,r) = varfield(i,j)
                          foundPt = .true.
                          exit
                       endif
                    enddo
                 enddo
                 kk = kk+1
              enddo
           endif
        enddo
     enddo
  endif

end subroutine AGRMET_fillgaps
