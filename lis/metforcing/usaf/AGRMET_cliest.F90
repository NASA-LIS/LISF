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
! !ROUTINE: AGRMET_cliest
! \label{AGRMET_cliest}
!
! !REVISION HISTORY:
!
!
!    04 dec 97  initial version................capt andrus/dnxm(agromet)
!    31 mar 99  changed array and loop sizes to hemisphere size.  
!               .......................................... mr moore/dnxm
!     7 oct 99  ported to ibm sp-2, updated prolog, incorporated
!               FORTRAN 90 features.................capt hidalgo/agrmet
!    11 Mar 10  Changed program names in messages to LIS.
!               ..............................Chris Franks/16WS/WXE/SEMS
!
! 29Jul2005 Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine AGRMET_cliest( n, cest, cliprc, clmult, quad9r)

! !USES: 
  use LIS_coreMod,  only   : LIS_rc
  use LIS_LMLCMod,  only   : LIS_LMLC

  implicit none
 ! !ARGUMENTS: 
  integer                       :: n 
  real                          :: cest(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                          :: cliprc(LIS_rc%lnc(n),LIS_rc%lnr(n))   
  real                          :: clmult  
  real                          :: quad9r
 
! 
! !DESCRIPTION:
!
!    to create a precipitation estimate based on climatological precip  
!    amounts.   
!
!    \textbf{Method}
!    
!    - for all land points  \newline
!      - check validity of 3-hour climo precip  \newline
!        - if valid, replace initialized climo estimate with climo  
!          precip amount multiplied by the mult factor read in from 
!          the control file. \newline
!        - if not valid, print out a message to that effect. \newline
!
! The arguments and variables are: 
! \begin{description}
!  \item[cest]      pure climatological precipitation estimate (mm/3hr)   
!  \item[cliprc]    3-hour climo precip used for estimate 
!  \item[clmult]    alternate monthly weighting factor
!  \item[hemi]      hemisphere (1=nh, 2=sh)
!  \item[i]         loop counter  
!  \item[imax]      number of gridpoints in east/west direction
!  \item[j]         loop counter
!  \item[jmax]      number of gridpoints in north/south direction
!  \item[quad9r]    value of 9999.0 for intializing arrays
! \end{description}
!EOP
  integer :: i    
  integer :: j

!     ------------------------------------------------------------------
!     executable code begins here ...
!     compute estimate by multiplying climatological 3-hour precip  
!     amount by the mult factor which was read in from the control file.
!     ------------------------------------------------------------------

  cest = quad9r
  do j = 1, LIS_rc%lnr(n)
     do i = 1, LIS_rc%lnc(n)

!         --------------------------------------------------------------
!         at land points, see if the climatological amount is valid. if 
!         it is not valid, generate a message to indicate there is a
!         a problem since the climo file should contain a value at every
!         land point.   
!         --------------------------------------------------------------
        
        if( LIS_LMLC(n)%landmask(i,j) .gt. 0 )then
           if( cliprc(i,j) .gt. -9990.0 )then   
              cest(i,j) = cliprc(i,j) * clmult  
           else
              cest(i,j) = quad9r
              !              write(LIS_logunit,6000) hemi, i, j  
           endif
        endif
     enddo
  enddo
  return

!     ------------------------------------------------------------------
!     format section
!     ------------------------------------------------------------------

6000 format (/, 1x, 55('-'), &
          /, 3x, 'program:  LIS'  &
          /, 3x, '  routine:  AGRMET_cliest',   &
          /, 3x, '  an invalid climo amount for land was found',&
          /, 3x, '  hemi, i, j are ', 3(i4), &
          /, 3x, '  value filled with 9999.0',&   
          /, 1x, 55('-'))
  
end subroutine AGRMET_cliest

