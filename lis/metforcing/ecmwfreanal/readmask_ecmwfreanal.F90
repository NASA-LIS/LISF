!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readmask_ecmwfreanal
!  \label{readmask_ecmwfreanal}
!
! !REVISION HISTORY:
!  16 Apr 2002: Urszula Jambor; Initial Code
!  24 Nov 2003: Sujay Kumar; Included in LIS
! 
! !INTERFACE:
subroutine readmask_ecmwfreanal(n)
! !USES: 
  use LIS_coreMod, only : LIS_rc  ! LIS non-model-specific 1-D variables
  use LIS_logMod, only : LIS_logunit
  use ecmwfreanal_forcingMod, only : ecmwfreanal_struc

  implicit none
! !ARGUMENTS:   
  integer, intent(in) :: n 

! !DESCRIPTION:
! Reads in land-sea mask for 0.5 degree Reanalysis ECMWF data set,
! changes grid from ECMWF convention to LIS convention by calling
! {\tt ecmwfreanalgrid\_2\_lisgrid}, and transforms grid array into 1D vector 
! for later use. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[ecmwfreanalgrid\_2\_lisgrid](\ref{ecmwfreanalgrid_2_lisgrid}) \newline
!    transform the ECMWF Reanalysis data to the LIS grid 
!  \end{description}
!EOP

  integer :: i, j, c
  integer :: istat
  integer :: lmask(ecmwfreanal_struc(n)%ncold,ecmwfreanal_struc(n)%nrold)
  real    :: rlmask(ecmwfreanal_struc(n)%ncold,ecmwfreanal_struc(n)%nrold)

  !----------------------------------------------------------------
  !Open land-sea mask file
  !----------------------------------------------------------------
  write(LIS_logunit,*) 'opening mask file ',ecmwfreanal_struc(n)%emaskfile
  open(88,file=ecmwfreanal_struc(n)%emaskfile, form='formatted', iostat=istat)
  if (istat/=0) then
     write(LIS_logunit,*) 'Problem opening file: ', trim(ecmwfreanal_struc(n)%emaskfile)
  else
     do i = 1,ecmwfreanal_struc(n)%nrold
        read(88,'(720(2x,i1))') (lmask(j,i),j=1,ecmwfreanal_struc(n)%ncold)
     end do
  end if
  close(88)

  !----------------------------------------------------------------
  !Change grid from ECMWF convention to GLDAS convention
  !----------------------------------------------------------------

  rlmask = real(lmask)
  call ecmwfreanalgrid_2_lisgrid(ecmwfreanal_struc(n)%ncold, ecmwfreanal_struc(n)%nrold, rlmask)
  lmask = int(rlmask)

  !----------------------------------------------------------------
  !Transfer 2D grid to 1D vector
  !----------------------------------------------------------------

  c = 0
  do i=1, ecmwfreanal_struc(n)%nrold
     do j=1, ecmwfreanal_struc(n)%ncold
        c = c + 1
        ecmwfreanal_struc(n)%remask1d(c) = lmask(j,i)
     enddo
  enddo

end subroutine readmask_ecmwfreanal

