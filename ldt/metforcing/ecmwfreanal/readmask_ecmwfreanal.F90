!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readmask_ecmwfreanal
!  \label{readmask_ecmwfreanal}
!
! !REVISION HISTORY:
!  16 Apr 2002: Urszula Jambor; Initial Code
!  24 Nov 2003: Sujay Kumar; Included in LDT
! 
! !INTERFACE:
subroutine readmask_ecmwfreanal(n)
! !USES: 
  use LDT_coreMod, only : LDT_rc  ! LDT non-model-specific 1-D variables
  use LDT_logMod,  only : LDT_logunit
  use ecmwfreanal_forcingMod, only : ecmwfreanal_struc

  implicit none
! !ARGUMENTS:   
  integer, intent(in) :: n 

! !DESCRIPTION:
! Reads in land-sea mask for 0.5 degree Reanalysis ECMWF data set,
! changes grid from ECMWF convention to LDT convention by calling
! {\tt ecmwfreanalgrid\_2\_ldtgrid}, and transforms grid array into 1D vector 
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
!  \item[ecmwfreanalgrid\_2\_ldtgrid](\ref{ecmwfreanalgrid_2_ldtgrid}) \newline
!    transform the ECMWF Reanalysis data to the LDT grid 
!  \end{description}
!EOP

  integer :: i, j, c
  integer :: istat
  integer :: lmask(ecmwfreanal_struc(n)%nc,ecmwfreanal_struc(n)%nr)
  real    :: rlmask(ecmwfreanal_struc(n)%nc,ecmwfreanal_struc(n)%nr)

  !----------------------------------------------------------------
  !Open land-sea mask file
  !----------------------------------------------------------------
  write(LDT_logunit,*) 'Opening mask file ',trim(ecmwfreanal_struc(n)%emaskfile)
  open(88,file=ecmwfreanal_struc(n)%emaskfile, form='formatted', iostat=istat)
  if (istat/=0) then
     write(LDT_logunit,*) 'Problem opening file: ', trim(ecmwfreanal_struc(n)%emaskfile)
  else
     do i = 1,ecmwfreanal_struc(n)%nr
        read(88,'(720(2x,i1))') (lmask(j,i),j=1,ecmwfreanal_struc(n)%nc)
     end do
  end if
  close(88)

  !----------------------------------------------------------------
  !Change grid from ECMWF convention to GLDAS convention
  !----------------------------------------------------------------

  rlmask = real(lmask)
  call ecmwfreanalgrid_2_ldtgrid(ecmwfreanal_struc(n)%nc, &
                                 ecmwfreanal_struc(n)%nr, rlmask)
  lmask = int(rlmask)

  !----------------------------------------------------------------
  !Transfer 2D grid to 1D vector
  !----------------------------------------------------------------

  c = 0
  do i=1, ecmwfreanal_struc(n)%nr
     do j=1, ecmwfreanal_struc(n)%nc
        c = c + 1
        ecmwfreanal_struc(n)%remask1d(c) = lmask(j,i)
     enddo
  enddo

end subroutine readmask_ecmwfreanal

