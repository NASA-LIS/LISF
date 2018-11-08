!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT)
!
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: interp_stndata
! \label{interp_stndata}
!
! !REVISION HISTORY: 
!   08 Dec 04 Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine interp_stndata(wt,undef,vari,varo,npts,nstns)

  use LDT_coreMod,  only : LDT_rc

  implicit none
! !ARGUMENTS: 

  integer :: nstns,npts
  real    :: undef
  real    :: wt(npts,nstns)
  real    :: vari(nstns)
  real    :: varo(npts)  

! !DESCRIPTION: 
!  This subroutine interpolates station data into a gridded set. 
!  Currently works only on a lat/lon grid
!
!  The arguments are:
!  \begin{description}
!  \item[nstns]
!    number of stations
!  \item[npts]
!    number of points in the output field
!  \item[vari]
!    input variable (array of observations from the stations)
!  \item[varo]
!    output variable (gridded data on the output field)
!  \item[wt]
!    interpolation weights 
!  \end{description}
!EOP

  integer:: i,j
  integer :: nundefs

  varo = 0.0
  nundefs = 0

  do i=1,npts
     do j=1,nstns
        if(vari(j).ne.undef) then 
           varo(i) = varo(i) + wt(i,j)*vari(j)
        else
           nundefs = nundefs + 1
        endif
!        if(i==1.and.j==1) print*, i,j,varo(i),wt(i,j),vari(j),npts
     enddo
     if(nundefs .eq. nstns) then !all stns have undef values. 
        varo(i) = LDT_rc%udef
     endif
     nundefs = 0 
!     print*, i,varo(i)
  enddo
end subroutine interp_stndata
