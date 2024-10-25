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
! !ROUTINE: LIS_readPEDecSpaceAttributes
! \label{LIS_readPEDecSpaceAttributes}
!
! !REVISION HISTORY:
!  12 Jan 2012: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine LIS_readPEDecSpaceAttributes(decattribfile, &
     nparam, vname,selectOpt,varmin,varmax)
! !USES: 
  use LIS_coreMod,   only : LIS_rc
  use LIS_logMod,    only : LIS_logunit, LIS_endrun, & 
       LIS_getNextUnitNumber, LIS_releaseUnitNumber
  
  implicit none
! !ARGUMENTS: 
  character(len=*)   :: decattribfile
  integer            :: nparam
  character(len=*)   :: vname(nparam)
  integer            :: selectOpt(nparam)
  real               :: varmin(nparam)
  real               :: varmax(nparam)
! 
! !DESCRIPTION: 
!  This routine reads the observation attributes for each
!  observation type. The format of the attributes file is: 
!
!  Variable name 
!  standard deviation  Max Value  Min Value
!
!  The arguments are: 
!  \begin{description}
!   \item[n]           index of nest
!   \item[vname]       name of the observation variables
!   \item[ssdev]       error rate in observation
!   \item[varmin]      minimum value of the variable
!   \item[varmax]      maximum value of the variable
!  \end{description}
!EOP
  integer            :: ftn 
  integer            :: i, count
  character*100      :: temp

  write(LIS_logunit,*) 'Opening attributes for observations ',&
       trim(decattribfile)
  ftn = LIS_getNextUnitNumber()
  open(ftn,file=trim(decattribfile),status='old')
  read(ftn,*)
  count = 0
  do i=1,nparam
     read(ftn,*) selectOpt(i),vname(i), varmin(i),varmax(i)
     count = count+1
     if(count.gt.nparam) then 
        write(LIS_logunit,*) 'The number of variables in the decision '
        write(LIS_logunit,*) 'space attributes file does not match the '
        write(LIS_logunit,*) 'total number of variables (',nparam,')'
        write(LIS_logunit,*) 'Program stopping....'
        call LIS_endrun()
     endif
  enddo
  call LIS_releaseUnitNumber(ftn)
end subroutine LIS_readPEDecSpaceAttributes
