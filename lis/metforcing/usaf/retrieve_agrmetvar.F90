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
! !ROUTINE: retrieve_agrmetvar
! \label{retrieve_agrmetvar}
!
! !REVISION HISTORY:
! 29Jul2005; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine retrieve_agrmetvar(name,nc,nr,varfield)
! !USES: 
  use LIS_logMod, only : LIS_logunit, LIS_endrun

  implicit none
! !ARGUMENTS: 
  character*100       :: name
  integer             :: nc,nr
  real, intent(inout) :: varfield(nc,nr)
!
! !DESCRIPTION:
!  Routine to retrieve variables in the polar stereographic projection, 
!  in a hemisphere
!  
!  The arguments are: 
!  \begin{description}
!  \item[name]
!    name of the file to be read from 
!  \item[nc]
!    number of columns (east-west dimension) of the data
!  \item[nr]
!    number of rows (north-south dimension) of the data
!  \item[varfield]
!    the retrieved data
!  \end{description}
!EOP
  logical :: file_exists

  inquire (file=name, exist=file_exists)
  if(file_exists) then 
     open(110,file=name,form='unformatted',access='direct',&
          recl=nc*nr*4)
     read(110,rec=1) varfield
     close(110)
  else
     write(LIS_logunit,*) 'Missing file ',name
     write(LIS_logunit,*) 'Stopping program '
     call LIS_endrun
  endif
end subroutine retrieve_agrmetvar
