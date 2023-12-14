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
! !ROUTINE: AGRMET_readterrain
!  \label{AGRMET_readterrain}
!
! !REVISION HISTORY:
!
!    1  nov 05  Sujay Kumar, Initial specification
!
! !INTERFACE:
subroutine AGRMET_readterrain(n)
! !USES: 
  use LIS_coreMod,       only : LIS_rc
  use LIS_fileIOMod,     only : LIS_putget
  use LIS_logMod,        only : LIS_logunit,LIS_abort, LIS_endrun
  use AGRMET_forcingMod, only : agrmet_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
! 
! !DESCRIPTION: 
!  This routine reads the AGRMET terrain data, in polar 
!  stereographic projection at 48km. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!   
!  The routines invoked are: 
!  \begin{description}
!   \item[get\_agrmetterrain\_filename](\ref{get_agrmetterrain_filename}) \newline
!    puts together the name of the AGRMET terrain file based on the
!    location of the root directory
!   \item[LIS\_putget](\ref{LIS_putget}) \newline
!    retrieves the terrain data
!   \item[lis\_abort](\ref{LIS_abort}) \newline
!    program aborts in case of error in file open or read
!   \end{description}
!
!EOP
  integer       :: hemi
  logical       :: exists
  character*100 :: name
  character*255 :: message(20)
  character*30  :: routine_name

  data routine_name     / 'AGRMET_readterrain' /
  
  do hemi=1,2
     call get_agrmetterrain_filename(name,agrmet_struc(n)%terrainfile,hemi)
     inquire( file = trim(name), exist = exists)
     if ( exists ) then
        !write(LIS_logunit,*)' '
        write(LIS_logunit,*)'[INFO] READING ', trim(name)
        call LIS_putget( agrmet_struc(n)%alt(:,:,hemi), 'r', name, routine_name, &
                     agrmet_struc(n)%imax, agrmet_struc(n)%jmax )
     else
        write(LIS_logunit,*)
        write(LIS_logunit,*) "*****************************************************"
        write(LIS_logunit,*) "[ERR] LIS: ERROR OPENING FILE:" 
        write(LIS_logunit,*) "[ERR] ", trim(name)
        write(LIS_logunit,*) "[ERR] FILE DOES NOT EXIST."
        write(LIS_logunit,*) "[ERR] LIS WILL ABORT."
        write(LIS_logunit,*) "*****************************************************"
        message    = ' '
        message(1) = 'program:  LIS'
        message(2) = '  routine: AGRMET_readterrain'
        message(3) = '  error opening file '//trim(name)
        message(4) = '  file does not exist'
        message(5) = '  this could indicate serious data discrepancies'
        call lis_abort( message )
        call LIS_endrun

     endif
  enddo

end subroutine AGRMET_readterrain

!BOP
! 
! !ROUTINE: get_agrmetterrain_filename
!  \label{get_agrmetterrain_filename}
! 
! !INTERFACE: 
subroutine get_agrmetterrain_filename(name, dir,hemi)

  implicit none
! !ARGUMENTS:   
  integer, intent(in) :: hemi
  character*100       :: name
  character*100        :: dir
! 
! !DESCRIPTION: 
!  This routines generates the name of the AGRMET terrain file, by 
!  appending the root directory and the hemisphere information. 
!  The name of the file is expected to be: 
!  <dir>/terrain\_<hh>, where hh is the 'nh' or 'sh', for
!  northern and southern hemisphere, respectively. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[hemi]
!    index of the hemisphere (1-NH, 2-SH)
!   \item[dir]
!    path to the directory containing the terrain file
!   \item[name]
!    created filename
!  \end{description}
!EOP
  character*100 :: temp
  character*1 :: fbase(100),fhemi(3)
  integer :: c,i
  write(UNIT=temp, fmt='(a100)') dir  
  read(UNIT=temp, fmt='(100a1)') (fbase(i), i=1,100)

  if(hemi==1) then 
     write(UNIT=temp, fmt='(a3)') '_nh'
     read(UNIT=temp, fmt='(3a1)') fhemi
  else
     write(UNIT=temp, fmt='(a3)') '_sh'
     read(UNIT=temp, fmt='(3a1)') fhemi
  endif
  c = 0
  do i = 1, 100
     if ( (fbase(i) == ' ') .and. (c == 0) ) c = i-1
  end do
  write(UNIT=temp, fmt='(100a1)') (fbase(i), i=1,c), &
       (fhemi(i), i=1,3)

  read(UNIT=temp, fmt='(a100)') name

end subroutine get_agrmetterrain_filename
