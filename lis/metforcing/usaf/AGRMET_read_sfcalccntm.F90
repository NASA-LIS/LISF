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
! !ROUTINE: AGRMET_read_sfcalccntm
! \label{AGRMET_read_sfcalccntm}
! 
! !REVISION HISTORY:
!
!    1  nov 05  Sujay Kumar, Initial specification
!   29  May 15  Ryan Ruhge, Updated for faster readin
!
! !INTERFACE:
subroutine AGRMET_read_sfcalccntm(n)
! !USES: 
  use LIS_coreMod,       only : LIS_rc, LIS_ews_halo_ind, LIS_ewe_halo_ind,&
                                LIS_nss_halo_ind, LIS_nse_halo_ind, LIS_localPet
  use LIS_logMod,        only : LIS_logunit, LIS_abort, LIS_endrun, &
       LIS_getNextUnitNumber, LIS_releaseUnitNumber
  use AGRMET_forcingMod, only : agrmet_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
! 
! !DESCRIPTION: 
!  This routine reads the spreading radii used in the barnes analyses
!  for AGRMET surface calculations. The data is in polar stereographic
!  projection at 48km resolution. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!   
!  The routines invoked are: 
!  \begin{description}
!   \item[lis\_abort](\ref{LIS_abort}) \newline
!    program aborts in case of error in file open or read
!   \end{description}
!
!EOP

  logical       :: exists
  character*255 :: message(20)
  integer       :: ftn
  real          :: data_in(LIS_rc%gnc(n), LIS_rc%gnr(n))
  integer       :: istat
  
  inquire( file = trim(agrmet_struc(n)%sfcntmfile), exist = exists)
  if ( exists ) then
     !write(LIS_logunit,*)' '
     write(LIS_logunit,*)'[INFO] READING ', trim(agrmet_struc(n)%sfcntmfile)
     ftn= LIS_getNextUnitNumber()
     open(ftn, file=trim(agrmet_struc(n)%sfcntmfile), access='direct',&
          status='old', form="unformatted", recl=LIS_rc%gnr(n)*LIS_rc%gnc(n)*4)
     
     read(ftn, rec=1, iostat=istat) data_in
     agrmet_struc(n)%irad(:,:) =  data_in(&
        LIS_ews_halo_ind(n,LIS_localPet+1):&
        LIS_ewe_halo_ind(n,LIS_localPet+1), &
        LIS_nss_halo_ind(n,LIS_localPet+1): &
        LIS_nse_halo_ind(n,LIS_localPet+1))
     call LIS_releaseUnitNumber(ftn)

  else
     write(LIS_logunit,*)
     write(LIS_logunit,*) "*****************************************************"
     write(LIS_logunit,*) "[ERR] LIS: ERROR OPENING FILE:" 
     write(LIS_logunit,*) "[ERR] ", trim(agrmet_struc(n)%sfcntmfile)
     write(LIS_logunit,*) "[ERR] FILE DOES NOT EXIST."
     write(LIS_logunit,*) "[ERR] LIS WILL ABORT."
     write(LIS_logunit,*) "*****************************************************"
     message    = ' '
     message(1) = 'program:  LIS'
     message(2) = '  routine: AGRMET_read_sfcalccntm'
     message(3) = '  error opening file '//trim(agrmet_struc(n)%sfcntmfile)
     message(4) = '  file does not exist'
     message(5) = '  this could indicate serious data discrepancies'
     call lis_abort( message )
     call LIS_endrun
  endif
end subroutine AGRMET_read_sfcalccntm
