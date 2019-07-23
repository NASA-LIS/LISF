!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: read_SACHTET356_gfrac
! \label{read_SACHTET356_gfrac}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar;  Initial Specification
!  12 Feb 2013: K. Arsenault; Modified for use with SACHTET greenness parameter
!
! !INTERFACE:
 subroutine read_SACHTET356_gfrac(n, array, maskarray)

! !USES:
  use ESMF
  use LDT_coreMod, only : LDT_rc, LDT_domain
  use LDT_logMod,  only : LDT_logunit, LDT_getNextUnitNumber, &
                          LDT_releaseUnitNumber, LDT_endrun
  use LDT_xmrg_reader
  use LDT_gfracMod

  implicit none

! !ARGUMENTS: 
  integer, intent(in)    :: n    ! nest index
  real,    intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real, optional, intent(inout) :: maskarray(LDT_rc%lnc(n),LDT_rc%lnr(n))

! !DESCRIPTION:
!  This subroutine retrieves the greenness fraction climatology for the 
!  specified month and returns the values in the latlon projection
!  
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[array]
!   output field with the retrieved greenness fraction
!  \item[maskarray]
!   optional input field of reading in the mask array
!  \end{description}
!
!EOP      
  integer :: ftn
  integer :: c,r
  real    :: xfactor
  logical :: file_exists
! _________________________________________________

  array = LDT_rc%udef

  inquire(file=trim(LDT_gfrac_struc(n)%gfracFile), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) "[ERR] Greenness map ",trim(LDT_gfrac_struc(n)%gfracFile)," not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif
  select case ( LDT_gfrac_struc(n)%gfrac_gridtransform )
    case( "none" )
     write(LDT_logunit,*) "[INFO] Reading SAC-HTET greenness file: ",&
           trim(LDT_gfrac_struc(n)%gfracFile)
    case default
      write(LDT_logunit,*) "[ERR] The spatial transform option selected for SAC-HTET "
      write(LDT_logunit,*) "  greenness file is not recognized nor recommended."
      write(LDT_logunit,*) "  Please select for now:  none "
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
  end select

  call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
           LDT_rc%gridDesc(n,:), LDT_gfrac_struc(n)%gfracfile, &
           LDT_rc%udef, array )


end subroutine read_SACHTET356_gfrac
