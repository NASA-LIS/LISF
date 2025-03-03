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
! !ROUTINE: read_RUC_NCEPLIS_slopetype
! \label{read_RUC_NCEPLIS_slopetype}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  17 Jan 2011: David Mocko, added slope type index
!
! !INTERFACE:
subroutine read_RUC_NCEPLIS_slopetype(n, array)

! !USES:
  use LDT_coreMod,        only : LDT_rc
  use LDT_logMod,         only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod,      only : readLISdata 
  use RUC_parmsMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  real, intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),1)

! !DESCRIPTION:
!  This subroutine retrieves static, slope type data and reprojects
!  it to the latlon projection. 
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[array]
!   output field with the retrieved slope type data
!  \end{description}
!EOP

  integer :: ftn
  integer :: c, r
  logical :: file_exists
! ____________________________

  array = LDT_rc%udef

  inquire(file=trim(RUC_struc(n)%slopetypefile), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) "SLOPETYPE map ",trim(RUC_struc(n)%slopetypefile)," not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif
  select case ( RUC_struc(n)%slopetype_gridtransform )
    case( "none", "neighbor", "mode" )  ! Discrete data type
      write(LDT_logunit,*) "MSG: Reading RUC_NCEPLIS slopetype file: ",&
            trim(RUC_struc(n)%slopetypefile)
  case default
     write(LDT_logunit,*) "ERR: Since the Slope type field involves discrete data values,"
     write(LDT_logunit,*) "     only 'neighbor' or 'mode' are currently supported spatial"
     write(LDT_logunit,*) "     transform types.  Please check your entries for this parameter."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  end select

  ftn = LDT_getNextUnitNumber()
  open(ftn, file=RUC_struc(n)%slopetypefile, access='direct',status='old', &
       form="unformatted", recl=4)

  call readLISdata(n, ftn, RUC_struc(n)%slopetype_proj, &
                   RUC_struc(n)%slopetype_gridtransform,&
                   RUC_struc(n)%slopetype_gridDesc(:), 1, array)   ! 1 indicates 2D layer
 
  do r = 1, LDT_rc%lnr(n)
     do c = 1, LDT_rc%lnc(n)
        if( array(c,r,1) > 9. ) array(c,r,1) = 9.  ! Set slopetype upper limit to 9  
     enddo
  enddo

  call LDT_releaseUnitNumber(ftn)
  write(LDT_logunit, *) "MSG: Done reading RUC_NCEPLIS slope type file"

end subroutine read_RUC_NCEPLIS_slopetype
