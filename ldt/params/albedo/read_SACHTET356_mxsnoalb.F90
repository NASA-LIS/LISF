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
! !ROUTINE: read_SACHTET356_mxsnoalb
!  \label{read_SACHTET356_mxsnoalb}

! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar;  Initial Specification
!  25 Nov 2013: KR Arsenault; Added reader for SAC-HTET model
!
! !INTERFACE:
subroutine read_SACHTET356_mxsnoalb(n,array) 

! !USES:
  use LDT_coreMod,       only : LDT_rc, LDT_domain
  use LDT_logMod,        only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_xmrg_reader
  use LDT_albedoMod
!EOP      

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n
  real, intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n))
!
! !DESCRIPTION:
!  This subroutine retrieves the maximum snow expected over 
!  deep snow and returns the values for SAC-HTET model.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[array]
!   output field with the retrieved max snow albedo
!  \end{description}
!
!EOP      
  integer  :: ftn
  integer  :: c,r
  real     :: xfactor
  logical  :: file_exists

! __________________________________

   array = LDT_rc%udef

   inquire(file=trim(LDT_albedo_struc(n)%mxsnoalbfile), exist=file_exists)
   if(.not.file_exists) then 
      write(LDT_logunit,*) "Max Snow Albedo map ",&
           trim(LDT_albedo_struc(n)%mxsnoalbfile)," not found."
      write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
   endif
   select case ( LDT_albedo_struc(n)%mxsnoalb_gridtransform)
    case( "none" )
      write(LDT_logunit,*) "[INFO] Reading SAC-HTET max snow albedo file: ",&
        trim(LDT_albedo_struc(n)%mxsnoalbfile)
    case default
      write(LDT_logunit,*) "[ERR] The spatial transform option selected for SAC-HTET "
      write(LDT_logunit,*) "   max snow albedo file is not recognized nor recommended."
      write(LDT_logunit,*) "   Please select for now:  none "
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   end select
  
   call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
        LDT_rc%gridDesc(n,:), LDT_albedo_struc(n)%mxsnoalbfile,&
        LDT_rc%udef, array )


end subroutine read_SACHTET356_mxsnoalb
