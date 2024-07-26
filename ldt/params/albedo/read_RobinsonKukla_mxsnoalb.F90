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
! !ROUTINE: read_RobinsonKukla_mxsnoalb
!  \label{read_RobinsonKukla_mxsnoalb}

! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine read_RobinsonKukla_mxsnoalb(n,array) 

! !USES:
  use LDT_coreMod,       only : LDT_rc, LDT_domain
  use LDT_logMod,        only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod,     only : readLISdata 
  use LDT_albedoMod
!EOP      
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  real, intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
!
! !DESCRIPTION:
!  This subroutine retrieves the maximum snow expected over 
!  deep snow and returns the values in a latlon projection
!
!  Ref: D. A. Robinson and G. Kukla, Maximum surface albedo of 
!  seasonally snow-covered lands in the Northern Hemisphere, 
!  Journal of Climate and Applied Meteorology, vol. 24, pp. 402-411, 1985.
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
  real     :: alb_gridDesc(20)
  logical  :: file_exists
  integer  :: c,r

  ftn = LDT_getNextUnitNumber()

  inquire(file=trim(LDT_albedo_struc(n)%mxsnoalbfile), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) "Max Snow Albedo map ",trim(LDT_albedo_struc(n)%mxsnoalbfile)," not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif
  select case ( LDT_albedo_struc(n)%mxsnoalb_gridtransform )
    case( "none", "neighbor", "average", "bilinear", "budget-bilinear" )
      write(LDT_logunit,*) "[INFO] Reading NCEP-LIS Max Snow Albedo "
    case default
      write(LDT_logunit,*) "[ERR] The spatial transform option selected for NCEP-LIS"
      write(LDT_logunit,*) "   max snow albedo file is not recognized nor recommended."
      write(LDT_logunit,*) "   Please select: "
      write(LDT_logunit,*) "  ==  none, neighbor, average, bilinear, budget-bilinear. "
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
  end select
  
  open(ftn, file=LDT_albedo_struc(n)%mxsnoalbfile,&
       access='direct',status='old', &
       form="unformatted", recl=4)

  call readLISdata( n, ftn, LDT_albedo_struc(n)%mxsnoalb_proj, &
       LDT_albedo_struc(n)%mxsnoalb_gridtransform, &
       LDT_albedo_struc(n)%mxsnoalb_gridDesc, 1, array)   ! 1 indicates 2D layer

  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        if(array(c,r,1).lt.0) then
           array(c,r,1) = LDT_rc%udef
        endif
     enddo
  enddo
  call LDT_releaseUnitNumber(ftn)

end subroutine read_RobinsonKukla_mxsnoalb
