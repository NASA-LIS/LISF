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
! !ROUTINE: read_Briegleb_albedo
!  \label{read_Briegleb_albedo}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine read_Briegleb_albedo(n, array)
! !USES:
  use LDT_coreMod,        only : LDT_rc
  use LDT_logMod,         only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_albedoMod
  use LDT_timeMgrMod,     only : LDT_calendar
  use LDT_fileIOMod,      only : readLISdata

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n
  real, intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
!
! !DESCRIPTION:
!  This subroutine retrieves the albedo climatology for the 
!  specified month and returns the values in a latlon projection
!  
!  Ref: Briegleb, B.P., P. Minnis, V. Ramanathan, and E. Harrison, 1986: 
!  Comparison of regional clear-sky albedos inferred from satellite 
!  observations and model computations. J. Clim. Appl. Meteor., 25, 214-226
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[q]
!   time index (month or quarter)
!  \item[array]
!   output field with the retrieved albedo
!  \end{description}
!
!EOP      

  integer :: ftn
  integer :: c,r
  logical :: file_exists

  ftn = LDT_getNextUnitNumber()

  inquire(file=trim(LDT_albedo_struc(n)%albFile), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) "Albedo map ",trim(LDT_albedo_struc(n)%albFile)," not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif
  select case ( LDT_albedo_struc(n)%alb_gridtransform )
    case( "none", "neighbor", "average", "bilinear", "budget-bilinear" )
!      write(LDT_logunit,*) "[INFO] Reading NCEP-LIS based Monthly Albedo "
    case default
      write(LDT_logunit,*) "[ERR] The spatial transform option selected for NCEP-LIS "
      write(LDT_logunit,*) "  monthly albedo files is not recognized nor recommended."
      write(LDT_logunit,*) "  Please select: "
      write(LDT_logunit,*) "  ==  none, neighbor, average, bilinear, budget-bilinear. "
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
  end select

  open(ftn, file=trim(LDT_albedo_struc(n)%albFile), &
       access='direct',status='old', &
       form="unformatted", recl=4)
  
  call readLISdata( n, ftn, LDT_albedo_struc(n)%alb_proj, &
       LDT_albedo_struc(n)%alb_gridtransform, &
       LDT_albedo_struc(n)%alb_gridDesc(:), 1, array)    ! 1 indicates 2D layer

  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        if(array(c,r,1).lt.0) then
           array(c,r,1) = LDT_rc%udef
        endif
     enddo
  enddo

  call LDT_releaseUnitNumber(ftn)


end subroutine read_Briegleb_albedo
