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
! !ROUTINE: read_OzdoganGutman_irrigfrac
!  \label{read_OzdoganGutman_irrigfrac}

! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  16 Nov 2012: K. Arsenault; Modified for irrigation maps
!
! !INTERFACE:
subroutine read_OzdoganGutman_irrigfrac(n,array) 

! !USES:
  use LDT_coreMod,    only : LDT_rc, LDT_domain
  use LDT_logMod,     only : LDT_logunit, LDT_getNextUnitNumber, &
           LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod,  only : readLISdata
  use LDT_irrigationMod
!EOP      
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  real, intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
!
! !DESCRIPTION:
!  This subroutine retrieves the irrigation fraction for
!  each gridcell and returns the values in a latlon projection
!
!  Ref: 
!  Ozdogan, M., and G. Gutman, 2008: A new methodology to map irrigated areas
!   using multi-temporal MODIS and ancillary data: An application example in the 
!   continental U.S.  Remote Sensing of Environment, 112, 3520-3537.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[array]
!   output field MODIS irrigation fraction
!  \end{description}
!
!EOP      
  integer :: ftn
  logical :: file_exists
  integer :: c,r
! __________________________________________________________

!- Set parameter grid array inputs:
   LDT_irrig_struc(n)%irrig_proj           = "latlon"
   LDT_irrig_struc(n)%irrig_gridDesc(1)  = 0.          ! Latlon
   LDT_irrig_struc(n)%irrig_gridDesc(2)  = 464 
   LDT_irrig_struc(n)%irrig_gridDesc(3)  = 224 
   LDT_irrig_struc(n)%irrig_gridDesc(4)  = 25.0625     ! LL lat 
   LDT_irrig_struc(n)%irrig_gridDesc(5)  = -124.9375   ! LL lon 
   LDT_irrig_struc(n)%irrig_gridDesc(6)  = 128
   LDT_irrig_struc(n)%irrig_gridDesc(7)  = 52.9375     ! UR lat
   LDT_irrig_struc(n)%irrig_gridDesc(8)  = -67.0625    ! UR lon
   LDT_irrig_struc(n)%irrig_gridDesc(9)  = 0.125    
   LDT_irrig_struc(n)%irrig_gridDesc(10) = 0.125    
   LDT_irrig_struc(n)%irrig_gridDesc(20) = 64

  inquire(file=trim(LDT_irrig_struc(n)%irrigfracfile), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) "Irrigation fraction map ",&
                           trim(LDT_irrig_struc(n)%irrigfracfile)," not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif
  select case ( LDT_irrig_struc(n)%irrigfrac_gridtransform )
   case( "none", "neighbor", "average" )  ! continuous data type
     write(LDT_logunit,*) "[INFO] Reading MODIS (e.g., CONUS) irrigation fraction file: ",&
           trim(LDT_irrig_struc(n)%irrigfracfile)
   case default
     write(LDT_logunit,*) "[ERR] Irrigation fraction maps can be modified currently with "
     write(LDT_logunit,*) "     'neighbor' or 'average' spatial transform types."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  end select
  
  ftn = LDT_getNextUnitNumber()
  open(ftn, file=LDT_irrig_struc(n)%irrigfracfile,access='direct',status='old', &
       form="unformatted", recl=4)

  call readLISdata( n, ftn, LDT_irrig_struc(n)%irrig_proj, &
       LDT_irrig_struc(n)%irrigfrac_gridtransform, &
       LDT_irrig_struc(n)%irrig_gridDesc(:), 1, array)   ! 1 indicates 2D layer
  
  do r = 1,LDT_rc%lnr(n)
     do c = 1,LDT_rc%lnc(n)
        if( array(c,r,1).lt.0 ) then
           array(c,r,1) = LDT_rc%udef
        endif
     enddo
  enddo

  call LDT_releaseUnitNumber(ftn)

end subroutine read_OzdoganGutman_irrigfrac
