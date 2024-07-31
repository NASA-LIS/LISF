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
! !ROUTINE: read_BrieglebNative_qtralbedo
!  \label{read_BrieglebNative_qtralbedo}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  20 Aug 2013: KR Arsenault; Modified to read in raw data file
!
! !INTERFACE:
subroutine read_BrieglebNative_qtralbedo(n, array)

! !USES:
  use LDT_coreMod, only : LDT_rc
  use LDT_logMod,  only : LDT_logunit, LDT_getNextUnitNumber, &
           LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod
  use LDT_albedoMod

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n
  real,    intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),4)

!
! !DESCRIPTION:
!  This subroutine retrieves the albedo climatology for the 
!  specified quarter and returns the values in a latlon projection
!  
! - ALBEDO: Quarterly global albedo (snow-free) data at 1 deg.
!   !NOTE: ALBEDO is quarterly and valid at end of Jan, Apr,
!          Jul, and Oct. 
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
!   time index ( quarter)
!  \item[array]
!   output field with the retrieved albedo
!  \end{description}
!
!EOP      
  integer :: ftn
  integer :: c,r,i,k
  integer :: qtr
  logical :: file_exists
  integer, parameter :: incols = 361
  integer, parameter :: inrows = 180
  real    :: readalb(incols, inrows, 4)  !** input

!- Used for downscaling/upscaling:
  integer   :: mi                       ! Total number of input param grid array points
  integer   :: mo                       ! Total number of output LIS grid array points
  real,    allocatable :: gi1(:)                ! input parameter 1d grid
  logical*1,allocatable:: li1(:)                ! input logical mask (to match gi)
  real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n)) ! output lis 1d grid
  logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n)) ! output logical mask (to match go)
  real      :: param_gridDesc(20)       ! Input parameter grid desc array
! ______________________________________________

   readalb = LDT_rc%udef
   array = LDT_rc%udef

!- Set parameter grid array inputs:
   param_gridDesc(1)  = 0.         ! Latlon
   param_gridDesc(2)  = incols     ! ncols
   param_gridDesc(3)  = inrows     ! nrows
   param_gridDesc(4)  = -89.50     ! LL lat
   param_gridDesc(5)  = -180.0     ! LL lon
   param_gridDesc(6)  = 128
   param_gridDesc(7)  = 89.50      ! UR lat
   param_gridDesc(8)  = 180.0      ! UR lon
   param_gridDesc(9)  = 1.0        ! dy
   param_gridDesc(10) = 1.0        ! dx
   param_gridDesc(20) = 64

   inquire(file=trim(LDT_albedo_struc(n)%albFile), exist=file_exists)
   if(.not.file_exists) then 
      write(LDT_logunit,*) 'Albedo map ',trim(LDT_albedo_struc(n)%albFile),' not found'
      write(LDT_logunit,*) 'Program stopping ...'
      call LDT_endrun
   endif

   select case ( LDT_albedo_struc(n)%alb_gridtransform )
     case( "none", "neighbor", "average", "bilinear", "budget-bilinear" )
       write(LDT_logunit,*) "[INFO] Reading NCEP Native Quarterly Albedo "
       if( LDT_albedo_struc(n)%alb_gridtransform == "average" ) then
          write(LDT_logunit,*) "[INFO] This 'Native' albedo map is at 1.0 deg resolution,"
          write(LDT_logunit,*) "     SO Please make sure that when selecting 'average' it is not"
          write(LDT_logunit,*) "     for LIS run domain resolutions less than 1.0 degree."
       endif
     case default
       write(LDT_logunit,*) "[ERR] The spatial transform option selected for NCEP 'Native'"
       write(LDT_logunit,*) "  quarterly albedo files (@1.0 deg resolution) is not recognized"
       write(LDT_logunit,*) "  or not recommended.  Please select: "
       write(LDT_logunit,*) "  ==  none, neighbor, average, bilinear, budget-bilinear. "
       write(LDT_logunit,*) "Program stopping ..."
       call LDT_endrun
   end select

   write(LDT_logunit,*)"[INFO] Reading NCEP Native albedo quarterly file."
   ftn = LDT_getNextUnitNumber()
   open(ftn, file=trim(LDT_albedo_struc(n)%albFile),form='formatted', status='old')

!- Read Albedo file (at 1.0 deg resolution) for 3-month quarters:
   do qtr=1,4
      read(ftn,'(361(f6.2))') ((readalb(c,r,qtr),c=1,361),r=1,180)
   end do

!- Enter spatial downscaling/upscaling options to bring the max snow alb
!  1 deg domain to the LIS-run domain ...
   mi = incols*inrows
   allocate( gi1(mi), li1(mi) )
   mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)

!- Loop over each quarterly season:
   do qtr = 1, 4

  !- Assign 2-D array to 1-D for aggregation routines:
     gi1 = LDT_rc%udef
     li1 = .false.
     lo1 = .false.
     i = 0
     do r = 1, inrows
        do c = 1, incols;  i = i + 1
           gi1(i) = readalb(c,r,qtr)
           if( gi1(i) > 0. )  li1(i) = .true.   ! Exclude ocean points
        enddo
     enddo

  !- Transform parameter from original grid to LIS output grid:
     call LDT_transform_paramgrid(n, LDT_albedo_struc(n)%alb_gridtransform, &
              param_gridDesc, mi, 1, gi1, li1, mo, go1, lo1 )


  !- Convert 1D to 2D grid output arrays:
     i = 0
     do r = 1, LDT_rc%lnr(n)
        do c = 1, LDT_rc%lnc(n); i = i + 1
           if( go1(i) <= 0.) then   ! Assign Ocean/water values as undefined
             array(c,r,qtr) = LDT_rc%udef
           else
             array(c,r,qtr) = go1(i)/100.
           end if
        enddo
     enddo
 
   end do  ! End seasonal quarter loop
   deallocate( li1, gi1 )

   call LDT_releaseUnitNumber(ftn)


end subroutine read_BrieglebNative_qtralbedo
