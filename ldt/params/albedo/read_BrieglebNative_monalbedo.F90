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
! !ROUTINE: read_BrieglebNative_monalbedo
!  \label{read_BrieglebNative_monalbedo}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  20 Aug 2013: KR Arsenault; Modified to read in raw data file
!
! !INTERFACE:
subroutine read_BrieglebNative_monalbedo(n, array)

! !USES:
  use LDT_coreMod,   only : LDT_rc
  use LDT_logMod,    only : LDT_logunit, LDT_getNextUnitNumber, &
          LDT_releaseUnitNumber, LDT_endrun
  use LDT_albedoMod
  use LDT_fileIOMod

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n
  real, intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n))
!
! !DESCRIPTION:
!  This subroutine retrieves the albedo climatology for the 
!  specified month and returns the values in a latlon projection
!  
!  albedo values represent:
!      0    : ocean/water flag 
!      1-69 : surface albedo (percent)
!      70   : glacial/perpetual snow (percent, with gfrac=1)
!
!  Albedo:  the original monthly albedo files covered only 55S to 75N, so
!   the global greenness fraction (including the filled in Antarctic) was
!   used to fill in the regions south of 55S and north of 75N as follows:
!   ---------  ------
!   greenness  albedo
!   ---------  ------
!   1 (land*)  70 (glacial)
!   0 (sea)     0 (sea)
!  *south of 55S and north of 75N, the only non-zero greenness values=1
!
!  Ref: Briegleb, B.P., P. Minnis, V. Ramanathan, and E. Harrison, 1986: 
!  Comparison of regional clear-sky albedos inferred from satellite 
!  observations and model computations. J. Clim. Appl. Meteor., 25, 214-226
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[array]
!   output field with the retrieved albedo
!  \end{description}
!
!EOP      
  integer :: ftn
  integer :: c,r,i
  logical :: file_exists
  integer, parameter :: incols = 2500
  integer, parameter :: inrows = 1250
  integer :: readalb(incols, inrows)  !** input

!- Used for downscaling/upscaling:
  integer   :: mi                       ! Total number of input param grid array points
  integer   :: mo                       ! Total number of output LIS grid array points
  real,    allocatable :: gi1(:)                ! input parameter 1d grid
  logical*1,allocatable:: li1(:)                ! input logical mask (to match gi)
  real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n)) ! output lis 1d grid
  logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n)) ! output logical mask (to match go)
  real      :: param_gridDesc(20)       ! Input parameter grid desc array
! __________________________________________________________________

   readalb = LDT_rc%udef
   array = LDT_rc%udef

!- Set parameter grid array inputs:
   param_gridDesc(1)  = 0.         ! Latlon
   param_gridDesc(2)  = incols     ! ncols
   param_gridDesc(3)  = inrows     ! nrows
   param_gridDesc(4)  = -89.928   ! LL lat
   param_gridDesc(5)  = -180.00    ! LL lon
   param_gridDesc(6)  = 128
   param_gridDesc(7)  = 89.928     ! UR lat
   param_gridDesc(8)  = 179.856    ! UR lon
   param_gridDesc(9)  = 0.144      ! dy
   param_gridDesc(10) = 0.144      ! dx
   param_gridDesc(20) = 64

   inquire(file=trim(LDT_albedo_struc(n)%albFile), exist=file_exists)
   if(.not.file_exists) then 
      write(LDT_logunit,*) "[ERR] NCEP Native Monthly Albedo map ",&
                           trim(LDT_albedo_struc(n)%albFile)," not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif

   select case ( LDT_albedo_struc(n)%alb_gridtransform )
     case( "none", "neighbor", "average", "bilinear", "budget-bilinear" )
       write(LDT_logunit,*) "[INFO] Reading NCEP Native Monthly Albedo "
     case default
       write(LDT_logunit,*) "[ERR] The spatial transform option selected for NCEP 'Native'"
       write(LDT_logunit,*) "   monthly albedo files (@0.144 deg resolution) is not recognized"
       write(LDT_logunit,*) "   nor recommended.  Please select: "
       write(LDT_logunit,*) "  ==  none, neighbor, average, bilinear, budget-bilinear. "
       write(LDT_logunit,*) "Program stopping ..."
       call LDT_endrun
   end select

   ftn = LDT_getNextUnitNumber()
   open(ftn, file=trim(LDT_albedo_struc(n)%albFile),form="formatted",status='old' )

!- Read file in:
   read(ftn,'(2500(i2))')((readalb(c,r),c=1,incols),r=1,inrows)

!- Enter spatial downscaling/upscaling options to bring the snow-free albedo
!  1 deg domain to the LIS-run domain ...
    mi = incols*inrows
    allocate( gi1(mi), li1(mi) )
    gi1 = LDT_rc%udef
    li1 = .false.
    mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
    lo1 = .false.

!- Assign 2-D array to 1-D for aggregation routines:
   i = 0
   do r = 1, inrows
      do c = 1, incols;  i = i + 1
         gi1(i) = readalb(c,r)
!         if( gi1(i) >= 0. )  li1(i) = .true.  ! Include ocean points
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
            array(c,r) = LDT_rc%udef
         else
            array(c,r) = go1(i)/100.
         end if
      enddo
   enddo
   deallocate( li1, gi1 )

  call LDT_releaseUnitNumber(ftn)

end subroutine read_BrieglebNative_monalbedo
