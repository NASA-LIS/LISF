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
! !ROUTINE: read_RUC_NCEPslopetype
! \label{read_RUC_NCEPslopetype}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  16 Aug 2013: KR Arsenault, added ability to read in native slopetype file
!
! !INTERFACE:
subroutine read_RUC_NCEPslopetype(n, array)

! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod
  use RUC_ParmsMod

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n
  real, intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n))

! !DESCRIPTION:
!  This subroutine retrieves static, slope type data and reprojects
!  it to the latlon projection. 
!
!  SLOPE CLASS      PERCENT SLOPE
!       1               0-8
!       2               8-30
!       3               > 30
!       4               0-30
!       5               0-8 & > 30
!       6               8-30 & > 30
!       7               0-8, 8-30, > 30
!      13               GLACIAL ICE
!       0               OCEAN/WATER
!
!  REFERENCE:
!   Staub, B. and C. Rosenzweig, 1987: Global digital datasets of soil
!   type, soil texture, surface slope, and other properties: documentation
!   of archived data tape.  NASA Technical Memorandum #100685.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[array]
!   output field with the retrieved slope type data
!  \end{description}
!EOP
  integer   :: ftn
  integer   :: r, c, i
  logical   :: file_exists
  integer, parameter :: incols = 360
  integer, parameter :: inrows = 180
  integer*2 :: readslope(incols, inrows)  !** input

!- Used for downscaling/upscaling:
  integer   :: mi                       ! Total number of input param grid array points
  integer   :: mo                       ! Total number of output LIS grid array points
  real,    allocatable :: gi1(:)        ! input parameter 1d grid
  logical*1,allocatable:: li1(:)        ! input logical mask (to match gi)
  real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n)) ! output lis 1d grid
  logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n)) ! output logical mask (to match go)
  real      :: param_gridDesc(20)       ! Input parameter grid desc array

! _________________________________________

   readslope = LDT_rc%udef
   array = LDT_rc%udef 

!- Set parameter grid array inputs:
   param_gridDesc(1)  = 0.         ! Latlon
   param_gridDesc(2)  = incols     ! ncols
   param_gridDesc(3)  = inrows     ! nrows
   param_gridDesc(4)  = -89.50     ! LL lat
   param_gridDesc(5)  = -179.50    ! LL lon
   param_gridDesc(7)  = 89.50      ! UR lat
   param_gridDesc(8)  = 179.50     ! UR lon
   param_gridDesc(9)  = 1.0        ! dy
   param_gridDesc(10) = 1.0        ! dx
   param_gridDesc(20) = 64

   inquire(file=trim(RUC_struc(n)%slopetypefile), exist=file_exists)
   if(.not.file_exists) then 
      write(LDT_logunit,*) "SLOPETYPE map ",trim(RUC_struc(n)%slopetypefile)," not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
   select case ( RUC_struc(n)%slopetype_gridtransform )
    case( "none", "neighbor", "mode" )  ! Discrete data type
      write(LDT_logunit,*) "[INFO] Reading NCEP Native slopetype file: ",&
            trim(RUC_struc(n)%slopetypefile)
      if( RUC_struc(n)%slopetype_gridtransform == "mode" ) then
         write(LDT_logunit,*) "[INFO] This 'Native' Slope type map is at 1.0 deg resolution,"
         write(LDT_logunit,*) "     SO Please make sure that when selecting 'mode' it is not"
         write(LDT_logunit,*) "     for LIS run domain resolutions less than 1.0 degree."
      endif
   case default
      write(LDT_logunit,*) "[ERR] Since the Slope type field involves discrete data values,"
      write(LDT_logunit,*) "  only 'neighbor' or 'mode' are currently supported spatial"
      write(LDT_logunit,*) "  transform types.  Please check your entries for this parameter."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   end select

   ftn = LDT_getNextUnitNumber()
   open(ftn, file=RUC_struc(n)%slopetypefile, form='formatted', status='old')

!- Read file in:
   do r = 1, inrows
      read(ftn,'(360i2)') (readslope(c, r), c=1, incols)
   end do

!- Enter spatial downscaling/upscaling options to bring the slopetype
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
         gi1(i) = readslope(c,r) 
         if( gi1(i) >= 0. )  li1(i) = .true.  ! Include ocean points
!         if( gi1(i) > 0. )  li1(i) = .true.   ! Exclude ocean points
         if( gi1(i) > 9.  )  gi1(i) = 9.      ! Set upper bound on slope type to be 9.
      enddo
   enddo

!- Transform parameter from original grid to LIS output grid:
   call LDT_transform_paramgrid(n, RUC_struc(n)%slopetype_gridtransform, &
            param_gridDesc, mi, 1, gi1, li1, mo, go1, lo1 )

!- Convert 1D to 2D grid output arrays:
   i = 0
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n); i = i + 1
!         if( go1(i) <= 0.) then   ! Assign Ocean/water values as undefined
!           array(c,r) = LDT_rc%udef
!         else
          array(c,r) = go1(i)
!         end if
      enddo
   enddo
   deallocate( li1, gi1 )

   call LDT_releaseUnitNumber(ftn)

   write(LDT_logunit, *) "[INFO] Done reading native slope type file"

end subroutine read_RUC_NCEPslopetype
