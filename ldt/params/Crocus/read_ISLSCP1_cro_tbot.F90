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
! !ROUTINE: read_ISLSCP1_cro_tbot
! \label{read_ISLSCP1_cro_tbot}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  25 Jun 2013: KR Arsenault; Modified to read in raw data file
!  23 Sep 2019: Mahdi Navari, modiffied for Crocus 
!
! !INTERFACE:
subroutine read_ISLSCP1_cro_tbot(n, array)

! !USES:
  use ESMF
  use LDT_coreMod,  only : LDT_rc, LDT_domain
  use LDT_logMod,   only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod
  use Crocus_parmsMod

  implicit none

! !ARGUMENTS: 
  integer,   intent(in) :: n
  real,   intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n))

! !DESCRIPTION:
!  This subroutine retrieves the bottom temperature climatology for the 
!  specified month and returns the values in the latlon projection
!  
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[array]
!   output field with the retrieved greenness fraction
!  \end{description}
!
!EOP      
  integer :: ftn
  integer :: i, c, r, nrec, nn
  logical :: file_exists
  integer :: rec_len, length
  integer, parameter :: rec_len_lat = 360
!  integer, parameter :: rec_len_lon = 180
  integer*2, dimension(360) :: soil_int
  real, dimension(360,180)  :: soilt

!- Used for downscaling/upscaling:
  integer   :: mi                       ! Total number of input param grid array points
  integer   :: mo                       ! Total number of output LIS grid array points
  real      :: param_gridDesc(20)       ! Input parameter grid desc array
  real,    allocatable :: gi1(:)                ! input parameter 1d grid
  logical*1,allocatable:: li1(:)                ! input logical mask (to match gi)
  real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n)) ! output lis 1d grid
  logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n)) ! output logical mask (to match go)

! _____________________________________________________________

!- Set parameter grid array inputs:
   param_gridDesc(1)  = 0.         ! Latlon
   param_gridDesc(2)  = 360.       ! ncols
   param_gridDesc(3)  = 180.       ! nrows
   param_gridDesc(4)  = -89.50     ! LL lat
   param_gridDesc(5)  = -179.50    ! LL lon
   param_gridDesc(7)  = 89.50      ! UR lat
   param_gridDesc(8)  = 179.50     ! UR lon
   param_gridDesc(9)  = 1.0        ! dy
   param_gridDesc(10) = 1.0        ! dx
   param_gridDesc(20) = 64

  inquire(file=trim(Crocus_struc(n)%tbotFile), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) "ISLSCP1 TBOT map ",trim(Crocus_struc(n)%tbotFile)," not found"
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif
  select case ( Crocus_struc(n)%tbot_gridtransform )
    case( "none", "neighbor", "average", "bilinear", "budget-bilinear" )  ! Discrete data type
     write(LDT_logunit,*) "[INFO] Reading ISCLSP1 (1.0 deg res) Native Bottom Temperature file: ",&
           trim(Crocus_struc(n)%tbotfile)
     if( Crocus_struc(n)%tbot_gridtransform == "average" ) then
        write(LDT_logunit,*) "[INFO] This 'Native' bottom temperature map is at 1.0 deg resolution,"
        write(LDT_logunit,*) "     SO Please make sure that when selecting 'average' it is not"
        write(LDT_logunit,*) "     for LIS run domain resolutions less than 1.0 degree."
     endif
   case default
     write(LDT_logunit,*) "[ERR] The ISLSCP1 bottom temperaure field is continuous and "
     write(LDT_logunit,*) "     current spatial transforms supported include:"
     write(LDT_logunit,*) "  == none, neighbor, average, bilinear, budget-bilinear "
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  end select


  ftn = LDT_getNextUnitNumber()

!- Record length for the data
!   ~ Each record has 360 data points, and each 2-byte is
!    used to represent a value for soil T.
  rec_len = rec_len_lat
  length  = rec_len/2
  length  = length*4

!- Open 60-minute deep soil temperature data file: 
  open (ftn, file=trim(Crocus_struc(n)%tbotFile), form="unformatted", status='old', &
        access='direct', convert='big_endian',recl=length )

! begin_record and end_record may be calculated as follows:
! nrec: from north to south (max 1080)
! irec: from west to east   (max 2160)

!- Read bottom soil depth temperature field:
   do nrec = 1, 180    ! y-pts (nr)
      read (ftn, rec=nrec) soil_int
      do nn = 1, 360   ! x-pts (nc)
         soilt(nn,nrec) = float(soil_int(nn))/100.
      end do
   end do

!- Enter spatial downscaling/upscaling options to bring the ISLSCP1
!   1 deg domain to the LIS-run domain ...
   mi = 360*180
   allocate( gi1(mi), li1(mi) )
   gi1 = LDT_rc%udef
   li1 = .false.
   mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
   lo1 = .false.

!- Assign 2-D array to 1-D for aggregation routines:
   i = 0
   do r = 1, 180
      do c = 1, 360;  i = i + 1
         gi1(i) = soilt(c,180-r+1)   ! Flip row-order
       ! Using temperature lower limit for constraint
         if( gi1(i) > 220. )  li1(i) = .true.
      enddo
   enddo

!   Note: The boundary soil temperature data starts at
!      latitude 89.5 and longitude -179.5. The file is in direct access
!      format.  Each direct access record contains data in one latitude
!      circle beginning at -179.5 degree longitude and ending
!      at +179.5 degree. The data are arranged to start at northernmost
!      latitude (north pole) and end at south pole.

!- Transform parameter from original grid to LIS output grid:
   call LDT_transform_paramgrid(n, Crocus_struc(n)%tbot_gridtransform, &
            param_gridDesc, mi, 1, gi1, li1, mo, go1, lo1 )
                                  
!- Convert 1D to 2D grid output arrays:
   i = 0
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n); i = i + 1
         if( go1(i) <= 220.) then   ! Using temperature lower limit for constraint
           array(c,r) = LDT_rc%udef
         else
           array(c,r) = go1(i)
         end if
      enddo
   enddo
   deallocate( li1, gi1 )

  call LDT_releaseUnitNumber(ftn)

end subroutine read_ISLSCP1_cro_tbot

