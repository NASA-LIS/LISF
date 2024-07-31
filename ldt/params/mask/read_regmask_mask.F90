!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: read_regmask_mask
! \label{read_regmask_mask}
!
! !REVISION HISTORY:
!  03 Sep  2012: KR Arsenault; Read in Regional Mask File
!  23 Mar  2016: KR Arsenault; Updated to support more GIS-original files
!  13 Aug  2018: David M. Mocko; read in a regular 1dg4r mask file
!
! !INTERFACE:
subroutine read_regmask_mask( n, regmask )

! !USES:
  use LDT_coreMod,    only : LDT_rc
  use LDT_logMod,     only : LDT_logunit, LDT_getNextUnitNumber, &
                             LDT_releaseUnitNumber, LDT_endrun
  use LDT_gridmappingMod
  use LDT_fileIOMod

  implicit none
! !ARGUMENTS:
  integer, intent(in)  :: n          ! nest index
  real,    intent(out) :: regmask(LDT_rc%lnc(n),LDT_rc%lnr(n),1)

! !DESCRIPTION:
!  This subroutine retrieves regional mask information for clipping
!   out designated area.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
 !   index of the nest
!  \item[regmask]
 !   output field with regional mask information
!  \end{description}
!EOP
   integer :: ftn
   integer :: i, t, c, r
   integer :: line, pnc, pnr
   logical :: file_exists

   integer :: mi                         ! Total number of input param grid array points
   integer :: mo                         ! Total number of output LIS grid array points
   real    :: param_gridDesc(20)         ! Input parameter grid desc array
   real    :: subparam_gridDesc(20)      ! Subsetted parameter grid desc array
   real,    allocatable  :: gi(:)        ! Input parameter 1d grid
   logical*1,allocatable :: li(:)        ! Input logical mask (to match gi)
   real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))            ! Output lis 1d grid
   logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))            ! Output logical mask (to match go)

   real, allocatable :: read_maskmap(:,:)  ! Read input parameter

!________________________________________________________________________

   regmask = 0.

!- Check file and open:
   inquire(file=trim(LDT_rc%regfile(n)), exist=file_exists)
   if(.not.file_exists) then
      write(LDT_logunit,*) "Regional mask map ",trim(LDT_rc%regfile(n))," not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif

!- Grab pnc, pnr values
   pnc = LDT_rc%reg_gridDesc(n,2)
   pnr = LDT_rc%reg_gridDesc(n,3)

   write(LDT_logunit,*) "[INFO] Reading Regional Mask File: ",&
        trim(LDT_rc%regfile(n))
   ftn = LDT_getNextUnitNumber()
   open(ftn,file=LDT_rc%regfile(n),form='unformatted',status='old',&
            access='direct',recl=4)

   allocate( read_maskmap(pnc, pnr) )
   read_maskmap = LDT_rc%udef
   i = 0
   do r = 1, pnr
      do c = 1, pnc;  i = i + 1
         read(ftn, rec=i) read_maskmap(c,r)
      enddo
   enddo

!- Read and aggregate data:
   select case ( LDT_rc%reg_gridtransform(n) )

    case( "none", "mode", "neighbor" )

    !- Transform parameter from original grid to LIS output grid:
       mi = LDT_rc%reg_gridDesc(n,2) * LDT_rc%reg_gridDesc(n,3)
       mo = LDT_rc%lnc(n) * LDT_rc%lnr(n)
       allocate( gi(mi), li(mi) )
       gi = LDT_rc%udef;  go1 = LDT_rc%udef
       li = .false.;  lo1 = .false.

    !- Assign 2-D array to 1-D for aggregation routines:
       i = 0
       do r = 1, pnr
          do c = 1, pnc;  i = i + 1
             gi(i) = read_maskmap(c,r)
             if( gi(i) .ne. LDT_rc%udef ) li(i) = .true.
          enddo
       enddo
       deallocate( read_maskmap )

       call LDT_transform_paramgrid(n, LDT_rc%reg_gridtransform(n), &
                LDT_rc%reg_gridDesc(n,:), mi, 1, gi, li, mo, go1, lo1 )

    !- Convert 1D map to 2D grid arrays:
       i = 0
       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n)
             i = i + 1
             regmask(c,r,1) = go1(i)
          enddo
       enddo

    case default
      write(LDT_logunit,*) &
      "[ERR] This spatial transform for regional mask not available ..."
      write(LDT_logunit,*)" Program stopping ..."
      call LDT_endrun
   end select
   deallocate( gi, li )

   call LDT_releaseUnitNumber(ftn)

   write(LDT_logunit,*) "[INFO] Done reading regional mask file"

 end subroutine read_regmask_mask

