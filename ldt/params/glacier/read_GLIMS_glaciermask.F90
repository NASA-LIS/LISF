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
! !ROUTINE: read_GLIMS_glaciermask
!  \label{read_GLIMS_glaciermask}
!
! !REVISION HISTORY:
!  30 Mar 2018: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine read_GLIMS_glaciermask(n, localmask )

! !USES:
  use LDT_coreMod,   only : LDT_rc, LDT_localPet
  use LDT_gridmappingMod
  use LDT_logMod
  use LDT_glacierMod
  use LDT_fileIOMod

#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  use netcdf
#endif

  implicit none

! !ARGUMENTS: 
  integer, intent(in)  :: n
  real,    intent(out) :: localmask(LDT_rc%lnc(n),LDT_rc%lnr(n),1)

! !DESCRIPTION:
!  This subroutine reads the landmask data and returns the 
!   mask and surface type arrays.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of nest
!   \item[localmask]
!    landmask for the region of interest
!   \end{description}
!
!EOP      
  integer, parameter :: nr = 18000
  integer :: ftn, ios1,maskid
  logical :: file_exists
  integer :: c, r, t, i, line
  integer :: glpnc, glpnr             ! Parameter (global) total columns and rows
  integer :: subpnc, subpnr           ! Parameter subsetted columns and rows
  real    :: subparam_gridDesc(20)    ! Input parameter grid desc array
  real    :: cnt_0mask_1lc
  real    :: cnt_1mask_0lc
  integer, allocatable :: lat_line(:,:)
  integer, allocatable :: lon_line(:,:)
  real,    allocatable :: read_inputparm(:,:)  ! Read input parameter
  integer :: mi                       ! Total number of input param grid array points
  integer :: mo                       ! Total number of output LIS grid array points
  real,    allocatable  :: gi1(:)      ! Input parameter 1d grid
  logical*1,allocatable :: li1(:)      ! Input logical mask (to match gi)

  real                 :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! Output lis 1d grid
  logical*1            :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! Output logical mask (to match go)
  real      :: param_gridDesc(20)       ! Input parameter grid desc array

!- Check for and open landmask file:
   inquire(file=trim(LDT_rc%glaciermask(n)), exist=file_exists)
   if( file_exists ) then 
      write(LDT_logunit,*)"[INFO] GLIMS mask -- Reading ",trim(LDT_rc%glaciermask(n))
   else
      write(LDT_logunit,*) "[ERR] GLIMS map: ",trim(LDT_rc%glaciermask(n))," does not exist."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif

! -------------------------------------------------------------------
!    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
! -------------------------------------------------------------------
!- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
!- Using Landcover grid description as default for now:
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, LDT_glacier_struc(n)%mask_proj, &
        LDT_glacier_struc(n)%mask_gridDesc(:), &
        glpnc, glpnr, subpnc, subpnr,  &
        subparam_gridDesc, lat_line, lon_line )
   
   allocate( read_inputparm(subpnc, subpnr) )
   read_inputparm = 0.
   
   ! -------------------------------------------------------------------

#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
   call LDT_verify(nf90_open(path=LDT_rc%glaciermask(n),&
        mode=nf90_nowrite,ncid=ftn),&
        'nf90_open file failed in read_GLIMS_glaciermask')
   
   call LDT_verify(nf90_inq_varid(ftn,"glaciermask",maskid),&
        'nf90_inq_varid failed in read_GLIMS_glaciermask')

   call LDT_verify(nf90_get_var(ftn,maskid,read_inputparm,&
        start=(/lon_line(1,1),nr-(lat_line(1,1)+subpnr)+1/),&
        count=(/subpnc,subpnr/)),&
        'nf90_get_var failed in read_GLIMS_glaciermask')
   call LDT_verify(nf90_close(ftn))

#endif


!- Enter spatial downscaling/upscaling options to bring the data
!  to the target domain ...
    mi = subpnc*subpnr
    allocate( gi1(mi), li1(mi) )
    gi1 = LDT_rc%udef
    li1 = .false.
    mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
    lo1 = .false.

!- Assign 2-D array to 1-D for aggregation routines:
   do r = 1, subpnr
      do c = 1, subpnc
         i = c+(r-1)*subpnc
         gi1(i) = read_inputparm(c,r)
         if( gi1(i) .ge. 0. )  li1(i) = .true. 
      enddo
   enddo
   param_gridDesc = 0

!- Set parameter grid array inputs:
   param_gridDesc(1)  = 0.         ! Latlon
   param_gridDesc(2)  = subpnc
   param_gridDesc(3)  = subpnr
   param_gridDesc(4)  = -89.995 + lat_line(1,1)*0.01   ! LL lat
   param_gridDesc(5)  = -179.995 + lon_line(1,1)*0.01    ! LL lon
   param_gridDesc(6)  = 128
   param_gridDesc(7)  = -89.995 + (lat_line(1,1) + subpnc)*0.01   ! UR lat
   param_gridDesc(8)  = 179.856 + (lon_line(1,1) + subpnr)*0.01   ! UR lon
   param_gridDesc(9)  = 0.01     ! dy
   param_gridDesc(10) = 0.01      ! dx
   param_gridDesc(20) = 64

   deallocate( lat_line, lon_line )

!- Transform parameter from original grid to LIS output grid:
   call LDT_transform_paramgrid(n, LDT_glacier_struc(n)%mask_gridtransform, &
        param_gridDesc, mi, 1, gi1, li1, mo, go1, lo1 )

!- Convert 1D to 2D grid output arrays:
   i = 0
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n); i = i + 1
         if( go1(i) <=  LDT_rc%gridcell_glacier_frac(n)) then   
           localmask(c,LDT_rc%lnr(n)-r+1,1) = LDT_rc%udef
        else
           localmask(c,LDT_rc%lnr(n)-r+1,1) = 1.0
        end if
     enddo
  enddo

  deallocate( li1, gi1 )

end subroutine read_GLIMS_glaciermask
