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
! !ROUTINE: read_STATSGOv1_soilfractions
! \label{read_STATSGOv1_soilfractions}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  03 Aug  2012: KR Arsenault; Expanded to all soil fractions and tiling
!
! !INTERFACE:
subroutine read_STATSGOv1_soilfractions(n, num_bins, soilsfgrd, &
                       sandave, clayave, siltave )
! !USES:
  use LDT_coreMod,     only : LDT_rc
  use LDT_logMod,      only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod
  use LDT_gridmappingMod
  use LDT_paramTileInputMod, only: param_2dbin_areacalc

  implicit none
! !ARGUMENTS: 
  integer, intent(in)   :: n          ! nest index
  integer, intent(in)   :: num_bins   ! number of bins for tiling
  real,    intent(out)  :: soilsfgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real,    intent(out)  :: sandave(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real,    intent(out)  :: clayave(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real,    intent(out)  :: siltave(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)

! !DESCRIPTION:
!  This subroutine retrieves STATSGO v1 sand, silt, and clay fraction data 
!   and reprojects it to the latlon projection. 
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[num_bins]
!   number of bins (or types/classes)
!  \item[soilsfgrd]
!   output field with merged soils fraction of grid
!  \item[sandave]
!   output field with average sand binned fraction 
!  \item[clayave]
!   output field with average clay binned fraction 
!  \item[siltave]
!   output field with average silt binned fraction 
!  \end{description}
!EOP
   integer :: ftn_sa, ftn_cl, ftn_si
   logical :: file_exists
   integer :: i, t, c, r, line
   integer :: glpnc, glpnr             ! Parameter (global) total columns and rows
   integer :: subpnc, subpnr           ! Parameter subsetted columns and rows
   integer :: mi                       ! Total number of input param grid array points
   integer :: mo                       ! Total number of output LIS grid array points
   real    :: subparam_gridDesc(20)    ! Input parameter grid desc array
   integer, allocatable :: lat_line(:,:), lon_line(:,:)
   real,    allocatable :: read_sand_sub(:,:)    ! Read input parameter
   real,    allocatable :: read_clay_sub(:,:)    ! Read input parameter
   real,    allocatable :: read_silt_sub(:,:)    ! Read input parameter
   integer, allocatable :: n11(:)                ! Maps each input grid point to output grid.
   real,    allocatable :: gi1(:), gi2(:), gi3(:) ! input parameter 1d grid
   logical*1,allocatable:: li1(:), li2(:), li3(:) ! input logical mask (to match gi)
   real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! output lis 1d grid
   real      :: go2(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! output lis 1d grid
   real      :: go3(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! output lis 1d grid
   logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! output logical mask (to match go)
   logical*1 :: lo2(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! output logical mask (to match go)
   logical*1 :: lo3(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! output logical mask (to match go)
!________________________________________________________________________

  sandave = 0.; clayave = 0.; siltave = 0.
  soilsfgrd = 0.

   write(LDT_logunit,*) '[INFO] Reading STATSGO v1 sand,clay and silt files: ',&
        trim(LDT_rc%safile(n)),", ",trim(LDT_rc%clfile(n)),",",trim(LDT_rc%sifile(n)) 

   call LDT_checkDomainExtents(n,LDT_rc%soil_gridDesc(n,:))

! -------------------------------------------------------------------
!   CONFIRM AND OPEN SOIL FRACTION FILES
! -------------------------------------------------------------------

!- Sand file:
   inquire(file=trim(LDT_rc%safile(n)), exist=file_exists)
   if(.not.file_exists) then 
      write(LDT_logunit,*) 'Sand map ',trim(LDT_rc%safile(n)),' not found'
      write(LDT_logunit,*) 'Program stopping ...'
      call LDT_endrun
   endif
   ftn_sa = LDT_getNextUnitNumber()
   open(ftn_sa,file=LDT_rc%safile(n),form='unformatted',status='old',&
               access='direct',recl=4)

!- Clay file:
   inquire(file=trim(LDT_rc%clfile(n)), exist=file_exists)
   if(.not.file_exists) then
      write(LDT_logunit,*) 'Clay map ',trim(LDT_rc%clfile(n)),' not found'
      write(LDT_logunit,*) 'Program stopping ...'
      call LDT_endrun
   endif
   ftn_cl = LDT_getNextUnitNumber()
   open(ftn_cl,file=LDT_rc%clfile(n),form='unformatted',status='old',&
               access='direct',recl=4)

!- Silt file (optional):
   inquire(file=trim(LDT_rc%sifile(n)), exist=file_exists)
   if(.not.file_exists) then
      write(LDT_logunit,*) 'Silt map ',trim(LDT_rc%sifile(n)),' not found'
      write(LDT_logunit,*) 'Program stopping ...'
      call LDT_endrun
   endif
   ftn_si = LDT_getNextUnitNumber()
   open(ftn_si,file=LDT_rc%sifile(n),form='unformatted',status='old',&
               access='direct',recl=4)

! -------------------------------------------------------------------
!   PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
! -------------------------------------------------------------------
!- Map Parameter Grid Info to LIS Target Grid/Projection Info --
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, LDT_rc%soils_proj, LDT_rc%soil_gridDesc(n,:), &
                  glpnc, glpnr, subpnc, subpnr,  &
                  subparam_gridDesc, lat_line, lon_line )

!- Initialize parameter read-in array:
   allocate( read_sand_sub(subpnc,subpnr) )
   allocate( read_clay_sub(subpnc,subpnr) )
   allocate( read_silt_sub(subpnc,subpnr) )
   read_sand_sub = LDT_rc%udef
   read_clay_sub = LDT_rc%udef
   read_silt_sub = LDT_rc%udef
! -------------------------------------------------------------------
!    READ IN SOIL PARAMETER FIELDS (NON-TILED/TILED OPTIONS)
! -------------------------------------------------------------------
   line = 0
   do r = 1, subpnr
      do c = 1, subpnc
         line = (lat_line(c,r)-1)*glpnc + lon_line(c,r)
         read(ftn_sa,rec=line) read_sand_sub(c,r)
         read(ftn_cl,rec=line) read_clay_sub(c,r)
         read(ftn_si,rec=line) read_silt_sub(c,r)
      enddo
   enddo

! -------------------------------------------------------------------
!    AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
! -------------------------------------------------------------------
    if( LDT_rc%soils_gridtransform(n) == "average" .or.  &
        LDT_rc%soils_gridtransform(n) == "tile" ) then
        mi = subpnc*subpnr
        mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
        allocate( gi1(mi), gi2(mi), gi3(mi) )
        allocate( li1(mi), li2(mi), li3(mi) )
        allocate( n11(mi) )
        gi1 = LDT_rc%udef; gi2 = LDT_rc%udef; gi3 = LDT_rc%udef
        li1 = .false.; li2 = .false.;  li3 = .false.
        lo1 = .false.; lo2 = .false.;  lo3 = .false.
     !- Assign 2-D array to 1-D for aggregation routines:
        i = 0
        do r = 1, subpnr
           do c = 1, subpnc;  i = i + 1
              gi1(i) = read_sand_sub(c,r)
              gi2(i) = read_clay_sub(c,r)
              gi3(i) = read_silt_sub(c,r)
              if( gi1(i) .ne. LDT_rc%udef )  li1(i) = .true.
              if( gi2(i) .ne. LDT_rc%udef )  li2(i) = .true.
              if( gi3(i) .ne. LDT_rc%udef )  li3(i) = .true.
           enddo
        enddo

     !- Create mapping between parameter domain and LIS grid domain:
        call upscaleByAveraging_input( subparam_gridDesc, LDT_rc%gridDesc(n,:), &
                                       mi, mo, n11 )
    end if ! End transform check

 ! ----------------------------------------------------------

 !- Aggregate (if needed) and assign to final output fields:
    select case ( LDT_rc%soils_gridtransform(n) )

   !- Files of SAME resolution:
      case( "none" ) 
         write(LDT_logunit,*) " No aggregation applied for soil fraction files ... "
         do r=1,LDT_rc%lnr(n)
            do c=1,LDT_rc%lnc(n)
               if(read_sand_sub(c,r) < 0.) read_sand_sub(c,r) = LDT_rc%udef
               if(read_clay_sub(c,r) < 0.) read_clay_sub(c,r) = LDT_rc%udef
               if(read_silt_sub(c,r) < 0.) read_silt_sub(c,r) = LDT_rc%udef
               do t = 1, num_bins
               !- Single soil fraction layer, write to first bin:
                  if( t == 1 ) then
                    sandave(c,r,t) = read_sand_sub(c,r)
                    clayave(c,r,t) = read_clay_sub(c,r)
                    siltave(c,r,t) = read_silt_sub(c,r)
                  endif
               end do
            enddo
         enddo
         soilsfgrd(:,:,:) = 1.0

   !- Average finer scale points to output grid:
      case( "average" )
         call upscaleByAveraging( mi, mo, LDT_rc%udef, n11,  &
                                  li1, gi1, lo1, go1 )
         call upscaleByAveraging( mi, mo, LDT_rc%udef, n11,  &
                                  li2, gi2, lo2, go2 )
         call upscaleByAveraging( mi, mo, LDT_rc%udef, n11,  &
                                  li3, gi3, lo3, go3 )
      !- Convert 1D to 3D grid output arrays:
         i = 0
         do r = 1, LDT_rc%lnr(n)
            do c = 1, LDT_rc%lnc(n)
               i = i + 1
               if( go1(i) < 0.) go1(i) = LDT_rc%udef
               if( go2(i) < 0.) go2(i) = LDT_rc%udef
               if( go3(i) < 0.) go3(i) = LDT_rc%udef
               do t = 1, num_bins
                  if( t == 1 ) then
                    sandave(c,r,t) = go1(i)
                    clayave(c,r,t) = go2(i)
                    siltave(c,r,t) = go3(i)
                  endif
               end do
            enddo
         enddo
         soilsfgrd(:,:,:) = 1.0

   !- Tiling: Two-field combination for bin-based grid fraction and sand,clay estimates:
      case( "tile" )

         call param_2dbin_areacalc( n, num_bins, mi, mo, n11, -0.99, gi1, gi2, &
                                    sandave, clayave, soilsfgrd )
         do r=1,LDT_rc%lnr(n)
            do c=1,LDT_rc%lnc(n)
               do t = 1, num_bins
               !- Calculate silt from sand + clay:
                  if( sandave(c,r,t) /= LDT_rc%udef .and. &
                     clayave(c,r,t) /= LDT_rc%udef ) then
                     siltave(c,r,t) = 1.0 - (sandave(c,r,t)+clayave(c,r,t))
                  else
                     siltave(c,r,t) = LDT_rc%udef
                  endif
               end do
            enddo
         enddo

      case default
         write(LDT_logunit,*) ' This spatial transform for soil fraction not available at this time ... '
         write(LDT_logunit,*) 'Program stopping ...'
         call LDT_endrun

    end select
! _____________________________________________

   deallocate( read_sand_sub, read_clay_sub, read_silt_sub )
   deallocate( li1, li2, li3 )
   deallocate( gi1, gi2, gi3 )
   deallocate( lat_line, lon_line, n11 )

   call LDT_releaseUnitNumber(ftn_sa)
   call LDT_releaseUnitNumber(ftn_cl)
   call LDT_releaseUnitNumber(ftn_si)

   write(LDT_logunit,*) '[INFO] Done reading STATSGO v1 soil fraction files'

end subroutine read_STATSGOv1_soilfractions
