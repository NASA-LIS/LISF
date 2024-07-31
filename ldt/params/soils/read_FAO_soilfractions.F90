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
! !ROUTINE: read_FAO_soilfractions
! \label{read_FAO_soilfractions}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  03 Aug  2012: KR Arsenault; Expanded to all soil fractions and tiling
!
! !INTERFACE:
subroutine read_FAO_soilfractions(n, num_bins, soilsfgrd, &
                                  sandave, clayave, siltave )
! !USES:
  use LDT_coreMod,     only : LDT_rc
  use LDT_logMod,      only : LDT_logunit, LDT_getNextUnitNumber, &
          LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod,   only : readLISdata
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
!  This subroutine retrieves FAO soil, silt, clayfraction data 
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
   integer :: subpnc, subpnr, glpnc, glpnr
   integer :: mi                            ! Total number of input param grid array points
   integer :: mo                            ! Total number of output LIS grid array points
   real    :: subparam_gridDesc(20)         ! Input parameter grid desc array
   integer, allocatable :: lat_line(:,:), lon_line(:,:)
   integer, allocatable :: n11(:)           ! Maps each input grid point to output grid.
   real,    allocatable :: gi1(:), gi2(:)   ! input parameter 1d grid
   logical*1,allocatable:: li1(:), li2(:)   ! input logical mask (to match gi)
   real,    allocatable :: read_sand_sub(:,:)    ! Read input parameter
   real,    allocatable :: read_clay_sub(:,:)    ! Read input parameter
   real    :: temp(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
!________________________________________________________________________

   write(LDT_logunit,*) "[INFO] Reading FAO sand, clay and silt files: ",&
        trim(LDT_rc%safile(n)),", ",trim(LDT_rc%clfile(n)),",",trim(LDT_rc%sifile(n)) 

! -------------------------------------------------------------------
!  CHECK FOR AND OPEN SOIL FRACTION FILES
! -------------------------------------------------------------------
!- Sand file:
   inquire(file=trim(LDT_rc%safile(n)), exist=file_exists)
   if(.not.file_exists) then 
      write(LDT_logunit,*) "Sand map ",trim(LDT_rc%safile(n))," not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
   ftn_sa = LDT_getNextUnitNumber()
   open(ftn_sa,file=LDT_rc%safile(n),form='unformatted',status='old',&
               access='direct',recl=4)

!- Clay file:
   inquire(file=trim(LDT_rc%clfile(n)), exist=file_exists)
   if(.not.file_exists) then
      write(LDT_logunit,*) "Clay map ",trim(LDT_rc%clfile(n))," not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
   ftn_cl = LDT_getNextUnitNumber()
   open(ftn_cl,file=LDT_rc%clfile(n),form='unformatted',status='old',&
               access='direct',recl=4)

!- Silt file (optional):
   inquire(file=trim(LDT_rc%sifile(n)), exist=file_exists)
   if(.not.file_exists) then
      write(LDT_logunit,*) "Silt map ",trim(LDT_rc%sifile(n))," not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
   ftn_si = LDT_getNextUnitNumber()
   open(ftn_si,file=LDT_rc%sifile(n),form='unformatted',status='old',&
               access='direct',recl=4)

   sandave  = 0.; clayave = 0.; siltave = 0.
   soilsfgrd = 0.

!- Account for single vs. tiled soil fraction layers:
   select case ( LDT_rc%soils_gridtransform(n) )

  !- Single layer transformation:
     case( "none", "average", "neighbor", "bilinear", "budget-bilinear" )

      ! Read in sand fraction field:
        temp = LDT_rc%udef
        call readLISdata(n, ftn_sa, LDT_rc%soils_proj, LDT_rc%soils_gridtransform(n), &
                    LDT_rc%soil_gridDesc(n,:), 1, temp )  ! 1 indicates 2D layer
        sandave(:,:,1) = temp(:,:,1)

      ! Read in clay fraction field:
        temp = LDT_rc%udef
        call readLISdata(n, ftn_cl, LDT_rc%soils_proj, LDT_rc%soils_gridtransform(n), &
                    LDT_rc%soil_gridDesc(n,:), 1, temp )  ! 1 indicates 2D layer
        clayave(:,:,1) = temp(:,:,1)

      ! Read in silt fraction field:
        temp = LDT_rc%udef
        call readLISdata(n, ftn_si, LDT_rc%soils_proj, LDT_rc%soils_gridtransform(n), &
                    LDT_rc%soil_gridDesc(n,:), 1, temp )  ! 1 indicates 2D layer
        siltave(:,:,1) = temp(:,:,1)

      ! Make sure and negative values are set to universal undefined value:
        do r=1,LDT_rc%lnr(n)
           do c=1,LDT_rc%lnc(n)
              if( sandave(c,r,1) < 0. ) sandave(c,r,1) = LDT_rc%udef
              if( clayave(c,r,1) < 0. ) clayave(c,r,1) = LDT_rc%udef
              if( siltave(c,r,1) < 0. ) siltave(c,r,1) = LDT_rc%udef
           enddo
        enddo
        soilsfgrd(:,:,1) = 1.0


  !- Tiling: Two-field combination for bin-based grid fraction and sand,clay estimates:
     case( "tile" )

         subparam_gridDesc = 0.
         call LDT_RunDomainPts( n, LDT_rc%soils_proj, LDT_rc%soil_gridDesc(n,:), &
                     glpnc, glpnr, subpnc, subpnr, subparam_gridDesc, lat_line, lon_line )

      !- Initialize parameter read-in array:
         allocate( read_sand_sub(subpnc,subpnr) )
         allocate( read_clay_sub(subpnc,subpnr) )
         read_sand_sub = LDT_rc%udef
         read_clay_sub = LDT_rc%udef
      !- Read-in sand and clay:
         line = 0
         do r = 1, subpnr
            do c = 1, subpnc
               line = (lat_line(c,r)-1)*glpnc + lon_line(c,r)
               read(ftn_sa,rec=line) read_sand_sub(c,r)
               read(ftn_cl,rec=line) read_clay_sub(c,r)
               if( read_sand_sub(c,r) < 0. ) read_sand_sub(c,r) = LDT_rc%udef
               if( read_clay_sub(c,r) < 0. ) read_clay_sub(c,r) = LDT_rc%udef
            enddo
         enddo
         mi = subpnc*subpnr
         mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
         allocate( gi1(mi), gi2(mi), li1(mi), li2(mi), n11(mi) )
         gi1 = LDT_rc%udef; gi2 = LDT_rc%udef; li1 = .false.; li2 = .false.
      !- Assign 2-D array to 1-D for aggregation routines:
         i = 0
         do r = 1, subpnr
            do c = 1, subpnc;  i = i + 1
               gi1(i) = read_sand_sub(c,r)
               gi2(i) = read_clay_sub(c,r)
               if( gi1(i) .ne. LDT_rc%udef )  li1(i) = .true.
               if( gi2(i) .ne. LDT_rc%udef )  li2(i) = .true.
            enddo
         enddo
      !- Create mapping between parameter domain and LIS grid domain:
         call upscaleByAveraging_input( subparam_gridDesc, LDT_rc%gridDesc(n,:), &
                                        mi, mo, n11 )

       ! Calculate tiled soil fractions:
         call param_2dbin_areacalc( n, num_bins, mi, mo, n11, LDT_rc%udef, gi1, gi2, &
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
         deallocate( gi1, gi2, li1, li2, n11, lat_line, lon_line )
         deallocate( read_sand_sub, read_clay_sub )

      case default
         write(LDT_logunit,*) "[INFO] This spatial transform, ",trim(LDT_rc%soils_gridtransform(n)),&
                              ", for soil fraction is not available at this time ... "
         write(LDT_logunit,*) "Program stopping ..."
         call LDT_endrun

    end select

! _____________________________________________

   call LDT_releaseUnitNumber(ftn_sa)
   call LDT_releaseUnitNumber(ftn_cl)
   call LDT_releaseUnitNumber(ftn_si)

   write(LDT_logunit,*) "[INFO] Done reading FAO soil fraction files."

end subroutine read_FAO_soilfractions
