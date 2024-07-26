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
! !ROUTINE: read_Special_soilfractions
! \label{read_Special_soilfractions}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  03 Aug  2012: KR Arsenault; Expanded to all soil fractions and tiling
!
! !INTERFACE:
subroutine read_Special_soilfractions(n, num_bins, soilsfgrd, &
                                  sandave, clayave, siltave )
! !USES:
  use LDT_coreMod,   only : LDT_rc
  use LDT_logMod,    only : LDT_logunit, LDT_verify, LDT_endrun
  use LDT_gridmappingMod
  use LDT_fileIOMod, only : LDT_transform_paramgrid
  use LDT_paramTileInputMod, only: param_2dbin_areacalc
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in)   :: n          ! nest index
  integer, intent(in)   :: num_bins   ! number of bins for tiling
  real,    intent(out)  :: soilsfgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real,    intent(out)  :: sandave(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real,    intent(out)  :: clayave(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real,    intent(out)  :: siltave(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
 
  real                  :: gravelave(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)

! !DESCRIPTION:
!  This subroutine retrieves Special soil, silt, clayfraction data 
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
  integer :: i, t, c, r, line
  integer :: ftn, ierr
  logical :: file_exists

! Netcdf file read-in entries:
  integer :: sandid, clayid, siltid, gravelid
  integer :: latid, lonid
  integer :: nrows, ncols
! 
  integer :: subpnc, subpnr, glpnc, glpnr  ! Number of rows, cols
  integer :: mi                            ! Total number of input param grid array points
  integer :: mo                            ! Total number of output LIS grid array points
  real    :: param_gridDesc(20)            ! Input parameter grid desc array
  real    :: subparam_gridDesc(20)         ! Subsetted input parameter grid desc array
  integer, allocatable :: lat_line(:,:), lon_line(:,:)
  integer, allocatable :: n11(:)           ! Maps each input grid point to output grid.
  real,    allocatable :: gisand(:), giclay(:), gisilt(:), gigrav(:)   ! input parameter 1d grid
  logical*1,allocatable:: li1(:), li2(:), li3(:), li4(:)   ! input logical mask (to match gi)
  real    :: gosand(LDT_rc%lnc(n)*LDT_rc%lnr(n),num_bins)  ! output parameter 1d grid
  real    :: goclay(LDT_rc%lnc(n)*LDT_rc%lnr(n),num_bins)   
  real    :: gosilt(LDT_rc%lnc(n)*LDT_rc%lnr(n),num_bins)  
  real    :: gograv(LDT_rc%lnc(n)*LDT_rc%lnr(n),num_bins) 
  logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n)), lo2(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! output logical mask (to match go)
  logical*1 :: lo3(LDT_rc%lnc(n)*LDT_rc%lnr(n)), lo4(LDT_rc%lnc(n)*LDT_rc%lnr(n))

  real, allocatable :: read_file(:,:)        ! Read in parameter
  real, allocatable :: yrev_file(:,:)        ! Y-reverse parameter
  real, allocatable :: read_sand_sub(:,:)    ! Read input parameter
  real, allocatable :: read_clay_sub(:,:)    ! Read input parameter
  real, allocatable :: read_silt_sub(:,:)    ! Read input parameter
  real, allocatable :: read_gravel_sub(:,:)  ! Read input parameter

!________________________________________________________________________

   write(LDT_logunit,*) "MSG: Reading Special sand, clay and silt files: ",&
      trim(LDT_rc%safile(n)),", ",trim(LDT_rc%clfile(n)),",",trim(LDT_rc%sifile(n)) 

! -------------------------------------------------------------------
!  CHECK FOR, OPEN AND READ SOIL FRACTION FILES
! -------------------------------------------------------------------

!=====  Sand file:  =====

   inquire(file=trim(LDT_rc%safile(n)), exist=file_exists)
   if(.not.file_exists) then 
      write(LDT_logunit,*) "Sand map ",trim(LDT_rc%safile(n))," not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
   ierr = nf90_open(path=trim(LDT_rc%safile(n)),mode=NF90_NOWRITE,ncid=ftn)
   call LDT_verify(ierr,'error opening special sand data')

   ierr = nf90_inq_varid(ftn,'sand_percent',sandid)
   call LDT_verify(ierr, 'nf90_inq_varid failed for texture in read_Special_soilfractions')

   ierr = nf90_inq_dimid(ftn,'lon',lonid)
   call LDT_verify(ierr,'nf90_inq_dimid failed for longitude in read_Special_soilfractions')

   ierr = nf90_inq_dimid(ftn,'lat',latid)
   call LDT_verify(ierr,'nf90_inq_dimid failed for latitude in read_Special_soilfractions')

   ierr = nf90_inquire_dimension(ftn,lonid,len=ncols)
   call LDT_verify(ierr,'nf90_inquire_dimension for longitude')

   ierr = nf90_inquire_dimension(ftn,latid,len=nrows)
   call LDT_verify(ierr,'nf90_inquire_dimension for latitude')

   allocate( read_file(ncols,nrows) )
   ierr = nf90_get_var(ftn, sandid, read_file)
   call LDT_verify(ierr, 'nf90_get_var failed for sand fraction')

   ierr = nf90_close(ftn)
   call LDT_verify(ierr, 'nf90_close failed in read_Special_soilfractions')
#endif
! ________________________

!- Set parameter grid array inputs:
   param_gridDesc(1)  = 0.    ! Latlon
   param_gridDesc(2)  = float(ncols)
   param_gridDesc(3)  = float(nrows)
   param_gridDesc(4)  = -60.0000   ! LL lat
   param_gridDesc(5)  = -180.0000  ! LL lon
   param_gridDesc(6)  = 128
   param_gridDesc(7)  = 90.0000    ! UR lat
   param_gridDesc(8)  = 180.0000   ! UR lon
   param_gridDesc(9)  = 0.009000900090009
   param_gridDesc(10) = 0.009000900090009
   param_gridDesc(20) = 64

!- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, LDT_rc%lc_proj, param_gridDesc(:), &
            glpnc, glpnr, subpnc, subpnr, subparam_gridDesc, lat_line, lon_line )

 ! Reverse-Y and Convert 8-bit unsigned integers:
   allocate( yrev_file(ncols,nrows) )
   yrev_file = LDT_rc%udef
   i = 0
   do r = nrows, 1, -1
      i = i + 1
      do c = 1, ncols
         yrev_file(c,i) = read_file(c,r)
      end do
   end do
   deallocate( read_file )
!- Subset parameter read-in array:
   allocate( read_sand_sub(subpnc, subpnr) )
   read_sand_sub = LDT_rc%udef
   do r = 1, subpnr
      do c = 1, subpnc
         read_sand_sub(c,r) = yrev_file(lon_line(c,r),lat_line(c,r))
      enddo
   enddo
   deallocate(yrev_file)

!=====  Clay file:  =====

   inquire(file=trim(LDT_rc%clfile(n)), exist=file_exists)
   if(.not.file_exists) then
      write(LDT_logunit,*) "Clay map ",trim(LDT_rc%clfile(n))," not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
   ierr = nf90_open(path=trim(LDT_rc%clfile(n)),mode=NF90_NOWRITE,ncid=ftn)
   call LDT_verify(ierr,'error opening special clay data')

   ierr = nf90_inq_varid(ftn,'clay_percent',clayid)
   call LDT_verify(ierr, 'nf90_inq_varid failed for clay fraction in read_Special_soilfractions')

   allocate( read_file(ncols,nrows) )
   ierr = nf90_get_var(ftn, clayid, read_file)
   call LDT_verify(ierr, 'nf90_get_var failed for clay fraction')

   ierr = nf90_close(ftn)
   call LDT_verify(ierr, 'nf90_close failed in read_Special_soilfractions')
#endif
 ! Reverse-Y and Convert 8-bit unsigned integers:
   allocate( yrev_file(ncols,nrows) )
   yrev_file = LDT_rc%udef
   i = 0
   do r = nrows, 1, -1
      i = i + 1
      do c = 1, ncols
         yrev_file(c,i) = read_file(c,r)
      end do
   end do
   deallocate( read_file )
!- Subset parameter read-in array:
   allocate( read_clay_sub(subpnc, subpnr) )
   read_clay_sub = LDT_rc%udef
   do r = 1, subpnr
      do c = 1, subpnc
         read_clay_sub(c,r) = yrev_file(lon_line(c,r),lat_line(c,r))
      enddo
   enddo
   deallocate(yrev_file)

!=====  Silt file:  =====

   inquire(file=trim(LDT_rc%sifile(n)), exist=file_exists)
   if(.not.file_exists) then
      write(LDT_logunit,*) "Silt map ",trim(LDT_rc%sifile(n))," not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
   ierr = nf90_open(path=trim(LDT_rc%sifile(n)),mode=NF90_NOWRITE,ncid=ftn)
   call LDT_verify(ierr,'error opening special silt data')

   ierr = nf90_inq_varid(ftn,'silt_percent',siltid)
   call LDT_verify(ierr, 'nf90_inq_varid failed for silt fraction in read_Special_soilfractions')

   allocate( read_file(ncols,nrows) )
   ierr = nf90_get_var(ftn, siltid, read_file)
   call LDT_verify(ierr, 'nf90_get_var failed for silt fraction')

   ierr = nf90_close(ftn)
   call LDT_verify(ierr, 'nf90_close failed in read_Special_soilfractions')
#endif

 ! Reverse-Y and Convert 8-bit unsigned integers:
   allocate( yrev_file(ncols,nrows) )
   yrev_file = LDT_rc%udef
   i = 0
   do r = nrows, 1, -1
      i = i + 1
      do c = 1, ncols
         yrev_file(c,i) = read_file(c,r)
      end do
   end do
   deallocate( read_file )
!- Subset parameter read-in array:
   allocate( read_silt_sub(subpnc, subpnr) )
   read_silt_sub = LDT_rc%udef
   do r = 1, subpnr
      do c = 1, subpnc
         read_silt_sub(c,r) = yrev_file(lon_line(c,r),lat_line(c,r))
      enddo
   enddo
   deallocate(yrev_file)

!=====  Gravel file:  =====

   inquire(file=trim(LDT_rc%gravelfile(n)), exist=file_exists)
   if(.not.file_exists) then
      write(LDT_logunit,*) "Silt map ",trim(LDT_rc%gravelfile(n))," not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
   ierr = nf90_open(path=trim(LDT_rc%gravelfile(n)),mode=NF90_NOWRITE,ncid=ftn)
   call LDT_verify(ierr,'error opening special gravel data')

   ierr = nf90_inq_varid(ftn,'gravel_percent',gravelid)
   call LDT_verify(ierr,'nf90_inq_varid failed for gravel fraction in read_Special_soilfractions')

   allocate( read_file(ncols,nrows) )
   ierr = nf90_get_var(ftn, gravelid, read_file)
   call LDT_verify(ierr, 'nf90_get_var failed for gravel fraction')

   ierr = nf90_close(ftn)
   call LDT_verify(ierr, 'nf90_close failed in read_Special_soilfractions')
#endif
 ! Reverse-Y and Convert 8-bit unsigned integers:
   allocate( yrev_file(ncols,nrows) )
   yrev_file = LDT_rc%udef
   i = 0
   do r = nrows, 1, -1
      i = i + 1
      do c = 1, ncols
         yrev_file(c,i) = read_file(c,r)
      end do
   end do
   deallocate( read_file )
!- Subset parameter read-in array:
   allocate( read_gravel_sub(subpnc, subpnr) )
   read_gravel_sub = LDT_rc%udef
   do r = 1, subpnr
      do c = 1, subpnc
         read_gravel_sub(c,r) = yrev_file(lon_line(c,r),lat_line(c,r))
      enddo
   enddo
   deallocate(yrev_file)

! ===================================================================

! -------------------------------------------------------------------
!   UPSCALING/DOWNSCALING GRIDS TO LIS OUTPUT GRID
! -------------------------------------------------------------------
   mi = nrows*ncols
   allocate( gisand(mi), li1(mi) )
   allocate( giclay(mi), li2(mi) )
   allocate( gisilt(mi), li3(mi) )
   allocate( gigrav(mi), li4(mi) )
   gisand = LDT_rc%udef;  giclay  = LDT_rc%udef
   gisilt = LDT_rc%udef;  gigrav  = LDT_rc%udef
   li1 = .false.; li2 = .false.; li3 = .false.; li4 = .false.

   mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
   lo1 = .false.; lo2 = .false.; lo3 = .false.; lo4 = .false.
   sandave  = 0.; clayave = 0.; siltave = 0.; gravelave = 0.
   soilsfgrd = 0.

!- Assign 2-D array to 1-D for aggregation routines:
   i = 0
   do r = 1, subpnr
      do c = 1, subpnc;  i = i + 1
         gisand(i) = read_sand_sub(c,r)
         giclay(i) = read_clay_sub(c,r)
         gisilt(i) = read_silt_sub(c,r)
         gigrav(i) = read_gravel_sub(c,r)
         if( gisand(i) < 0. .or. gisand(i) > 100.0 ) gisand(i) = LDT_rc%udef
         if( giclay(i) < 0. .or. giclay(i) > 100.0 ) giclay(i) = LDT_rc%udef
         if( gisilt(i) < 0. .or. gisilt(i) > 100.0 ) gisilt(i) = LDT_rc%udef
         if( gigrav(i) < 0. .or. gigrav(i) > 100.0 ) gigrav(i) = LDT_rc%udef

         if( gisand(i) .ne. LDT_rc%udef ) gisand(i) = gisand(i)/100.
         if( giclay(i) .ne. LDT_rc%udef ) giclay(i) = giclay(i)/100.
         if( gisilt(i) .ne. LDT_rc%udef ) gisilt(i) = gisilt(i)/100.
         if( gigrav(i) .ne. LDT_rc%udef ) gigrav(i) = gigrav(i)/100.

         if( gisand(i) .ne. LDT_rc%udef )  li1(i) = .true.
         if( giclay(i) .ne. LDT_rc%udef )  li2(i) = .true.
         if( gisilt(i) .ne. LDT_rc%udef )  li3(i) = .true.
         if( gigrav(i) .ne. LDT_rc%udef )  li4(i) = .true.
      enddo
   enddo
   deallocate( read_sand_sub, read_clay_sub, read_silt_sub, read_gravel_sub )


!- Account for single vs. tiled soil fraction layers:
   select case ( LDT_rc%soils_gridtransform(n) )

  !- Single layer transformation:
     case( "none", "average", "neighbor", "bilinear", "budget-bilinear" )

    !- Transform parameter from original grid to LIS output grid:
       call LDT_transform_paramgrid(n, LDT_rc%soils_gridtransform(n), &
                subparam_gridDesc, mi, 1, gisand, li1, mo, gosand(:,1), lo1 )

       call LDT_transform_paramgrid(n, LDT_rc%soils_gridtransform(n), &
                subparam_gridDesc, mi, 1, giclay, li2, mo, goclay(:,1), lo2 )

       call LDT_transform_paramgrid(n, LDT_rc%soils_gridtransform(n), &
                subparam_gridDesc, mi, 1, gisilt, li3, mo, gosilt(:,1), lo3 )

       call LDT_transform_paramgrid(n, LDT_rc%soils_gridtransform(n), &
                subparam_gridDesc, mi, 1, gigrav, li4, mo, gograv(:,1), lo4 )

       deallocate( gisand, giclay, gisilt, gigrav )
       deallocate( li1, li2, li3, li4 )

    !- Convert 1D soil fractions to 2D grid arrays:
       i = 0
       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n)
             i = i + 1
             sandave(c,r,1) = gosand(i,1)
             clayave(c,r,1) = goclay(i,1)
             siltave(c,r,1) = gosilt(i,1)
             gravelave(c,r,1) = gograv(i,1)
           ! Make sure all negative values are set to universal undefined value:
              if( sandave(c,r,1) < 0. ) sandave(c,r,1) = LDT_rc%udef
              if( clayave(c,r,1) < 0. ) clayave(c,r,1) = LDT_rc%udef
              if( siltave(c,r,1) < 0. ) siltave(c,r,1) = LDT_rc%udef
              if( gravelave(c,r,1) < 0. ) gravelave(c,r,1) = LDT_rc%udef
          enddo
       enddo
       soilsfgrd(:,:,:) = 1.0


  !- Tiling: Two-field combination for bin-based grid fraction and sand,clay estimates:
     case( "tile" )

        allocate( n11(mi) )
     !- Create mapping between parameter domain and LIS grid domain:
        call upscaleByAveraging_input( subparam_gridDesc, LDT_rc%gridDesc(n,:), &
                                       mi, mo, n11 )

      ! Calculate tiled soil fractions:
        call param_2dbin_areacalc( n, num_bins, mi, mo, n11,    &
                                   LDT_rc%udef, gisand, giclay, &
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
        deallocate( gisand, giclay, gisilt, gigrav )
        deallocate( li1, li2, li3, li4 )
        deallocate( n11 )

     case default
        write(LDT_logunit,*) "ERR MSG: This spatial transform, ",&
             trim(LDT_rc%soils_gridtransform(n)),&
             ", for soil fraction is not available at this time ... "
        write(LDT_logunit,*) "Program stopping ..."
        call LDT_endrun

    end select

! _____________________________________________


   write(LDT_logunit,*) "MSG: Done reading Special soil fraction files."

end subroutine read_Special_soilfractions
