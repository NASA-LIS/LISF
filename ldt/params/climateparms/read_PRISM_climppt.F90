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
! !ROUTINE: read_PRISM_climppt
! \label{read_PRISM_climppt}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  10 Oct 2012: K. Arsenault; Added PRISM climatology datasets
!
! !INTERFACE:
 subroutine read_PRISM_climppt(n, out_ncols, out_nrows, &
                 out_gridDesc, ppt_data)

! !USES:
  use ESMF
  use LDT_coreMod,   only : LDT_rc, LDT_domain
  use LDT_logMod,    only : LDT_logunit, LDT_getNextUnitNumber, &
               LDT_releaseUnitNumber, LDT_endrun
  use LDT_climateParmsMod

  implicit none

! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: out_ncols
  integer, intent(in)    :: out_nrows
  real,    intent(in)    :: out_gridDesc(20)
  real,    intent(inout) :: ppt_data(out_ncols,out_nrows)

! !DESCRIPTION:
!  This subroutine retrieves the PRISM precipitation data climatology.
!  
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[ppt_data]
!   output field with PRISM precipitation climatogy data
!  \end{description}
!
!EOP      
  integer   :: ftn
  logical   :: file_exists
  integer   :: c, r, i, iret 
  integer   :: ncols, nrows
  real      :: xllcorner, yllcorner
  real      :: cellsize
  real      :: nodata_value
  integer   :: mi                        ! Total number of input param grid array points
  integer   :: mo                        ! Total number of output LIS grid array points
  real      :: param_gridDesc(20)
  real,     allocatable :: gi1(:)        ! input parameter 1d grid
  logical*1,allocatable :: li1(:)        ! input logical mask (to match gi)
  real      :: go1(out_ncols*out_nrows)  ! output lis 1d grid
  logical*1 :: lo1(out_ncols*out_nrows)  ! output logical mask (to match go)

!- Grid transform arrays:
  integer, allocatable     :: n11(:)     ! Map array for aggregating methods
  integer, allocatable     :: n113(:)    ! Map array for nearest neighbor interp
  integer, allocatable     :: n111(:)    ! Map array for bilinear interp
  integer, allocatable     :: n121(:)
  integer, allocatable     :: n211(:)
  integer, allocatable     :: n221(:)
  real, allocatable        :: w111(:),w121(:)
  real, allocatable        :: w211(:),w221(:)
  integer, allocatable     :: n112(:,:)  ! Map array for budget/conservative interp
  integer, allocatable     :: n122(:,:)
  integer, allocatable     :: n212(:,:)
  integer, allocatable     :: n222(:,:)
  real, allocatable        :: w112(:,:),w122(:,:)
  real, allocatable        :: w212(:,:),w222(:,:)

  character(20)     :: temp
  real, allocatable :: read_input(:,:)
  real, parameter   :: scale_factor = 0.01

! __________________________________________________________________________________________

   ppt_data = LDT_rc%udef

   inquire(file=trim(LDT_climate_struc(n)%climpptfile), exist=file_exists)
   if(.not. file_exists) then 
      write(LDT_logunit,*) "[ERR] PRISM Clim PPT map ",&
            trim(LDT_climate_struc(n)%climpptfile)," not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
   select case ( LDT_climate_struc(n)%clim_gridtransform )
    case( "none", "neighbor", "average", "bilinear", "budget-bilinear" )
    case default
      write(LDT_logunit,*) "[ERR] The spatial transform option selected for PRISM Climatological" 
      write(LDT_logunit,*) "     precipitation file is not recognized nor recommended."
      write(LDT_logunit,*) "     Please select: "
      write(LDT_logunit,*) "  ==  none, neighbor, average, bilinear, budget-bilinear "
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   end select

   ftn = LDT_getNextUnitNumber()
   open(ftn, file=trim(LDT_climate_struc(n)%climpptfile),form="formatted",status="old" )
  
   read(ftn,'(a5,1x,i4)'   )  temp, ncols
   read(ftn,'(a5,1x,i4)'   )  temp, nrows
   read(ftn,'(a9,1x,f12.7)')  temp, xllcorner
   read(ftn,'(a9,1x,f12.7)')  temp, yllcorner
   read(ftn,'(a9,1x,f10.8)')  temp, cellsize
   read(ftn,'(a12,1x,f5.0)')  temp, nodata_value

!- Read original dataset:
   allocate( read_input(ncols,nrows) )
   do r = nrows, 1, -1
      read(ftn,*) ( read_input(c,r), c=1,ncols )
   end do

! -------------------------------------------------------------------
!     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
! -------------------------------------------------------------------

!- Set parameter grid array inputs:
   param_gridDesc(1)  = 0.    ! Latlon
   param_gridDesc(2)  = float(ncols)
   param_gridDesc(3)  = float(nrows)
   param_gridDesc(4)  =  24.062500     ! Edges of gridcell (ArcGIS)
   param_gridDesc(5)  = -125.020833    ! Edges of gridcell
   param_gridDesc(7)  =  49.9375
   param_gridDesc(8)  = -66.479166
   param_gridDesc(9)  = 0.008333
   param_gridDesc(10) = 0.008333
   param_gridDesc(20) = 64

!   call LDT_checkDomainExtents(n, param_gridDesc(:))

   mi = ncols * nrows
   mo = out_ncols*out_nrows
   allocate( gi1(mi), li1(mi) )
   gi1 = LDT_rc%udef
   li1 = .false.
   lo1 = .false.

!- Assign 2-D array to 1-D for aggregation routines:
   i = 0
   do r = 1, nrows
      do c = 1, ncols;  i = i + 1
         gi1(i) = read_input(c,r) * scale_factor
         if( gi1(i) .ne. LDT_rc%udef ) then
           li1(i) = .true.
         endif
      enddo
   enddo
   deallocate( read_input )

! -------------------------------------------------------------------
!     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
! -------------------------------------------------------------------

!- Select grid spatial transform option:
   select case ( LDT_climate_struc(n)%clim_gridtransform )

  !- Aggregate by calculating average of each output gridcell:
     case ( "average" )
       allocate( n11(mi) )
       write(LDT_logunit,*)"[INFO] Regridding: Applying average to input parameter"
       call upscaleByAveraging_input( param_gridDesc, out_gridDesc, mi, mo, n11 )
       call upscaleByAveraging( mi, mo, LDT_rc%udef, n11,li1, gi1, lo1(:), go1(:) )

  !- Select neighboring point:
     case ( "neighbor" )
       allocate( n113(mo) )
       write(LDT_logunit,*)"[INFO] Regridding: Applying nearest neighbor to input parameter"
       call neighbor_interp_input( n, param_gridDesc, n113 )
       call neighbor_interp( out_gridDesc, li1, gi1, lo1(:), go1(:), mi, mo, &
                             LDT_domain(n)%lat, LDT_domain(n)%lon, n113, LDT_rc%udef, iret )
       deallocate( n113 )

  !- Bilinear interpolation:
     case ( "bilinear" )
       allocate( n111(mo), n121(mo), n211(mo), n221(mo) )
       allocate( w111(mo), w121(mo), w211(mo), w221(mo) )
       write(LDT_logunit,*)"[INFO] Regridding: Applying bilinear interp to input parameter"
       call bilinear_interp_input( n, param_gridDesc, &
            n111, n121, n211, n221, w111, w121, w211, w221)
       call bilinear_interp(out_gridDesc, li1, gi1, lo1(:), go1(:), &
            mi, mo, LDT_domain(n)%lat, LDT_domain(n)%lon, &
            w111, w121, w211, w221, n111, n121, n211, n221,&
            LDT_rc%udef, iret)
       deallocate( n111, n121, n211, n221, w111, w121, w211, w221 )

  !- Budget-bilinear interpolation:
     case ( "budget-bilinear" )
       allocate( n112(mo,25), n122(mo,25), n212(mo,25), n222(mo,25) )
       allocate( w112(mo,25), w122(mo,25), w212(mo,25), w222(mo,25) )
       write(LDT_logunit,*)"[INFO] Regridding: Applying budget-bilinear interp to input parameter"
       call conserv_interp_input( n, param_gridDesc, &
            n112, n122, n212, n222, w112, w122, w212, w222 )
       call conserv_interp( out_gridDesc, li1, gi1, lo1(:), go1(:), &
            mi, mo, LDT_domain(n)%lat, LDT_domain(n)%lon, &
            w112, w122, w212, w222, n112, n122, n212, n222, &
            LDT_rc%udef, iret)
       deallocate( n112, n122, n212, n222, w112, w122, w212, w222 )

  !- When no transform is performed:
     case ( "none" )
        write(LDT_logunit,*) "[INFO] No aggregation applied for parameter file ... "
        go1(:) = gi1(:)

    case default
        write(*,*) "[ERR] This spatial transformation option ("&
                    //trim(LDT_climate_struc(n)%clim_gridtransform)//") "
        write(*,*) "  is not currently supported."
        write(*,*) "Program stopping ...."
        call LDT_endrun
   end select

!- Convert 1D to 2D grid output arrays:
   i = 0
   do r = 1, out_nrows
      do c = 1, out_ncols;  i = i + 1
         if( go1(i) < 0.) then
            ppt_data(c,r) = LDT_rc%udef
         else
            ppt_data(c,r) = go1(i)
         endif
      enddo
   enddo
   deallocate( gi1, li1 )

   call LDT_releaseUnitNumber(ftn)

end subroutine read_PRISM_climppt
