!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: read_WorldClim_climppt
! \label{read_WorldClim_climppt}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  10 May 2013: K. Arsenault; Added climatology datasets
!
! !INTERFACE:
 subroutine read_WorldClim_climppt(n, out_ncols, out_nrows, &
                 out_gridDesc, ppt_data)

! !USES:
  use ESMF
  use LDT_coreMod,  only : LDT_rc, LDT_domain
  use LDT_logMod,   only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use fbil_module
  use LDT_climateParmsMod

  implicit none

! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: out_ncols
  integer, intent(in)    :: out_nrows
  real,    intent(in)    :: out_gridDesc(20)
  real,    intent(inout) :: ppt_data(out_ncols,out_nrows)

! !DESCRIPTION:
!  This subroutine retrieves the WorldClim precipitation data climatology.
!  
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[ppt_data]
!   output field with WorldClim precipitation climatogy data
!  \end{description}
!
!EOP      
  integer   :: iret
  logical   :: file_exists
  integer   :: ferror_bil
  integer, save :: file_status
  integer   :: c, r, i
  integer   :: mi                        ! Total number of input param grid array points
  integer   :: mo                        ! Total number of output LIS grid array points
  integer, allocatable :: n11(:)         ! Map array for aggregating methods
  integer, pointer     :: n113(:)        ! Map array for nearest neighbor interp
  real,    allocatable :: gi1(:)         ! input parameter 1d grid
  logical*1,allocatable:: li1(:)         ! input logical mask (to match gi)
  real      :: go1(out_ncols*out_nrows)  ! output lis 1d grid
  logical*1 :: lo1(out_ncols*out_nrows)  ! output logical mask (to match go)
  real      :: param_gridDesc(20)
!  real      :: subparam_gridDesc(20)

!- BIL Format file read paramters:
  integer, parameter :: xd=43200, yd=18000       ! Dimension of original data

  type(charN)                    :: bil_filename
  type(File_bil_header), pointer :: hdrInfPtr
  type(File_bil_header), target  :: hdrInfW
  real*4, pointer, dimension(:),save  :: real4ptr1dL  ! 4-byte real array

  character(len=200) :: last_file_read = "none"

! __________________________________________________________________________________________

!- Set parameter grid array inputs:
   param_gridDesc(1)  = 0.    ! Latlon
   param_gridDesc(2)  = float(xd)
   param_gridDesc(3)  = float(yd)
   param_gridDesc(4)  = -59.99166667   ! Edges of gridcell (ArcGIS)
   param_gridDesc(5)  =-179.99583333   ! Edges of gridcell 
   param_gridDesc(7)  =  89.99583333
   param_gridDesc(8)  = 179.99583333
   param_gridDesc(9)  =   0.008333333
   param_gridDesc(10) =   0.008333333
   param_gridDesc(20) = 64
!   print *, "param griddesc: ", param_gridDesc(2:10)
!   print *, "output griddec: ",  out_gridDesc(2:10)

   ferror_bil = 1    ! 1 - indicates bil file read success
   ppt_data = LDT_rc%udef

!- Check if file has already been read in:
   if( last_file_read /= LDT_climate_struc(n)%climpptfile ) then
     last_file_read = LDT_climate_struc(n)%climpptfile

     inquire(file=trim(LDT_climate_struc(n)%climpptfile)//".bil", exist=file_exists)
     if(.not. file_exists) then 
        write(LDT_logunit,*) "WorldClim Clim PPT map ",trim(LDT_climate_struc(n)%climpptfile)," not found."
        write(LDT_logunit,*) "Program stopping ..."
        call LDT_endrun
     endif

  !- Assign filename pointer and initialize header file arrays:
     allocate(bil_filename%str)
     bil_filename%str = LDT_climate_struc(n)%climpptfile

     hdrInfPtr => hdrInfW
     hdrInfW%refCoordsIsNull = .true.
     nullify(hdrInfW%refCoords)     ! This must be done just prior to the read call
     hdrInfW%abortFileReadIfHeaderDataRegionDoesntSubsumeAnalysisRegion = .false.

!     print *, "Reading in ClimPPT data:: ", trim(bil_filename%str)

  !- Read BIL file header information:
     call read_bil_or_float_headerFile(bil_filename, hdrInfW)
 
     if( fileReadFailedBasedOnValXremainingForVarY(&
         int4VarY=hdrInfPtr%ncols, int4ValX=gInitValNcols) .neqv. .true. ) then

     ! Perform actual data file read into the provided memory location
     !  Stage data into a 1d memory per size & order requirements in its header:

       if ( .not. associated(real4ptr1dL) ) then
          allocate( real4ptr1dL(hdrInfPtr%ncols*hdrInfPtr%nrows), STAT=file_status)
       endif

       if( file_status == gALLOCATE_SUCCESS ) then
          if( interpret_bilHdrInfoAndPopulateMemLoc( &
               trim(bil_filename%str)//'.bil', &       ! WorldClim filename
               hdrInfPtr, real4ptr1d=real4ptr1dL ) &
            .eqv. .false.) then
            write(LDT_logunit,*) " ... Read file failed at populate mem location"
            ferror_bil = 0
          endif
       else
          write(LDT_logunit,*) " ... Read file failed at allocate mem location"
          ferror_bil = 0
       endif
    else
       write(LDT_logunit,*) " ... Read file failed at header file read"
       ferror_bil = 0
    endif

  !- Flip (y-reverse) for writing out:
  !  (flipping data from 90N->90S to 90S->90N LIS standard)
     if( flipY(real4ptr1d=real4ptr1dL, &
         ncols=hdrInfPtr%ncols, nrows=hdrInfPtr%nrows) .eqv. .false.) then
         write(LDT_logunit,*) " ... Read file failed at y-reverse."
         ferror_bil = 0
     endif
     deallocate(bil_filename%str)
   endif  ! Check if file has already been read-in


!- Data read in successfully:
   if( ferror_bil == 1 ) then

! -------------------------------------------------------------------
!     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
! -------------------------------------------------------------------
     mi = xd * yd
     mo = out_ncols*out_nrows

     select case( LDT_climate_struc(n)%clim_gridtransform )

   !- Create mapping between parameter domain and LIS grid domain:
      case( "average", "neighbor" )

        if( LDT_climate_struc(n)%clim_gridtransform == "average" ) then
           allocate( n11(mi) )
           write(LDT_logunit,*) "#Regridding: Applying average to input parameter"
           call upscaleByAveraging_input( param_gridDesc, &
                            out_gridDesc(:), mi, mo, n11 )

        elseif( LDT_climate_struc(n)%clim_gridtransform == "neighbor" ) then
           allocate( n113(mo) )
           write(LDT_logunit,*) "#Regridding: Applying nearest neighbor to input parameter"
           call neighbor_interp_input( n, param_gridDesc, n113 )
        endif

        allocate( gi1(mi), li1(mi) )
        gi1 = LDT_rc%udef
        li1 = .false.
        lo1 = .false.

     !- Assign 2-D array to 1-D for aggregation routines:
!        print *, " Setting World Clim Data Array for Averaging (large loop)  ..."
        do i = 1, mi
           gi1(i) = real4ptr1dL(i)
           if( gi1(i) .ne. LDT_rc%udef ) then
             li1(i) = .true.
           endif
        enddo

     !- Average finer scale points to output grid:
        if( LDT_climate_struc(n)%clim_gridtransform == "average" ) then
          call upscaleByAveraging( mi, mo, LDT_rc%udef, n11,  &
                                 li1, gi1, lo1, go1 )
          deallocate( n11 )

     !- Use nearest neighbor to locate finer scale points to output grid:
        elseif( LDT_climate_struc(n)%clim_gridtransform == "neighbor" ) then
           call neighbor_interp( out_gridDesc, li1, gi1, lo1(:), go1(:), mi, mo, &
                                 LDT_domain(n)%lat, LDT_domain(n)%lon, &
                                 n113, LDT_rc%udef, iret )
           deallocate( n113 )
        endif
        deallocate( li1, gi1 )

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

   != Other spatial transforms ...
      case default
        write(LDT_logunit,*) "- Other spatial transforms are not yet supported for"
        write(LDT_logunit,*) "- WorldClim Clim PPT maps ... please select 'average' for now."
        write(LDT_logunit,*) "- Program stopping ..."
        call LDT_endrun

     end select  ! end grid transform condition

   elseif( ferror_bil == 0 ) then
      write(LDT_logunit,*) "Missing WorldCLIM PPT data ", trim(LDT_climate_struc(n)%climpptfile)

   endif

end subroutine read_WorldClim_climppt
