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
! !ROUTINE: read_maskfile
!  \label{read_maskfile}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  01 June 2012: KR Arsenault; Restructured to simply read in a mask file
!  01 Sept 2013: KR Arsenault; Added ability to upscale mask readin
!  06 Oct  2017: Eric Kemp; Added special netCDF readers for AFWA.
!
! !INTERFACE:
subroutine read_maskfile(n, vegtype, fgrd, localmask )

! !USES:
  use LDT_coreMod,   only : LDT_rc, LDT_localPet
  use LDT_logMod,    only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_gridmappingMod
  use LDT_fileIOMod, only : LDT_transform_paramgrid

  implicit none

! !ARGUMENTS: 
  integer, intent(in)  :: n
  real,    intent(in)  :: vegtype(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real,    intent(in)  :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%nt)
  real,    intent(out) :: localmask(LDT_rc%lnc(n),LDT_rc%lnr(n))

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
  integer :: ftn, ios1
  logical :: file_exists
  integer :: file_size 
  integer :: rec_length
  integer :: c, r, t, i, j, line, iret
  integer :: glpnc, glpnr             ! Parameter (global) total columns and rows
  integer :: subpnc, subpnr           ! Parameter subsetted columns and rows
  real    :: subparam_gridDesc(20)    ! Input parameter grid desc array
  real    :: cnt_0mask_1lc
  real    :: cnt_1mask_0lc
  integer :: mi                       ! Total number of input param grid array points
  integer :: mo                       ! Total number of output LIS grid array points
  real,    allocatable  :: gi(:)      ! Input parameter 1d grid
  logical*1,allocatable :: li(:)      ! Input logical mask (to match gi)

  real                 :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! Output lis 1d grid
  logical*1            :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! Output logical mask (to match go)
  integer, allocatable :: lat_line(:,:)
  integer, allocatable :: lon_line(:,:)
  real,    allocatable :: read_inputparm(:,:)  ! Read input parameter
  integer, allocatable :: read_inputparm_int(:,:)  ! Read input parameter
  real,    allocatable :: read_parmarray(:)    ! Read input parameter
!_________________________________________________________________________________

   LDT_rc%nmaskpts = 0.
   localmask = 0.

   ! EMK...Special handling for UKMO_IGBP_Native_PFT land mask, which is in
   ! netCDF format.
   if( LDT_rc%mask_source(n) == "UKMO_IGBP_Native_PFT" ) then
      call read_maskfile_UKMO_IGBP_Native_PFT(n, vegtype, fgrd, localmask)
      return
   end if
   ! EMK...Special handling for UKMO_CAP_Netcdf land mask, which is in
   ! netCDF format.
   if( LDT_rc%mask_source(n) == "UKMO_CAP_Netcdf" ) then
      call read_maskfile_UKMO_CAP_Netcdf(n, vegtype, fgrd, localmask)
      return
   end if

!- Check for and open landmask file:
   inquire(file=trim(LDT_rc%mfile(n)), exist=file_exists)
   if( file_exists ) then 
      write(LDT_logunit,*) "[INFO] Reading landmask ",trim(LDT_rc%mfile(n)), & 
                           " (",LDT_localPet,")"
   else
      write(LDT_logunit,*) "[ERR] Landmask map: ",trim(LDT_rc%mfile(n))," does not exist."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif

   if(LDT_rc%mask_source(n) =="MCD12Q1") then
      call read_maskfile_MCD12Q1(n,vegtype, fgrd, localmask)
      return
   endif
   
! -------------------------------------------------------------------
!    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
! -------------------------------------------------------------------
!- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
!- Using Landcover grid description as default for now:
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, LDT_rc%mask_proj, LDT_rc%mask_gridDesc(n,:), &
                          glpnc, glpnr, subpnc, subpnr,  &
                          subparam_gridDesc, lat_line, lon_line )

!- Check parameter resolution against run domain:
   if( (subparam_gridDesc(9)-(LDT_rc%gridDesc(n,9)/LDT_rc%lis_map_resfactor(n)))>0.001 ) then
     if( LDT_rc%mask_gridtransform(n) == "none" ) then
       write(LDT_logunit,*) "[ERR] The landmask file resolution is set as less than"
       write(LDT_logunit,*) "  the LIS run domain resolution, but NO upscaling method selected."
       write(LDT_logunit,*) "  Please either read in a mask file with the same resolution as"
       write(LDT_logunit,*) "  the run domain or select: 'mode' or 'neighbor' (mode recommended)"
       write(LDT_logunit,*) " Stopping ..."
       call LDT_endrun
     endif
   endif

! -------------------------------------------------------------------

   rec_length = 0
   inquire( file=trim(LDT_rc%mfile(n)), size=file_size )
   if( abs(file_size) > 2000000000 ) then
      write(LDT_logunit,*) "[INFO] Reading a mask file > 2GB ..."
      rec_length = glpnc * 4
   else
      rec_length = 4
   endif

!== Open land/water mask file:
   ftn = LDT_getNextUnitNumber()
   open(ftn, file=trim(LDT_rc%mfile(n)),form='unformatted', &
        recl=rec_length, access='direct', iostat=ios1)

! == (1) READ IN MASK PARAMETER DATA: == 

!++ CLSM F2.5 ++
 ! Read in global mask field for current CLSM F2.5 needs:
 !   (Currently important for reading in other CLSM F2.5 parameters)

!   allocate(LDT_rc%global_mask(subpnc, subpnr))  ! Read in parallel section of file 
   allocate(LDT_rc%global_mask(glpnc, glpnr))    ! Read in full global mask file
   LDT_rc%global_mask = LDT_rc%udef

   if( LDT_rc%lsm == "CLSMF2.5") then
     line = 0
     do r = 1, glpnr
        do c = 1, glpnc
           line = line + 1
           read(ftn,rec=line)  LDT_rc%global_mask(c,r)
        enddo
     enddo

!- Generate total number of accounted mask points:
     do r = 1, glpnr
        do c = 1, glpnc
           if( LDT_rc%global_mask(c,r) >= 1 ) then
              LDT_rc%nmaskpts(n) = LDT_rc%nmaskpts(n) + 1
           endif
        end do
     end do

   endif
 !  ( to be removed later after full CLSM-F2.5 preprocesser
 !    is implemented into LDT)

! ++++++++++++++

 ! Parallel process set up:
   allocate( read_inputparm(subpnc, subpnr) )
   allocate( read_inputparm_int(subpnc, subpnr) )
   allocate( read_parmarray(glpnc) )
   read_inputparm = 0.
   read_inputparm_int = 0
   read_parmarray = 0.
   line = 0
   do r = 1, subpnr
    ! Read in global array:
      if( rec_length == glpnc*4 ) then
        read(ftn,rec=lat_line(1,r)) read_parmarray(:)
        do c = 1, subpnc
         ! Handle MODIS 44W landmask separately:
           if( LDT_rc%mask_source(n) == "MOD44W" ) then
             if( read_parmarray(lon_line(c,r)) > 0. ) then
               read_inputparm(c,r) = 0.
             elseif( read_parmarray(lon_line(c,r)) == 0. ) then
               read_inputparm(c,r) = 1.
             endif
           else  ! all other mask type files
             read_inputparm(c,r) = read_parmarray(lon_line(c,r))
           endif
        end do

    ! Read subsetted array:
      else
        if( LDT_rc%mask_source(n) == "HYMAP" ) then
          do c = 1, subpnc
             line = (lat_line(c,r)-1)*glpnc + lon_line(c,r)
             read(ftn,rec=line) read_inputparm_int(c,r)
             read_inputparm(c,r) = float(read_inputparm_int(c,r))
             if( read_inputparm(c,r) ==LDT_rc%udef.or. read_inputparm(c,r) ==0.) then
               read_inputparm(c,r) = 0.
             else
               read_inputparm(c,r) = 1.
             endif
          enddo
        else
          do c = 1, subpnc
             line = (lat_line(c,r)-1)*glpnc + lon_line(c,r) 
             read(ftn,rec=line) read_inputparm(c,r)
          enddo
        endif
      endif
   enddo
   deallocate( lat_line, lon_line )

! == (2) Downscale/Aggregate mask points, if option turned on == 

   mi = subpnr*subpnc
   mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
   allocate( gi(mi), li(mi) )
   gi = 0
   li = .false.
   lo1 = .false.

   select case ( LDT_rc%mask_gridtransform(n) )

   ! 'mode' -- Use for aggregating
   ! 'neighbor' -- Can be used for aggregating or downscaling
 
     case ( "mode", "neighbor" )  ! Used for aggregating 

       write(LDT_logunit,*) "[INFO] Readin Landmask file is being aggregated with"
       write(LDT_logunit,*) "   the '",trim(LDT_rc%mask_gridtransform(n)),"' option."

!       LDT_rc%global_mask = read_inputparm

    !- Assign 2-D array to 1-D for aggregation routines:
       i = 0
       do r = 1, subpnr
          do c = 1, subpnc;  i = i + 1
             gi(i) = read_inputparm(c,r)
             if( gi(i) .ne. LDT_rc%udef ) li(i) = .true.
          enddo
       enddo

    !- Transform parameter from original grid to LIS output grid:
       call LDT_transform_paramgrid(n, LDT_rc%mask_gridtransform(n), &
                subparam_gridDesc, mi, 1, gi, li, mo, go1, lo1 )
 
    !- Convert 1D to 2D grid output arrays:
       i = 0
       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n); i = i + 1
             if( go1(i) < 0.) then
               localmask(c,r) = LDT_rc%udef
             else
               localmask(c,r) = go1(i)
             end if
          enddo
       enddo

   case default ! 'none' -- no aggregation/downscaling of mask
      write(LDT_logunit,*) "[INFO] No upscaling/downscaling of Landmask file is performed."
      localmask = read_inputparm
!      LDT_rc%global_mask = localmask
   end select

   deallocate( read_inputparm )
   deallocate( read_parmarray )

! == (3) CHECK MASK FOR CONSISTENCY AGAINST VEGTYPE/TILE MAP: == 

   if( LDT_rc%inc_water_pts ) then

   !- NON-tiled option:
      select case( LDT_rc%lc_gridtransform(n) ) 

        case( "none", "mode" )

          cnt_0mask_1lc = 0.
          cnt_1mask_0lc = 0.
          do r = 1, LDT_rc%lnr(n)
             do c = 1, LDT_rc%lnc(n)
                if( vegtype(c,r) == LDT_rc%waterclass .and. &     ! vegtype=water,mask=land
                   localmask(c,r) /= 0 ) then
                  cnt_1mask_0lc = cnt_1mask_0lc + 1

                elseif( vegtype(c,r) /= LDT_rc%waterclass .and. & ! vegtype=land, mask=water
                        localmask(c,r) == 0 ) then
                  cnt_0mask_1lc = cnt_0mask_1lc + 1
                endif
             end do
          end do

     !- TILED option:
        case( "tile" )

          cnt_0mask_1lc = 0.
          cnt_1mask_0lc = 0.
          do r=1,LDT_rc%lnr(n)
             do c=1,LDT_rc%lnc(n)

             ! fgrd=water, mask=land
               if( fgrd(c,r,LDT_rc%waterclass) >= LDT_rc%gridcell_water_frac(n) .and. &  
                   localmask(c,r) == 1 ) then
                  cnt_1mask_0lc = cnt_1mask_0lc + 1
               endif
             ! fgrd=land, mask=water
               if( fgrd(c,r,LDT_rc%waterclass) < LDT_rc%gridcell_water_frac(n) .and. & 
                       localmask(c,r) == 0 ) then
                  cnt_0mask_1lc = cnt_0mask_1lc + 1
               endif

             end do
          end do
      end select

      if( cnt_1mask_0lc > 0. ) then
          write(LDT_logunit,*)" MISMATCH: Exists between the mask and landcover maps ..."
          write(LDT_logunit,*)" ... mismatch count: ", cnt_1mask_0lc
          write(LDT_logunit,*)" To make sure that landmask and landcover fields agree, select:"
          write(LDT_logunit,*)"  'neighbor' in the 'Landcover fill option:' entry ... "
      endif

!- Do not include water points in tiled land cover map:
   else
      cnt_0mask_1lc = 0.
      cnt_1mask_0lc = 0.
   !- TILED option:
      if ( LDT_rc%lc_gridtransform(n) == "tile" ) then
         do r=1,LDT_rc%lnr(n)
            do c=1,LDT_rc%lnc(n)
             ! fgrd=water, mask=land
!               if( sum(fgrd(c,r,1:LDT_rc%nt)) <= 0. .and. &
               if( sum(fgrd(c,r,1:LDT_rc%nt)) <= LDT_rc%gridcell_water_frac(n) .and. &
                   localmask(c,r) == 1 ) then
                  cnt_1mask_0lc = cnt_1mask_0lc + 1
               endif
            end do
         end do
      endif

      if( cnt_1mask_0lc > 0. ) then
          write(LDT_logunit,*)"[WARN] MISMATCH: Exists between the mask and landcover maps ..."
          write(LDT_logunit,*)" ... mismatch count: ", cnt_1mask_0lc
          write(LDT_logunit,*)" To make sure that landmask and landcover fields agree, select:"
          write(LDT_logunit,*)"  'neighbor' in the 'Landcover fill option:' entry ... "
!          call LDT_endrun
      endif

    endif !end inc_water_pts condition

#if 0
  !- Tiled option:
     if( LDT_rc%lc_gridtransform(n) == "tile" ) then
       do r = 1, LDT_rc%lnr(n)
          do c = 1, LDT_rc%lnc(n)
          !- If water is the dominant tile type:
             if( fgrd(c,r,LDT_rc%waterclass) >= LDT_rc%gridcell_water_frac(n) ) then  
                if( localmask(c,r) == 1. ) localmask(c,r) = 0.   ! Make sure about consistency
             end if
          end do
       end do
     end if
#endif

#if 0
!- Generate total number of accounted mask points:
! - Only used with CLSM F2.5 ...
! - Only include below if make global mask for parallel runs ...
   do r = 1,LDT_rc%lnr(n)
      do c = 1,LDT_rc%lnc(n)
         if( localmask(c,r) >= 1 ) then
            LDT_rc%nmaskpts(n) = LDT_rc%nmaskpts(n) + 1
         endif
      end do 
   end do
#endif

   call LDT_releaseUnitNumber(ftn)

end subroutine read_maskfile

