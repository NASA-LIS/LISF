!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: MIRCA2000_crops_module
!  \label{MIRCA2000_crops_module}
!
! !REVISION HISTORY:
!  21  May 2019: H Beaudoing; Implemented MIRCA2000 crop map features
!                             The routine is based on Sarith Mahanama's 
!                             ReadProcess_MIRCA 
!                             section of preprocess_irrig.F90
!
! !INTERFACE:

module MIRCA2000_crops_module
! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_logunit, LDT_verify, LDT_endrun
  use LDT_gridmappingMod
  use LDT_fileIOMod
  use LDT_LSMCropModifier_Mod

 implicit none
 PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: read_MIRCA2000_croptype

CONTAINS

!BOP
!
! !ROUTINE: read_MIRCA2000_croptype
!  \label{read_MIRCA2000_croptype}
!
! !REVISION HISTORY:
!  21  May 2019: H Beaudoing; Implemented MIRCA2000 crop map features
!  21  Feb 2020: H Beaudoing; added crop calendar support
!
! !INTERFACE:
subroutine read_MIRCA2000_croptype(n,num_types,fgrd)
                                   

! !USES:
  use LDT_coreMod,     only : LDT_rc
  use LDT_logMod,      only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_verify, LDT_endrun
  use LDT_gridmappingMod    
  use LDT_fileIOMod
  use LDT_paramTileInputMod, only: param_index_fgrdcalc
  use LDT_LSMCropModifier_Mod

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: num_types
  real, intent(inout) :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_types)

!
! !DESCRIPTION:
!  This subroutine reads the MIRCA2000 Portmann et al. (2010) crop data 
!  and returns the distribution of vegetation in each grid cell, in a lat/lon
!  projection.  
!  The crop fraction inputs are avaiable at monthly. Currently, the monthly
!  information is condensed to annual values in the reading routine.
!
!  Crop planting and harvesting dates are returned by default using the 
!  monthly crop fractions for two seasons.  If user opted to use crop calendar
!  dates in ldt.config, LDT_LSMCrop_struc%planting and harvesting values are
!  updated.
!
!  Number of crops for MIRCA2000 can be specified in ldt.config 
!  (via LDT_rc%numcrop) to either 26 or 52 where:
!  LDT_rc%numcrop = 26 => irrigated crop tiles 
!  LDT_rc%numcrop = 52 => explicit irrigated + rainfed crop tiles
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!     index of nest
!   \item[num_types]
!     number of landcover and/or just crop types
!   \item[fgrd]
!     fraction of grid covered by each vegetation and crop type
!   \end{description}
!EOP      

   integer, parameter :: IN_NCROPS = 26 ! MIRCA crop types

! local variables
   character(len=140) :: tempfilei,tempfiler
   logical :: file_exists
   integer :: i, ii, j, t, c, r, mc
   integer :: err

   real      :: irrigated(LDT_rc%lnc(n),LDT_rc%lnr(n))
   real      :: rainfed(LDT_rc%lnc(n),LDT_rc%lnr(n))
   real      :: cropdom(LDT_rc%lnc(n),LDT_rc%lnr(n))
   real :: croptype_frac(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%numcrop(n)) !26or52
   real :: crop_array(LDT_rc%numcrop(n)) !26or52

   character(len=20), allocatable, dimension(:) :: croptype_array
!     fraction of cropland covered by irrigated crops and rainfed crops
   real    :: irrig_frac(LDT_rc%lnc(n),LDT_rc%lnr(n),IN_NCROPS) 
   real    :: rainfed_frac(LDT_rc%lnc(n),LDT_rc%lnr(n),IN_NCROPS)
   logical :: cropcaflag
!     planting/harvesting dates per season for irrigated and rainfed crops
   real, allocatable, dimension (:,:,:) :: irrigplant, irrigharvest, &
                                           rainfedplant, rainfedharvest

!==========================================================
   cropcaflag = LDT_LSMCrop_struc(n)%cropcalendar
   
   allocate (irrigplant(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_LSMCrop_struc(n)%multicroppingmax))
   allocate (irrigharvest(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_LSMCrop_struc(n)%multicroppingmax))
   allocate (rainfedplant(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_LSMCrop_struc(n)%multicroppingmax))
   allocate (rainfedharvest(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_LSMCrop_struc(n)%multicroppingmax))

   croptype_frac = 0.0
   irrig_frac = 0.0
   rainfed_frac = 0.0
   cropdom  = 0.0
   fgrd     = 0.0
   crop_array = 0.0

!- Check if MIRCA2000 Crop map directory exists:
   tempfilei = trim(LDT_LSMCrop_struc(n)%croptfile)//"01_irrigated_12.flt"
   inquire(file=tempfilei, exist=file_exists)
   if(.not.file_exists) then
      write(LDT_logunit,*) "[ERR] The MIRCA2000 crop type map diretory: ",&
                            trim(LDT_LSMCrop_struc(n)%croptfile)," is not correct. "
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif
   write(LDT_logunit,*)"[INFO] Reading MIRCA2000 crop file: ",trim(LDT_LSMCrop_struc(n)%croptfile)

!- Retrieve the croptype array for user-selected classification:
   allocate( croptype_array(LDT_rc%numcrop(n)), STAT=err )
   croptype_array = "none"
   call readcropinventory( n, LDT_rc%crop_classification(n), &
                           LDT_rc%numcrop(n), croptype_array ) 

!- Open each Crop-type file:'monthly_growing_area_grids/crop_*'
   do i = 1, IN_NCROPS

    tempfilei = trim(LDT_LSMCrop_struc(n)%croptfile)//trim(croptype_array(i))//"_irrigated_12.flt"
    tempfiler = trim(LDT_LSMCrop_struc(n)%croptfile)//trim(croptype_array(i))//"_rainfed_12.flt"

    write(LDT_logunit,*)"  Opening irrigated crop type file: ",trim(tempfilei)
    call readMIRCA2000Cropfiles( n, tempfilei, cropcaflag, irrigated,  &
                                 irrigplant, irrigharvest )
    write(LDT_logunit,*)"  Opening rainfed crop type file: ",trim(tempfiler)
    call readMIRCA2000Cropfiles( n, tempfiler, cropcaflag, rainfed, &
                                 rainfedplant, rainfedharvest )

    if ( LDT_rc%numcrop(n) .eq. 26 ) then    ! irrigated crops
      ! croptype_frac array (lnc,lnr,26)
      do r = 1, LDT_rc%lnr(n)
        do c = 1, LDT_rc%lnc(n)
           if ( irrigated(c,r).ne.LDT_rc%udef ) then
            croptype_frac(c,r,i) = irrigated(c,r) + rainfed(c,r)
            irrig_frac(c,r,i) = irrigated(c,r)
            rainfed_frac(c,r,i) = rainfed(c,r)
            if ( cropcaflag ) then  
             do mc = 1, LDT_LSMCrop_struc(n)%multicroppingmax 
               LDT_LSMCrop_struc(n)%plantday%value4d(c,r,i,mc) = irrigplant(c,r,mc)
               LDT_LSMCrop_struc(n)%harvestday%value4d(c,r,i,mc) = irrigharvest(c,r,mc)
             enddo
            endif
           else
            croptype_frac(c,r,i) = LDT_rc%udef
            irrig_frac(c,r,i) = LDT_rc%udef
            rainfed_frac(c,r,i) = LDT_rc%udef
            write(LDT_logunit,*) "[WARN] The MIRCA2000 irrigated and rainfed ",&
                                 " both should be missing at "c,r,irrigated(c,r),rainfed(c,r)  
            ! save rainfed plant/harvest days in case of later use
            if ( cropcaflag ) then  
             if ( rainfed(c,r).ne.LDT_rc%udef ) then
               do mc = 1, LDT_LSMCrop_struc(n)%multicroppingmax 
                 LDT_LSMCrop_struc(n)%plantday%value4d(c,r,i,mc) = rainfedplant(c,r,mc)
                 LDT_LSMCrop_struc(n)%harvestday%value4d(c,r,i,mc) = rainfedharvest(c,r,mc)
               enddo
             else
                LDT_LSMCrop_struc(n)%plantday%value4d(c,r,i,mc) = LDT_rc%udef
                LDT_LSMCrop_struc(n)%harvestday%value4d(c,r,i,mc) = LDT_rc%udef
             endif
            endif
           endif
        enddo
      enddo
     else if ( LDT_rc%numcrop(n) .eq. 52 ) then   ! irrigated & rainfed crops
     ! croptype_frac array (lnc,lnr,52)
      do r = 1, LDT_rc%lnr(n)
        do c = 1, LDT_rc%lnc(n)
           if ( irrigated(c,r).ne.LDT_rc%udef .and. &
                rainfed(c,r).ne.LDT_rc%udef) then
            if ( irrigated(c,r).ne.LDT_rc%udef ) then
             croptype_frac(c,r,i) = irrigated(c,r)
             irrig_frac(c,r,i) = irrigated(c,r)
             if ( cropcaflag ) then  
               do mc = 1, LDT_LSMCrop_struc(n)%multicroppingmax 
                 LDT_LSMCrop_struc(n)%plantday%value4d(c,r,i,mc) = irrigplant(c,r,mc)
                 LDT_LSMCrop_struc(n)%harvestday%value4d(c,r,i,mc) = irrigharvest(c,r,mc)
               enddo
             endif
            elseif( rainfed(c,r).ne.LDT_rc%udef ) then
             croptype_frac(c,r,i+26) = rainfed(c,r) 
             rainfed_frac(c,r,i) = rainfed(c,r) 
             if ( cropcaflag ) then  
               do mc = 1, LDT_LSMCrop_struc(n)%multicroppingmax 
                 LDT_LSMCrop_struc(n)%plantday%value4d(c,r,i+26,mc) = rainfedplant(c,r,mc)
                 LDT_LSMCrop_struc(n)%harvestday%value4d(c,r,i+26,mc) = rainfedharvest(c,r,mc)
               enddo
             endif
            endif
           else
            croptype_frac(c,r,i) = LDT_rc%udef
            croptype_frac(c,r,i+26) = LDT_rc%udef
            irrig_frac(c,r,i) = LDT_rc%udef
            rainfed_frac(c,r,i) = LDT_rc%udef
            if ( cropcaflag ) then  
              do mc = 1, LDT_LSMCrop_struc(n)%multicroppingmax 
                LDT_LSMCrop_struc(n)%plantday%value4d(c,r,i,mc) = LDT_rc%udef
                LDT_LSMCrop_struc(n)%plantday%value4d(c,r,i+26,mc) = LDT_rc%udef
                LDT_LSMCrop_struc(n)%harvestday%value4d(c,r,i,mc) = LDT_rc%udef
                LDT_LSMCrop_struc(n)%harvestday%value4d(c,r,i+26,mc) = LDT_rc%udef
              enddo
            endif
           endif
        enddo
      enddo
     else
      write(LDT_logunit,*) "number of crop type is not supported for MIRCA: ",&
                           LDT_rc%numcrop(n)
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
     endif  ! LDT_rc%numcrop(n)

   end do    ! i: End Crop Type File Read

!- Write final gridcell fractions:
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n)
!- Find dominant crop type out of inventory list:
         crop_array(:) = croptype_frac(c,r,:)
         cropdom(c,r) = maxloc(crop_array,1)

         if( maxval(croptype_frac(c,r,:)) == LDT_rc%udef ) then
            fgrd(c,r,:) = LDT_rc%udef
            LDT_LSMCrop_struc(n)%irrigcrop%value(c,r,:) = LDT_rc%udef
            LDT_LSMCrop_struc(n)%rainfedcrop%value(c,r,:) = LDT_rc%udef
         else
           do i = 1, LDT_rc%numcrop(n)
            fgrd(c,r,i) = croptype_frac(c,r,i)
            if ( LDT_rc%numcrop(n) .eq. 26 ) then 
              LDT_LSMCrop_struc(n)%irrigcrop%value(c,r,i) = irrig_frac(c,r,i)
              LDT_LSMCrop_struc(n)%rainfedcrop%value(c,r,i) = rainfed_frac(c,r,i)
            else if ( LDT_rc%numcrop(n) .eq. 52 ) then 
              if ( i .le. 26 ) then 
                LDT_LSMCrop_struc(n)%irrigcrop%value(c,r,i) = irrig_frac(c,r,i)
              else
                LDT_LSMCrop_struc(n)%rainfedcrop%value(c,r,i-26) = rainfed_frac(c,r,i-26)
              endif
            endif  ! LDT_rc%numcrop(n)
           enddo
         endif
         ! FINAL CROPCALENDAR VARIABLES  remove 998 and 999....
         if ( cropcaflag ) then  
           do mc = 1, LDT_LSMCrop_struc(n)%multicroppingmax 
             do i = 1, LDT_rc%numcrop(n)
              if ( LDT_LSMCrop_struc(n)%plantday%value4d(c,r,i,mc) == 998. ) then
                LDT_LSMCrop_struc(n)%plantday%value4d(c,r,i,mc) = LDT_rc%udef
              elseif ( LDT_LSMCrop_struc(n)%plantday%value4d(c,r,i,mc) == 999. ) then
              endif
              if ( LDT_LSMCrop_struc(n)%harvestday%value4d(c,r,i,mc) == 998. ) then
                LDT_LSMCrop_struc(n)%harvestday%value4d(c,r,i,mc) = LDT_rc%udef
              elseif ( LDT_LSMCrop_struc(n)%harvestday%value4d(c,r,i,mc) == 999. ) then
              endif
             enddo
           enddo
         endif
      enddo
   enddo
   deallocate( croptype_array )
  write(LDT_logunit,*) "[INFO] Done reading 'MIRCA2000' Crop map type. "

!- CROP CALENDAR ----
   if(allocated(irrigplant))  deallocate (irrigplant)
   if(allocated(irrigharvest))  deallocate (irrigharvest)
   if(allocated(rainfedplant))  deallocate (rainfedplant)
   if(allocated(rainfedharvest))  deallocate (rainfedharvest)


end subroutine read_MIRCA2000_croptype


 subroutine readMIRCA2000Cropfiles( n, tempfile, cropcaflag, croptypefrac, &
                                    plantday, harvestday )
                                  
! !USES:
    use LDT_coreMod
    use LDT_gridmappingMod
    use LDT_paramTileInputMod
    use LDT_LSMCropModifier_Mod
    use LDT_fileIOMod
    use LDT_logMod

    implicit none
! !ARGUMENTS:
   integer,       intent(in) :: n
   character(140),intent(in) :: tempfile
   logical,       intent(in) :: cropcaflag
   real,         intent(out) :: croptypefrac(LDT_rc%lnc(n),LDT_rc%lnr(n))
   real,intent(out) :: plantday(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_LSMCrop_struc(n)%multicroppingmax) ! lon x lat x seasons
   real,intent(out) :: harvestday(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_LSMCrop_struc(n)%multicroppingmax) ! lon x lat x seasons

! MIRCA2000 monthly_growing_area_grids crop area map (5 arcmin):
   integer, parameter :: IN_cols = 4320
   integer, parameter :: IN_rows = 2160
   real,    parameter :: IN_xres = 1.0/12.0 ! or 360deg/4320 col points
   real,    parameter :: IN_yres = 1.0/12.0 ! or 180deg/2160 row points
   real, parameter    :: PI = 3.1415926535898
   real, parameter    :: radius = 6371.0 ! [km]

   integer :: glpnc, glpnr             ! Parameter (global) total columns and rows
   integer :: subpnc, subpnr           ! Parameter subsetted columns and rows
   real    :: param_gridDesc(20)       ! Input parameter grid desc fgrd
   real    :: subparam_gridDesc(20)    ! Input parameter grid desc array
   integer, allocatable  :: lat_line(:,:), lon_line(:,:)

   logical :: file_exists
   integer :: ftn, ierr, nrec, ios 
   integer :: i, j, t, c, r, m, nc, nr
   integer :: mi                       ! Total number of input param grid array points
   integer :: mo                       ! Total number of output LIS grid array points
   real,    allocatable  :: gi1(:)     ! Input parameter 1d grid
   logical*1,allocatable :: li1(:)     ! Input logical mask (to match gi)
   real      :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))            ! Output lis 1d grid
   logical*1 :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))            ! Output logical mask (to match go)

   real, allocatable :: var_read(:)           ! Read in parameter
   real, allocatable :: var_in(:,:,:)         ! Read in parameter
   real, allocatable :: read_crop_ann(:,:)    ! Read input annual mean

   real              :: latc, lonc, area, d2r
! crop calendar varibles
   integer :: mc   ! crop seasons-- hard-coded to two for now
   integer :: mci
   real,    allocatable  :: gi1p(:)   ! Input parameter 1d grid
   logical*1,allocatable :: li1p(:)   ! Input logical mask (to match gi)
   real      :: var_out(LDT_rc%lnc(n),LDT_rc%lnr(n),12)  ! 12-mon croptype frac
   real      :: go1p(LDT_rc%lnc(n)*LDT_rc%lnr(n))        ! Output lis 1d grid
   logical*1 :: lo1p(LDT_rc%lnc(n)*LDT_rc%lnr(n))        ! Output logical mask (to match go)

!_____________________________________________________________________________
!- Set parameter grid array inputs:
   LDT_LSMCrop_struc(n)%crop_proj = "latlon"
   param_gridDesc(1)  = 0.          ! Latlon
   param_gridDesc(2)  = IN_cols
   param_gridDesc(3)  = IN_rows
   param_gridDesc(4)  = -90.0  + (IN_yres/2) ! LL lat
   param_gridDesc(5)  = -180.0 + (IN_xres/2) ! LL lon
   param_gridDesc(6)  = 128
   param_gridDesc(7)  =  90.0 - (IN_yres/2)  ! UR lat
   param_gridDesc(8)  = 180.0 - (IN_xres/2)  ! UR lon
   param_gridDesc(9)  = IN_yres     ! dy: 0.083333
   param_gridDesc(10) = IN_xres     ! dx: 0.083333
   param_gridDesc(20) = 64

!- Map Parameter Grid Info to LIS Target Grid/Projection/Subset Info:
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, LDT_LSMCrop_struc(n)%crop_proj, &
        param_gridDesc, glpnc, glpnr, subpnc, subpnr, &
        subparam_gridDesc, lat_line, lon_line )


!- Perform check on grid and projection choices selected:
   if( LDT_rc%lis_map_proj(n) == "latlon"  ) then 
     if( param_gridDesc(10) .gt. (LDT_rc%gridDesc(n,9)/LDT_rc%lis_map_resfactor(n))) then
        write(LDT_logunit,*) "[ERR] 'Native' MIRCA2000 crop files have been selected, "
        write(LDT_logunit,*) "    with a resolution (0.0833deg), but the LIS run domain resolution"
        write(LDT_logunit,*) "    selected is finer than that. Downscaling crop type information is "
        write(LDT_logunit,*) "    currently not supported in LDT."
        write(LDT_logunit,*) "Program stopping ..."
        call LDT_endrun
     endif
   else
    write(LDT_logunit,*) "[ERR] projection other than latlon is not supported "
    write(LDT_logunit,*) " due to only UpscaleByAveraging transformation "
    write(LDT_logunit,*) " is appropriate for MIRCA2000 crop map "
    write(LDT_logunit,*) "Program stopping ..."
    call LDT_endrun
   endif

   if( allocated(var_read) ) then
       deallocate(var_read)
   endif
   mc = LDT_LSMCrop_struc(n)%multicroppingmax
   allocate( var_read(IN_cols*12) )
   allocate( var_in(subpnc,subpnr,12) )
!
!  READ IN MIRCA2000 CROP AREA FIELDS 
!
   ftn = LDT_getNextUnitNumber()
   open(ftn, file=trim(tempfile),status='old',form='unformatted',&
        access='direct',convert='little_endian',recl=4*IN_cols*12,iostat=ios)
!- Reverse-Y and read in 12-month values and subset
   nrec = 0
   do nr=subpnr, 1, -1
    nrec = IN_rows - lat_line(1,nr) + 1
    read(ftn,rec=nrec) var_read
    do nc = 1, subpnc
      do m = 1, 12
        var_in(nc,nr,m) = var_read((lon_line(nc,1)-1)*12 + m)
      end do
    end do
   end do

   deallocate( var_read )
   call LDT_releaseUnitNumber(ftn)

   write(LDT_logunit,*) "[INFO] Done reading MIRCA2000 crop type "


! -------------------------------------------------------------------
!     Compute annual maximum of the 12 month data
! -------------------------------------------------------------------
    allocate(read_crop_ann(subpnc,subpnr))
     do nr = 1, subpnr
       do nc = 1, subpnc
         read_crop_ann(nc,nr) = maxval(var_in(nc,nr,:))
         ! assign missing to non-existing crop type
         if( read_crop_ann(nc,nr) .le. 0. ) read_crop_ann(nc,nr) = LDT_rc%udef
       end do
     end do


! -------------------------------------------------------------------
!     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
! -------------------------------------------------------------------
    mi = subpnc*subpnr
    allocate( gi1(mi), li1(mi) )
    gi1 = LDT_rc%udef;  li1 = .false.
    mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
    lo1 = .false.;  go1 = 0.

 !- Assign 2-D array to 1-D for aggregation routines:
    i = 0
    do nr = 1, subpnr
       do nc = 1, subpnc
          i = i + 1
          gi1(i) = read_crop_ann(nc,nr)
          if( gi1(i) < 0. ) gi1(i) = LDT_rc%udef
          if( gi1(i) .ne. LDT_rc%udef )  li1(i) = .true.
       enddo
    enddo
    deallocate( read_crop_ann )

 !- Upscale routine:
    select case ( LDT_LSMCrop_struc(n)%crop_gridtransform )

   !- Single layer transformation:
      case( "none", "average", "neighbor", "bilinear", "budget-bilinear" )

       !- Transform parameter from original grid to LIS output grid:
          call LDT_transform_paramgrid(n, LDT_LSMCrop_struc(n)%crop_gridtransform, &
                   subparam_gridDesc, mi, 1, gi1, li1, mo, go1, lo1 )

       !- Convert 1D to 2D grid arrays and units from area (Ha) to fraction
          d2r = PI/180.
          j = 0
          do r = 1, LDT_rc%lnr(n)
             latc = (r-1) * LDT_rc%gridDesc(n,9) + LDT_rc%gridDesc(n,4)
             area = (sin(d2r*(latc+0.5*LDT_rc%gridDesc(n,9))) - sin(d2r*(latc-0.5*LDT_rc%gridDesc(n,9))))*(LDT_rc%gridDesc(n,10)*d2r)
             area = 100. * area * radius * radius ! in hectares as in MIRCA
             do c = 1, LDT_rc%lnc(n)
                j = j + 1
                if ( go1(j) < 0. ) then
                 croptypefrac(c,r) = 0.0
                else
                 croptypefrac(c,r) = go1(j)/area
                endif
              ! Make sure all negative values are set to universal undefined value:
                if( croptypefrac(c,r) < 0.  ) croptypefrac(c,r) = LDT_rc%udef
                if( croptypefrac(c,r) > 100.) croptypefrac(c,r) = LDT_rc%udef
             enddo
          enddo
! HKB: save monthly fraction for CROPCALENDAR in_plantday/in_harvestday 
!      upscale/downscale using the same options as in croptypefrac
      !- Assign 2-D array to 1-D for aggregation routines:
      allocate( gi1p(mi), li1p(mi) )
      do m = 1, 12
        lo1p = .false.;  go1p = 0.
        i = 0
        do nr = 1, subpnr
           do nc = 1, subpnc
              i = i + 1
              gi1p(i) = var_in(nc,nr,m)
              if( gi1p(i) < 0. ) gi1p(i) = LDT_rc%udef
              if( gi1p(i) .ne. LDT_rc%udef )  li1p(i) = .true.
           enddo
        enddo
        !- Transform parameter from original grid to LIS output grid:
        call LDT_transform_paramgrid(n, LDT_LSMCrop_struc(n)%crop_gridtransform, &
                  subparam_gridDesc, mi, 1, gi1p, li1p, mo, go1p, lo1p )

        !- Convert 1D to 2D grid arrays and units from area (Ha) to fraction
        j = 0
        do r = 1, LDT_rc%lnr(n)
          latc = (r-1) * LDT_rc%gridDesc(n,9) + LDT_rc%gridDesc(n,4)
          area = (sin(d2r*(latc+0.5*LDT_rc%gridDesc(n,9))) - sin(d2r*(latc-0.5*LDT_rc%gridDesc(n,9))))*(LDT_rc%gridDesc(n,10)*d2r)
          area = 100. * area * radius * radius ! in hectares as in MIRCA
          do c = 1, LDT_rc%lnc(n)
             j = j + 1
             if ( go1p(j) < 0. ) then
               var_out(c,r,m) = LDT_rc%udef
             else
               var_out(c,r,m) = go1p(j)/area
             endif
             if ( var_out(c,r,m) < 0. ) var_out(c,r,m) = LDT_rc%udef
             if ( var_out(c,r,m) > 100. ) var_out(c,r,m) = LDT_rc%udef
          enddo  ! c 
        enddo  ! r 
      enddo  ! m
      deallocate(gi1p, li1p)

     case default
       write(*,*)"[ERR] Other spatial grid transformations are not currently supported"
       write(*,*)"   for the tiled 'MIRCA2000' crop type maps. Please select either:"
       write(*,*)"  -- neighbor, bilinear, budget-bilinear (to downscale)"
       write(*,*)"  -- average (to upscale) or none"
       write(*,*)" Stopping program ..."
       call LDT_endrun

    end select  ! End vegtype/cnt aggregation method
! -------------------------------------------------------------------
!HKB: compute crop planting/harvesting days here with "var_out"
!     (lnc,lnr,12), per crop with 2 seasons
! -------------------------------------------------------------------
    if ( cropcaflag ) then
     call getplantharvest2mc(mc,LDT_rc%lnc(n),LDT_rc%lnr(n),var_out, &
                            plantday,harvestday)
    else
     plantday = LDT_rc%udef
     harvestday = LDT_rc%udef
    endif
    deallocate( gi1, li1 )
    deallocate( lat_line, lon_line )
    deallocate( var_in )


 end subroutine readMIRCA2000Cropfiles
 
 subroutine getplantharvest2mc(mc,nx,ny,var_out,plantday,harvestday)
                                 
! This subroutine retrieves Crop planting/Harvesting days based on the 12-mon
! crop fraction in the MIRCA data. Two seasons of multiple cropping is 
! accounted for in this version.
! This routine provides initial estimate at native resolution.  The dates will
! be aggregated to model resolution then needs further adjustment when 
! irrigation fraction and crop fractions are merged.
!
! Code is adopted from preprocess_irrig.F90, originally written by 
! Sarith Mahanama
!
  implicit none
  integer,       intent(in) :: mc ! crop seasons
  integer,       intent(in) :: nx,ny
  real, dimension(nx,ny,12),intent(in) :: var_out
  real, dimension(nx,ny,mc),intent(out) :: plantday 
  real, dimension(nx,ny,mc),intent(out) :: harvestday

  integer :: r,c,m
  integer :: nc, day1, dayL, day1_2, dayL_2
  integer, dimension (12) :: fmonth, fmonth2, fmonth3
  integer, dimension (12)  :: DOY_MidMonth, DOY_BegMonth, DOY_EndMonth
  integer, allocatable , dimension (:) :: crop_mons
  logical, dimension (4)  :: found = .false.

  data DOY_BegMonth / 1, 32, 60,  91, 121, 152, 182, 213, 244, 274, 305, 335/
  data DOY_MidMonth /15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349/
  data DOY_EndMonth /31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 366/

  plantday = 998.
  harvestday = 998.

  do r = 1, ny
    do c = 1, nx
     if ( maxval(var_out(c,r,:)) > 0.0 ) then
       forall (m=1:12) fmonth3(m) = m
       fmonth  = 0
       fmonth2 = 0
       nc = count (var_out (c,r,:) > 0.)
       if(nc > 0) then

         if(nc == 12) then
         ! year around

          day1 = 1
          dayL = 366
          day1_2 = 998
          dayL_2 = 998

         ! print '(a20, 2i4)', 'INPUT nc=12 : at ',c,r
         ! print *, '..................................'
         else

          fmonth  = 0
          fmonth2 = 0
          day1 = 998
          dayL = 998
          day1_2 = 998
          dayL_2 = 998

          forall (m=1:12) fmonth(m) = ceiling (var_out (c,r,m))
         ! print '(a15, 12i2, a4, 2i4)', 'INPUT fmonth : ',  fmonth, ' at ',c,r
         ! print *, '..................................'

          fmonth2(1) = 1
          do m = 2,12
             if(fmonth(m) == fmonth(m-1)) then
                fmonth2 (m) = fmonth2(m-1)
             else
                fmonth2 (m) = fmonth2(m-1) + 1
             endif
          end do

          if(maxval (fmonth2) > 3) then
          ! This crop grows in 2 seasons
          ! ............................
            allocate (crop_mons (1:nc))
            crop_mons = pack(fmonth3, mask = (fmonth > 0))
            found = .false.

            if(fmonth(1) == 1) then
                if(fmonth(12) == 0) then
                   ! Season begins on Jan 1
                   day1   = DOY_BegMonth(crop_mons(1))
                   found(1) = .true.
                   do m = 1, nc-1
                      if((crop_mons(m+1) - crop_mons(m)) > 1) then
                         dayL   = DOY_EndMonth(crop_mons(m))
                         day1_2   = DOY_BegMonth(crop_mons(m+1))
                         dayL_2   = DOY_EndMonth(crop_mons(nc))
                         found(2) = .true.
                         found(3) = .true.
                         found(4) = .true.
                         exit
                      endif
                   enddo

                else

                   ! season one begins in the fall
                   do m = 1, nc-1
                      if((crop_mons(m+1) - crop_mons(m)) > 1) then
                         if(.not.found(2)) then
                            dayL   = DOY_EndMonth(crop_mons(m))
                            day1_2 = DOY_BegMonth(crop_mons(m+1))
                            found(2) = .true.
                            found(3) = .true.
                         elseif (.not.found(4)) then
                            found(4) = .true.
                            found(1) = .true.
                            dayL_2 = DOY_EndMonth(crop_mons(m))
                            day1   = DOY_BegMonth(crop_mons(m+1))
                         endif
                      endif
                   end do
                endif
            else

                ! season 1 brings in the spring
                day1   = DOY_BegMonth(crop_mons(1))
                found(1) = .true.
                do m = 1, nc-1
                   if((crop_mons(m+1) - crop_mons(m)) > 1) then
                      dayL   = DOY_EndMonth(crop_mons(m))
                      day1_2   = DOY_BegMonth(crop_mons(m+1))
                      dayL_2   = DOY_EndMonth(crop_mons(nc))
                      found(2) = .true.
                      found(3) = .true.
                      found(4) = .true.
                      exit
                   endif
                end do

            endif
            deallocate (crop_mons)

          else
          ! Single crop season
          ! ..................

             if((fmonth(1) == 0).and.(fmonth(12) == 0)) then
                day1 = DOY_BegMonth (maxloc(fmonth, 1))
                dayL = DOY_EndMonth (maxloc(fmonth2, 1)-1)
             else
                if((fmonth(1) == 1).and.(fmonth(12) == 1)) then
                   day1 = DOY_BegMonth (maxloc(fmonth2, 1,mask=(fmonth2 > 2)))
                   dayL = DOY_EndMonth (maxloc(fmonth2, 1,mask=(fmonth2 == 2))-1)
                endif
                if((fmonth(1) == 0).and.(fmonth(12) == 1)) then
                   day1 = DOY_BegMonth (maxloc(fmonth2, 1))
                   dayL = DOY_EndMonth (12)
                endif
                if((fmonth(1) == 1).and.(fmonth(12) == 0)) then
                   day1 = DOY_BegMonth (1)
                   dayL = DOY_EndMonth (maxloc(fmonth2, 1)-1)
                endif
             endif

          endif   ! maxval(fmonth2) > 3
         endif   ! nc ==12

       plantday(c,r,1) = float(day1)
       plantday(c,r,2) = float(day1_2)
       harvestday(c,r,1) = float(dayL)
       harvestday(c,r,2) = float(dayL_2)
       !print  '(a18, 4i4)', 'Plant & Harvest : ', day1, dayL, day1_2, dayL_2

       endif  ! nc > 0
     endif  ! var_in > 0
    end do ! c
  end do ! r

 end subroutine getplantharvest2mc

end module MIRCA2000_crops_module
