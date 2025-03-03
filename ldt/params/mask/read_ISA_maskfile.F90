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
! !ROUTINE: read_ISA_maskfile
!  \label{read_ISA_maskfile}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  01 June 2012: KR Arsenault; Restructured to simply read in a mask file
!  10 July 2013: KR Arsenault; Added for ISA landcover case
!
! !INTERFACE:
subroutine read_ISA_maskfile(n, fgrd, localmask )

! !USES:
  use LDT_coreMod,   only : LDT_rc, LDT_localPet
  use LDT_logMod,    only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_gridmappingMod

  implicit none

! !ARGUMENTS: 
  integer, intent(in)  :: n
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
  integer :: c, r, t, i, line
  integer :: glpnc, glpnr             ! Parameter (global) total columns and rows
  integer :: subpnc, subpnr           ! Parameter subsetted columns and rows
  real    :: subparam_gridDesc(20)    ! Input parameter grid desc array
  real    :: cnt_0mask_1lc
  real    :: cnt_1mask_0lc

  integer, allocatable :: lat_line(:,:)
  integer, allocatable :: lon_line(:,:)
  real,    allocatable :: read_inputparm(:,:)  ! Read input parameter
!_________________________________________________________________________________

   LDT_rc%nmaskpts = 0.
   localmask = 0.

!- Check for and open landmask file:
   inquire(file=trim(LDT_rc%mfile(n)), exist=file_exists)
   if( file_exists ) then 
      write(LDT_logunit,*)"[INFO] ISA mask -- Reading ",trim(LDT_rc%mfile(n)), & 
                          " (",LDT_localPet,")"
   else
      write(LDT_logunit,*) "[ERR] ISA Landmask map: ",trim(LDT_rc%mfile(n))," does not exist."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
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

   allocate( read_inputparm(subpnc, subpnr) )
   read_inputparm = 0.

! -------------------------------------------------------------------

!- Open land/water mask file:
   ftn = LDT_getNextUnitNumber()
   open(ftn, file=trim(LDT_rc%mfile(n)),form='unformatted', recl=4, &
           access='direct', convert='little_endian', iostat=ios1)

 ! Global mask field for current CLSM F2.5 needs:
   allocate(LDT_rc%global_mask(subpnc,subpnr))  !  Temporary ...
   LDT_rc%global_mask = LDT_rc%udef

! == (1) READ IN MASK PARAMETER DATA: == 
   line = 0
   do r = 1, subpnr
      do c = 1, subpnc
         line = (lat_line(c,r)-1)*glpnc + lon_line(c,r) 
         read(ftn,rec=line) read_inputparm(c,r)
      enddo
   enddo
   localmask = read_inputparm
   LDT_rc%global_mask = localmask

! == (2) CHECK MASK FOR CONSISTENCY AGAINST VEGTYPE/TILE MAP: == 

   if( LDT_rc%inc_water_pts ) then

   !- NON-tiled option:
     select case( LDT_rc%lc_gridtransform(n) )

    !- NON-tiled option:
       case( "none", "mode" )
    
         write(LDT_logunit,*)" ERR:  ISA landmask is not available for 2D "
         write(LDT_logunit,*)"       ISA landcover grid transforms (e.g., mode)."
         write(LDT_logunit,*)" Stopping ..."
         call LDT_endrun

#if 0
         cnt_0mask_1lc = 0.
         cnt_1mask_0lc = 0.
         do r=1,LDT_rc%lnr(n)
            do c=1,LDT_rc%lnc(n)
               if( vegtype(c,r) == LDT_rc%waterclass .and. &     ! vegtype=water,mask=land
                   localmask(c,r) /= 0 ) then
                  cnt_1mask_0lc = cnt_1mask_0lc + 1

               elseif( vegtype(c,r) /= LDT_rc%waterclass .and. & ! vegtype=land, mask=water
                       localmask(c,r) == 0 ) then
                  cnt_0mask_1lc = cnt_0mask_1lc + 1
               endif
            end do
         end do
#endif

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
             ! fgrd=land, mask=water
               elseif( fgrd(c,r,LDT_rc%waterclass) < LDT_rc%gridcell_water_frac(n) .and. & 
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
!          write(LDT_logunit,*)" Program stopping ..."
!          call LDT_endrun
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
               if( sum(fgrd(c,r,1:LDT_rc%nt)) <= 0. .and. &
                   localmask(c,r) == 1 ) then
                  cnt_1mask_0lc = cnt_1mask_0lc + 1
!                  write(450,*) "(fgrd=water,mask=land):",cnt_1mask_0lc, &
!                               sum(fgrd(c,r,1:LDT_rc%nt)), localmask(c,r)
               endif
            end do
         end do
      endif

      if( cnt_1mask_0lc > 0. ) then
          write(LDT_logunit,*)" MISMATCH: Exists between the mask and landcover maps ..."
          write(LDT_logunit,*)" ... mismatch count: ", cnt_1mask_0lc
          write(LDT_logunit,*)" To make sure that landmask and landcover fields agree, select:"
          write(LDT_logunit,*)"  'neighbor' in the 'Landcover fill option:' entry ... "
!          write(LDT_logunit,*)" Program stopping ..."
!          call LDT_endrun
      endif

    endif !end inc_water_pts condition

!- Generate total number of accounted mask points:
   do r = 1,LDT_rc%lnr(n)
      do c = 1,LDT_rc%lnc(n)
         if( localmask(c,r) >= 1 ) then
            LDT_rc%nmaskpts(n) = LDT_rc%nmaskpts(n) + 1
         endif
      end do 
   end do

   call LDT_releaseUnitNumber(ftn)

end subroutine read_ISA_maskfile
