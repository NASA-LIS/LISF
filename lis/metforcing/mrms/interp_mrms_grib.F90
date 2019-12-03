!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! Jonathan Case; 13 Feb 2015: revised stg4 code for MRMS QPE.
! Jessica Erlingis; 5 September 2017: revised for MRMS operational products
! 
! !ROUTINE: interp_mrms_grib
! \label{interp_mrms_grib}
!
! !INTERFACE: 
subroutine interp_mrms_grib (n, findex, nmrms, ppt_field, lb, lis_gds, &
                        nc, nr, varfield )

! !USES:  
  use mrms_grib_forcingMod, only : mrms_grib_struc
  use LIS_coreMod, only: LIS_rc, LIS_domain

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n 
  integer, intent(in)    :: findex
  integer                :: nc      ! Number of columns (in the E-W dimension) in the LIS grid
  integer                :: nr      ! Number of rows (in the N-S dimension) in the LIS grid
  integer                :: nmrms   ! Number of points in original MRMS grid
  real                   :: lis_gds(50)
  real                   :: ppt_field(nmrms)
  logical*1              :: lb(nmrms)
  real, dimension(nc,nr) :: varfield

! !DESCRIPTION:
!   This subroutine interpolates a given MRMS field 
!   to the LIS grid. 
!  The arguments are: 
!  \begin{description}
! \item[n]
!  index of the nest
! \item[nmrms]
!  number of elements in the input grid
! \item[ppt\_field]
!  input data array to be interpolated
! \item[lb]
!  input bitmap
! \item[lis\_gds]
!  array description of the LIS grid
! \item[nc]
!  number of columns (in the east-west dimension) in the LIS grid
! \item[nr]
!  number of rows (in the north-south dimension) in the LIS grid
! \item[varfield]
!  output interpolated field
!  \end{description} 
!
!  The routines invoked are: 
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp})\\
!    Spatially interpolate the forcing data using bilinear interpolation, or
!  \item[conserv\_interp](\ref{conserv_interp})\\
!    spatially interpolate the forcing data using conservative (budget bil.) interpolation
! \end{description}
!EOP

  integer :: iret
  integer :: mo
  integer :: count1, i, j

  real    :: lis1d(nc*nr)   ! LIS Grid  
  logical*1 :: lo(nc*nr)    ! LIS Grid

! J.Case -- local variables
  real         :: gridDesci(LIS_rc%nnest, 50)
  character*50 :: dom
  real         :: res

!=== End variable declarations

!--------------------------------------------------------------------
! NOTE:: Recommended to use budget bilinear for precip forcing fields
!--------------------------------------------------------------------
  mrms_grib_struc(n)%mi = nmrms
  mo = nc * nr   ! LIS Grid

  lo = .true.    !-Initialize output bitmap
  !lb = .true.   ! Set in read_mrms

! define MRMS grid attributes
  dom=LIS_rc%lis_map_proj
  res=LIS_rc%gridDesc(n,9)
  gridDesci(n,:) = 0.0
  gridDesci(n,1) = 0                       ! Projection type (Lat/Lon)
  gridDesci(n,2) = mrms_grib_struc(n)%ncol ! X-dir amount of points
  gridDesci(n,3) = mrms_grib_struc(n)%nrow ! y-dir amount of points
  gridDesci(n,4) =   54.995                ! Starting latitude point
  gridDesci(n,5) = -129.995                ! Starting longitude point
  gridDesci(n,6) = 128                     ! (not used)
  gridDesci(n,7) =   20.005                ! Ending latitude point
  gridDesci(n,8) =  -60.005                ! Ending longitude point 
  gridDesci(n,9) =  0.01                   ! spatial resolution in W-E dirn (deg)
  gridDesci(n,10) = 0.01                   ! spatial resolution in S-N dirn (deg)
  gridDesci(n,20) = 64                     ! N-S ordering (number divisible by 32; same as in NLDAS2)

!-- Interpolate to LIS grid

! J.Case (2/13/2015) -- Use interpolation only if LIS resolution is < MRMS resolution.
!                       Otherwise, use upscaling by default.
  if ( (((dom.eq."latlon").or.(dom.eq."gaussian")) .and. (res.eq.0.01)) .or. & ! JE add equal case
       (((dom.eq."mercator").or.(dom.eq."lambert").or.                       &
         (dom.eq."polar").or.(dom.eq."UTM")) .and. (res.eq.1.0)) ) then

         call neighbor_interp(lis_gds, lb, ppt_field, lo, &
            lis1d,mrms_grib_struc(n)%mi, mo, &
            LIS_domain(n)%lat, LIS_domain(n)%lon, &
            mrms_grib_struc(n)%n113,LIS_rc%udef,iret)

  elseif ( (((dom.eq."latlon").or.(dom.eq."gaussian")) .and. (res.lt.0.01)) .or. & ! JE LE
       (((dom.eq."mercator").or.(dom.eq."lambert").or.                       &
         (dom.eq."polar").or.(dom.eq."UTM")) .and. (res.lt.1.0)) ) then

    if (trim(LIS_rc%met_interp(findex)) .eq. "bilinear") then 
! == BILINEAR INTERPOLATION ==== 
       call bilinear_interp( lis_gds,  lb, ppt_field, lo, &
            lis1d, mrms_grib_struc(n)%mi, mo, &
            LIS_domain(n)%lat, LIS_domain(n)%lon, &          !JE switch from mrms_grib_struc(n)%rl*1 to LIS_domain
            mrms_grib_struc(n)%w111, mrms_grib_struc(n)%w121, &
            mrms_grib_struc(n)%w211, mrms_grib_struc(n)%w221, &
            mrms_grib_struc(n)%n111, mrms_grib_struc(n)%n121, &
            mrms_grib_struc(n)%n211, mrms_grib_struc(n)%n221, LIS_rc%udef, iret)
    elseif (trim(LIS_rc%met_interp(findex)) .eq. "budget-bilinear") then 
! == BUDGET BILINEAR INTERPOLATION ==== 
       call conserv_interp( lis_gds,  lb, ppt_field,  lo, &
            lis1d, mrms_grib_struc(n)%mi, mo, & 
            LIS_domain(n)%lat, LIS_domain(n)%lon,&  !JE switch from mrms_grib_struc(n)%rl*2 to LIS_domain
            mrms_grib_struc(n)%w112, mrms_grib_struc(n)%w122,&
            mrms_grib_struc(n)%w212, mrms_grib_struc(n)%w222,&
            mrms_grib_struc(n)%n112, mrms_grib_struc(n)%n122,&
            mrms_grib_struc(n)%n212, mrms_grib_struc(n)%n222,LIS_rc%udef,iret)
    endif

  else !! use upscaling
! === UPSCALING ==== 
     call upscaleByAveraging( nmrms, mo, LIS_rc%udef, mrms_grib_struc(n)%n11, &
          lb, ppt_field, lo, lis1d)

  endif !! dom/res check
! J.Case (end interpolation mods)

!--------------------------------------------------------------------    
! Create 2D array (on LIS LIS_domain) for main program. 
!--------------------------------------------------------------------    
   count1 = 0
   do j = 1, nr
      do i = 1, nc    
         varfield(i,j) = lis1d(i+count1)
      enddo
      count1 = count1 + nc
   enddo


end subroutine interp_mrms_grib
