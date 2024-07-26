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
! !ROUTINE: interp_AWRAL
! \label{interp_AWRAL}
!
! !INTERFACE: 
subroutine interp_AWRAL (n, findex, nAWRAL, ppt_field, lb, lis_gds, &
                        nc, nr, varfield )

! !USES:  
  use AWRAL_forcingMod, only : AWRAL_struc
  use LIS_coreMod

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n 
  integer, intent(in)    :: findex
  integer                :: nc      ! Number of columns (in the E-W dimension) in the LIS grid
  integer                :: nr      ! Number of rows (in the N-S dimension) in the LIS grid
  integer                :: nAWRAL   ! Number of points in original STAGE IV grid
  real                   :: lis_gds(50)
  real                   :: ppt_field(nAWRAL)
  logical*1              :: lb(nAWRAL)
  real, dimension(nc,nr) :: varfield

! !DESCRIPTION:
!   This subroutine interpolates a given STAGE4 field 
!   to the LIS grid. 
!  The arguments are: 
!  \begin{description}
! \item[n]
!  index of the nest
! \item[nAWRAL]
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
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    Spatially interpolate the forcing data using bilinear interpolation, or
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative (budget bil.) interpolation
! \end{description}
!EOP

  integer :: iret
  integer :: mo
  integer :: count1, i, j

  real    :: lis1d(nc*nr)   ! LIS Grid  
  logical*1 :: lo(nc*nr)    ! LIS Grid

!=== End variable declarations

!--------------------------------------------------------------------
! NOTE:: Recommended to use budget bilinear for precip forcing fields
!--------------------------------------------------------------------
  AWRAL_struc(n)%mi = nAWRAL
  mo = nc * nr   ! LIS Grid

  lo = .true.    !-Initialize output bitmap

  if(AWRAL_struc(n)%interp_flag) then 

    if (trim(LIS_rc%met_interp(findex)) .eq. "bilinear") then 
       ! == BILINEAR INTERPOLATION ==== 
       call bilinear_interp( lis_gds,  lb, ppt_field, lo, &
            lis1d, AWRAL_struc(n)%mi, mo, &
            LIS_domain(n)%lat, LIS_domain(n)%lon,&
            AWRAL_struc(n)%w111, AWRAL_struc(n)%w121, &
            AWRAL_struc(n)%w211, AWRAL_struc(n)%w221, &
            AWRAL_struc(n)%n111, AWRAL_struc(n)%n121, &
            AWRAL_struc(n)%n211, AWRAL_struc(n)%n221, LIS_rc%udef, iret)

    elseif (trim(LIS_rc%met_interp(findex)) .eq. "budget-bilinear") then 
       ! == BUDGET BILINEAR INTERPOLATION ==== 
       call conserv_interp( lis_gds,  lb, ppt_field,  lo, &
            lis1d, AWRAL_struc(n)%mi, mo, & 
            LIS_domain(n)%lat, LIS_domain(n)%lon,&
            AWRAL_struc(n)%w112, AWRAL_struc(n)%w122,&
            AWRAL_struc(n)%w212, AWRAL_struc(n)%w222,&
            AWRAL_struc(n)%n112, AWRAL_struc(n)%n122,&
            AWRAL_struc(n)%n212, AWRAL_struc(n)%n222,LIS_rc%udef,iret)
    endif
 else
    call upscaleByAveraging(AWRAL_struc(n)%mi, mo, LIS_rc%udef, &
         AWRAL_struc(n)%n111, lb, ppt_field, lo, lis1d)
 endif
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


end subroutine interp_AWRAL

