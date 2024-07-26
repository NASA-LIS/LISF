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
! !ROUTINE: interp_stg4
! \label{interp_stg4}
!
! !INTERFACE: 
subroutine interp_stg4 (n, findex, nstg4, ppt_field, lb, lis_gds, &
                        nc, nr, varfield )

! !USES:  
  use stg4_forcingMod, only : stg4_struc
  use LIS_coreMod, only: LIS_rc, LIS_domain

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n 
  integer, intent(in)    :: findex
  integer                :: nc      ! Number of columns (in the E-W dimension) in the LIS grid
  integer                :: nr      ! Number of rows (in the N-S dimension) in the LIS grid
  integer                :: nstg4   ! Number of points in original STAGE IV grid
  real                   :: lis_gds(50)
  real                   :: ppt_field(nstg4)
  logical*1              :: lb(nstg4)
  real, dimension(nc,nr) :: varfield

! !DESCRIPTION:
!   This subroutine interpolates a given STAGE4 field 
!   to the LIS grid. 
!  The arguments are: 
!  \begin{description}
! \item[n]
!  index of the nest
! \item[nstg4]
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
  stg4_struc(n)%mi = nstg4
  mo = nc * nr   ! LIS Grid

  lo = .true.    !-Initialize output bitmap
  !lb = .true.   ! Set in read_stg4


!-- Interpolate to LIS grid

    if (trim(LIS_rc%met_interp(findex)) .eq. "bilinear") then 
! == BILINEAR INTERPOLATION ==== 
       call bilinear_interp( lis_gds,  lb, ppt_field, lo, &
            lis1d, stg4_struc(n)%mi, mo, &
            LIS_domain(n)%lat, LIS_domain(n)%lon,&
            stg4_struc(n)%w111, stg4_struc(n)%w121, &
            stg4_struc(n)%w211, stg4_struc(n)%w221, &
            stg4_struc(n)%n111, stg4_struc(n)%n121, &
            stg4_struc(n)%n211, stg4_struc(n)%n221, LIS_rc%udef, iret)

    elseif (trim(LIS_rc%met_interp(findex)) .eq. "budget-bilinear") then 
! == BUDGET BILINEAR INTERPOLATION ==== 
       call conserv_interp( lis_gds,  lb, ppt_field,  lo, &
            lis1d, stg4_struc(n)%mi, mo, & 
            LIS_domain(n)%lat, LIS_domain(n)%lon,&
            stg4_struc(n)%w112, stg4_struc(n)%w122,&
            stg4_struc(n)%w212, stg4_struc(n)%w222,&
            stg4_struc(n)%n112, stg4_struc(n)%n122,&
            stg4_struc(n)%n212, stg4_struc(n)%n222,LIS_rc%udef,iret)
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


end subroutine interp_stg4

