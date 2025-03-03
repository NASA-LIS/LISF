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
! !ROUTINE: reproject_RFE2Daily
! \label{reproject_RFE2Daily}
!
! !REVISION HISTORY:
!  30 May 2010; Soni Yatheendradas, Initial LIS version for FEWSNET
!  20 Mar 2013; KR Arsenault, Cleaned up code and documentation
!
! !INTERFACE:
subroutine reproject_RFE2Daily(n,findex,num_pts,f1d,lb1d, &
                               lis_gds, nc,nr, varfield)
! !USES:
  use LIS_coreMod,          only : LIS_rc, LIS_domain
  use LIS_logMod,           only : LIS_logunit
  use RFE2Daily_forcingMod, only : RFE2Daily_struc

  implicit none

! !ARGUMENTS:
  integer, intent(in)    :: n
  integer, intent(in)    :: findex
  integer, intent(in)    :: num_pts
  real, intent(in)       :: f1d(num_pts)
  logical*1, intent(in)  :: lb1d(num_pts)
  real, intent(in)       :: lis_gds(50)
  integer, intent(in)    :: nc
  integer, intent(in)    :: nr
  real, intent(inout)    :: varfield(nc,nr)
!
! !DESCRIPTION:
!  This subroutine interpolates a given RFE2Daily field
!   to the LIS grid.
!
!  The arguments are:
!  \begin{description}
! \item[n]
!  index of the nest
!  \item[findex]
!   index of the supplemental forcing source
! \item[num\_pts]
!   number of elements in the input grid
! \item[f1d]
!   1-d input data array to be interpolated
! \item[lb1d]
!   input 1-d bitmap
! \item[lis\_gds]
!   array description of the LIS grid
! \item[nc]
!   number of columns (in the east-west dimension) in the LIS grid
! \item[nr]
!   number of rows (in the north-south dimension) in the LIS grid
! \item[varfield]
!   output interpolated field
!  \end{description}
!  
!  The routines invoked are:
!  \begin{description}
!  \item[upscaleByAveraging](\ref{upscaleByAveraging}) \newline
!    upscales scalar data from a finer grid to a coarser grid 
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
!  \item[neighbor\_interp](\ref{neighbor_interp}) \newline
!    spatially interpolate the forcing data using nearest neighbor interpolation
! \end{description}
!EOP

!==== Local Variables=======================
  integer   :: c,r,t
  integer   :: iret,mo
  integer   :: mi
  integer   :: count1,i,j
  logical*1 :: lo(nc*nr)
  real, dimension(nc*nr) :: lis1d

!=== End Variable Definition =======================

   lo = .true. ! SY: Not required
   mo = nc*nr

!- Upscale finer scale forcing field to coarser scale run domain:
  select case( LIS_rc%met_interp(findex) )

    case( "average" )   ! Upscaling 
      mi = RFE2Daily_struc(n)%mi
      call upscaleByAveraging(mi, mo, LIS_rc%udef, &
           RFE2Daily_struc(n)%n111, &
           lb1d, f1d, lo, lis1d)
      
!- Downscale (or interpolate) or at same resolution as run domain:
    case( "bilinear" )   
      call bilinear_interp(lis_gds,lb1d,f1d,lo,lis1d,&
             RFE2Daily_struc(n)%mi,mo,&
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             RFE2Daily_struc(n)%w111,  RFE2Daily_struc(n)%w121,&
             RFE2Daily_struc(n)%w211,  RFE2Daily_struc(n)%w221,&
             RFE2Daily_struc(n)%n111,  RFE2Daily_struc(n)%n121,&
             RFE2Daily_struc(n)%n211,  RFE2Daily_struc(n)%n221,&
             LIS_rc%udef,iret)

    case( "budget-bilinear" )   
      call conserv_interp(lis_gds,lb1d,f1d,lo,lis1d,&
          RFE2Daily_struc(n)%mi,mo,&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          RFE2Daily_struc(n)%w112, RFE2Daily_struc(n)%w122,&
          RFE2Daily_struc(n)%w212, RFE2Daily_struc(n)%w222,&
          RFE2Daily_struc(n)%n112, RFE2Daily_struc(n)%n122,&
          RFE2Daily_struc(n)%n212, RFE2Daily_struc(n)%n222,&
          LIS_rc%udef,iret)

    case( "neighbor" )   
      call neighbor_interp(lis_gds,lb1d,f1d,lo,lis1d,&
          RFE2Daily_struc(n)%mi,mo,&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          RFE2Daily_struc(n)%n113,LIS_rc%udef,iret)

  end select

!-----------------------------------------------------------------------
! Create 2D array for main program. 
!-----------------------------------------------------------------------
  count1 = 0
  do j = 1, nr
     do i = 1, nc
        varfield(i,j) = lis1d(i+count1)
     enddo
     count1 = count1 + nc
  enddo

end subroutine reproject_RFE2Daily
