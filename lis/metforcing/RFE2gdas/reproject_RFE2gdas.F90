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
! !ROUTINE: reproject_RFE2gdas
! \label{reproject_RFE2gdas}
!
! !INTERFACE:
subroutine reproject_RFE2gdas(n,findex,month,nRFE2gdas,f1d,f2d,lb1d,lb2d, &
     lis_gds, nc,nr, varfield)
! !USES:
  use LIS_coreMod
  use LIS_spatialDownscalingMod
  use LIS_logMod
  use RFE2gdas_forcingMod

  implicit none

! !ARGUMENTS:
  integer, intent(in)    :: n
  integer, intent(in)    :: findex
  integer, intent(in)    :: month
  integer, intent(in)    :: nRFE2gdas
  real, intent(in)       :: f1d(nRFE2gdas)
  real, intent(in)       :: f2d(NINT(RFE2gdas_struc(n)%gridDesci(2)),&
                                NINT(RFE2gdas_struc(n)%gridDesci(3)))
  logical*1, intent(in)  :: lb1d(nRFE2gdas)
  logical*1, intent(in)  :: lb2d(NINT(RFE2gdas_struc(n)%gridDesci(2)),&
                                 NINT(RFE2gdas_struc(n)%gridDesci(3)))
  real,  intent(in)      :: lis_gds(50)
  integer, intent(in)    :: nc
  integer, intent(in)    :: nr
  real, intent(inout)    :: varfield(nc,nr)
!
! !DESCRIPTION:
!   This subroutine interpolates a given RFE2gdas field
!   to the LIS grid.
!  The arguments are:
!  \begin{description}
! \item[n]
!  index of the nest
!  \item[findex]
!    index of the supplemental forcing source
! \item[nRFE2gdas]
!  number of elements in the input grid
! \item[f1d]
!  1-d input data array to be interpolated
! \item[f2d]
!  2-d input data array to be interpolated
! \item[lb1d]
!  input 1-d bitmap
! \item[lb2d]
!  input 2-d bitmap
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
!
!==== Local Variables=======================

  integer   :: iret,mo
  integer   :: count1,i,j
  real, dimension(nc*nr) :: lis1d
  real, dimension(nc,nr) :: lis2d
  logical*1 :: lo(nc*nr)
  logical*1 :: lo2(nc,nr)
  integer   :: c,r,t
  real, allocatable :: chirpsprec_1d(:)
  logical*1, allocatable :: lb(:)

!=== End Variable Definition =======================

   lo = .true.
   mo = nc*nr

  if(LIS_rc%pcp_downscale(findex).ne.0) then 
    !input_data becomes the ratio field. 
     call LIS_generatePcpClimoRatioField(n,findex,"RFE2gdas",&
          month, RFE2gdas_struc(n)%mi, &
          f1d, lb1d)     
  endif

!- Upscale finer scale forcing field to coarser scale run domain:
  select case( LIS_rc%met_interp(findex) )

    case( "average" )   ! Upscaling 

#if 0
! Soni's former way:
      call upscaleByAveraging(LIS_rc%lnc(n), LIS_rc%lnr(n), &
         NINT(RFE2gdas_struc(n)%gridDesci(2)), &
         NINT(RFE2gdas_struc(n)%gridDesci(3)),&
         RFE2gdas_struc(n)%stc, RFE2gdas_struc(n)%str, &
         RFE2gdas_struc(n)%enc, RFE2gdas_struc(n)%enr, LIS_rc%udef,&
         lb2d, f2d, lo2,lis2d)
#endif

      call upscaleByAveraging( RFE2gdas_struc(n)%mi, &
           mo, LIS_rc%udef, &
           RFE2gdas_struc(n)%n111, &
           lb1d, f1d, lo, lis1d)


!- Downscale (or interpolate) or at same resolution as run domain:
    case( "bilinear" )
      call bilinear_interp(lis_gds,lb1d,f1d,lo,lis1d,&
          RFE2gdas_struc(n)%mi,mo,&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          RFE2gdas_struc(n)%w111, RFE2gdas_struc(n)%w121,&
          RFE2gdas_struc(n)%w211,RFE2gdas_struc(n)%w221,&
          RFE2gdas_struc(n)%n111,RFE2gdas_struc(n)%n121,&
          RFE2gdas_struc(n)%n211,RFE2gdas_struc(n)%n221,&
          LIS_rc%udef,iret)

    case( "budget-bilinear" )
      call conserv_interp(lis_gds,lb1d,f1d,lo,lis1d,&
          RFE2gdas_struc(n)%mi,mo,&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          RFE2gdas_struc(n)%w112,RFE2gdas_struc(n)%w122,&
          RFE2gdas_struc(n)%w212,RFE2gdas_struc(n)%w222,&
          RFE2gdas_struc(n)%n112,RFE2gdas_struc(n)%n122,&
          RFE2gdas_struc(n)%n212,RFE2gdas_struc(n)%n222,&
          LIS_rc%udef,iret)

    case( "neighbor" )
      call neighbor_interp(lis_gds,lb1d,f1d,lo,lis1d,&
          RFE2gdas_struc(n)%mi,mo,&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          RFE2gdas_struc(n)%n113,LIS_rc%udef,iret)
  end select

  if(LIS_rc%pcp_downscale(findex).ne.0) then 

     call LIS_pcpClimoDownscaling(n,findex,month,&
          nc*nr, lis1d, lo)     
  endif

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

end subroutine reproject_RFE2gdas
