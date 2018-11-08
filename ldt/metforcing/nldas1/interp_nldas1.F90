!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: interp_nldas1
! \label{interp_nldas1}
!
! !INTERFACE:
subroutine interp_nldas1(n, findex, pcp_flag, nldas1,f,lb,ldt_gds,nc,nr, &
     varfield)
! !USES:
  use LDT_coreMod,        only : LDT_rc, LDT_domain
  use nldas1_forcingMod,  only : nldas1_struc
  
  implicit none

! !ARGUMENTS:   
  integer, intent(in)   :: n 
  integer, intent(in)   :: findex
  logical, intent(in)   :: pcp_flag
  integer, intent(in)   :: nldas1
  real, intent(in)      :: f(nldas1)
  logical*1, intent(in) :: lb(nldas1)
  real, intent(in)      :: ldt_gds(20)
  integer, intent(in)   :: nc
  integer, intent(in)   :: nr
  real, intent(inout)   :: varfield(nc,nr)
!
! !DESCRIPTION:
!   This subroutine interpolates a given NLDAS1 field 
!   to the LDT grid. 
!  The arguments are: 
!  \begin{description}
! \item[n]
!  index of the nest
! \item[kpds]
!  grib decoding array
! \item[ngdas]
!  number of elements in the input grid
! \item[f]
!  input data array to be interpolated
! \item[lb]
!  input bitmap
! \item[ldt\_gds]
!  array description of the LDT grid
! \item[nc]
!  number of columns (in the east-west dimension) in the LDT grid
! \item[nr]
!  number of rows (in the north-south dimension) in the LDT grid
! \item[varfield]
!  output interpolated field
!  \end{description} 
! 
!
!  The routines invoked are: 
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
!  \item[neighbor\_interp](\ref{neighbor_interp}) \newline
!    spatially interpolate the forcing data using neighbor interpolation
! \end{description}
!EOP
  integer :: iret
  integer :: mo
  integer :: count1,i,j
  real, dimension(nc*nr) :: ldt1d
  logical*1 :: lo(nc*nr)

!=== End variable declarations

!-----------------------------------------------------------------------
! Setting interpolation options (ip=0,bilinear)
! (km=1, one parameter, ibi=1,use undefined bitmap
! (needed for soil moisture and temperature only)
! Use budget bilinear (ip=3) for precip forcing fields
!-----------------------------------------------------------------------
  mo = nc*nr
!-----------------------------------------------------------------------
! Initialize output bitmap. Important for soil moisture and temp.
!-----------------------------------------------------------------------
  lo = .true.

!-----------------------------------------------------------------------  
! Interpolate to LDT grid
!-----------------------------------------------------------------------
  select case( LDT_rc%met_gridtransform(findex) )
   case( "bilinear" )
     call bilinear_interp(ldt_gds,lb,f,lo,ldt1d,nldas1_struc(n)%mi,mo,&
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          nldas1_struc(n)%w111, nldas1_struc(n)%w121,&
          nldas1_struc(n)%w211,nldas1_struc(n)%w221,&
          nldas1_struc(n)%n111,nldas1_struc(n)%n121,&
          nldas1_struc(n)%n211,nldas1_struc(n)%n221,LDT_rc%udef,iret)
   case( "budget-bilinear" )
     if (pcp_flag) then     
        call conserv_interp(ldt_gds,lb,f,lo,ldt1d,nldas1_struc(n)%mi,mo,& 
             LDT_domain(n)%lat, LDT_domain(n)%lon,&
             nldas1_struc(n)%w112,nldas1_struc(n)%w122,&
             nldas1_struc(n)%w212,nldas1_struc(n)%w222,&
             nldas1_struc(n)%n112,nldas1_struc(n)%n122,&
             nldas1_struc(n)%n212,nldas1_struc(n)%n222,LDT_rc%udef,iret)
     else 
        call bilinear_interp(ldt_gds,lb,f,lo,ldt1d,nldas1_struc(n)%mi,mo,&
             LDT_domain(n)%lat, LDT_domain(n)%lon,&
             nldas1_struc(n)%w111,nldas1_struc(n)%w121,&
             nldas1_struc(n)%w211,nldas1_struc(n)%w221,&
             nldas1_struc(n)%n111,nldas1_struc(n)%n121,&
             nldas1_struc(n)%n211,nldas1_struc(n)%n221,LDT_rc%udef,iret)
     endif
   case( "neighbor" )
        call neighbor_interp(ldt_gds,lb,f,lo,ldt1d,nldas1_struc(n)%mi,mo,&
             LDT_domain(n)%lat, LDT_domain(n)%lon,&
             nldas1_struc(n)%n113,LDT_rc%udef,iret)
  end select

!-----------------------------------------------------------------------    
! Create 2D array for main program. Also define a "soil" mask
! due to different geography between LDAS & LDAS. For LDAS land 
! points not included in LDAS geography dataset only.
!-----------------------------------------------------------------------    
  count1 = 0
  do j = 1, nr
     do i = 1, nc
        varfield(i,j) = ldt1d(i+count1)
     enddo
     count1 = count1 + nc
  enddo

end subroutine interp_nldas1
