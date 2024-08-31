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
! !ROUTINE: interp_TRMM3B42V7.F90
! \label{interp_TRMM3B42V7}
!
! !REVISION HISTORY:
!  3 Jul 2013: Soni Yatheendradas: Modified to add findex argument towards
!              spatial interpolation choice
!
! !INTERFACE: 

subroutine interp_TRMM3B42V7(n, nTRMM3B42V7, f, lb, ldt_gds, nc, nr, varfield, &
                             findex) 
! !USES:
  use TRMM3B42V7_forcingMod, only : TRMM3B42V7_struc
  use LDT_coreMod,           only : LDT_rc , LDT_domain
  use LDT_logMod,            only : LDT_logunit, LDT_endrun 

  implicit none

! !ARGUMENTS:
  integer, intent(in)    :: n
  integer                :: nc
  integer                :: nr
  integer                :: nTRMM3B42V7
  real                   :: ldt_gds(20)
  real                   :: f(nTRMM3B42V7)
  logical*1              :: lb(nTRMM3B42V7)
  real, dimension(nc,nr) :: varfield
  integer, intent(in)    :: findex 

!
! !DESCRIPTION:
!   This subroutine interpolates a given 3B42 V7 field
!   to the LDT grid.
!  The arguments are:
!  \begin{description}
! \item[n]
!  index of the nest
! \item[nTRMM3B42V7]
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
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
! \end{description}
!EOP


  integer :: iret
  integer :: mo
  integer :: count,i,j

  real, dimension(nc*nr) :: ldt1d

  logical*1 :: geogmask(nc,nr)
  logical*1 :: lo(nc*nr)

  !=== End variable declarations
  !--------------------------------------------------------------------
  ! Setting interpolation options (ip=0,bilinear)
  ! (km=1, one parameter, ibi=1,use undefined bitmap
  ! (needed for soil moisture and temperature only)
  ! Use budget bilinear (ip=3) for precip forcing fields
  !--------------------------------------------------------------------
  mo = nc*nr

  !--------------------------------------------------------------------  
  ! Initialize output bitmap. Important for soil moisture and temp.
  !--------------------------------------------------------------------  
  lb = .true.
  lo = .true.


  TRMM3B42V7_struc(n)%mi = nTRMM3B42V7

  select case( LDT_rc%met_gridtransform(findex) )

   case( "budget-bilinear" )

     call conserv_interp(ldt_gds,lb,f,lo,ldt1d,TRMM3B42V7_struc(n)%mi,mo,&
        LDT_domain(n)%lat, LDT_domain(n)%lon, &
        TRMM3B42V7_struc(n)%w112,&
        TRMM3B42V7_struc(n)%w122,TRMM3B42V7_struc(n)%w212,&
        TRMM3B42V7_struc(n)%w222,&
        TRMM3B42V7_struc(n)%n112,TRMM3B42V7_struc(n)%n122,&
        TRMM3B42V7_struc(n)%n212,&
        TRMM3B42V7_struc(n)%n222,-9999.0,iret)

   case( "neighbor" )  ! SY
     call neighbor_interp(ldt_gds,lb,f,lo,ldt1d,TRMM3B42V7_struc(n)%mi,mo,&
        LDT_domain(n)%lat, LDT_domain(n)%lon, &
        TRMM3B42V7_struc(n)%n112, LDT_rc%udef,iret)

   case default
     write(LDT_logunit,*) 'This interpolation not defined for TRMM data'
     write(LDT_logunit,*) 'Program stopping ... '
     call LDT_endrun()

  end select

  !--------------------------------------------------------------------    
  ! Create 2D array for main program. Also define a "soil" mask
  ! due to different geography between GDAS & LDAS. For LDAS land 
  ! points not included in GDAS geography dataset only.
  !--------------------------------------------------------------------    
  count = 0
  do j = 1, nr
     do i = 1, nc
        varfield(i,j) = ldt1d(i+count)
        geogmask(i,j) = lo(i+count)
     enddo
     count = count + nc
  enddo

end subroutine interp_TRMM3B42V7

