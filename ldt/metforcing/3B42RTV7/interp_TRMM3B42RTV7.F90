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
! !ROUTINE: interp_TRMM3B42RTV7
! \label{interp_TRMM3B42RTV7}
!
! !REVISION HISTORY: 
!  06 Jan 2015: KR Arsenault; Added support for latest V7 RT dataset
!
! !INTERFACE: 
subroutine interp_TRMM3B42RTV7(n, nx, ny, finput, ldt_gds, &
                               nc, nr, varfield, findex)

! !USES:
  use TRMM3B42RTV7_forcingMod, only : TRMM3B42RTV7_struc
  use LDT_coreMod, only : LDT_rc, LDT_domain
  use LDT_logMod,  only : LDT_logunit, LDT_endrun 

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n 
  integer, intent(in) :: nx
  integer, intent(in) :: ny
  integer, intent(in) :: nc
  integer, intent(in) :: nr
  real,    intent(in) :: ldt_gds(20)
  real,    intent(in) :: finput(nx,ny)
  real,    intent(inout) :: varfield(nc,nr)
  integer, intent(in) :: findex 
!
! !DESCRIPTION:
!   This subroutine interpolates a given TRMM 3B42RT V7 field
!   to the LDT grid.
!
!  The arguments are:
!  \begin{description}
! \item[n]
!  index of the nest
! \item[nx]
!  number of columns (in the east-west dimension) in the TRMM 3B42RT V7 grid
! \item[ny]
!  number of rows (in the north-south dimension) in the TRMM 3B42RTV7 grid
! \item[finput]
!  input data array to be interpolated
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
!  The routines invoked are:
!  \begin{description}
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
! \end{description}
!EOP
  integer    :: count1,i,j,v
  integer    :: mo
  integer    :: iret

  real       :: f(nx*ny)      ! Native Input grid
  logical*1  :: lb(nx*ny)
  real       :: ldt1d(nc*nr)  ! LDT-output grid
  logical*1  :: lo(nc*nr)

!=== End variable declarations

  TRMM3B42RTV7_struc(n)%mi = nx * ny
  mo = nc*nr

  v = 0
  Do j=1, ny
     Do i=1, nx
        v = v + 1
        f(v) = finput(i, j) 
     End Do
  End Do
  
!--------------------------------------------------------------------  
! Initialize output bitmap
!--------------------------------------------------------------------  
  lb = .true.
  lo = .true.

  where ( f < 0 )
     lb = .false.
  endwhere

! Interpolate:
  select case( LDT_rc%met_gridtransform(findex) )

    case( "budget-bilinear" )
      call conserv_interp(ldt_gds,lb,f,lo,ldt1d,&
           TRMM3B42RTV7_struc(n)%mi,mo,&
           LDT_domain(n)%lat, LDT_domain(n)%lon, &
           TRMM3B42RTV7_struc(n)%w112,TRMM3B42RTV7_struc(n)%w122,&
           TRMM3B42RTV7_struc(n)%w212,TRMM3B42RTV7_struc(n)%w222,&
           TRMM3B42RTV7_struc(n)%n112,TRMM3B42RTV7_struc(n)%n122,&
           TRMM3B42RTV7_struc(n)%n212,TRMM3B42RTV7_struc(n)%n222,&
           LDT_rc%udef,iret)

    case( "neighbor" )
     ! SY: begin for nearest neighbor
       call neighbor_interp(ldt_gds,lb,f,lo,ldt1d,&
            TRMM3B42RTV7_struc(n)%mi,mo,&
            LDT_domain(n)%lat, LDT_domain(n)%lon, &
            TRMM3B42RTV7_struc(n)%n113, LDT_rc%udef,iret)

    case default
      write(LDT_logunit,*) "This interpolation not defined for TRMM data "
      write(LDT_logunit,*) "Program stopping ... "
      call LDT_endrun()
  end select 

!--------------------------------------------------------------------    
! Create 2D array writing out final precip grid:
!--------------------------------------------------------------------    
  count1 = 0
  do j = 1, nr
     do i = 1, nc
        varfield(i,j) = ldt1d(i+count1)
     enddo
     count1 = count1 + nc
  enddo

end subroutine interp_TRMM3B42RTV7

