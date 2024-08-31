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
! !ROUTINE: interp_stg2
! \label{interp_stg2}
!
! !INTERFACE: 
subroutine interp_stg2 (n, findex, nstg2, ppt_field, lb, ldt_gds, &
                        nc, nr, varfield )

! !USES:  
  use LDT_coreMod, only: LDT_rc, LDT_domain
  use stg2_forcingMod, only : stg2_struc

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n 
  integer, intent(in)    :: findex
  integer                :: nc      ! Number of columns (in the E-W dimension) in the LDT grid
  integer                :: nr      ! Number of rows (in the N-S dimension) in the LDT grid
  integer                :: nstg2   ! Number of points in original STAGE II grid
  real                   :: ldt_gds(20)
  real                   :: ppt_field(nstg2)
  logical*1              :: lb(nstg2)
  real, dimension(nc,nr) :: varfield

! !DESCRIPTION:
!   This subroutine interpolates a given STG2 field 
!   to the LDT grid. 
!  The arguments are: 
!  \begin{description}
! \item[n]
!  index of the nest
! \item[nstg2]
!  number of elements in the input grid
! \item[ppt\_field]
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

  real    :: ldt1d(nc*nr)   ! LDT Grid  
  logical*1 :: lo(nc*nr)    ! LDT Grid

!=== End variable declarations

!--------------------------------------------------------------------
! NOTE:: Recommended to use budget bilinear for precip forcing fields
!--------------------------------------------------------------------
  stg2_struc(n)%mi = nstg2
  mo = nc * nr   ! LDT Grid
  lo = .true.    !-Initialize output bitmap
!  lb = .true.   ! lb is set in read_stg2

!- Assign undefined value for STAGE II
!        do i = 1, stg2_struc(n)%mi
!           if (ppt_field(i) < 0 ) then
!              udef = ppt_field(i)
!              print *, "udef:: ",udef
!              exit
!           end if
!        end do

!-- Interpolate to LDT grid

    select case( LDT_rc%met_gridtransform(findex) )

     case( "bilinear" )
       call bilinear_interp( ldt_gds, lb, ppt_field, lo, &
            ldt1d, stg2_struc(n)%mi, mo, &
            LDT_domain(n)%lat, LDT_domain(n)%lon,&
            stg2_struc(n)%w111, stg2_struc(n)%w121, &
            stg2_struc(n)%w211, stg2_struc(n)%w221, &
            stg2_struc(n)%n111, stg2_struc(n)%n121, &
            stg2_struc(n)%n211, stg2_struc(n)%n221, LDT_rc%udef, iret)

     case( "budget-bilinear" )
       call conserv_interp( ldt_gds, lb, ppt_field, lo, &
            ldt1d, stg2_struc(n)%mi, mo, & 
            LDT_domain(n)%lat, LDT_domain(n)%lon,&
            stg2_struc(n)%w112, stg2_struc(n)%w122,&
            stg2_struc(n)%w212, stg2_struc(n)%w222,&
            stg2_struc(n)%n112, stg2_struc(n)%n122,&
            stg2_struc(n)%n212, stg2_struc(n)%n222,LDT_rc%udef,iret)

    end select

!--------------------------------------------------------------------    
! Create 2D array (on LDT_domain) for main program. 
!--------------------------------------------------------------------    
   count1 = 0
   do j = 1, nr
      do i = 1, nc    
         varfield(i,j) = ldt1d(i+count1)
      enddo
      count1 = count1 + nc
   enddo


end subroutine interp_stg2

