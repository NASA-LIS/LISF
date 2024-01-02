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
! !ROUTINE: gdas_reset_interp_input
!  \label{gdas_reset_interp_input}
!
! !REVISION HISTORY:
!  01 Feb 2016: James Geiger: Initial specification
!
! !INTERFACE:
subroutine gdas_reset_interp_input(n, findex, gridDesci)
! !USES:

   use LIS_coreMod,        only : LIS_rc, LIS_isatAfinerResolution
   use LIS_logMod,         only : LIS_logunit, LIS_endrun
   use gdas_forcingMod,    only : gdas_struc

   implicit none
! !ARGUMENTS: 
   integer, intent(in) :: n
   integer, intent(in) :: findex
   real, intent(in)    :: gridDesci(50)

! !DESCRIPTION:
! Resets the neighbours and weights arrays used for spatially
! interpolating the GDAS forcing data to the LIS running domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    forcing index
!  \item[gridDesci]
!    array of magic numbers describing the GDAS forcing domain
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[upscaleByAveraging\_input](\ref{upscaleByAveraging_input}) \newline
!    computes the neighbors for upscaling by averaging
!   \item[LIS\_isatAfinerResolution](\ref{LIS_isatAfinerResolution}) \newline
!    determines whether LIS' running domain is at a finer resolution
!    than the given resolution
!  \end{description}
!EOP

   integer :: rc

   deallocate(gdas_struc(n)%n111, stat=rc)
   deallocate(gdas_struc(n)%n121, stat=rc)
   deallocate(gdas_struc(n)%n211, stat=rc)
   deallocate(gdas_struc(n)%n221, stat=rc)
   deallocate(gdas_struc(n)%w111, stat=rc)
   deallocate(gdas_struc(n)%w121, stat=rc)
   deallocate(gdas_struc(n)%w211, stat=rc)
   deallocate(gdas_struc(n)%w221, stat=rc)

   deallocate(gdas_struc(n)%n112, stat=rc)
   deallocate(gdas_struc(n)%n122, stat=rc)
   deallocate(gdas_struc(n)%n212, stat=rc)
   deallocate(gdas_struc(n)%n222, stat=rc)
   deallocate(gdas_struc(n)%w112, stat=rc)
   deallocate(gdas_struc(n)%w122, stat=rc)
   deallocate(gdas_struc(n)%w212, stat=rc)
   deallocate(gdas_struc(n)%w222, stat=rc)

   if ( LIS_isatAfinerResolution(n,gridDesci(9)) ) then

      gdas_struc(n)%met_interp = LIS_rc%met_interp(findex)

      write(LIS_logunit,*) 'MSG: The GDAS forcing resolution is coarser ' // &
                           'than the running domain.'
      write(LIS_logunit,*) '     Interpolating with the ' // &
                           trim(gdas_struc(n)%met_interp) // ' method.'

      if ( gdas_struc(n)%met_interp == "bilinear" ) then

         allocate(gdas_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(gdas_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(gdas_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(gdas_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(gdas_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(gdas_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(gdas_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(gdas_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

         call bilinear_interp_input(n,gridDesci,        &
                 gdas_struc(n)%n111,gdas_struc(n)%n121, &
                 gdas_struc(n)%n211,gdas_struc(n)%n221, &
                 gdas_struc(n)%w111,gdas_struc(n)%w121, &
                 gdas_struc(n)%w211,gdas_struc(n)%w221)

      elseif ( gdas_struc(n)%met_interp == "budget-bilinear" ) then

         allocate(gdas_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(gdas_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(gdas_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(gdas_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(gdas_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(gdas_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(gdas_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(gdas_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

         call bilinear_interp_input(n,gridDesci,        &
                 gdas_struc(n)%n111,gdas_struc(n)%n121, &
                 gdas_struc(n)%n211,gdas_struc(n)%n221, &
                 gdas_struc(n)%w111,gdas_struc(n)%w121, &
                 gdas_struc(n)%w211,gdas_struc(n)%w221)

         allocate(gdas_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(gdas_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(gdas_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(gdas_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(gdas_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(gdas_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(gdas_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(gdas_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

         call conserv_interp_input(n,gridDesci,         &
                 gdas_struc(n)%n112,gdas_struc(n)%n122, &
                 gdas_struc(n)%n212,gdas_struc(n)%n222, &
                 gdas_struc(n)%w112,gdas_struc(n)%w122, &
                 gdas_struc(n)%w212,gdas_struc(n)%w222)
      endif
   else
      gdas_struc(n)%met_interp = LIS_rc%met_upscale(findex)

      write(LIS_logunit,*) 'MSG: The GDAS forcing resolution is finer ' // &
                           'than the running domain.'
      write(LIS_logunit,*) '     Upscaling with the ' // &
                           trim(gdas_struc(n)%met_interp) // ' method.'
   
   
      select case( gdas_struc(n)%met_interp )
      case( "average" )
         allocate(gdas_struc(n)%n111(gdas_struc(n)%mi))
   
         call upscaleByAveraging_input(gridDesci,                   &
                                       LIS_rc%gridDesc(n,:),        &
                                       gdas_struc(n)%mi,            &
                                       LIS_rc%lnc(n)*LIS_rc%lnr(n), &
                                       gdas_struc(n)%n111)
      case default
         write(LIS_logunit,*) 'The specified spatial interpolation option '
         write(LIS_logunit,*) 'is not supported for GDAS.'
         write(LIS_logunit,*) 'LIS is stopping.'
         call LIS_endrun()
      end select
   endif
end subroutine gdas_reset_interp_input
