!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! !ROUTINE: read_nam242_elev
! \label{read_nam242_elev}
!
! !REVISION HISTORY:
!     Sep 2012: NOHRSC/NOAA: Initial specification
!
! !INTERFACE:
subroutine read_nam242_elev(n, findex, nam242elev, elevdiff)

! !USES:
  use LDT_coreMod,        only : LDT_rc, LDT_domain
  use LDT_metforcingMod,  only : LDT_forc
  use LDT_logMod,         only : LDT_logunit, LDT_getNextUnitNumber, &
                                 LDT_releaseUnitNumber, LDT_endrun,  &
                                 LDT_verify
  use LDT_fileIOMod
  use nam242_forcingMod,  only : nam242_struc
#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
  real, intent(inout) :: nam242elev(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real, intent(inout) :: elevdiff(LDT_rc%met_nc(findex), LDT_rc%met_nr(findex))
!
! !DESCRIPTION:
!
!  Opens, reads, and interpolates NAM model elevation to the LDT
!  grid. The data will be used to perform any topographical 
!  adjustments to the forcing.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[LDT\_readData](\ref{LDT_readData}) \newline
!    Abstract method to read the elevation of the forcing
!    data in the same map projection used in LDT.
!  \end{description}
!EOP
  integer   :: ndata, outpts, ftn, iret, igrib, c, r, icount
  real      :: nam242elevin(LDT_rc%met_nc(findex)*LDT_rc%met_nr(findex))
  real      :: elev_regrid(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  logical*1 :: lb(LDT_rc%met_nc(findex)*LDT_rc%met_nr(findex))
  logical*1 :: lb_regrid(LDT_rc%lnc(n)*LDT_rc%lnr(n))

#if(defined USE_GRIBAPI) 
  nam242elev = LDT_rc%udef
  elevdiff   = LDT_rc%udef

  ndata = LDT_rc%met_nc(findex)*LDT_rc%met_nr(findex)
  nam242elevin = 0.0
  lb = .true.
  lb_regrid = .true.

  write(LDT_logunit,*) 'Reading the NAM elevation ',trim(nam242_struc(n)%elevfile)
     
  call grib_open_file(ftn,trim(nam242_struc(n)%elevfile),'r',iret)
  call LDT_verify(iret,'error grib_open_file in read_nam242_elev')

  call grib_new_from_file(ftn,igrib,iret)
  call LDT_verify(iret,'error in grib_new_from_file in read_nam242_elev')

  call grib_get(igrib,'values',nam242elevin,iret)
  call LDT_verify(iret, 'error in grib_get:values in read_nam242_elev')
     
  !- Interp elevation field to output field:
  outpts = LDT_rc%lnr(n)*LDT_rc%lnc(n)
  call LDT_transform_paramgrid(n, LDT_rc%met_gridtransform_parms(findex),   &
           LDT_rc%met_gridDesc(findex,:), ndata, 1, nam242elevin, lb, &
           outpts, elev_regrid, lb_regrid)

  !- Convert 1D to 2D elevation file:
  icount = 0
  do r = 1, LDT_rc%lnr(n)
     do c = 1, LDT_rc%lnc(n)
        icount = icount + 1
        nam242elev(c,r) = elev_regrid(icount)
     end do
  end do

  call grib_release(igrib,iret)
  call LDT_verify(iret,'error in grib_release in read_nam242_elev')

  call grib_close_file(ftn)
#else
  write(LDT_logunit,*) 'ERR: read_nam242_elev: ' //              &
                       'NAM242 support requires the GRIB_API library. ', &
                       'Please recompile LDT.'

  call LDT_endrun
#endif
end subroutine read_nam242_elev
