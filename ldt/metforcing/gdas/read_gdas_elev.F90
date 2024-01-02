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
! !ROUTINE: read_gdas_elev
! \label{read_gdas_elev}
!
! !REVISION HISTORY:
!
!  17Dec2004; Sujay Kumar; Initial Specificaton
!  09Aug2013; KR Arsenault; Modified for use in LDT
!
! !INTERFACE:
subroutine read_gdas_elev(n, findex, gdaselev, elevdiff)

! !USES:
  use LDT_coreMod,       only : LDT_rc, LDT_domain
  use LDT_metforcingMod, only : LDT_forc
  use LDT_logMod,        only : LDT_logunit, LDT_verify, LDT_endrun
  use LDT_fileIOMod
  use gdas_forcingMod,   only : gdas_struc
#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!- Terrain height will be set to run domain:
  real, intent(inout) :: gdaselev(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real, intent(inout) :: elevdiff(LDT_rc%met_nc(findex), LDT_rc%met_nr(findex))
!
! !DESCRIPTION:
!
!  Opens, reads, and interpolates GDAS model elevation to the LDT
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
!  \item[LDT_transform_paramgrid](\ref{LDT_transform_paramgrid}) \newline
!    Routine to read the elevation of the forcing
!    data in the same map projection used in LDT.
!  \end{description}
!EOP
  integer     :: c, r
  integer     :: i, j, t, iret ! Loop indices and error flags
  integer     :: icount
  integer     :: ftn, ndata, outpts, index
  integer     :: igrib
  logical     :: file_exists            ! Check file status 
  real        :: missingValue
  integer     :: pds5_val, pds7_val
  real        :: gdaselevin(LDT_rc%met_nc(findex)*LDT_rc%met_nr(findex))
  logical*1   :: lb(LDT_rc%met_nc(findex)*LDT_rc%met_nr(findex))
  real        :: elev_regrid(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  logical*1   :: lb_regrid(LDT_rc%lnc(n)*LDT_rc%lnr(n))
! ________________________________________________________________

  gdaselev = LDT_rc%udef
  elevdiff = LDT_rc%udef

  select case( LDT_rc%metforc_parms(findex) )

    case( "GDAS_T126" ) ! period 1980-2000
       gdas_struc(n)%elevfile = trim(LDT_rc%gdasT126elevfile(n))
    case( "GDAS_T170" ) ! period 2000-2002
       gdas_struc(n)%elevfile = trim(LDT_rc%gdast170elevfile(n))
    case( "GDAS_T254" ) ! period 2002-2005
       gdas_struc(n)%elevfile = trim(LDT_rc%gdast254elevfile(n))
    case( "GDAS_T382" ) ! period 2005-2010
       gdas_struc(n)%elevfile = trim(LDT_rc%gdast382elevfile(n))
    case( "GDAS_T574" ) ! period 2010-present
       gdas_struc(n)%elevfile = trim(LDT_rc%gdast574elevfile(n))
    case( "GDAS_T1534" ) ! period 2015-present
       gdas_struc(n)%elevfile = trim(LDT_rc%gdast1534elevfile(n))
    case default
       write(LDT_logunit,*) "ERR: GDAS resolution - elev grid not available"
       call LDT_endrun
  end select

!- Check initially if file exists:
  inquire( file=gdas_struc(n)%elevfile, exist=file_exists )   ! Check if file exists
  if (.not. file_exists)  then
     write(LDT_logunit,*) "GDAS elevation file missing: ",trim(gdas_struc(n)%elevfile)
  endif
  write(LDT_logunit,*) "Opening/reading GDAS elev file:: ",trim(gdas_struc(n)%elevfile)

#if(defined USE_GRIBAPI) 

  ndata = LDT_rc%met_nc(findex)*LDT_rc%met_nr(findex)

  gdaselevin = 0.0
  lb = .true.
  lb_regrid = .true.

  call grib_open_file(ftn,trim(gdas_struc(n)%elevfile),'r',iret)
  call LDT_verify(iret,'error grib_open_file in read_gdas_elev')

  call grib_new_from_file(ftn,igrib,iret)
  call LDT_verify(iret,'error in grib_new_from_file in read_gdas_elev')

  call grib_get(igrib,'indicatorOfParameter',pds5_val,iret)
  call LDT_verify(iret, 'error in grib_get: indicatorOfParameter in read_gdas_elev')

  call grib_get(igrib,'level',pds7_val,iret)
  call LDT_verify(iret, 'error in grib_get: level in read_gdas_elev')

  if( pds5_val == 8 .and. pds7_val == 0 ) then
     call grib_get(igrib,'values',gdaselevin,iret)
     call LDT_verify(iret, 'error in grib_get:values in read_gdas_elev')

     call grib_get(igrib,'missingValue',missingValue,iret)
     call LDT_verify(iret, 'error in grib_get:missingValue in read_gdas_elev')

  !- Interp elevation field to output field:
     outpts = LDT_rc%lnr(n)*LDT_rc%lnc(n)
     call LDT_transform_paramgrid(n, LDT_rc%met_gridtransform_parms(findex), &
              LDT_rc%met_gridDesc(findex,:), ndata, 1, gdaselevin, lb, &
              outpts, elev_regrid, lb_regrid )

     lb = .false.
     do t=1,ndata
        if(gdaselevin(t).ne.missingvalue) lb(t) = .true.
     end do

  !- Convert 1D to 2D elevation file:
     icount = 0
     do r = 1, LDT_rc%lnr(n)
        do c = 1, LDT_rc%lnc(n)
           icount = icount + 1
           gdaselev(c,r) = elev_regrid(icount)
        end do
     end do

#if 0
     if( LDT_rc%metforc_parms(findex) == "GDAS_T1534" ) then
       open( 500, file="gdaselev_t1534_globe3KM.1gd4r", &
             form="unformatted", access="direct", recl=LDT_rc%lnr(n)*LDT_rc%lnc(n)*4 )
       write(500, rec=1) gdaselev
     endif
     stop
#endif

     call grib_release(igrib,iret)
     call LDT_verify(iret,'error in grib_release in read_gdas_elev')

  else
     write(LDT_logunit,*) 'Could not retrieve entries in file: ',trim(gdas_struc(n)%elevfile)
     return
  endif

  call grib_close_file(ftn)

#endif

end subroutine read_gdas_elev
