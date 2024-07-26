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
! !ROUTINE: read_ecmwf_elev
! \label{read_ecmwf_elev}
!
! !REVISION HISTORY:
!
!  17Dec2004; Sujay Kumar; Initial Specificaton
!  09Jun2014; KR Arsenault; Modified for use in LDT
!
! !INTERFACE:
subroutine read_ecmwf_elev(n, findex, ecmwfelev, elevdiff)

! !USES:
  use LDT_coreMod,       only : LDT_rc, LDT_domain
  use LDT_metforcingMod, only : LDT_forc
  use LDT_logMod,        only : LDT_logunit, LDT_verify, LDT_endrun
  use LDT_fileIOMod
  use ecmwf_forcingMod,  only : ecmwf_struc
#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!- Terrain height will be set to run domain:
  real, intent(inout) :: ecmwfelev(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real, intent(inout) :: elevdiff(LDT_rc%met_nc(findex), LDT_rc%met_nr(findex))

! !DESCRIPTION:
!
!  Opens, reads, and interpolates ECMWF model elevation to the LDT
!  grid. The data will be used to perform any topographical 
!  adjustments to the forcing.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP
  integer   :: c, r
  integer   :: i, j, t, iret     ! Loop indices and error flags
  integer   :: icount
  integer   :: ftn, ndata, outpts, index
  integer   :: igrib
  real      :: missingValue
  integer   :: pds5_val, pds7_val
  real      :: ecmwfelevin(LDT_rc%met_nc(findex)*LDT_rc%met_nr(findex))
  logical*1 :: lb(LDT_rc%met_nc(findex)*LDT_rc%met_nr(findex))
  real      :: elev_regrid(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  logical*1 :: lb_regrid(LDT_rc%lnc(n)*LDT_rc%lnr(n))

  logical   :: file_exists        ! Check file status 
  character(100) :: ecmwf_filename

! ________________________________________________________________

  ecmwfelev = LDT_rc%udef
  elevdiff = LDT_rc%udef

!  select case( LDT_rc%metforc(findex) ) 
  select case( LDT_rc%metforc_parms(findex) ) 
    case( "ECMWF_S23R4" ) ! Change = 0; period 2001-2002
       ecmwf_filename = trim(LDT_rc%elevfileifs23r4(n))
    case( "ECMWF_S25R1" ) ! Change = 1; period 04/2002-02/2006
       ecmwf_filename = trim(LDT_rc%elevfileifs25r1(n))
    case( "ECMWF_S30R1" ) ! Change = 2; period 02/2006-06/2008
       ecmwf_filename = trim(LDT_rc%elevfileifs30r1(n))
    case( "ECMWF_S33R1" ) ! Change = 3; period 06/2008-03/2009
       ecmwf_filename = trim(LDT_rc%elevfileifs33r1(n))
    case( "ECMWF_S35R2" ) ! Change = 4; period 03/2009-10/2009
       ecmwf_filename = trim(LDT_rc%elevfileifs35r2(n))
    case( "ECMWF_S35R3" ) ! Change = 5; period 10/2009-01/2010
       ecmwf_filename = trim(LDT_rc%elevfileifs35r3(n))
    case( "ECMWF_S36R1" ) ! Change = 6; period 01/2010-05/2011
       ecmwf_filename = trim(LDT_rc%elevfileifs36r1(n))
    case( "ECMWF_S37R2" ) ! Change = 7; period 05/2011 onwards
       ecmwf_filename = trim(LDT_rc%elevfileifs37r2(n))
    case default
       write(LDT_logunit,*) "ERR: ECMWF elevation file not available"
       call LDT_endrun
  end select

!- Check initially if file exists:
  inquire( file=ecmwf_filename, exist=file_exists )   ! Check if file exists
  if (.not. file_exists)  then
     write(LDT_logunit,*) "ECMWF elevation file missing: ",trim(ecmwf_filename)
  endif
  write(LDT_logunit,*) "Opening/reading ECMWF elev file:: ",trim(ecmwf_filename)


#if(defined USE_GRIBAPI) 

  ndata = LDT_rc%met_nc(findex)*LDT_rc%met_nr(findex)
  ecmwfelevin = 0.0
  lb = .true.
  lb_regrid = .true.

  call grib_open_file(ftn,trim(ecmwf_filename),'r',iret)
  call LDT_verify(iret,'error grib_open_file in read_ecmwf_elev')

  call grib_new_from_file(ftn,igrib,iret)
  call LDT_verify(iret,'error in grib_new_from_file in read_ecmwf_elev')

  call grib_get(igrib,'indicatorOfParameter',pds5_val,iret)
  call LDT_verify(iret, 'error in grib_get: indicatorOfParameter in read_ecmwf_elev')

  call grib_get(igrib,'level',pds7_val,iret)
  call LDT_verify(iret, 'error in grib_get: level in read_ecmwf_elev')

  if( pds5_val == 129 .and. pds7_val == 0 ) then  ! 129 = height field
     call grib_get(igrib,'values', ecmwfelevin, iret)
     call LDT_verify(iret, 'error in grib_get:values in read_ecmwf_elev')

     call grib_get(igrib,'missingValue',missingValue,iret)
     call LDT_verify(iret, 'error in grib_get:missingValue in read_ecmwf_elev')

  !- Interp elevation field to output field:
     outpts = LDT_rc%lnr(n)*LDT_rc%lnc(n)
     call LDT_transform_paramgrid(n, LDT_rc%met_gridtransform_parms(findex), &
              LDT_rc%met_gridDesc(findex,:), ndata, 1, ecmwfelevin, lb, &
              outpts, elev_regrid, lb_regrid )

     lb = .false.
     do t=1,ndata
        if( ecmwfelevin(t).ne.missingvalue ) lb(t) = .true.
     end do
  !- Convert 1D to 2D elevation file:
     icount = 0
     do r = 1, LDT_rc%lnr(n)
        do c = 1, LDT_rc%lnc(n)
           icount = icount + 1
         ! Native ECMWF "terrain" field is in geopotential (m2/s2):
           ecmwfelev(c,r) = elev_regrid(icount) / 9.8  ! Converted to height (meters)
        end do
     end do

     call grib_release(igrib,iret)
     call LDT_verify(iret,'error in grib_release in read_ecmwf_elev')

  else
     write(LDT_logunit,*)"Could not retrieve entries in file: ",trim(ecmwf_filename)
     return
  endif

  call grib_close_file(ftn)

#endif

end subroutine read_ecmwf_elev
