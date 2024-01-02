!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"

module clm2_accumulMod

#if (defined ACCUM)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Subroutines [accini], [accslf], and [accext] can be used to 
! accumulate specified fields over user-defined intervals. This 
! interval and the type of accumulation is unique to each field
! processed. Subroutine [accini] defines the fields to be processed
! and the type of accumulation. Subroutine [accslf] does the actual
! accumulation. Subroutine [accext] extracts the current value of
! a specified field.
!
! Fields are accumulated for each patch by calls to subroutine [accslf].
! To accumulate a field, it must first be defined in subroutine [accini] 
! and then accumulated by calls to [accslf]. 
!
! Method: 
! 
! Author: Sam Levis
! 
!-----------------------------------------------------------------------
! $Id: accumulMod.F90,v 1.6 2004/11/24 22:57:07 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2_shr_const_mod, only: SHR_CONST_CDAY
  use clm2_varcon         
  use clm_varmap         
  use infnan
  implicit none

  private

!
! Public interfaces
!
  public :: accini        ! Initialization
  public :: accext        ! Extract
  public :: accslf        ! Single level field
  public :: accumvar_ini  ! Allocate dynamic memory
!
! Public data
!
  integer, public :: naccflds              !actual number of accumulated fields
  real(r8), public, allocatable :: accval(:,:) !accumulated field

!
! Data private to this module
!
  integer, parameter :: maxaccflds =  10   !maximum number of accumulated fields 
  character(len=40) :: accdes(maxaccflds)  !field description
  character(len= 8) :: accnam(maxaccflds)  !field name
  character(len= 8) :: accuni(maxaccflds)  !field units
  character(len= 8) :: acctyp(maxaccflds)  !field type: tavg, runm, runa, inst
  character(len= 8) :: ntavg = 'timeavg '  !time average of field
  character(len= 8) :: nrunm = 'runmean '  !running mean of field
  character(len= 8) :: nruna = 'runaccum'  !running accumulation of field
  character(len= 8) :: nins  = 'instant '  !instantaneous field value
  integer  :: accper(maxaccflds)           !field accumulation period
  real(r8) :: accbeg(maxaccflds)           !field initial or reset value

!=======================================================================
contains
!=======================================================================

  subroutine accini

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize field list for accslf
!
! Method: 
! Four types of accumulations are possible:
! o average over time interval
! o running mean over time interval
! o running accumulation over time interval
! o instantaneous at time of call
!
! Time average fields are only valid at the end of the averaging interval.
! Running means are valid once the length of the simulation exceeds the
! averaging interval. Accumulated fields are continuously accumulated.
! The trigger value -99999. resets the accumulation to zero.
!
! This subroutine sets:
! o number        of accumulated fields: naccflds
! o name          of accumulated fields: accnam
! o units         of accumulated fields: accuni
! o type          of accumulated fields: acctyp
! o description   of accumulated fields: accdes
! o period        of accumulation      : accper (- = days. + = iterations)
! o initial value of accumulated fields: accbeg
! 
! Author: Sam Levis
! 
!-----------------------------------------------------------------------

    use LIS_precisionMod	
    use clm_varctl  , only : nsrest            
    use LIS_timeMgrMod, only : get_step_size
    use spmdMod     , only : LIS_masterproc
    implicit none

! ------------------------ local variables ------------------------
    integer i,k                !do loop indices
    integer dtime              !time step size  
! -----------------------------------------------------------------
    print*,'not supposed to be called..'
#if 0 
    if (LIS_masterproc) then
       write(6,*) 'Initializing variables for time accumulation .....'
       write(6,'(72a1)') ("-",i=1,60)
    endif

! -----------------------------------------------------------------
! Begin list of user-defined fields to accumulate
! -----------------------------------------------------------------

    naccflds = 0

! The following is an example of a time-average field. 
! The accumulation period is set to 30 days. The initial value is zero.

    naccflds = naccflds + 1
    accnam(naccflds) = 'TDA     '
    accuni(naccflds) = 'K       '
    acctyp(naccflds) = ntavg
    accdes(naccflds) = '30-day average 2-m temperature'
    accper(naccflds) = -30
    accbeg(naccflds) = 0.   

! The following are examples running means. 
! The accumulation period is set to 10 days for a 10-day running mean. 

    naccflds = naccflds + 1
    accnam(naccflds) = 'T10     '
    accuni(naccflds) = 'K       '
    acctyp(naccflds) = nrunm
    accdes(naccflds) = '10-day running mean of 2-m temperature'
    accper(naccflds) = -10       
    accbeg(naccflds) = tfrz + 20. 

    naccflds = naccflds + 1
    accnam(naccflds) = 'FNPSN10 '
    accuni(naccflds) = 'UMOL/M2S'
    acctyp(naccflds) = nrunm
    accdes(naccflds) = '10-day running mean net cpy photosynth'
    accper(naccflds) = -10
    accbeg(naccflds) = 0.

    naccflds = naccflds + 1
    accnam(naccflds) = 'PREC365 '
    accuni(naccflds) = 'MM H2O/S'
    acctyp(naccflds) = nrunm
    accdes(naccflds) = '365-day running mean tot. precip.'
    accper(naccflds) = -365
    accbeg(naccflds) = 0.

! The following are examples of accumulated fields. 
! These types of fields are accumulated until a trigger value resets
! the accumulation to zero (see subroutine accslf). 
! Hence, [accper] is not valid.

    naccflds = naccflds + 1
    accnam(naccflds) = 'AGDD0   '
    accuni(naccflds) = 'K       '
    acctyp(naccflds) = nruna
    accdes(naccflds) = 'growing degree-days base 0C'
    accper(naccflds) = bigint    !not used
    accbeg(naccflds) = 0.

    naccflds = naccflds + 1
    accnam(naccflds) = 'AGDD5   '
    accuni(naccflds) = 'K       '
    acctyp(naccflds) = nruna
    accdes(naccflds) = 'growing degree-days base -5C'
    accper(naccflds) = bigint    !not used
    accbeg(naccflds) = 0.

    naccflds = naccflds + 1
    accnam(naccflds) = 'AGDDTW  '
    accuni(naccflds) = 'K       '
    acctyp(naccflds) = nruna
    accdes(naccflds) = 'growing degree-days base twmax'
    accper(naccflds) = bigint    !not used
    accbeg(naccflds) = 0.

    naccflds = naccflds + 1
    accnam(naccflds) = 'AGDD    '
    accuni(naccflds) = 'K       '
    acctyp(naccflds) = nruna
    accdes(naccflds) = 'growing degree-days base 5C'
    accper(naccflds) = bigint    !not used
    accbeg(naccflds) = 0.

! -----------------------------------------------------------------
! End list of user-defined fields to accumulate
! -----------------------------------------------------------------

    if (naccflds > maxaccflds) then
       write (6,*) 'ACCINI error: user-defined accumulation fields ', &
            'equal to ',naccflds,' exceeds maxaccflds '
       call endrun
    end if

! convert [accper] from days to iterations

    dtime = get_step_size()
    where (accper(:) < 0)
       accper(:) = -accper(:) * nint(SHR_CONST_CDAY) / dtime
    endwhere

! initialize all fields for an initial run

    if (nsrest == 0) then
       do i = 1, naccflds
          do k = begpatch, endpatch
             accval(k,i) = accbeg(i) 
          end do
       end do
    endif

! echo fields

    if (LIS_masterproc) then
       write(6,*)
       write(6,*) 'Accumulated fields'
       write(6,1002) 
       write(6,'(72a1)') ("_",i=1,71)
       do i = 1, naccflds
          if (accper(i) /= bigint) then
             write (6,1003) i,accnam(i),accuni(i),acctyp(i),accper(i),accbeg(i),accdes(i)
          else
             write (6,1004) i,accnam(i),accuni(i),acctyp(i),accbeg(i),accdes(i)
          endif
       end do
       write(6,'(72a1)') ("_",i=1,71)
       write(6,*)
       write(6,'(72a1)') ("-",i=1,60)
       write(6,*) 'Successfully initialized variables for accumulation'
       write(6,*)
    endif

1002 format(' No',' Name    ',' Units   ',' Type    ','Period',' Inival',' Description')
1003 format((1x,i2),(1x,a8),(1x,a8),(1x,a8), (1x,i5),(1x,f4.0),(1x,a40))
1004 format((1x,i2),(1x,a8),(1x,a8),(1x,a8),'  N.A.',(1x,f4.0),(1x,a40))
#endif
    return
  end subroutine accini

!=======================================================================

  subroutine accext (name, field, nstep)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Extract accumulated field 
!
! Method: 
! This routine extracts the field [name] from the [accslf] array.
! It extracts the current value except if the field type is a
! time average. In this case, an absurd value is assigned to 
! indicate the time average is not yet valid.
! 
! Author: Sam Levis
! 
!-----------------------------------------------------------------------

! ------------------------ arguments ------------------------------
    character(len=*), intent(in) :: name             !field name
    real(r8), intent(out):: field(begpatch:endpatch) !field values for current time step
    integer , intent(in) :: nstep                    !timestep index
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
    integer i          !do loop index
    integer n          !field index
    integer k          !patch index
! -----------------------------------------------------------------

! find field index. return if "name" is not on list

    n = 0
    do i = 1, naccflds
       if (name == accnam(i)) n = i
    end do
    if (n == 0) then
       write(6,*) 'ACCEXT error: field name ',name,' not found'
       return
    endif

!$OMP PARALLEL DO PRIVATE (K)
    do k = begpatch,endpatch
       if (acctyp(n) == ntavg .and. mod(nstep,accper(n)) /= 0) then
          field(k) = 1.e36     !assign absurd value when avg not ready
       else    
          field(k) = accval(k,n)
       endif
    end do

    return
  end subroutine accext

!=======================================================================

  subroutine accslf (name, field, nstep)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! accumulate single level field over specified time interval
!
! Method: 
! This subroutine accumulates the appropriate field in the array [accval].
! [field] is a "small" vector of length [litvec], which corresponds
! to locations [numpatchi] to [numpatchf] in the "big" vector of length [numpatch]
! Save accumulated field values in the array [accval].
! 
! Author: Sam Levis
! 
!-----------------------------------------------------------------------

! ------------------------ arguments ------------------------------
    character(len=*), intent(in) :: name             !field name
    real(r8), intent(in) :: field(begpatch:endpatch) !field values for current time step
    integer , intent(in) :: nstep                    !time step index
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
    integer i                   !do loop index
    integer n                   !field index
    integer k                   !patch index
    integer accper2(maxaccflds) !temporary replacement of accper
! -----------------------------------------------------------------

! find field index. return if "name" is not on list

    n = 0
    do i = 1, naccflds
       if (name == accnam(i)) n = i
    end do
    if (n == 0) then
       write(6,*) 'ACCSLF error: field name ',name,' not found'
       return
    endif

!$OMP PARALLEL DO PRIVATE (K)
    do k = begpatch,endpatch

! reset accumulated field value if necessary
! running mean and instantaneous fields never reset, but...

       if (acctyp(n) == ntavg) then          !time average field
          if ((mod(nstep,accper(n)) == 1) .and. (nstep /= 0)) then
             accval(k,n) = 0._r4             !every acc. period
          endif
       else if (acctyp(n) == nruna) then  !running accumulation field
          if (nint(field(k)) == -99999) then
             accval(k,n) = 0._r4        !reset at trigger -99999 and
          endif
       else if (acctyp(n) == nrunm) then
          accper2(n) = min (nstep,accper(n)) !reset accper until > nstep
       else if (acctyp(n) == nins) then
          accval(k,n) = 1.e36                !reset only as check
       else
          write(6,*) 'ACCSLF error: incorrect or no field type was '
          write(6,*) ' specified for field ',name
          accval(k,n) = 1.e36
       end if
       
! accumulate field

       if (acctyp(n) == ntavg) then
          accval(k,n) = accval(k,n) + field(k)
          if (mod(nstep,accper(n)) == 0) then !normalize at end of acc. prd
             accval(k,n) = accval(k,n) / accper(n)
          endif
       else if (acctyp(n) == nrunm) then
          accval(k,n) = ( (accper2(n)-1) * accval(k,n) + field(k) ) / accper2(n)
       else if (acctyp(n) == nruna) then      !running accumulation field
          accval(k,n) = min(max(accval(k,n) + field(k), 0.), 99999.)
       else if (acctyp(n) == nins) then       !instantaneous field value
          accval(k,n) = field(k)
       end if
    end do
       
    return
  end subroutine accslf

!=======================================================================

  subroutine accumvar_ini

! Allocate and initialize dynamic memoray 

    allocate (accval(begpatch:endpatch,maxaccflds))
    accval(:,:) = 1.e36

  end subroutine accumvar_ini

!=======================================================================

#endif

end module clm2_accumulMod
