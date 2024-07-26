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
! Macros for tracing - Requires ESMF 7_1_0+
#ifdef ESMF_TRACE
#define TRACE_ENTER(region) call ESMF_TraceRegionEnter(region)
#define TRACE_EXIT(region) call ESMF_TraceRegionExit(region)
#else
#define TRACE_ENTER(region)
#define TRACE_EXIT(region)
#endif
module LIS_vegDataMod
!BOP
!
! !MODULE: LIS_vegDataMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read various vegetation
!  datasets.
!
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the
!  greenness fraction, LAI, SAI and roughness data. Both real-time
!  and climatological datasets are supported. The climatological
!  data is expected to be provided in the parameter input file (from LPT)
!
! !REVISION HISTORY:
!
!  08 Aug 2005: Sujay Kumar; Initial implementation
!  10 Aug 2010: Jonathan Case: Modified for SPORT daily GVF data.
!
  use ESMF
  use LIS_precisionMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_histDataMod
  use LIS_fileIOMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none

  PRIVATE

!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LIS_greenness_setup    !allocates memory for required structures
  public :: LIS_read_greenness     !reads the greenness data
  public :: LIS_diagnosegfrac      !diagnoses greenness data for history output
  public :: LIS_greenness_finalize !cleanup allocated structures
  public :: LIS_greenness_reset    !resets alarms and data structures
  public :: LIS_read_shdmax        !read the shdmax map
  public :: LIS_read_shdmin        !read the shdmin map

  public :: LIS_roughness_setup   !allocates memory for required structures
  public :: LIS_read_roughness    !reads the roughness data
  public :: LIS_diagnoseroughness !maps the roughness data to history writer
  public :: LIS_roughness_finalize !cleanup allocated structures
  public :: LIS_roughness_reset    !resets datastructures

  public :: LIS_lai_setup ! allocates memory for required structures
  public :: LIS_read_lai  ! reads the LAI data
  public :: LIS_lai_finalize ! cleanup allocated structures
  public :: LIS_diagnoseLAI !maps LAI data to the history writer
  public :: LIS_lai_reset  !resets datastructures
  public :: LIS_read_laimax        !read the laimax map
  public :: LIS_read_laimin        !read the laimin map

  public :: LIS_sai_setup ! allocates memory for required structures
  public :: LIS_read_sai  ! reads the SAI data
  public :: LIS_diagnoseSAI !maps the SAI data to the LIS history writer
  public :: LIS_sai_finalize ! cleanup allocated structures
  public :: LIS_sai_reset  !resets the sai datastructures and alarams

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LIS_gfrac !data structure containing greenness fraction data.
  public :: LIS_roughness !data structure containing roughness data.
  public :: LIS_lai ! data structure containing LAI data.
  public :: LIS_sai ! data structure containing SAI data.
!EOP
  type, public :: gfrac_type_dec
     integer       :: realtimemode
     character(len=LIS_CONST_PATH_LEN) :: gfracfile
     real          :: gfracInterval
     character*50  :: gfracIntervalType
     logical       :: firstInstance
     real, allocatable :: greenness(:)
     real, allocatable :: vegmp1(:)
     real, allocatable :: vegmp2(:)
     real*8           :: time1, time2
     integer, allocatable :: n111(:), n121(:)
     integer, allocatable :: n211(:), n221(:)
     real, allocatable    :: w111(:), w121(:)
     real, allocatable    :: w211(:), w221(:)
  end type gfrac_type_dec

  type(gfrac_type_dec), allocatable :: LIS_gfrac(:)

  type, public :: roughness_type_dec
     character(len=LIS_CONST_PATH_LEN) :: roughnessfile
     character*50  :: roughnessIntervalType
     real          :: roughnessInterval
     logical       :: firstInstance
     real, allocatable :: roughness(:)
     real, allocatable :: z0v1(:)
     real, allocatable :: z0v2(:)
     integer, allocatable :: n111(:), n121(:)
     integer, allocatable :: n211(:), n221(:)
     real, allocatable    :: w111(:), w121(:)
     real, allocatable    :: w211(:), w221(:)
  end type roughness_type_dec

  type(roughness_type_dec), allocatable :: LIS_roughness(:)

  type, public :: lai_type_dec
     character(len=LIS_CONST_PATH_LEN) :: laifile
     character*50  :: laiIntervalType
     real          :: laiInterval
     logical       :: firstInstance
     type(ESMF_time)   :: time1
     type(ESMF_time)   :: time2
     real(r8), allocatable :: tlai(:)
     real(r8), allocatable :: lai1(:)
     real(r8), allocatable :: lai2(:)
     integer, allocatable :: n111(:), n121(:)
     integer, allocatable :: n211(:), n221(:)
     real, allocatable    :: w111(:), w121(:)
     real, allocatable    :: w211(:), w221(:)
  end type lai_type_dec

  type(lai_type_dec), allocatable :: LIS_lai(:)

  type, public ::  sai_type_dec
     character(len=LIS_CONST_PATH_LEN) :: saifile
     character*50  :: saiIntervalType
     real          :: saiInterval
     real(r8), allocatable :: tsai(:)
     real(r8), allocatable :: sai1(:)
     real(r8), allocatable :: sai2(:)
  end type sai_type_dec

  type(sai_type_dec), allocatable :: LIS_sai(:)

contains

!BOP
!
! !ROUTINE: LIS_greenness_setup
! \label{LIS_greenness_setup}
!
! !INTERFACE:
  subroutine LIS_greenness_setup
! !USES:

! !DESCRIPTION:
!
! Allocates memory for data structures for reading
! the greenness fraction datasets. This routine also
! reads the greenness datasets as initial values.
!
!  The routines invoked are:
!  \begin{description}
!   \item[gfracsetup](\ref{gfracsetup}) \newline
!    calls the registry to invoke the gfrac setup methods.
!   \item[LIS\_registerAlarm](\ref{LIS_registerAlarm}) \newline
!    registers the alarm for reading greenness datasets
!   \item[LIS\_computeTemporalWeights](\ref{LIS_computeTemporalWeights})
!      \newline
!    computes the interpolation weights
!   \item[read\_gfracclimo](\ref{read_gfracclimo}) \newline
!    reads the greenness climatology data from the LIS parameter data file.
!   \item[readgfrac](\ref{readgfrac}) \newline
!    invokes the method from the registry to read the real-time (time-varying)
!    greenness data.
!  \end{description}
!
!EOP
    implicit none
    integer   :: n
    integer   :: ndoms
    integer   :: rc,ios,nid
    integer :: i
    real :: wt1, wt2
    real, allocatable :: value1(:) ! temporary value holder for t1
    real, allocatable :: value2(:) ! temporary value holder for t2
    integer   :: t1, t2
    logical   :: file_exists

    TRACE_ENTER("green_setup")
    ndoms = 0
    do n=1,LIS_rc%nnest
       if(LIS_rc%usegreennessmap(n).ne."none") then
          ndoms = ndoms+1
       endif
    enddo

    if(ndoms.gt.0) then
       allocate(LIS_gfrac(LIS_rc%nnest))
       do n=1,LIS_rc%nnest
          LIS_gfrac(n)%firstInstance = .true.
          LIS_gfrac(n)%realtimemode = 0
       enddo

       do n=1,LIS_rc%nnest
          allocate(LIS_gfrac(n)%greenness(LIS_rc%ntiles(n)))
          allocate(LIS_gfrac(n)%vegmp1(LIS_rc%ntiles(n)))
          allocate(LIS_gfrac(n)%vegmp2(LIS_rc%ntiles(n)))
          LIS_gfrac(n)%vegmp1 = 0.0
          LIS_gfrac(n)%vegmp2 = 0.0
          LIS_gfrac(n)%greenness = 0.0

          if(LIS_rc%usegreennessmap(n).eq."LDT") then

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

             inquire(file=LIS_rc%paramfile(n), exist=file_exists)
             if(file_exists) then

                ios = nf90_open(path=trim(LIS_rc%paramfile(n)),&
                     mode=NF90_NOWRITE,ncid=nid)
                call LIS_verify(ios,'Error in nf90_open in read_gfracclimo')

                ios = nf90_get_att(nid, NF90_GLOBAL, &
                     'GREENNESS_DATA_INTERVAL', &
                     LIS_gfrac(n)%gfracIntervalType)
                call LIS_verify(ios, &
                     'Error in nf90_get_att in read_gfracclimo')

                ios = nf90_close(nid)
                call LIS_verify(ios,'Error in nf90_close in read_gfracclimo')
             else
                write(LIS_logunit,*) '[ERR] ',LIS_rc%paramfile(n), &
                     ' does not exist'
                write(LIS_logunit,*) '[ERR] program stopping ...'
                call LIS_endrun
             endif
#endif
             if(LIS_gfrac(n)%gfracIntervalType.eq."monthly") then
                LIS_gfrac(n)%gfracInterval = 2592000
             endif
!The intervaltype and interval is set in the plugin, now
!register the alarm.
             call LIS_registerAlarm("LIS gfrac read alarm",LIS_rc%ts, &
                  LIS_gfrac(n)%gfracInterval,&
                  intervalType=LIS_gfrac(n)%gfracIntervalType)

             call LIS_computeTemporalWeights(LIS_rc,&
                  LIS_gfrac(n)%gfracIntervalType, &
                  t1,t2,wt1,wt2,"LIS gfrac read alarm")

             allocate(value1(LIS_rc%ntiles(n)))
             allocate(value2(LIS_rc%ntiles(n)))

             value1 = LIS_rc%udef
             value2 = LIS_rc%udef

             call read_gfracclimo(n,t1,value1)
             call read_gfracclimo(n,t2,value2)

             do i=1,LIS_rc%ntiles(n)
                LIS_gfrac(n)%vegmp1(i) = value1(i)
                LIS_gfrac(n)%vegmp2(i) = value2(i)
             enddo
             deallocate(value1)
             deallocate(value2)
          else
             call gfracsetup(trim(LIS_rc%usegreennessmap(n))//char(0),n)

             call readgfrac(trim(LIS_rc%usegreennessmap(n))//char(0),&
                  n,wt1,wt2,LIS_gfrac(n)%vegmp1,LIS_gfrac(n)%vegmp2)

          endif
          do i=1,LIS_rc%ntiles(n)
             LIS_gfrac(n)%greenness(i) = (wt1*LIS_gfrac(n)%vegmp1(i))+&
                  (wt2*LIS_gfrac(n)%vegmp2(i))
          enddo
       end do
    endif
    TRACE_EXIT("green_setup")

  end subroutine LIS_greenness_setup

!BOP
!
! !ROUTINE: LIS_read_greenness
! \label{LIS_read_greenness}
!
! !INTERFACE:
  subroutine LIS_read_greenness(n)
! !USES:

    implicit none
! !ARGUMENTS:
    integer, intent(in) :: n
!
! !DESCRIPTION:
!
!  Reads the greenness fraction climalotogy and temporally interpolates
!  it to the current day.
!
!  The arguments are:
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!   \item[LIS\_computeTemporalWeights](\ref{LIS_computeTemporalWeights})
!     \newline
!    computes the interpolation weights
!   \item[read\_gfracclimo](\ref{read_gfracclimo}) \newline
!    reads the greenness climatology data from the LIS parameter data file.
!   \item[readgfrac](\ref{readgfrac}) \newline
!    invokes the method from the registry to read the real-time (time-varying)
!    greenness data.
!  \end{description}
!
!EOP

    logical :: gfracAlarmCheck
    integer :: i
    integer         :: t1, t2
    real            :: wt1, wt2
    real, allocatable    :: value1(:) ! temporary value holder for t1
    real, allocatable    :: value2(:) ! temporary value holder for t2

    TRACE_ENTER("green_read")
    allocate(value1(LIS_rc%ntiles(n)))
    allocate(value2(LIS_rc%ntiles(n)))

    if(LIS_rc%usegreennessmap(n).ne."none") then
       if(LIS_rc%usegreennessmap(n).eq."LDT") then
          gfracAlarmCheck = LIS_isAlarmRinging(LIS_rc,&
               "LIS gfrac read alarm",&
               LIS_gfrac(n)%gfracIntervalType)
!move the tindex to computetemporalweights?
          call LIS_computeTemporalWeights(LIS_rc,&
               LIS_gfrac(n)%gfracIntervalType, &
               t1,t2,wt1,wt2, "LIS gfrac read alarm")
          if(gfracAlarmCheck) then

             value1 = LIS_rc%udef
             value2 = LIS_rc%udef

             call read_gfracclimo(n,t1,value1)
             call read_gfracclimo(n,t2,value2)

             do i=1,LIS_rc%ntiles(n)
                LIS_gfrac(n)%vegmp1(i) = value1(i)
                LIS_gfrac(n)%vegmp2(i) = value2(i)
             enddo
          endif
       else

          call readgfrac(trim(LIS_rc%usegreennessmap(n))//char(0),&
               n,wt1,wt2,LIS_gfrac(n)%vegmp1,LIS_gfrac(n)%vegmp2)
       endif

       do i=1,LIS_rc%ntiles(n)
          LIS_gfrac(n)%greenness(i) = (wt1*LIS_gfrac(n)%vegmp1(i))+&
               (wt2*LIS_gfrac(n)%vegmp2(i))
       end do
    endif
    deallocate(value1)
    deallocate(value2)
    TRACE_EXIT("green_read")

  end subroutine LIS_read_greenness

!BOP
!
! !ROUTINE: LIS_greenness_finalize
! \label{LIS_greenness_finalize}
!
! !INTERFACE:
  subroutine LIS_greenness_finalize
! !USES:

!
! !DESCRIPTION:
!
! Deallocates objects created in this module
!
!EOP
    implicit none
    integer :: n
    integer :: ndoms

    ndoms = 0
    do n=1,LIS_rc%nnest
       if(LIS_rc%usegreennessmap(n).ne."none") ndoms = ndoms+1
    enddo

    if(ndoms.gt.0) then
       do n=1,LIS_rc%nnest
          deallocate(LIS_gfrac(n)%greenness)
          deallocate(LIS_gfrac(n)%vegmp1)
          deallocate(LIS_gfrac(n)%vegmp2)
       enddo
       deallocate(LIS_gfrac)
    endif
  end subroutine LIS_greenness_finalize

!BOP
!
! !ROUTINE: LIS_diagnosegfrac
! \label{LIS_diagnosegfrac}
!
! !INTERFACE:
  subroutine LIS_diagnosegfrac(n)
! !USES:

! !ARGUMENTS:
    implicit none
    integer, intent(in)   :: n

! !DESCRIPTION:
!  This routine maps the greenness data to the LIS history writer.
!
!  The arguments are:
!  \begin{description}
!  \item[n] index of the nest \newline
!  \end{description}
!
!  The routines called are:
!  \begin{description}
!  \item[LIS\_diagnoseOutputVar](\ref{LIS_diagnoseSurfaceOutputVar}) \newline
!    maps the greenness data to the LIS history writer
!  \end{description}
!EOP
    integer :: t
    real, pointer :: temp(:)

    TRACE_ENTER("green_diag")
    allocate(temp(LIS_rc%ntiles(n)))

    if ((LIS_rc%usegreennessmap(n).ne."none").and.    &
! For the following list of LSMs, the LIS_MOC_GREENNESS attribute is set
! within the src/surfacemodels/land/[lsm_name]/[lsm_name]_main.F90 files.
        (LIS_rc%lsm.ne."CLSM F2.5")          .and.    &
        (LIS_rc%lsm.ne."HySSIB")             .and.    &
        (LIS_rc%lsm.ne."Noah.2.7.1")         .and.    &
        (LIS_rc%lsm.ne."Noah.3.2")           .and.    &
        (LIS_rc%lsm.ne."Noah.3.3")           .and.    &
        (LIS_rc%lsm.ne."Noah.3.6")           .and.    &
        (LIS_rc%lsm.ne."Noah.3.9")           .and.    &
        (LIS_rc%lsm.ne."Noah-MP.3.6")        .and.    &
        (LIS_rc%lsm.ne."Noah-MP.4.0.1")) then
       temp = LIS_rc%udef
       do t=1,LIS_rc%ntiles(n)
          if(LIS_domain(n)%tile(t)%index.ne.-1) then
             temp(t) = LIS_gfrac(n)%greenness(t)
          endif
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_GREENNESS,vlevel=1,&
               value=temp(t),unit="-", &
               direction="-")
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_GREENNESS,vlevel=1,&
               value=temp(t)*100.0,unit="%", &
               direction="-")
       enddo
    endif
    deallocate(temp)
    TRACE_EXIT("green_diag")

  end subroutine LIS_diagnosegfrac

!BOP
!
! !ROUTINE: LIS_greenness_reset
! \label{LIS_greenness_reset}
!
! !INTERFACE:
  subroutine LIS_greenness_reset
! !USES:

! !DESCRIPTION:
!
! Resets the data structures for reading
! the greenness fraction datasets
!
!  The routines invoked are:
!  \begin{description}
!   \item[gfracsetup](\ref{gfracsetup}) \newline
!    calls the registry to invoke the gfrac data reading methods.
!   \item[LIS\_computeTemporalWeights](\ref{LIS_computeTemporalWeights})
!       \newline
!    computes the interpolation weights
!   \item[read\_gfracclimo](\ref{read_gfracclimo}) \newline
!    reads the greenness climatology data from the LIS parameter data file.
!   \item[readgfrac](\ref{readgfrac}) \newline
!    invokes the method from the registry to read the real-time (time-varying)
!    greenness data.
!  \end{description}
!
!EOP
    implicit none
    integer   :: n
    integer   :: ndoms
    integer   :: rc
    integer :: i
    real :: wt1, wt2
    real, allocatable :: value1(:) ! temporary value holder for t1
    real, allocatable :: value2(:) ! temporary value holder for t2
    integer       :: t1, t2

    TRACE_ENTER("green_reset")
    ndoms = 0
    do n=1,LIS_rc%nnest
       if(LIS_rc%usegreennessmap(n).ne."none") then
          ndoms = ndoms+1
       endif
    enddo

    if(ndoms.gt.0) then

       do n=1,LIS_rc%nnest
          LIS_gfrac(n)%firstInstance = .true.
          LIS_gfrac(n)%vegmp1 = 0.0
          LIS_gfrac(n)%vegmp2 = 0.0
          LIS_gfrac(n)%greenness = 0.0

          if(LIS_rc%usegreennessmap(n).eq."LDT") then

             call LIS_computeTemporalWeights(LIS_rc,&
                  LIS_gfrac(n)%gfracIntervalType, &
                  t1,t2,wt1,wt2,"LIS gfrac read alarm")

             allocate(value1(LIS_rc%ntiles(n)))
             allocate(value2(LIS_rc%ntiles(n)))

             value1 = LIS_rc%udef
             value2 = LIS_rc%udef

             call read_gfracclimo(n,t1,value1)
             call read_gfracclimo(n,t2,value2)

             do i=1,LIS_rc%ntiles(n)
                LIS_gfrac(n)%vegmp1(i) = value1(i)
                LIS_gfrac(n)%vegmp2(i) = value2(i)
             enddo
             deallocate(value1)
             deallocate(value2)
          else
             call gfracsetup(trim(LIS_rc%usegreennessmap(n))//char(0),n)

             call readgfrac(trim(LIS_rc%usegreennessmap(n))//char(0),&
                  n,wt1,wt2,LIS_gfrac(n)%vegmp1,LIS_gfrac(n)%vegmp2)

          endif

          do i=1,LIS_rc%ntiles(n)
             LIS_gfrac(n)%greenness(i) = (wt1*LIS_gfrac(n)%vegmp1(i))+&
                  (wt2*LIS_gfrac(n)%vegmp2(i))
          enddo

       end do
    endif
    TRACE_EXIT("green_reset")
  end subroutine LIS_greenness_reset


!BOP
!
! !ROUTINE: read_gfracclimo
!  \label{read_gfracclimo}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
  subroutine read_gfracclimo(n,time,array)
! !USES:

    implicit none
! !ARGUMENTS:
    integer, intent(in)         :: n
    integer, intent(in)         :: time
    real, intent(inout)         :: array(LIS_rc%ntiles(n))
! !DESCRIPTION:
!  This subroutine reads the greenness data climatology from the LIS
!  parameter data file
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[time]
!    time index (in months or quarters) of the data to be read
!   \item[array]
!    array containing the greenness values
!   \end{description}
!
!EOP

    integer :: ios1
    integer :: ios,nid,gfracid,ncId, nrId,mid
    integer :: nc,nr,t,months
    integer :: mo
    real, allocatable :: gfrac(:,:,:)
    real    :: localgfrac(LIS_rc%lnc(n),LIS_rc%lnr(n))
    logical :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    mo = time

    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then

       write(LIS_logunit,*)'[INFO] Reading greenness map for month ', mo,&
            'from ',trim(LIS_rc%paramfile(n))

       ios = nf90_open(path=trim(LIS_rc%paramfile(n)),&
            mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in read_gfracclimo')

       ios = nf90_inq_dimid(nid,"east_west",ncId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_gfracclimo')

       ios = nf90_inq_dimid(nid,"north_south",nrId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_gfracclimo')

       ios = nf90_inq_dimid(nid,"month",mId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_gfracclimo')

       ios = nf90_inquire_dimension(nid,ncId, len=nc)
       call LIS_verify(ios, &
            'Error in nf90_inquire_dimension in read_gfracclimo')

       ios = nf90_inquire_dimension(nid,nrId, len=nr)
       call LIS_verify(ios, &
            'Error in nf90_inquire_dimension in read_gfracclimo')

       ios = nf90_inquire_dimension(nid,mId, len=months)
       call LIS_verify(ios, &
            'Error in nf90_inquire_dimension in read_gfracclimo')

       ios = nf90_inq_varid(nid,'GREENNESS',gfracid)
       call LIS_verify(ios,'GREENNESS field not found in the LIS param file')

       ios = nf90_get_var(nid,gfracid,localgfrac,&
             start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
             LIS_nss_halo_ind(n,LIS_localPet+1),mo/),&
             count=(/LIS_rc%lnc(n),LIS_rc%lnr(n),1/))
       call LIS_verify(ios,'Error in nf90_get_var in read_gfracclimo')

       ios = nf90_close(nid)
       call LIS_verify(ios,'Error in nf90_close in read_gfracclimo')

       do t=1,LIS_rc%ntiles(n)
          array(t) = localgfrac(LIS_domain(n)%tile(t)%col,&
               LIS_domain(n)%tile(t)%row)
       enddo
    else
       write(LIS_logunit,*) '[ERR] gfrac map: ',LIS_rc%paramfile(n), &
            ' does not exist'
       write(LIS_logunit,*) '[ERR] program stopping ...'
       call LIS_endrun
    endif
#endif
  end subroutine read_gfracclimo

!BOP
!
! !ROUTINE: LIS_read_shdmin
! \label{LIS_read_shdmin}
!
! !INTERFACE:
  subroutine LIS_read_shdmin(n,array)
! !USES:

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n
  real,    intent(inout) :: array(LIS_rc%lnc(n),LIS_rc%lnr(n))
!
! !DESCRIPTION:
!  Reads the static albedo upper bound over deep snow
!  for each domain.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[array]
!    shdmin for the region of interest
!   \end{description}
!
!EOP

  integer :: ios1
  integer :: ios,nid,shdminid,ncId, nrId
  integer :: nc,nr,c,r
  logical :: file_exists

  TRACE_ENTER("green_readmin")
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  inquire(file=LIS_rc%paramfile(n), exist=file_exists)
  if(file_exists) then

     write(LIS_logunit,*)'[INFO] Reading SHDMIN map from ',&
          trim(LIS_rc%paramfile(n))

     ios = nf90_open(path=trim(LIS_rc%paramfile(n)),&
          mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error in nf90_open in LIS_read_shdmin')

     ios = nf90_inq_dimid(nid,"east_west",ncId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in LIS_read_shdmin')

     ios = nf90_inq_dimid(nid,"north_south",nrId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in LIS_read_shdmin')

     ios = nf90_inquire_dimension(nid,ncId, len=nc)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in LIS_read_shdmin')

     ios = nf90_inquire_dimension(nid,nrId, len=nr)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in LIS_read_shdmin')

     ios = nf90_inq_varid(nid,'SHDMIN',shdminid)
     call LIS_verify(ios,'SHDMIN field not found in the LIS param file')

     ios = nf90_get_var(nid,shdminid,array,&
          start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
          LIS_nss_halo_ind(n,LIS_localPet+1)/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/))
     call LIS_verify(ios,'Error in nf90_get_var in LIS_read_shdmin')

     ios = nf90_close(nid)
     call LIS_verify(ios,'Error in nf90_close in LIS_read_shdmin')

  else
     write(LIS_logunit,*) '[ERR] SHDMIN map: ',&
          LIS_rc%paramfile(n), ' does not exist'
     write(LIS_logunit,*) '[ERR] program stopping ...'
     call LIS_endrun
  endif

#endif
  TRACE_EXIT("green_readmin")

  end subroutine LIS_read_shdmin

!BOP
!
! !ROUTINE: LIS_read_shdmax
! \label{LIS_read_shdmax}
!
! !INTERFACE:
  subroutine LIS_read_shdmax(n,array)
! !USES:

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n
  real,    intent(inout) :: array(LIS_rc%lnc(n),LIS_rc%lnr(n))
!
! !DESCRIPTION:
!  Reads the static albedo upper bound over deep snow
!  for each domain.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[array]
!    shdmax for the region of interest
!   \end{description}
!
!EOP

  integer :: ios1
  integer :: ios,nid,shdmaxid,ncId, nrId
  integer :: nc,nr,c,r
  logical :: file_exists

  TRACE_ENTER("green_readmax")
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  inquire(file=LIS_rc%paramfile(n), exist=file_exists)
  if(file_exists) then

     write(LIS_logunit,*)'[INFO] Reading SHDMAX map from ',&
          trim(LIS_rc%paramfile(n))

     ios = nf90_open(path=trim(LIS_rc%paramfile(n)),&
          mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error in nf90_open in LIS_read_shdmax')

     ios = nf90_inq_dimid(nid,"east_west",ncId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in LIS_read_shdmax')

     ios = nf90_inq_dimid(nid,"north_south",nrId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in LIS_read_shdmax')

     ios = nf90_inquire_dimension(nid,ncId, len=nc)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in LIS_read_shdmax')

     ios = nf90_inquire_dimension(nid,nrId, len=nr)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in LIS_read_shdmax')

     ios = nf90_inq_varid(nid,'SHDMAX',shdmaxid)
     call LIS_verify(ios,'SHDMAX field not found in the LIS param file')

     ios = nf90_get_var(nid,shdmaxid,array,&
          start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
          LIS_nss_halo_ind(n,LIS_localPet+1)/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/))
     call LIS_verify(ios,'Error in nf90_get_var in LIS_read_shdmax')

     ios = nf90_close(nid)
     call LIS_verify(ios,'Error in nf90_close in LIS_read_shdmax')

  else
     write(LIS_logunit,*) '[ERR] SHDMAX map: ',&
          LIS_rc%paramfile(n), ' does not exist'
     write(LIS_logunit,*) '[ERR] program stopping ...'
     call LIS_endrun
  endif

#endif
  TRACE_EXIT("green_readmax")

  end subroutine LIS_read_shdmax

!BOP
!
! !ROUTINE: LIS_read_laimin
! \label{LIS_read_laimin}
!
! !INTERFACE:
  subroutine LIS_read_laimin(n,array)
! !USES:

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n
  real,    intent(inout) :: array(LIS_rc%lnc(n),LIS_rc%lnr(n))
!
! !DESCRIPTION:
!  Reads the static albedo upper bound over deep snow
!  for each domain.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[array]
!    laimin for the region of interest
!   \end{description}
!
!EOP

  integer :: ios1
  integer :: ios,nid,laiminid,ncId, nrId
  integer :: nc,nr,c,r
  logical :: file_exists

  TRACE_ENTER("lai_readmin")
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  inquire(file=LIS_rc%paramfile(n), exist=file_exists)
  if(file_exists) then

     write(LIS_logunit,*)'[INFO] Reading LAIMIN map from ',&
          trim(LIS_rc%paramfile(n))

     ios = nf90_open(path=trim(LIS_rc%paramfile(n)),&
          mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error in nf90_open in LIS_read_laimin')

     ios = nf90_inq_dimid(nid,"east_west",ncId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in LIS_read_laimin')

     ios = nf90_inq_dimid(nid,"north_south",nrId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in LIS_read_laimin')

     ios = nf90_inquire_dimension(nid,ncId, len=nc)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in LIS_read_laimin')

     ios = nf90_inquire_dimension(nid,nrId, len=nr)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in LIS_read_laimin')

     ios = nf90_inq_varid(nid,'LAIMIN',laiminid)
     call LIS_verify(ios,'LAIMIN field not found in the LIS param file')

     ios = nf90_get_var(nid,laiminid,array,&
          start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
          LIS_nss_halo_ind(n,LIS_localPet+1)/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/))
     call LIS_verify(ios,'Error in nf90_get_var in LIS_read_laimin')

     ios = nf90_close(nid)
     call LIS_verify(ios,'Error in nf90_close in LIS_read_laimin')

  else
     write(LIS_logunit,*) '[ERR] LAIMIN map: ',&
          LIS_rc%paramfile(n), ' does not exist'
     write(LIS_logunit,*) '[ERR] program stopping ...'
     call LIS_endrun
  endif

#endif
  TRACE_EXIT("lai_readmin")

  end subroutine LIS_read_laimin

!BOP
!
! !ROUTINE: LIS_read_laimax
! \label{LIS_read_laimax}
!
! !INTERFACE:
  subroutine LIS_read_laimax(n,array)
! !USES:

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n
  real,    intent(inout) :: array(LIS_rc%lnc(n),LIS_rc%lnr(n))
!
! !DESCRIPTION:
!  Reads the static albedo upper bound over deep snow
!  for each domain.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[array]
!    laimax for the region of interest
!   \end{description}
!
!EOP

  integer :: ios1
  integer :: ios,nid,laimaxid,ncId, nrId
  integer :: nc,nr,c,r
  logical :: file_exists

  TRACE_ENTER("lai_readmax")
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  inquire(file=LIS_rc%paramfile(n), exist=file_exists)
  if(file_exists) then

     write(LIS_logunit,*)'[INFO] Reading LAIMAX map from ',&
          trim(LIS_rc%paramfile(n))

     ios = nf90_open(path=trim(LIS_rc%paramfile(n)),&
          mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error in nf90_open in LIS_read_laimax')

     ios = nf90_inq_dimid(nid,"east_west",ncId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in LIS_read_laimax')

     ios = nf90_inq_dimid(nid,"north_south",nrId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in LIS_read_laimax')

     ios = nf90_inquire_dimension(nid,ncId, len=nc)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in LIS_read_laimax')

     ios = nf90_inquire_dimension(nid,nrId, len=nr)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in LIS_read_laimax')

     ios = nf90_inq_varid(nid,'LAIMAX',laimaxid)
     call LIS_verify(ios,'LAIMAX field not found in the LIS param file')

     ios = nf90_get_var(nid,laimaxid,array,&
          start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
          LIS_nss_halo_ind(n,LIS_localPet+1)/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/))          
     call LIS_verify(ios,'Error in nf90_get_var in LIS_read_laimax')

     ios = nf90_close(nid)
     call LIS_verify(ios,'Error in nf90_close in LIS_read_laimax')

  else
     write(LIS_logunit,*) '[ERR] LAIMAX map: ',&
          LIS_rc%paramfile(n), ' does not exist'
     write(LIS_logunit,*) '[ERR] program stopping ...'
     call LIS_endrun
  endif

#endif
  TRACE_EXIT("lai_readmax")

  end subroutine LIS_read_laimax

!BOP
!
! !ROUTINE: LIS_roughness_setup
! \label{LIS_roughness_setup}
!
! !INTERFACE:
  subroutine LIS_roughness_setup
! !USES:

! !DESCRIPTION:
!
! Allocates memory for data structures for reading
! the roughness fraction datasets
!
!  The routines invoked are:
!  \begin{description}
!   \item[gfracsetup](\ref{gfracsetup}) \newline
!    calls the registry to invoke the gfrac setup methods.
!  \end{description}
!
!EOP
    implicit none
    integer   :: n
    integer   :: ndoms
    integer   :: rc,ios,nid
    integer :: i
    real :: wt1, wt2
    real, allocatable :: value1(:) ! temporary value holder for t1
    real, allocatable :: value2(:) ! temporary value holder for t2
    integer   :: t1, t2
    logical   :: file_exists

    TRACE_ENTER("rough_setup")
    ndoms = 0
    do n=1,LIS_rc%nnest
       if(LIS_rc%useroughnessmap(n).ne."none") then
          ndoms = ndoms+1
       endif
    enddo

    if(ndoms.gt.0) then
       allocate(LIS_roughness(LIS_rc%nnest))

       do n=1,LIS_rc%nnest
          LIS_roughness(n)%firstInstance = .true.
       enddo

       do n=1,LIS_rc%nnest
          allocate(LIS_roughness(n)%roughness(LIS_rc%ntiles(n)))
          allocate(LIS_roughness(n)%z0v1(LIS_rc%ntiles(n)))
          allocate(LIS_roughness(n)%z0v2(LIS_rc%ntiles(n)))
          LIS_roughness(n)%z0v1 = 0.0
          LIS_roughness(n)%z0v2 = 0.0
          LIS_roughness(n)%roughness = 0.0

          if(LIS_rc%useroughnessmap(n).eq."LDT") then

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

             inquire(file=LIS_rc%paramfile(n), exist=file_exists)
             if(file_exists) then

                ios = nf90_open(path=trim(LIS_rc%paramfile(n)),&
                     mode=NF90_NOWRITE,ncid=nid)
                call LIS_verify(ios, &
                     'Error in nf90_open in read_roughnessclimo')

                ios = nf90_get_att(nid, NF90_GLOBAL, &
                     'ROUGHNESS_DATA_INTERVAL', &
                     LIS_roughness(n)%roughnessIntervalType)
                call LIS_verify(ios,&
                     'Error in nf90_get_att in read_roughnessclimo')

                ios = nf90_close(nid)
                call LIS_verify(ios,&
                     'Error in nf90_close in read_roughnessclimo')
             else
                write(LIS_logunit,*) '[ERR] ',LIS_rc%paramfile(n), &
                     ' does not exist'
                write(LIS_logunit,*) '[ERR] program stopping ...'
                call LIS_endrun
             endif
#endif
             if(LIS_roughness(n)%roughnessIntervalType.eq."monthly") then
                LIS_roughness(n)%roughnessInterval = 2592000
             endif
!The intervaltype and interval is set in the plugin, now
!register the alarm.
             call LIS_registerAlarm("LIS roughness read alarm",LIS_rc%ts, &
                  LIS_roughness(n)%roughnessInterval,&
                  intervalType=LIS_roughness(n)%roughnessIntervalType)

             call LIS_computeTemporalWeights(LIS_rc,&
                  LIS_roughness(n)%roughnessIntervalType, &
                  t1,t2,wt1,wt2,"LIS roughness read alarm")

             allocate(value1(LIS_rc%ntiles(n)))
             allocate(value2(LIS_rc%ntiles(n)))

             value1 = LIS_rc%udef
             value2 = LIS_rc%udef

             call read_roughnessclimo(n,t1,value1)
             call read_roughnessclimo(n,t2,value2)

             do i=1,LIS_rc%ntiles(n)
                LIS_roughness(n)%z0v1(i) = value1(i)
                LIS_roughness(n)%z0v2(i) = value2(i)
             enddo
          else
             call roughnesssetup(trim(LIS_rc%useroughnessmap(n))//char(0),n)
             call readroughness(trim(LIS_rc%useroughnessmap(n))//char(0), &
                  n, wt1, wt2, LIS_roughness(n)%z0v1, LIS_roughness(n)%z0v2)
          endif

          do i=1,LIS_rc%ntiles(n)
             LIS_roughness(n)%roughness(i) = (wt1*LIS_roughness(n)%z0v1(i))+&
                  (wt2*LIS_roughness(n)%z0v2(i))
          enddo
       end do
    endif
    TRACE_EXIT("rough_setup")

  end subroutine LIS_roughness_setup

!BOP
!
! !ROUTINE: LIS_read_roughness
! \label{LIS_read_roughness}
!
! !INTERFACE:
  subroutine LIS_read_roughness(n)
! !USES:
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_isAlarmRinging, LIS_computeTemporalWeights

    implicit none
! !ARGUMENTS:
    integer, intent(in) :: n
!
! !DESCRIPTION:
!
!  Reads the roughness fraction climalotogy and temporally interpolates
!  it to the current day.
!
!  The arguments are:
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!   \item[LIS\_computeTemporalWeights](\ref{LIS_computeTemporalWeights})
!      \newline
!    computes the interpolation weights
!   \item[readroughness](\ref{readroughness}) \newline
!    invokes the generic method in the registry to read the
!    roughness climatology data
!  \end{description}
!
!EOP

    logical :: roughnessAlarmCheck
    integer :: i
    integer         :: t1, t2
    real            :: wt1, wt2
    real, allocatable   :: value1(:) ! temporary value holder for t1
    real, allocatable   :: value2(:) ! temporary value holder for t2

    TRACE_ENTER("rough_setup")
    if(LIS_rc%useroughnessmap(n).ne."none") then
       if(LIS_rc%useroughnessmap(n).eq."LDT") then
          roughnessAlarmCheck = LIS_isAlarmRinging(LIS_rc,&
               "LIS roughness read alarm",&
               LIS_roughness(n)%roughnessIntervalType)

          call LIS_computeTemporalWeights(LIS_rc,&
               LIS_roughness(n)%roughnessIntervalType, &
               t1,t2,wt1,wt2)
          if(roughnessAlarmCheck) then

             value1 = LIS_rc%udef
             value2 = LIS_rc%udef

             call read_roughnessclimo(n,t1,value1)
             call read_roughnessclimo(n,t2,value2)

             do i=1,LIS_rc%ntiles(n)
                LIS_roughness(n)%z0v1(i) = value1(i)
                LIS_roughness(n)%z0v2(i) = value2(i)
             enddo
          endif
       else
          call readroughness(trim(LIS_rc%useroughnessmap(n))//char(0), &
               n, wt1, wt2, LIS_roughness(n)%z0v1, LIS_roughness(n)%z0v2)
       endif

       do i=1,LIS_rc%ntiles(n)
          LIS_roughness(n)%roughness(i) = (wt1*LIS_roughness(n)%z0v1(i))+&
               (wt2*LIS_roughness(n)%z0v2(i))
       end do
    endif
    TRACE_EXIT("rough_setup")
  end subroutine LIS_read_roughness

!BOP
!
! !ROUTINE: LIS_roughness_finalize
! \label{LIS_roughness_finalize}
!
! !INTERFACE:
  subroutine LIS_roughness_finalize
! !USES:
    use LIS_coreMod, only : LIS_rc
!
! !DESCRIPTION:
!
! Deallocates objects created in this module
!
!EOP
    implicit none
    integer :: n
    integer :: ndoms

    ndoms = 0
    do n=1,LIS_rc%nnest
       if(LIS_rc%useroughnessmap(n).ne."none") ndoms = ndoms+1
    enddo

    if(ndoms.gt.0) then
       do n=1,LIS_rc%nnest
          deallocate(LIS_roughness(n)%roughness)
          deallocate(LIS_roughness(n)%z0v1)
          deallocate(LIS_roughness(n)%z0v2)
       enddo
       deallocate(LIS_roughness)
    endif
  end subroutine LIS_roughness_finalize

!BOP
!
! !ROUTINE: LIS_diagnoseroughness
! \label{LIS_diagnoseroughness}
!
! !INTERFACE:
  subroutine LIS_diagnoseroughness(n)
! !USES:
    use LIS_coreMod,     only : LIS_rc, LIS_domain
    use LIS_histDataMod, only : LIS_diagnoseSurfaceOutputVar, LIS_MOC_ROUGHNESS
! !ARGUMENTS:
    implicit none
    integer, intent(in)   :: n

! !DESCRIPTION:
!  This routine maps the roughness data to the LIS history writer.
!
!  The arguments are:
!  \begin{description}
!  \item[n] index of the nest \newline
!  \end{description}
!
!  The routines called are:
!  \begin{description}
!  \item[LIS\_diagnoseOutputVar](\ref{LIS_diagnoseSurfaceOutputVar}) \newline
!    generic routine to map a single variable to the LIS
!    history writer
!  \end{description}
!EOP
    integer :: t
    real, pointer :: temp(:)

    TRACE_ENTER("rough_diag")
    allocate(temp(LIS_rc%ntiles(n)))

    if(LIS_rc%useroughnessmap(n).ne."none") then
       temp = LIS_rc%udef
       do t=1,LIS_rc%ntiles(n)
          if(LIS_domain(n)%tile(t)%index.ne.-1) then
             temp(t) = LIS_roughness(n)%roughness(t)
          endif
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ROUGHNESS,vlevel=1,&
               value=temp(t),unit="m", &
               direction="-")
       enddo
    endif
    deallocate(temp)
    TRACE_EXIT("rough_diag")

  end subroutine LIS_diagnoseroughness

!BOP
!
! !ROUTINE: LIS_roughness_reset
! \label{LIS_roughness_reset}
!
! !INTERFACE:
  subroutine LIS_roughness_reset
! !USES:
    use LIS_coreMod,    only : LIS_rc, LIS_config
    use LIS_logMod,     only : LIS_verify
    use LIS_timeMgrMod, only : LIS_computeTemporalWeights
! !DESCRIPTION:
!
! Resets the data structures for reading
! the roughness fraction datasets
!
!  The routines invoked are:
!  \begin{description}
!   \item[roughnesssetup](\ref{roughnesssetup}) \newline
!    calls the registry to invoke the roughness data reading methods.
!  \end{description}
!
!EOP
    implicit none
    integer   :: n
    integer   :: ndoms
    integer   :: rc
    integer :: i
    real :: wt1, wt2
    real, allocatable :: value1(:) ! temporary value holder for t1
    real, allocatable :: value2(:) ! temporary value holder for t2
    integer       :: t1, t2

    TRACE_ENTER("rough_reset")
    ndoms = 0
    do n=1,LIS_rc%nnest
       if(LIS_rc%useroughnessmap(n).ne."none") then
          ndoms = ndoms+1
       endif
    enddo

    if(ndoms.gt.0) then

       do n=1,LIS_rc%nnest

          LIS_roughness(n)%z0v1 = 0.0
          LIS_roughness(n)%z0v2 = 0.0
          LIS_roughness(n)%roughness = 0.0

          if(LIS_rc%useroughnessmap(n).ne."LDT") then
             call roughnesssetup(trim(LIS_rc%useroughnessmap(n))//char(0),n)
          else
 !            LIS_roughness(n)%roughnessIntervalType = "monthly"
          endif

          !Read the data for the first time
          call LIS_computeTemporalWeights(LIS_rc,&
               LIS_roughness(n)%roughnessIntervalType, &
               t1,t2,wt1,wt2)

          allocate(value1(LIS_rc%ntiles(n)))
          allocate(value2(LIS_rc%ntiles(n)))

          if(LIS_rc%useroughnessmap(n).eq."LDT") then
             call read_roughnessclimo(n,t1,value1)
             call read_roughnessclimo(n,t2,value2)
          else
             !EMK Fixed argument list
             !call readroughness(trim(LIS_rc%useroughnessmap(n))//char(0),&
             !     n,t1,value1)
             !call readroughness(trim(LIS_rc%useroughnessmap(n))//char(0),&
             !     n,t2,value2)
             call readroughness(trim(LIS_rc%useroughnessmap(n))//char(0), &
                  n, wt1, wt2, value1, value2)
          endif

          do i=1,LIS_rc%ntiles(n)
             LIS_roughness(n)%z0v1(i) = value1(i)
             LIS_roughness(n)%z0v2(i) = value2(i)
          enddo
          deallocate(value1)
          deallocate(value2)

          do i=1,LIS_rc%ntiles(n)
             LIS_roughness(n)%roughness(i) = (wt1*LIS_roughness(n)%z0v1(i))+&
                  (wt2*LIS_roughness(n)%z0v2(i))
          enddo
       end do
    endif
    TRACE_EXIT("rough_reset")
  end subroutine LIS_roughness_reset

!BOP
!
! !ROUTINE: read_roughnessclimo
!  \label{read_roughnessclimo}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
  subroutine read_roughnessclimo(n,time,array)
! !USES:

    use LIS_coreMod,        only : LIS_rc, LIS_domain, LIS_localPet,&
         LIS_ews_ind, LIS_ewe_ind,&
         LIS_nss_ind, LIS_nse_ind, LIS_ews_halo_ind,LIS_ewe_halo_ind, &
         LIS_nss_halo_ind, LIS_nse_halo_ind
    use LIS_logMod,         only : LIS_logunit, LIS_getNextUnitNumber, &
         LIS_releaseUnitNumber, LIS_endrun, LIS_verify

    implicit none
! !ARGUMENTS:
    integer, intent(in)         :: n
    integer, intent(in)         :: time
    real, intent(inout)         :: array(LIS_rc%ntiles(n))
! !DESCRIPTION:
!  This subroutine reads the greenness data climatology
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \end{description}
!
!EOP

    integer :: ios1
    integer :: ios,nid,roughnessid,ncId, nrId,mid
    integer :: nc,nr,t,months
    integer :: mo
    real, allocatable :: roughness(:,:,:)
    real    :: localroughness(LIS_rc%lnc(n),LIS_rc%lnr(n))
    logical :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    mo = time

    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then

       write(LIS_logunit,*)'[INFO] Reading roughness map for month ', mo, &
            'from ',trim(LIS_rc%paramfile(n))

       ios = nf90_open(path=trim(LIS_rc%paramfile(n)),&
            mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in read_roughnessclimo')

       ios = nf90_inq_dimid(nid,"east_west",ncId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_roughnessclimo')

       ios = nf90_inq_dimid(nid,"north_south",nrId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_roughnessclimo')

       ios = nf90_inq_dimid(nid,"month",mId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_roughnessclimo')

       ios = nf90_inquire_dimension(nid,ncId, len=nc)
       call LIS_verify(ios, &
            'Error in nf90_inquire_dimension in read_roughnessclimo')

       ios = nf90_inquire_dimension(nid,nrId, len=nr)
       call LIS_verify(ios, &
            'Error in nf90_inquire_dimension in read_roughnessclimo')

       ios = nf90_inquire_dimension(nid,mId, len=months)
       call LIS_verify(ios, &
            'Error in nf90_inquire_dimension in read_roughnessclimo')

       ios = nf90_inq_varid(nid,'ROUGHNESS',roughnessid)
       call LIS_verify(ios,'ROUGHNESS field not found in the LIS param file')

       ios = nf90_get_var(nid,roughnessid,localroughness,&
            start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1),mo/),&
            count=(/LIS_rc%lnc(n),LIS_rc%lnr(n),1/))
       call LIS_verify(ios,'Error in nf90_get_var in read_roughnessclimo')

       ios = nf90_close(nid)
       call LIS_verify(ios,'Error in nf90_close in read_roughnessclimo')

       do t=1,LIS_rc%ntiles(n)
          array(t) = localroughness(LIS_domain(n)%tile(t)%col,&
               LIS_domain(n)%tile(t)%row)
       enddo
    else
       write(LIS_logunit,*) '[ERR] roughness map: ',LIS_rc%paramfile(n), &
            ' does not exist'
       write(LIS_logunit,*) '[ERR] program stopping ...'
       call LIS_endrun
    endif
#endif
  end subroutine read_roughnessclimo

#if 0
!BOP
!
! !ROUTINE: LIS_lai_setup
! \label{LIS_lai_setup_old}
!
! !INTERFACE:
  subroutine LIS_lai_setup
! !USES:
    use LIS_coreMod,    only : LIS_rc, LIS_Config, LIS_domain
    use LIS_timeMgrMod, only : LIS_computeTemporalWeights, &
         LIS_registerAlarm, LIS_calendar
    use LIS_logMod,     only : LIS_verify, LIS_logunit

! !DESCRIPTION:
!
! Allocates memory and other structures for reading
! LAI datasets
!
!  The routines invoked are:
!  \begin{description}
!   \item[LIS\_registerAlarm](\ref{LIS_registerAlarm}) \newline
!    registers the alarm for reading LAI datasets
!   \item[LIS\_computeTemporalWeights](\ref{LIS_computeTemporalWeights})
!     \newline
!    computes the interpolation weights
!   \item[read\_laiclimo](\ref{read_laiclimo}) \newline
!    reads the climatological lai data
!   \item[laisetup](\ref{laisetup}) \newline
!    initializes the realtime lai reader
!   \item[readlai](\ref{readlai}) \newline
!    reads the realtime lai data
!  \end{description}
!EOP
    implicit none
    integer :: n, i
    integer :: rc
    integer :: ndoms
    real, allocatable :: value1(:) ! temporary value holder for mo1
    real, allocatable :: value2(:) ! temporary value holder for mo2
    real          :: wt1,wt2
    integer       :: t1,t2, t_delta, t_delta_interval
    type(ESMF_Time)   :: time
    type(ESMF_TimeInterval) :: deltaT, deltaTinterval
    logical:: forward_search
    integer :: status

    TRACE_ENTER("lai_setup")
    ndoms = 0
    do n=1,LIS_rc%nnest
       if(LIS_rc%uselaimap(n).ne."none") then
          ndoms = ndoms+1
       endif
    enddo

    if(ndoms.gt.0) then !at least one nest requires lai to be supplied
       allocate(LIS_lai(LIS_rc%nnest))

       do n=1,LIS_rc%nnest
          LIS_lai(n)%firstInstance = .true.
       enddo
       !initialize variables/alarms
       do n=1,LIS_rc%nnest
          allocate(LIS_lai(n)%lai1(LIS_rc%ntiles(n)))
          allocate(LIS_lai(n)%lai2(LIS_rc%ntiles(n)))
          allocate(LIS_lai(n)%tlai(LIS_rc%ntiles(n)))

          LIS_lai(n)%lai1 = 0
          LIS_lai(n)%lai2 = 0
          LIS_lai(n)%tlai = 0

          !set alarms
          select case (LIS_rc%uselaimap(n))
          case ("LDT")
#if (defined USE_NETCDF3 || defined USE_NETCDF4)

             inquire(file=LIS_rc%paramfile(n), exist=file_exists)
             if(file_exists) then

                ios = nf90_open(path=trim(LIS_rc%paramfile(n)),&
                     mode=NF90_NOWRITE,ncid=nid)
                call LIS_verify(ios,'Error in nf90_open in read_laiclimo')

                ios = nf90_get_att(nid, NF90_GLOBAL, 'LAISAI_DATA_INTERVAL', &
                     LIS_lai(n)%laiIntervalType)
                call LIS_verify(ios,'Error in nf90_get_att in read_laiclimo')
                ios = nf90_close(nid)
                call LIS_verify(ios,'Error in nf90_close in read_laiclimo')
             else
                write(LIS_logunit,*) '[ERR] ',LIS_rc%paramfile(n), &
                     ' does not exist'
                write(LIS_logunit,*) '[ERR] program stopping ...'
                call LIS_endrun
             endif
#endif
             if(LIS_lai(n)%laiIntervalType.eq."monthly") then
                LIS_lai(n)%laiInterval = 2592000  ! 30 days
             endif

             call LIS_registerAlarm("LIS LAI climo read alarm", LIS_rc%ts, &
                  LIS_lai(n)%laiInterval, &
                  intervalType = LIS_lai(n)%laiIntervalType)
             call LIS_computeTemporalWeights(LIS_rc,&
                  LIS_lai(n)%laiIntervalType, &
                  t1,t2,wt1,wt2)
             allocate(value1(LIS_rc%ntiles(n)))
             allocate(value2(LIS_rc%ntiles(n)))
             call read_laiclimo(n,t1,value1)
             call read_laiclimo(n,t2,value2)
             do i=1,LIS_rc%ntiles(n)
                LIS_lai(n)%lai1(i) = value1(i)
                LIS_lai(n)%lai2(i) = value2(i)
             enddo
             deallocate(value1)
             deallocate(value2)

          case default ! is some kind of real-time
!             !presume for now that files no more frequent than daily
!             call LIS_registerAlarm("LIS LAI real-time read alarm", &
!                  LIS_rc%ts, &
!                  86400)
             call laisetup(trim(LIS_rc%uselaimap(n))//char(0),n)

             ! EMK Fixed argument list.
             !call readlai(trim(LIS_rc%uselaimap(n))//char(0), n,&
             !     wt1, wt2, LIS_lai(n)%lai1, LIS_lai(n)%time1)
             call readlai(trim(LIS_rc%uselaimap(n))//char(0), n,&
                  wt1, wt2, LIS_lai(n)%lai1, LIS_lai(n)%lai2)

          end select
          do i=1,LIS_rc%ntiles(n)
             LIS_lai(n)%tlai(i)=wt1*LIS_lai(n)%lai1(i)+wt2*LIS_lai(n)%lai2(i)
          enddo
       end do
#if 0
!read in initial values
       do n=1,LIS_rc%nnest
          select case (LIS_rc%uselaimap(n))
          case ("LDT")
             call LIS_computeTemporalWeights(LIS_rc,&
                  LIS_lai(n)%laiIntervalType, &
                  t1,t2,wt1,wt2)
             allocate(value1(LIS_rc%ntiles(n)))
             allocate(value2(LIS_rc%ntiles(n)))
             call read_laiclimo(n,t1,value1)
             call read_laiclimo(n,t2,value2)
             do i=1,LIS_rc%ntiles(n)
                LIS_lai(n)%lai1(i) = value1(i)
                LIS_lai(n)%lai2(i) = value2(i)
             enddo
             deallocate(value1)
             deallocate(value2)

             do i=1,LIS_rc%ntiles(n)
                LIS_lai(n)%tlai(i) = wt1* LIS_lai(n)%lai1(i)+ &
                     wt2*LIS_lai(n)%lai2(i)
             end do
          case default ! is some kind of real-time
             call laisetup(trim(LIS_rc%uselaimap(n))//char(0),n)

             forward_search=.false.
             call readlai(trim(LIS_rc%uselaimap(n))//char(0),n,&
                  forward_search,LIS_lai(n)%lai1,LIS_lai(n)%time1)

             forward_search=.true.
             call readlai(trim(LIS_rc%uselaimap(n))//char(0),n,&
                  forward_search,LIS_lai(n)%lai2,LIS_lai(n)%time2)

             call ESMF_TimeSet(time, yy=LIS_rc%yr, &
                  mm=LIS_rc%mo, &
                  dd=LIS_rc%da, &
                  h =LIS_rc%hr, &
                  m =LIS_rc%mn, &
                  s =LIS_rc%ss, &
                  calendar = LIS_calendar, &
                  rc = status)
             deltaT=time-LIS_lai(n)%time1
             deltaTinterval=LIS_lai(n)%time2-LIS_lai(n)%time1

             call ESMF_TimeIntervalGet(deltaT, s=t_delta)
             call ESMF_TimeIntervalGet(deltaTinterval, s=t_delta_interval)

             !get current values from linear interpolation
             wt1=real(t_delta_interval-t_delta)/real(t_delta_interval)
             wt2=real(t_delta)/real(t_delta_interval)
             do i=1,LIS_rc%ntiles(n)
                LIS_lai(n)%tlai(i)= &
                     wt1*LIS_lai(n)%lai1(i)+wt2*LIS_lai(n)%lai2(i)
             enddo

             if (time >= LIS_lai(n)%time2) then
                LIS_lai(n)%lai1=LIS_lai(n)%lai2
                LIS_lai(n)%time1=LIS_lai(n)%time2

                forward_search=.true.
                call readlai(trim(LIS_rc%uselaimap(n))//char(0),n,&
                     forward_search,LIS_lai(n)%lai2,LIS_lai(n)%time2)
             end if
          end select
       enddo
    endif
#endif
    TRACE_EXIT("lai_setup")
  end subroutine LIS_lai_setup

!BOP
!
! !ROUTINE: LIS_read_lai
! \label{LIS_read_lai_old}
!
! !INTERFACE:
  subroutine LIS_read_lai(n)
! !USES:
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_isAlarmRinging,LIS_computeTemporalWeights, &
         LIS_calendar
    use LIS_logMod,     only : LIS_verify, LIS_logunit

    implicit none
! !ARGUMENTS:
    integer, intent(in) :: n
!
! !DESCRIPTION:
!
!  Reads the LAI climatology and temporally interpolates it to the
!  current day.
!
!  The arguments are:
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!   \item[readlai](\ref{readlai}) \newline
!    invokes the generic method in the registry to read the
!    LAI climatology data
!  \end{description}
!
!EOP
    real, allocatable :: value1(:) ! temporary value holder for mo1
    real, allocatable :: value2(:) ! temporary value holder for mo2
    real, allocatable :: value(:)  !
    real              :: wt1,wt2
    integer           :: t1, t2, t_delta, t_delta_interval
    type(ESMF_Time)   :: time1, time2, time
    type(ESMF_TimeInterval) :: deltaT, deltaTinterval
    integer           :: i
    logical           :: laiAlarmCheckclimo,laiAlarmCheckrealtime
    logical           :: forward_search
    logical           :: file_exists
    integer  :: status

    TRACE_ENTER("lai_read")
    select case (LIS_rc%uselaimap(n))
    case("none")
       !nothing
    case ("LDT")
       laiAlarmCheckclimo = LIS_isAlarmRinging(LIS_rc,&
            "LIS LAI climo read alarm",&
            LIS_lai(n)%laiIntervalType)

       call LIS_computeTemporalWeights(LIS_rc,&
            LIS_lai(n)%laiIntervalType, t1,t2,wt1,wt2)

       if(laiAlarmCheckclimo) then
          allocate(value1(LIS_rc%ntiles(n)))
          allocate(value2(LIS_rc%ntiles(n)))

          call read_LAIclimo(n,t1,value1)
          call read_LAIclimo(n,t2,value2)

          do i=1,LIS_rc%ntiles(n)
             LIS_lai(n)%lai1(i) = value1(i)
             LIS_lai(n)%lai2(i) = value2(i)
          enddo
          deallocate(value1)
          deallocate(value2)

       endif
       do i=1,LIS_rc%ntiles(n)
          LIS_lai(n)%tlai(i) = wt1* LIS_lai(n)%lai1(i)+ &
               wt2*LIS_lai(n)%lai2(i)
       end do
    case default ! is some kind of real-time
       call readlai(trim(LIS_rc%uselaimap(n))//char(0), n, &
            wt1, wt2, LIS_lai(n)%lai1, LIS_lai(n)%lai2)
       do i=1,LIS_rc%ntiles(n)
          LIS_lai(n)%tlai(i) = wt1* LIS_lai(n)%lai1(i)+ &
               wt2*LIS_lai(n)%lai2(i)
       enddo
    end select

#if 0
          call ESMF_TimeSet(time, yy=LIS_rc%yr, &
               mm=LIS_rc%mo, &
               dd=LIS_rc%da, &
               h =LIS_rc%hr, &
               m =LIS_rc%mn, &
               s =LIS_rc%ss, &
               calendar = LIS_calendar, &
               rc = status)
          deltaT=time-LIS_lai(n)%time1
          deltaTinterval=LIS_lai(n)%time2-LIS_lai(n)%time1

          call ESMF_TimeIntervalGet(deltaT, s=t_delta)
          call ESMF_TimeIntervalGet(deltaTinterval, s=t_delta_interval)

          !get current values from linear interpolation
          wt1=real(t_delta_interval-t_delta)/real(t_delta_interval)
          wt2=real(t_delta)/real(t_delta_interval)
          do i=1,LIS_rc%ntiles(n)
             LIS_lai(n)%tlai(i)=wt1*LIS_lai(n)%lai1(i)+wt2*LIS_lai(n)%lai2(i)
          enddo

          if (time >= LIS_lai(n)%time2) then
             LIS_lai(n)%lai1=LIS_lai(n)%lai2
             LIS_lai(n)%time1=LIS_lai(n)%time2

             forward_search=.true.
             call readlai(trim(LIS_rc%uselaimap(n))//char(0),n,&
                  forward_search,LIS_lai(n)%lai2,LIS_lai(n)%time2)
          end if
#endif
!       endif
    TRACE_EXIT("lai_read")

  end subroutine LIS_read_lai

#endif

!BOP
!
! !ROUTINE: LIS_lai_setup
! \label{LIS_lai_setup}
!
! !INTERFACE:
  subroutine LIS_lai_setup
! !USES:
    use LIS_coreMod,    only : LIS_rc, LIS_Config, LIS_domain
    use LIS_timeMgrMod, only : LIS_computeTemporalWeights, &
         LIS_registerAlarm
    use LIS_logMod,     only : LIS_verify, LIS_logunit

! !DESCRIPTION:
!
! Allocates memory and other structures for reading
! LAI datasets
!
!  The routines invoked are:
!  \begin{description}
!   \item[laisetup](\ref{laisetup}) \newline
!    initializes the realtime lai reader
!   \item[LIS\_registerAlarm](\ref{LIS_registerAlarm}) \newline
!    registers the alarm for reading LAI datasets
!   \item[LIS\_computeTemporalWeights](\ref{LIS_computeTemporalWeights})
!      \newline
!    computes the interpolation weights
!   \item[read\_LAIclimo](\ref{read_laiclimo}) \newline
!    reads the climatological lai data
!   \item[readlai](\ref{readlai}) \newline
!    reads the realtime lai data
!  \end{description}
!EOP
    implicit none
    integer :: n, i
    integer :: rc
    integer :: ndoms
    real, allocatable :: value1(:) ! temporary value holder for mo1
    real, allocatable :: value2(:) ! temporary value holder for mo2
    real          :: wt1,wt2
    integer       :: t1,t2

    TRACE_ENTER("lai_setup")
    ndoms = 0
    do n=1,LIS_rc%nnest
       if(LIS_rc%uselaimap(n).ne."none") then
          ndoms = ndoms+1
       endif
    enddo

    if(ndoms.gt.0) then
       allocate(LIS_lai(LIS_rc%nnest))
       do n=1,LIS_rc%nnest
          LIS_lai(n)%firstInstance = .true.
       enddo

       do n=1,LIS_rc%nnest
          allocate(LIS_lai(n)%lai1(LIS_rc%ntiles(n)))
          allocate(LIS_lai(n)%lai2(LIS_rc%ntiles(n)))
          allocate(LIS_lai(n)%tlai(LIS_rc%ntiles(n)))
          LIS_lai(n)%lai1 = 0
          LIS_lai(n)%lai2 = 0
          LIS_lai(n)%tlai = 0

          if(LIS_rc%uselaimap(n).ne."none".and.&
               LIS_rc%uselaimap(n).ne."LDT") then
             call laisetup(trim(LIS_rc%uselaimap(n))//char(0),n)
          else
             LIS_lai(n)%laiInterval = 2592000
             LIS_lai(n)%laiIntervalType = "monthly"
          endif

          call LIS_registerAlarm("LIS LAI read alarm",&
               LIS_rc%ts, LIS_lai(n)%laiInterval, &
               intervalType = LIS_lai(n)%laiIntervalType)

          call LIS_computeTemporalWeights(LIS_rc,&
               LIS_lai(n)%laiIntervalType, t1,t2,wt1,wt2)

          if(LIS_rc%uselaimap(n).eq."LDT") then
             allocate(value1(LIS_rc%ntiles(n)))
             allocate(value2(LIS_rc%ntiles(n)))

             call read_LAIclimo(n,t1,value1)
             call read_LAIclimo(n,t2,value2)

             do i=1,LIS_rc%ntiles(n)
                LIS_lai(n)%lai1(i) = value1(i)
                LIS_lai(n)%lai2(i) = value2(i)
             enddo
             deallocate(value1)
             deallocate(value2)

          else
             call readlai(trim(LIS_rc%uselaimap(n))//char(0), n, &
                  wt1, wt2, LIS_lai(n)%lai1, LIS_lai(n)%lai2)
          endif

          do i=1,LIS_rc%ntiles(n)
             LIS_lai(n)%tlai(i) = wt1* LIS_lai(n)%lai1(i)+ &
                  wt2*LIS_lai(n)%lai2(i)
#if 0
!SVK : Should be cleaned up -- hardcoding checks for LSM is not clean
!  YDT: 10/18/2011 added a check for minimal lai values
!    Otherwise noah will crash in CANRES()
!             if ( LIS_rc%gfracsrc(n) .gt. 0 .and. LIS_gfrac(n)%greenness(i) .gt. 0) &
!                LIS_lai(n)%tlai(i) = max(LIS_lai(n)%tlai(i), 0.1)
!  YDT: 10/20/2011 more "elegant" solution: compute greenness from lai if lai is MODIS-RT (src=3).
             if ( LIS_rc%uselaimap(n) .eq. "MODIS" ) then !??
                LIS_gfrac(n)%greenness(i) = 1.0 - exp(-0.52 * LIS_lai(n)%tlai(i) )
                ! for Urban LC, Noah will assign gfrac=0.05. So need to set a min lai to not crash
                ! CNARES()
                if (LIS_domain(n)%tile(i)%vegt .eq. LIS_rc%urbanclass .and. LIS_rc%lsm .eq. "NOAH32" ) &
                   LIS_lai(n)%tlai(i) = max(LIS_lai(n)%tlai(i), 0.1)
             end if
#endif
          end do
       enddo
    endif
    TRACE_EXIT("lai_setup")

  end subroutine LIS_lai_setup

!BOP
!
! !ROUTINE: LIS_read_lai
! \label{LIS_read_lai}
!
! !INTERFACE:
  subroutine LIS_read_lai(n)
! !USES:
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_isAlarmRinging,LIS_computeTemporalWeights
    use LIS_logMod,     only : LIS_verify, LIS_logunit

    implicit none
! !ARGUMENTS:
    integer, intent(in) :: n
!
! !DESCRIPTION:
!
!  Reads the LAI climatology and temporally interpolates it to the
!  current day.
!
!  The arguments are:
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!   \item[readlai](\ref{readlai}) \newline
!    invokes the generic method in the registry to read the
!    LAI climatology data
!  \end{description}
!
!EOP
    real, allocatable :: value1(:) ! temporary value holder for mo1
    real, allocatable :: value2(:) ! temporary value holder for mo2
    real              :: wt1,wt2
    integer           :: t1, t2
    integer           :: i
    logical           :: laiAlarmCheck

    TRACE_ENTER("lai_read")
    if(LIS_rc%uselaimap(n).ne."none") then
       laiAlarmCheck = LIS_isAlarmRinging(LIS_rc,&
            "LIS LAI read alarm",&
            LIS_lai(n)%laiIntervalType)
       call LIS_computeTemporalWeights(LIS_rc,&
            LIS_lai(n)%laiIntervalType, t1,t2,wt1,wt2)

       if(laiAlarmCheck) then

          if(LIS_rc%uselaimap(n).eq."LDT") then
             allocate(value1(LIS_rc%ntiles(n)))
             allocate(value2(LIS_rc%ntiles(n)))

             call read_LAIclimo(n,t1,value1)
             call read_LAIclimo(n,t2,value2)

             do i=1,LIS_rc%ntiles(n)
                LIS_lai(n)%lai1(i) = value1(i)
                LIS_lai(n)%lai2(i) = value2(i)
             enddo
             deallocate(value1)
             deallocate(value2)
          else
             call readlai(trim(LIS_rc%uselaimap(n))//char(0), n, &
                  wt1, wt2, LIS_lai(n)%lai1, LIS_lai(n)%lai2)
          endif


       endif
       do i=1,LIS_rc%ntiles(n)
          LIS_lai(n)%tlai(i) = wt1* LIS_lai(n)%lai1(i)+ &
               wt2*LIS_lai(n)%lai2(i)
#if 0
!  YDT: 10/18/2011 added a check for minimal lai values
!    Otherwise noah will crash in CANRES()
!             if ( LIS_rc%gfracsrc(n) .gt. 0 .and. LIS_gfrac(n)%greenness(i) .gt. 0) &
!                LIS_lai(n)%tlai(i) = max(LIS_lai(n)%tlai(i), 0.1)
!  YDT: 10/20/2011 more "elegant" solution: compute greenness from lai if lai is MODIS-RT (src=3).
             if ( LIS_rc%uselaimap(n) .eq. "MODIS") then !??
                LIS_gfrac(n)%greenness(i) = 1.0 - exp(-0.52 * LIS_lai(n)%tlai(i) )
                ! for Urban LC, Noah will assign gfrac=0.05. So need to set a min lai to not crash
                ! CNARES()
                if (LIS_domain(n)%tile(i)%vegt .eq. LIS_rc%urbanclass ) &
                   LIS_lai(n)%tlai(i) = max(LIS_lai(n)%tlai(i), 0.1)
             end if
#endif
       end do
    endif
    TRACE_EXIT("lai_read")
  end subroutine LIS_read_lai

!BOP
!
! !ROUTINE: LIS_lai_finalize
! \label{LIS_lai_finalize}
!
! !INTERFACE:
  subroutine LIS_lai_finalize
! !USES:
    use LIS_coreMod, only : LIS_rc
!
! !DESCRIPTION:
!
! Deallocates objects created in this module
!
!EOP
    implicit none

    integer :: n
    integer :: ndoms

    ndoms = 0
    do n=1,LIS_rc%nnest
       if(LIS_rc%uselaimap(n).ne."none") then
          ndoms = ndoms+1
       endif
    enddo

    if(ndoms.gt.0) then
       do n=1,LIS_rc%nnest
          deallocate(LIS_lai(n)%lai1)
          deallocate(LIS_lai(n)%lai2)
          deallocate(LIS_lai(n)%tlai)
       enddo
       deallocate(LIS_lai)
    endif

  end subroutine LIS_lai_finalize

!BOP
!
! !ROUTINE: LIS_diagnoseLAI
! \label{LIS_diagnoseLAI}
!
! !INTERFACE:
  subroutine LIS_diagnoseLAI(n)
! !USES:
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_histDataMod, only : LIS_diagnoseSurfaceOutputVar, LIS_MOC_LAI

! !ARGUMENTS:
    implicit none
    integer, intent(in)   :: n

! !DESCRIPTION:
!  This routine writes the LIS LAI to the parameter
!  output file.
!
!  The arguments are:
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[ftn] file unit number to be used \newline
!  \end{description}
!
!  The routines called are:
!  \begin{description}
!  \item[LIS\_diagnoseOutputVar] (\ref{LIS_diagnoseSurfaceOutputVar})
!     \newline
!   This routine maps a variable to the history writing routines
!  \end{description}
!EOP
    real    :: temp(LIS_rc%ntiles(n))
    integer :: t

    TRACE_ENTER("lai_diag")
    if(LIS_rc%uselaimap(n).ne."none") then
       temp = LIS_rc%udef
       do t=1,LIS_rc%ntiles(n)
          if(LIS_domain(n)%tile(t)%index.ne.-1) then
             temp(t) = LIS_lai(n)%tlai(t)
          endif
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_LAI, vlevel=1,&
               value=temp(t),unit="-", &
               direction="-")
       enddo
    endif
    TRACE_EXIT("lai_diag")
  end subroutine LIS_diagnoseLAI

!BOP
!
! !ROUTINE: LIS_lai_reset
! \label{LIS_lai_reset}
!
! !INTERFACE:
  subroutine LIS_lai_reset
! !USES:
    use LIS_coreMod,    only : LIS_rc, LIS_Config, LIS_domain
    use LIS_timeMgrMod, only : LIS_computeTemporalWeights, &
         LIS_registerAlarm, LIS_calendar
    use LIS_logMod,     only : LIS_verify, LIS_logunit
! !DESCRIPTION:
!
! Resets data structures for reading
! LAI datasets
!
!  The routines invoked are:
!  \begin{description}
!   \item[laisetup](\ref{laisetup}) \newline
!    initializes the realtime lai reader
!   \item[LIS\_registerAlarm](\ref{LIS_registerAlarm}) \newline
!    registers the alarm for reading LAI datasets
!   \item[LIS\_computeTemporalWeights](\ref{LIS_computeTemporalWeights})
!      \newline
!    computes the interpolation weights
!   \item[read\_LAIclimo](\ref{read_laiclimo}) \newline
!    reads the climatological lai data
!   \item[readlai](\ref{readlai}) \newline
!    reads the realtime lai data
!  \end{description}
!EOP
    implicit none
    integer :: n, i
    integer :: rc,ios,nid
    logical :: file_exists
    integer :: ndoms
    real, allocatable :: value1(:) ! temporary value holder for mo1
    real, allocatable :: value2(:) ! temporary value holder for mo2
    real          :: wt1,wt2
    integer       :: t1,t2, t_delta, t_delta_interval
    type(ESMF_Time)   :: time
    type(ESMF_TimeInterval) :: deltaT, deltaTinterval
    logical:: forward_search
    integer :: status

    TRACE_ENTER("lai_reset")
    ndoms = 0
    do n=1,LIS_rc%nnest
       if(LIS_rc%uselaimap(n).ne."none") then
          ndoms = ndoms+1
       endif
    enddo

    if(ndoms.gt.0) then

       do n=1,LIS_rc%nnest

          LIS_lai(n)%lai1 = 0
          LIS_lai(n)%lai2 = 0
          LIS_lai(n)%tlai = 0

          if(LIS_rc%uselaimap(n).ne."LDT") then
             call laisetup(trim(LIS_rc%uselaimap(n))//char(0),n)
          else
#if (defined USE_NETCDF3 || defined USE_NETCDF4)

             inquire(file=LIS_rc%paramfile(n), exist=file_exists)
             if(file_exists) then

                ios = nf90_open(path=trim(LIS_rc%paramfile(n)),&
                     mode=NF90_NOWRITE,ncid=nid)
                call LIS_verify(ios,'Error in nf90_open in read_laiclimo')

                ios = nf90_get_att(nid, NF90_GLOBAL, 'LAISAI_DATA_INTERVAL', &
                     LIS_lai(n)%laiIntervalType)
                call LIS_verify(ios,'Error in nf90_get_att in read_laiclimo')
                ios = nf90_close(nid)
                call LIS_verify(ios,'Error in nf90_close in read_laiclimo')
             else
                write(LIS_logunit,*) '[ERR] ',LIS_rc%paramfile(n), &
                     ' does not exist'
                write(LIS_logunit,*) '[ERR] program stopping ...'
                call LIS_endrun
             endif
#endif
             if(LIS_lai(n)%laiIntervalType.eq."monthly") then
                LIS_lai(n)%laiInterval = 2592000  ! 30 days
             endif
          endif

          call LIS_registerAlarm("LIS LAI read alarm",&
               LIS_rc%ts, LIS_lai(n)%laiInterval, &
               intervalType = LIS_lai(n)%laiIntervalType)

          call LIS_computeTemporalWeights(LIS_rc,&
               LIS_lai(n)%laiIntervalType, t1,t2,wt1,wt2)

          allocate(value1(LIS_rc%ntiles(n)))
          allocate(value2(LIS_rc%ntiles(n)))

          if(LIS_rc%uselaimap(n).eq."LDT") then
             call read_LAIclimo(n,t1,value1)
             call read_LAIclimo(n,t2,value2)
          else
             !EMK Corrected function call.
             !call readlai(trim(LIS_rc%uselaimap(n))//char(0),n,t1,value1)
             !call readlai(trim(LIS_rc%uselaimap(n))//char(0),n,t2,value2)
             call readlai(trim(LIS_rc%uselaimap(n))//char(0), &
                  n, wt1, wt2, value1, value2)
          endif

          do i=1,LIS_rc%ntiles(n)
             LIS_lai(n)%lai1(i) = value1(i)
             LIS_lai(n)%lai2(i) = value2(i)
          enddo
          deallocate(value1)
          deallocate(value2)

          do i=1,LIS_rc%ntiles(n)
             LIS_lai(n)%tlai(i) = wt1* LIS_lai(n)%lai1(i)+ &
                  wt2*LIS_lai(n)%lai2(i)
          end do
       enddo
    endif
    TRACE_EXIT("lai_reset")
  end subroutine LIS_lai_reset

!BOP
!
! !ROUTINE: LIS_sai_setup
! \label{LIS_sai_setup}
!
! !INTERFACE:
  subroutine LIS_sai_setup
! !USES:
    use LIS_coreMod,    only : LIS_rc, LIS_config
    use LIS_timeMgrMod, only : LIS_computeTemporalWeights, &
         LIS_registerAlarm
    use LIS_logMod,     only : LIS_verify
! !DESCRIPTION:
!
! Allocates memory and other structures for reading
! SAI datasets
!
!  The routines invoked are:
!  \begin{description}
!   \item[saisetup](\ref{saisetup}) \newline
!    initializes the realtime sai reader
!   \item[LIS\_registerAlarm](\ref{LIS_registerAlarm}) \newline
!    registers the alarm for reading SAI datasets
!   \item[LIS\_computeTemporalWeights](\ref{LIS_computeTemporalWeights})
!      \newline
!    computes the interpolation weights
!   \item[read\_saiclimo](\ref{read_saiclimo}) \newline
!    reads the climatological sai data
!   \item[readsai](\ref{readsai}) \newline
!    reads the realtime sai data
!  \end{description}
!EOP
    implicit none
    integer :: n, i
    integer :: rc,ios,nid
    logical :: file_exists
    integer :: ndoms
    real, allocatable :: value1(:) ! temporary value holder for mo1
    real, allocatable :: value2(:) ! temporary value holder for mo2
    real              :: wt1,wt2
    integer           :: t1,t2

    TRACE_ENTER("sai_setup")
    ndoms = 0
    do n=1,LIS_rc%nnest
       if(LIS_rc%usesaimap(n).ne."none") then
          ndoms = ndoms+1
       endif
    enddo

    if(ndoms.gt.0) then
       allocate(LIS_sai(LIS_rc%nnest))

       do n=1,LIS_rc%nnest
          allocate(LIS_sai(n)%sai1(LIS_rc%ntiles(n)))
          allocate(LIS_sai(n)%sai2(LIS_rc%ntiles(n)))
          allocate(LIS_sai(n)%tsai(LIS_rc%ntiles(n)))
          LIS_sai(n)%sai1 = 0
          LIS_sai(n)%sai2 = 0
          LIS_sai(n)%tsai = 0

          if(LIS_rc%usesaimap(n).ne."LDT") then
             call saisetup(trim(LIS_rc%usesaimap(n))//char(0),n)
          else
#if (defined USE_NETCDF3 || defined USE_NETCDF4)

             inquire(file=LIS_rc%paramfile(n), exist=file_exists)
             if(file_exists) then

                ios = nf90_open(path=trim(LIS_rc%paramfile(n)),&
                     mode=NF90_NOWRITE,ncid=nid)
                call LIS_verify(ios,'Error in nf90_open in read_saiclimo')

                ios = nf90_get_att(nid, NF90_GLOBAL, 'LAISAI_DATA_INTERVAL', &
                     LIS_sai(n)%saiIntervalType)
                call LIS_verify(ios,'Error in nf90_get_att in read_saiclimo')

                ios = nf90_close(nid)
                call LIS_verify(ios,'Error in nf90_close in read_saiclimo')
             else
                write(LIS_logunit,*) '[ERR] ',LIS_rc%paramfile(n), &
                     ' does not exist'
                write(LIS_logunit,*) '[ERR] program stopping ...'
                call LIS_endrun
             endif
#endif
             if(LIS_sai(n)%saiIntervalType.eq."monthly") then
                LIS_sai(n)%saiInterval = 2592000  ! 30 days
             endif
          endif

          call LIS_registerAlarm("LIS SAI read alarm",&
               LIS_rc%ts, LIS_sai(n)%saiInterval, &
               intervalType = LIS_sai(n)%saiIntervalType)

          call LIS_computeTemporalWeights(LIS_rc,&
               LIS_sai(n)%saiIntervalType, t1,t2,wt1,wt2)

          allocate(value1(LIS_rc%ntiles(n)))
          allocate(value2(LIS_rc%ntiles(n)))

          if(LIS_rc%usesaimap(n).eq."LDT") then
             call read_saiclimo(n,t1,value1)
             call read_saiclimo(n,t2,value2)
          else
             call readsai(trim(LIS_rc%usesaimap(n))//char(0),n,t1,value1)
             call readsai(trim(LIS_rc%usesaimap(n))//char(0),n,t2,value2)
          endif

          do i=1,LIS_rc%ntiles(n)
             LIS_sai(n)%sai1(i) = value1(i)
             LIS_sai(n)%sai2(i) = value2(i)
          enddo
          deallocate(value1)
          deallocate(value2)

          do i=1,LIS_rc%ntiles(n)
             LIS_sai(n)%tsai(i) = wt1* LIS_sai(n)%sai1(i)+ &
                  wt2*LIS_sai(n)%sai2(i)
          end do
       enddo
    endif
    TRACE_EXIT("sai_setup")

  end subroutine LIS_sai_setup

!BOP
!
! !ROUTINE: LIS_read_sai
! \label{LIS_read_sai}
!
! !INTERFACE:
  subroutine LIS_read_sai(n)
! !USES:
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_computeTemporalWeights, &
         LIS_isAlarmRinging

    implicit none
! !ARGUMENTS:
    integer, intent(in) :: n
!
! !DESCRIPTION:
!
!  Reads the SAI climatology and temporally interpolates it to the
!  current day.
!
!  The arguments are:
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!   \item[readsai](\ref{readsai}) \newline
!    invokes the generic method in the registry to read the
!    SAI climatology data
!  \end{description}
!
!EOP
    real, allocatable :: value1(:) ! temporary value holder for mo1
    real, allocatable :: value2(:) ! temporary value holder for mo2
    real              :: wt1,wt2
    integer           :: t1,t2
    integer           :: i
    logical           :: saiAlarmCheck

    TRACE_ENTER("sai_read")
    if(LIS_rc%usesaimap(n).ne."none") then
       saiAlarmCheck = LIS_isAlarmRinging(LIS_rc,&
            "LIS SAI read alarm",&
            LIS_sai(n)%saiIntervalType)

       call LIS_computeTemporalWeights(LIS_rc,&
            LIS_sai(n)%saiIntervalType, t1,t2,wt1,wt2)

       if(saiAlarmCheck) then
          allocate(value1(LIS_rc%ntiles(n)))
          allocate(value2(LIS_rc%ntiles(n)))

          if(LIS_rc%usesaimap(n).eq."LDT") then
             call read_SAIclimo(n,t1,value1)
             call read_SAIclimo(n,t2,value2)
          else
             call readsai(trim(LIS_rc%usesaimap(n))//char(0),n,t1,value1)
             call readsai(trim(LIS_rc%usesaimap(n))//char(0),n,t2,value2)
          endif

          do i=1,LIS_rc%ntiles(n)
             LIS_sai(n)%sai1(i) = value1(i)
             LIS_sai(n)%sai2(i) = value2(i)
          enddo
          deallocate(value1)
          deallocate(value2)

       endif
       do i=1,LIS_rc%ntiles(n)
          LIS_sai(n)%tsai(i) = wt1* LIS_sai(n)%sai1(i)+ &
               wt2*LIS_sai(n)%sai2(i)
       end do
    endif
    TRACE_EXIT("sai_read")
  end subroutine LIS_read_sai

!BOP
!
! !ROUTINE: LIS_sai_finalize
! \label{LIS_sai_finalize}
!
! !INTERFACE:
  subroutine LIS_sai_finalize
! !USES:
    use LIS_coreMod, only : LIS_rc
!
! !DESCRIPTION:
!
! Deallocates objects created in this module
!
!EOP
    implicit none

    integer :: n
    integer :: ndoms

    ndoms = 0
    do n=1,LIS_rc%nnest
       if(LIS_rc%usesaimap(n).ne."none") then
          ndoms = ndoms+1
       endif
    enddo

    if(ndoms.gt.0) then
       do n=1,LIS_rc%nnest
          deallocate(LIS_sai(n)%sai1)
          deallocate(LIS_sai(n)%sai2)
          deallocate(LIS_sai(n)%tsai)
       enddo
       deallocate(LIS_sai)
    endif
  end subroutine LIS_sai_finalize

!BOP
!
! !ROUTINE: LIS_diagnoseSAI
! \label{LIS_diagnoseSAI}
!
! !INTERFACE:
  subroutine LIS_diagnoseSAI(n)
! !USES:
    use LIS_coreMod,     only : LIS_rc, LIS_domain
    use LIS_histDataMod, only : LIS_diagnoseSurfaceOutputVar, LIS_MOC_SAI

! !ARGUMENTS:
    implicit none
    integer, intent(in)   :: n

! !DESCRIPTION:
!  This routine diagnoses the LIS SAI to the parameter
!  output file.
!
!  The arguments are:
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[ftn] file unit number to be used \newline
!  \end{description}
!
!  The routines called are:
!  \begin{description}
!  \item[LIS\_diagnoseOutputVar] (\ref{LIS_diagnoseSurfaceOutputVar})
!     \newline
!   This routine maps a variable to the history writing routines
!  \end{description}
!EOP
    real    :: temp(LIS_rc%ntiles(n))
    integer :: t
    TRACE_ENTER("sai_diag")
    if(LIS_rc%usesaimap(n).ne."none") then
       temp = LIS_rc%udef
       do t=1,LIS_rc%ntiles(n)
          if(LIS_domain(n)%tile(t)%index.ne.-1) then
             temp(t) = LIS_sai(n)%tsai(t)
          endif
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SAI,vlevel=1,&
               value=temp(t),unit="-", &
               direction="-")
       enddo
    endif
    TRACE_EXIT("sai_diag")

  end subroutine LIS_diagnoseSAI

!BOP
!
! !ROUTINE: LIS_sai_reset
! \label{LIS_sai_reset}
!
! !INTERFACE:
  subroutine LIS_sai_reset
! !USES:
    use LIS_coreMod,    only : LIS_rc
    use LIS_timeMgrMod, only : LIS_computeTemporalWeights, &
         LIS_Clock
    use LIS_logMod,     only : LIS_verify
! !DESCRIPTION:
!
! resets the structures that hold
! SAI datasets
!
!  The routines invoked are:
!  \begin{description}
!   \item[saisetup](\ref{saisetup}) \newline
!    initializes the realtime sai reader
!   \item[LIS\_computeTemporalWeights](\ref{LIS_computeTemporalWeights})
!     \newline
!    computes the interpolation weights
!   \item[read\_saiclimo](\ref{read_saiclimo}) \newline
!    reads the climatological sai data
!   \item[readsai](\ref{readsai}) \newline
!    reads the realtime sai data
!  \end{description}
!EOP
    implicit none
    integer :: n, i
    integer :: rc,ios,nid
    logical :: file_exists
    integer :: ndoms
    real, allocatable :: value1(:) ! temporary value holder for mo1
    real, allocatable :: value2(:) ! temporary value holder for mo2
    real              :: wt1,wt2
    integer           :: t1,t2

    TRACE_ENTER("sai_reset")
    ndoms = 0
    do n=1,LIS_rc%nnest
       if(LIS_rc%usesaimap(n).ne."none") then
          ndoms = ndoms+1
       endif
    enddo

    if(ndoms.gt.0) then
       do n=1,LIS_rc%nnest

          LIS_sai(n)%sai1 = 0
          LIS_sai(n)%sai2 = 0
          LIS_sai(n)%tsai = 0

          if(LIS_rc%usesaimap(n).ne."LDT") then
             call saisetup(trim(LIS_rc%usesaimap(n))//char(0),n)
          else
#if (defined USE_NETCDF3 || defined USE_NETCDF4)

             inquire(file=LIS_rc%paramfile(n), exist=file_exists)
             if(file_exists) then

                ios = nf90_open(path=trim(LIS_rc%paramfile(n)),&
                     mode=NF90_NOWRITE,ncid=nid)
                call LIS_verify(ios,'Error in nf90_open in read_saiclimo')

                ios = nf90_get_att(nid, NF90_GLOBAL, 'LAISAI_DATA_INTERVAL', &
                     LIS_sai(n)%saiIntervalType)
                call LIS_verify(ios,'Error in nf90_get_att in read_saiclimo')

                ios = nf90_close(nid)
                call LIS_verify(ios,'Error in nf90_close in read_saiclimo')
             else
                write(LIS_logunit,*) '[ERR] ',LIS_rc%paramfile(n), &
                     ' does not exist'
                write(LIS_logunit,*) '[ERR] program stopping ...'
                call LIS_endrun
             endif
#endif
             if(LIS_sai(n)%saiIntervalType.eq."monthly") then
                LIS_sai(n)%saiInterval = 2592000  ! 30 days
             endif

          endif

          call LIS_computeTemporalWeights(LIS_rc,&
               LIS_sai(n)%saiIntervalType, t1,t2,wt1,wt2)

          allocate(value1(LIS_rc%ntiles(n)))
          allocate(value2(LIS_rc%ntiles(n)))

          if(LIS_rc%usesaimap(n).eq."LDT") then
             call read_saiclimo(n,t1,value1)
             call read_saiclimo(n,t2,value2)
          else
             call readsai(trim(LIS_rc%usesaimap(n))//char(0),n,t1,value1)
             call readsai(trim(LIS_rc%usesaimap(n))//char(0),n,t2,value2)
          endif

          do i=1,LIS_rc%ntiles(n)
             LIS_sai(n)%sai1(i) = value1(i)
             LIS_sai(n)%sai2(i) = value2(i)
          enddo
          deallocate(value1)
          deallocate(value2)

          do i=1,LIS_rc%ntiles(n)
             LIS_sai(n)%tsai(i) = wt1* LIS_sai(n)%sai1(i)+ &
                  wt2*LIS_sai(n)%sai2(i)
          end do
       enddo
    endif
    TRACE_EXIT("sai_reset")

  end subroutine LIS_sai_reset

!BOP
!
! !ROUTINE: read_laiclimo
!  \label{read_laiclimo}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
  subroutine read_laiclimo(n,time,array)
! !USES:

    use LIS_coreMod,        only : LIS_rc, LIS_domain, LIS_localPet,&
         LIS_ews_ind, LIS_ewe_ind,&
         LIS_nss_ind, LIS_nse_ind, LIS_ews_halo_ind,LIS_ewe_halo_ind, &
         LIS_nss_halo_ind, LIS_nse_halo_ind
    use LIS_logMod,         only : LIS_logunit, LIS_getNextUnitNumber, &
         LIS_releaseUnitNumber, LIS_endrun, LIS_verify

    implicit none
! !ARGUMENTS:
    integer, intent(in)         :: n
    integer, intent(in)         :: time
    real, intent(inout)         :: array(LIS_rc%ntiles(n))
! !DESCRIPTION:
!  This subroutine reads the lai data climatology
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \end{description}
!
!EOP

    integer :: ios1
    integer :: ios,nid,laiid,ncId, nrId,mid
    integer :: nc,nr,t,months
    integer :: mo
    real, allocatable :: lai(:,:,:)
    real    :: locallai(LIS_rc%lnc(n),LIS_rc%lnr(n))
    logical :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    mo = time

    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then

       write(LIS_logunit,*)'[INFO] Reading LAI map for month ', mo, &
            'from ',trim(LIS_rc%paramfile(n))

       ios = nf90_open(path=trim(LIS_rc%paramfile(n)),&
            mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in read_laiclimo')

       ios = nf90_inq_dimid(nid,"east_west",ncId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_laiclimo')

       ios = nf90_inq_dimid(nid,"north_south",nrId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_laiclimo')

       ios = nf90_inq_dimid(nid,"month",mId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_laiclimo')

       ios = nf90_inquire_dimension(nid,ncId, len=nc)
       call LIS_verify(ios,'Error in nf90_inquire_dimension in read_laiclimo')

       ios = nf90_inquire_dimension(nid,nrId, len=nr)
       call LIS_verify(ios,'Error in nf90_inquire_dimension in read_laiclimo')

       ios = nf90_inquire_dimension(nid,mId, len=months)
       call LIS_verify(ios,'Error in nf90_inquire_dimension in read_laiclimo')

       ios = nf90_inq_varid(nid,'LAI',laiid)
       call LIS_verify(ios,'LAI field not found in the LIS param file')

       ios = nf90_get_var(nid,laiid,locallai,&
            start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1),mo/),&
            count=(/LIS_rc%lnc(n),LIS_rc%lnr(n),1/))
       call LIS_verify(ios,'Error in nf90_get_var in read_laiclimo')

       ios = nf90_close(nid)
       call LIS_verify(ios,'Error in nf90_close in read_laiclimo')

       do t=1,LIS_rc%ntiles(n)
          array(t) = locallai(LIS_domain(n)%tile(t)%col,&
               LIS_domain(n)%tile(t)%row)
       enddo
    else
       write(LIS_logunit,*) '[ERR] lai map: ',LIS_rc%paramfile(n), &
            ' does not exist'
       write(LIS_logunit,*) '[ERR] program stopping ...'
       call LIS_endrun
    endif
#endif
  end subroutine read_laiclimo

!BOP
!
! !ROUTINE: read_saiclimo
!  \label{read_saiclimo}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
  subroutine read_saiclimo(n,time,array)
! !USES:

    use LIS_coreMod,        only : LIS_rc, LIS_domain, LIS_localPet,&
         LIS_ews_ind, LIS_ewe_ind,&
         LIS_nss_ind, LIS_nse_ind, LIS_ews_halo_ind,LIS_ewe_halo_ind, &
         LIS_nss_halo_ind, LIS_nse_halo_ind
    use LIS_logMod,         only : LIS_logunit, LIS_getNextUnitNumber, &
         LIS_releaseUnitNumber, LIS_endrun, LIS_verify

    implicit none
! !ARGUMENTS:
    integer, intent(in)         :: n
    integer, intent(in)         :: time
    real, intent(inout)         :: array(LIS_rc%ntiles(n))
! !DESCRIPTION:
!  This subroutine reads the greenness data climatology
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \end{description}
!
!EOP

    integer :: ios1
    integer :: ios,nid,saiid,ncId, nrId,mid
    integer :: nc,nr,t,months
    integer :: mo
    real, allocatable :: sai(:,:,:)
    real    :: localsai(LIS_rc%lnc(n),LIS_rc%lnr(n))
    logical :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    mo = time

    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then

       write(LIS_logunit,*)'[INFO] Reading SAI map for month ', mo, &
            'from ',trim(LIS_rc%paramfile(n))

       ios = nf90_open(path=trim(LIS_rc%paramfile(n)),&
            mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in read_saiclimo')

       ios = nf90_inq_dimid(nid,"east_west",ncId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_saiclimo')

       ios = nf90_inq_dimid(nid,"north_south",nrId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_saiclimo')

       ios = nf90_inq_dimid(nid,"month",mId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_saiclimo')

       ios = nf90_inquire_dimension(nid,ncId, len=nc)
       call LIS_verify(ios,'Error in nf90_inquire_dimension in read_saiclimo')

       ios = nf90_inquire_dimension(nid,nrId, len=nr)
       call LIS_verify(ios,'Error in nf90_inquire_dimension in read_saiclimo')

       ios = nf90_inquire_dimension(nid,mId, len=months)
       call LIS_verify(ios,'Error in nf90_inquire_dimension in read_saiclimo')

       ios = nf90_inq_varid(nid,'SAI',saiid)
       call LIS_verify(ios,'SAI field not found in the LIS param file')

       ios = nf90_get_var(nid,saiid,localsai,&
            start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1),mo/),&
            count=(/LIS_rc%lnc(n),LIS_rc%lnr(n),1/))
       call LIS_verify(ios,'Error in nf90_get_var in read_saiclimo')

       ios = nf90_close(nid)
       call LIS_verify(ios,'Error in nf90_close in read_saiclimo')

       do t=1,LIS_rc%ntiles(n)
          array(t) = localsai(LIS_domain(n)%tile(t)%col,&
               LIS_domain(n)%tile(t)%row)
       enddo
    else
       write(LIS_logunit,*) '[ERR] sai map: ',LIS_rc%paramfile(n), &
            ' does not exist'
       write(LIS_logunit,*) '[ERR] program stopping ...'
       call LIS_endrun
    endif
#endif
  end subroutine read_saiclimo

end module LIS_vegDataMod
