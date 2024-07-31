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
module LIS_albedoMod

!BOP
!
! !MODULE: LIS_albedoMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read various sources of
!  albedo parameter data.
!  \subsubsection{Overview}
!
!  This routines in this module reads two sources of albedo:
!   \begin{itemize}
!     \item{Albedo climatology} \newline
!     \item{Static, maximum albedo expected over deep snow} \newline
!   \end{itemize}
!  This module provides routines to read the albedo climatology (monthly)
!  data and allows the users to specify the frequency of climatology
!  (in months). The climatological data is temporally interpolated
!  between months to the current simulation date.
!
! !REVISION HISTORY:
!
!  8 Aug 2005: Sujay Kumar; Initial implementation
!
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_fileIOMod
  use LIS_timeMgrMod
  use LIS_histDataMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LIS_albedo_setup       !allocates memory for required variables
  public :: LIS_read_albedo        !reads climatological albedo
  public :: LIS_diagnosealbedo     !maps variables for history output
  public :: LIS_albedo_finalize    !cleanup allocated structures
  public :: LIS_albedo_reset       !resets the alarms and the data holders

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LIS_alb      !data structure containing albedo data

!EOP

  type, public :: alb_type_dec
     character(len=LIS_CONST_PATH_LEN) :: albfile
     logical       :: firstInstance
     real, allocatable :: albsf1(:)
     real, allocatable :: albsf2(:)
     real, allocatable :: albsf(:)
     real, allocatable :: mxsnalb(:)
     real, allocatable :: albedo(:)
     character*50  :: albIntervalType
     real          :: albInterval
     real*8        :: alarmTime
  end type alb_type_dec

  type(alb_type_dec), allocatable :: LIS_alb(:)



contains

!BOP
!
! !ROUTINE: LIS_albedo_setup
! \label{LIS_albedo_setup}
!
! !INTERFACE:
  subroutine LIS_albedo_setup
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
!
! !DESCRIPTION:
!
! Allocates memory for data structures for reading
! in albedo datasets
!
! The routines invoked are:
! \begin{description}
!  \item[LIS\_read\_mxsnalb](\ref{LIS_read_mxsnalb}) \newline
!    method to read the max snow albedo
!   \item[readalbedo](\ref{readalbedo}) \newline
!    invokes the generic method in the registry to read the
!    albedo climatology data
!   \item[albedosetup](\ref{albedosetup}) \newline
!    calls the registry to invoke the albedo setup method.
! \end{description}
!EOP
    implicit none
    integer         :: n
    integer         :: i
    integer         :: ndoms
    integer         :: t1, t2
    real            :: wt1, wt2
    integer         :: t,c,r
    real, allocatable :: value1(:,:) ! Temporary value holder for QQ1
    real, allocatable :: value2(:,:) ! Temporary value holder for QQ2
    integer         :: interval
    integer         :: rc,ios,nid
    logical         :: file_exists
!------------------------------------------------------------------------------
! If albedo datasets are to be used, the routine allocates memory for
! each domain. Further, the routine sets the alarms for reading the
! climatology datasets
!------------------------------------------------------------------------------
    TRACE_ENTER("alb_setup")
    allocate(lis_alb(LIS_rc%nnest))

    ndoms = 0
    do n=1,LIS_rc%nnest
       if(LIS_rc%usealbedomap(n).ne."none") ndoms = ndoms+1
    enddo

    if(ndoms.gt.0) then
       do n=1,LIS_rc%nnest
          LIS_alb(n)%firstInstance = .true.
       enddo
       do n=1,LIS_rc%nnest
#if (defined USE_NETCDF3 || defined USE_NETCDF4)

          inquire(file=LIS_rc%paramfile(n), exist=file_exists)
          if(file_exists) then
             ios = nf90_open(path=trim(LIS_rc%paramfile(n)),&
                  mode=NF90_NOWRITE,ncid=nid)
             call LIS_verify(ios,'Error in nf90_open in read_albclimo')

             ios = nf90_get_att(nid, NF90_GLOBAL, 'ALBEDO_DATA_INTERVAL', &
                  LIS_alb(n)%albIntervalType)
             call LIS_verify(ios,'Error in nf90_get_att in read_albclimo')

             ios = nf90_close(nid)
             call LIS_verify(ios,'Error in nf90_close in read_albclimo')
          else
             write(LIS_logunit,*) '[ERR] ',LIS_rc%paramfile(n), &
                  ' does not exist'
             write(LIS_logunit,*) '[ERR] program stopping ...'
             call LIS_endrun
          endif
#endif

          if(LIS_alb(n)%albIntervalType.eq."monthly") then
             LIS_alb(n)%albInterval = 2592000
          elseif(LIS_alb(n)%albIntervalType.eq."quarterly") then
             LIS_alb(n)%albInterval = 10368000
          endif
       enddo

       do n=1,LIS_rc%nnest

          allocate(lis_alb(n)%albsf1(LIS_rc%ntiles(n)))
          allocate(lis_alb(n)%albsf2(LIS_rc%ntiles(n)))
          allocate(lis_alb(n)%albsf(LIS_rc%ntiles(n)))
          allocate(lis_alb(n)%albedo(LIS_rc%ntiles(n)))

          lis_alb(n)%albsf1 = 0.0
          lis_alb(n)%albsf2 = 0.0
          lis_alb(n)%albsf = 0.0
          lis_alb(n)%albedo = 0.0

          if(LIS_rc%usealbedomap(n).ne."LDT") then
             call albedosetup(trim(LIS_rc%usealbedomap(n))//char(0),n)
          endif

          call LIS_registerAlarm("LIS albedo read alarm",&
               LIS_rc%ts, interval=LIS_alb(n)%albInterval,&
               intervalType = LIS_alb(n)%albIntervalType)

          call LIS_computeTemporalWeights(LIS_rc,&
               LIS_alb(n)%albIntervalType, t1,t2,wt1,wt2)

          allocate(value1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          allocate(value2(LIS_rc%lnc(n),LIS_rc%lnr(n)))

          if(LIS_rc%usealbedomap(n).eq."LDT") then
             call read_albedoclimo(n,t1,value1)
             call read_albedoclimo(n,t2,value2)


             do t=1,LIS_rc%ntiles(n)
                lis_alb(n)%albsf1(t) = value1(LIS_domain(n)%tile(t)%col,&
                     LIS_domain(n)%tile(t)%row)
                lis_alb(n)%albsf2(t) = value2(LIS_domain(n)%tile(t)%col,&
                     LIS_domain(n)%tile(t)%row)
             enddo
             deallocate(value1)
             deallocate(value2)
          else
             call albedosetup(trim(LIS_rc%usealbedomap(n))//char(0),n)
             call readalbedo(trim(LIS_rc%usealbedomap(n))//char(0), &
                  n, wt1, wt1, LIS_alb(n)%albsf1, LIS_alb(n)%albsf2)

          endif
!-------------------------------------------------------------------------
! Interpolate the albedo fractions to daily values
!-------------------------------------------------------------------------
          do i=1,LIS_rc%ntiles(n)
             if (lis_alb(n)%albsf1(i) .ne. -9999.000) then
                lis_alb(n)%albsf(i) = wt1 * lis_alb(n)%albsf1(i) + &
                     wt2 * lis_alb(n)%albsf2(i)
             endif
          end do
       enddo
    endif

    ndoms = 0
    do n=1,LIS_rc%nnest
       if(LIS_rc%usemxsnalbmap(n).eq."LDT") ndoms = ndoms+1
    enddo

    if(ndoms.gt.0) then
       do n=1,LIS_rc%nnest
          allocate(lis_alb(n)%mxsnalb(LIS_rc%ntiles(n)))
          lis_alb(n)%mxsnalb = 0.0
          call LIS_read_mxsnalb(n)
       enddo
    endif
    TRACE_EXIT("alb_setup")

  end subroutine LIS_albedo_setup

!BOP
!
! !ROUTINE: LIS_read_mxsnalb
! \label{LIS_read_mxsnalb}
!
! !INTERFACE:
  subroutine LIS_read_mxsnalb(n)
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n

!
! !DESCRIPTION:
!  Reads the static albedo upper bound over deep snow
!  for each domain.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[localmask]
!    mxsnalb for the region of interest
!   \end{description}
!
!EOP

  integer :: ios1
  integer :: ios,nid,mxsnalid,ncId, nrId
  integer :: nc,nr,t
  real    :: tmpalb(LIS_rc%lnc(n), LIS_rc%lnr(n))
  logical :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  inquire(file=LIS_rc%paramfile(n), exist=file_exists)
  if(file_exists) then

     write(LIS_logunit,*)'[INFO] Reading max snow albedo map from ',&
          trim(LIS_rc%paramfile(n))

     ios = nf90_open(path=LIS_rc%paramfile(n),&
          mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error in nf90_open in read_mxsnalb')

     ios = nf90_inq_dimid(nid,"east_west",ncId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in read_mxsnalb')

     ios = nf90_inq_dimid(nid,"north_south",nrId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in read_mxsnalb')

     ios = nf90_inquire_dimension(nid,ncId, len=nc)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in read_mxsnalb')

     ios = nf90_inquire_dimension(nid,nrId, len=nr)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in read_mxsnalb')

     ios = nf90_inq_varid(nid,'MXSNALBEDO',mxsnalid)
     call LIS_verify(ios,'MXSNALBEDO field not found in the LIS param file')

     ios = nf90_get_var(nid,mxsnalid,tmpalb,&
          start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
          LIS_nss_halo_ind(n,LIS_localPet+1)/),&
          count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/))
     call LIS_verify(ios,'Error in nf90_get_var in read_mxsnalb')

     ios = nf90_close(nid)
     call LIS_verify(ios,'Error in nf90_close in read_mxsnalb')

  else
     write(LIS_logunit,*) '[INFO] max snow albedo map: ', &
          LIS_rc%paramfile(n), ' does not exist'
     write(LIS_logunit,*) '[INFO] program stopping ...'
     call LIS_endrun
  endif

  do t=1,LIS_rc%ntiles(n)
     lis_alb(n)%mxsnalb(t) = tmpalb(LIS_domain(n)%tile(t)%col,&
          LIS_domain(n)%tile(t)%row)
  enddo
#endif

  end subroutine LIS_read_mxsnalb

!BOP
!
! !ROUTINE: LIS_read_albedo
! \label{LIS_read_albedo}
!
! !INTERFACE:
  subroutine LIS_read_albedo(n)
! !USES:

! !ARGUMENTS:
    integer, intent(in) :: n
!
! !DESCRIPTION:
!  Reads the snow free albedo climatology data and
!  temporally interpolates it to the current day
!
!  The arguments are:
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!  \end{description}
!
! The routines invoked are:
! \begin{description}
!   \item[LIS\_computeTemporalweights](\ref{LIS_computeTemporalWeights}) \newline
!    computes the temporal interpolation weights
!   \item[readalbedo](\ref{readalbedo}) \newline
!    invokes the generic method in the registry to read the
!    albedo climatology data
! \end{description}
!EOP
    real, allocatable :: value1(:,:) ! Temporary value holder for QQ1
    real, allocatable :: value2(:,:) ! Temporary value holder for QQ2

    integer         :: i, t
    integer         :: t1, t2
    real            :: wt1,wt2
    logical         :: albAlarmCheck

    TRACE_ENTER("alb_read")
    if(LIS_rc%usealbedomap(n).ne."none") then
       albAlarmCheck = LIS_isAlarmRinging(LIS_rc,&
            "LIS albedo read alarm",&
            LIS_alb(n)%albIntervalType)

       call LIS_computeTemporalWeights(LIS_rc,&
            LIS_alb(n)%albIntervalType, t1,t2,wt1,wt2)
       if(albAlarmCheck) then

          if(LIS_rc%usealbedomap(n).eq."LDT") then
             allocate(value1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
             allocate(value2(LIS_rc%lnc(n),LIS_rc%lnr(n)))

             call read_albedoclimo(n,t1,value1)
             call read_albedoclimo(n,t2,value2)

             do t=1,LIS_rc%ntiles(n)
                lis_alb(n)%albsf1(t) = value1(LIS_domain(n)%tile(t)%col,&
                     LIS_domain(n)%tile(t)%row)
                lis_alb(n)%albsf2(t) = value2(LIS_domain(n)%tile(t)%col,&
                     LIS_domain(n)%tile(t)%row)
             enddo
             deallocate(value1)
             deallocate(value2)

          else
             call readalbedo(trim(LIS_rc%usealbedomap(n))//char(0), &
                  n, wt1, wt2, LIS_alb(n)%albsf1, LIS_alb(n)%albsf2)
          endif

       endif
!-------------------------------------------------------------------------
! Interpolate the albedo fractions to daily values
!-------------------------------------------------------------------------
       do i=1,LIS_rc%ntiles(n)
          if (lis_alb(n)%albsf1(i) .ne. -9999.000) then
             lis_alb(n)%albsf(i) = wt1 * lis_alb(n)%albsf1(i) + &
                                   wt2 * lis_alb(n)%albsf2(i)
!-------------------------------------------------------------------------
! set albedo same as the snow free one, to be modified in agrmet forcing
!-------------------------------------------------------------------------
!             lis_alb(n)%albedo = lis_alb(n)%albsf(i)

          endif
       end do
    endif
  TRACE_EXIT("alb_read")
  end subroutine LIS_read_albedo

!BOP
! !ROUTINE: LIS_albedo_finalize
! \label{LIS_albedo_finalize}
!
! !INTERFACE: albedo_finalize
  subroutine LIS_albedo_finalize
!
! !USES:

    implicit none
!
! !DESCRIPTION:
!
! Deallocates objects created in this module
!
!EOP
    integer :: n
    integer :: ndoms

    ndoms = 0
    do n=1,LIS_rc%nnest
       if(LIS_rc%usealbedomap(n).ne."none") ndoms = ndoms+1
    enddo

    if(ndoms.gt.0) then
       do n=1,LIS_rc%nnest
          deallocate(lis_alb(n)%albsf1)
          deallocate(lis_alb(n)%albsf2)
          deallocate(lis_alb(n)%albsf)
          deallocate(lis_alb(n)%albedo)
       enddo
    endif

    ndoms = 0
    do n=1,LIS_rc%nnest
       if(LIS_rc%usemxsnalbmap(n).eq."LDT") ndoms = ndoms+1
    enddo

    if(ndoms.gt.0) then
       do n=1,LIS_rc%nnest
          deallocate(lis_alb(n)%mxsnalb)
       enddo
    endif

    deallocate(lis_alb)

  end subroutine LIS_albedo_finalize

!BOP
!
! !ROUTINE: LIS_diagnosealbedo
! \label{LIS_diagnosealbedo}
!
! !INTERFACE:
  subroutine LIS_diagnosealbedo(n)
! !USES:


! !ARGUMENTS:
    implicit none
    integer, intent(in)   :: n

! !DESCRIPTION:
!  This routine diagnoses the LIS albedo for history output.
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
    real    :: temp1, temp2
    integer :: t,c,r,gid

    TRACE_ENTER("alb_diag")
    if(LIS_rc%usealbedomap(n).ne."none") then
       do t=1,LIS_rc%ntiles(n)
          temp1 = LIS_alb(n)%albsf(t)
          if(temp1.ne.-9999.0) then
             temp2 = temp1*100.0
          else
             temp2 = -9999.0
          endif

          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNFRALBEDO, vlevel=1,&
               value=temp1,unit="-",direction="-")
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SNFRALBEDO, vlevel=1,&
               value=temp2,unit="%",direction="-")
       enddo
    endif

    if(LIS_rc%usemxsnalbmap(n).eq."LDT") then
       do t=1,LIS_rc%ntiles(n)
          temp1 = lis_alb(n)%mxsnalb(t)
          if(temp1.ne.-9999.0) then
             temp2 = temp1*100.0
          else
             temp2 = -9999.0
          endif

          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_MXSNALBEDO, vlevel=1,&
               value=temp1,unit="-",direction="-")
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_MXSNALBEDO, vlevel=1,&
               value=temp2,unit="%",direction="-")
       enddo
    endif
  TRACE_EXIT("alb_diag")
  end subroutine LIS_diagnosealbedo

!BOP
!
! !ROUTINE: LIS_albedo_reset
! \label{LIS_albedo_reset}
!
! !INTERFACE:
  subroutine LIS_albedo_reset
! !USES:

!
! !DESCRIPTION:
!
! Resets data structures for reading
! in albedo datasets
!
! The routines invoked are:
! \begin{description}
!  \item[LIS\_read\_mxsnalb](\ref{LIS_read_mxsnalb}) \newline
!    method to read the max snow albedo
!   \item[albedosetup](\ref{albedosetup}) \newline
!    calls the registry to invoke the albedo setup method.
!   \item[readalbedo](\ref{readalbedo}) \newline
!    invokes the generic method in the registry to read the
!    albedo climatology data
! \end{description}
!EOP
    implicit none
    integer         :: n
    integer         :: i
    integer         :: ndoms
    integer         :: t1, t2,t
    real            :: wt1, wt2
    integer         :: c,r
    real, allocatable :: value1(:,:) ! Temporary value holder for QQ1
    real, allocatable :: value2(:,:) ! Temporary value holder for QQ2
    integer         :: rc
!------------------------------------------------------------------------------
! If albedo datasets are to be used, the routine allocates memory for
! each domain. Further, the routine sets the alarms for reading the
! climatology datasets
!------------------------------------------------------------------------------
    TRACE_ENTER("alb_reset")
    ndoms = 0 
    do n=1,LIS_rc%nnest
       if(LIS_rc%usealbedomap(n).ne."none") ndoms = ndoms+1
    enddo

    if(ndoms.gt.0) then
       do n=1,LIS_rc%nnest

          lis_alb(n)%albsf1 = 0.0
          lis_alb(n)%albsf2 = 0.0
          lis_alb(n)%albsf = 0.0
          lis_alb(n)%albedo = 0.0

          if(LIS_rc%usealbedomap(n).ne."LDT") then
             call albedosetup(trim(LIS_rc%usealbedomap(n))//char(0),n)
          endif
          call LIS_computeTemporalWeights(LIS_rc,&
               LIS_alb(n)%albIntervalType, t1,t2,wt1,wt2)

          allocate(value1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
          allocate(value2(LIS_rc%lnc(n),LIS_rc%lnr(n)))

          if(LIS_rc%usealbedomap(n).eq."LDT") then
             call read_albedoclimo(n,t1,value1)
             call read_albedoclimo(n,t2,value2)
          else
             !EMK Corrected subroutine call.
             !call readalbedo(trim(LIS_rc%usealbedomap(n))//char(0),&
             !     n,t1,value1)
             !call readalbedo(trim(LIS_rc%usealbedomap(n))//char(0),&
             !     n,t2,value2)
             call readalbedo(trim(LIS_rc%usealbedomap(n))//char(0), &
                  n, wt1, wt2, value1, value2)
          endif


          do t=1,LIS_rc%ntiles(n)
             lis_alb(n)%albsf1(t) = value1(LIS_domain(n)%tile(t)%col,&
                  LIS_domain(n)%tile(t)%row)
             lis_alb(n)%albsf2(t) = value2(LIS_domain(n)%tile(t)%col,&
                  LIS_domain(n)%tile(t)%row)
          enddo

          deallocate(value1)
          deallocate(value2)
!-------------------------------------------------------------------------
! Interpolate the albedo fractions to daily values
!-------------------------------------------------------------------------
          do i=1,LIS_rc%ntiles(n)
             if (lis_alb(n)%albsf1(i) .ne. -9999.000) then
                lis_alb(n)%albsf(i) = wt1 * lis_alb(n)%albsf1(i) + &
                     wt2 * lis_alb(n)%albsf2(i)
             endif
          end do
       enddo
    endif

    ndoms = 0
    do n=1,LIS_rc%nnest
       if(LIS_rc%usemxsnalbmap(n).ne."none") ndoms = ndoms+1
    enddo

    if(ndoms.gt.0) then
       do n=1,LIS_rc%nnest
          lis_alb(n)%mxsnalb = 0.0
          call LIS_read_mxsnalb(n)
       enddo
    endif
    TRACE_EXIT("alb_reset")

  end subroutine LIS_albedo_reset

!BOP
!
! !ROUTINE: read_albedoclimo
!  \label{read_albedoclimo}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
  subroutine read_albedoclimo(n,time,array)
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

    implicit none
! !ARGUMENTS:
    integer, intent(in) :: n
    integer             :: time
    real, intent(inout) :: array(LIS_rc%lnc(n),LIS_rc%lnr(n))
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
    integer :: status
    integer :: ios1
    integer :: ios,nid,albedoid,ncId, nrId,mid
    integer :: nc,nr,t,months
    integer :: year, month, q
    real, allocatable :: albedo(:,:,:)
    logical :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    q = time

    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then

       if(LIS_alb(n)%albIntervalType.eq."quarterly") then
          write(LIS_logunit,*)'[INFO] Reading albedo map for quarter ', q
       elseif(LIS_alb(n)%albIntervaltype.eq."monthly") then
          write(LIS_logunit,*)'[INFO] Reading albedo map for month ', q
       endif

       ios = nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in read_albedoclimo')

       ios = nf90_inq_dimid(nid,"east_west",ncId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_albedoclimo')

       ios = nf90_inq_dimid(nid,"north_south",nrId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_albedoclimo')

       if(LIS_alb(n)%albIntervalType.eq."quarterly") then
          ios = nf90_inq_dimid(nid,"quarter",mId)
          call LIS_verify(ios,'Error in nf90_inq_dimid in read_albedoclimo')
       elseif(LIS_alb(n)%albIntervalType.eq."monthly") then
          ios = nf90_inq_dimid(nid,"month",mId)
          call LIS_verify(ios,'Error in nf90_inq_dimid in read_albedoclimo')
       endif

       ios = nf90_inquire_dimension(nid,ncId, len=nc)
       call LIS_verify(ios,&
            'Error in nf90_inquire_dimension in read_albedoclimo')

       ios = nf90_inquire_dimension(nid,nrId, len=nr)
       call LIS_verify(ios,&
            'Error in nf90_inquire_dimension in read_albedoclimo')

       ios = nf90_inquire_dimension(nid,mId, len=months)
       call LIS_verify(ios,&
            'Error in nf90_inquire_dimension in read_albedoclimo')

       ios = nf90_inq_varid(nid,'ALBEDO',albedoid)
       call LIS_verify(ios,'ALBEDO field not found in the LIS param file')

       ios = nf90_get_var(nid,albedoid,array,&
            start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
            LIS_nss_halo_ind(n,LIS_localPet+1),q/),&
            count=(/LIS_rc%lnc(n),LIS_rc%lnr(n),1/))
       call LIS_verify(ios,'Error in nf90_get_var in read_albedoclimo')

       ios = nf90_close(nid)
       call LIS_verify(ios,'Error in nf90_close in read_albedoclimo')

    else
       write(LIS_logunit,*) '[ERR] albedo map: ',LIS_rc%paramfile(n), &
            ' does not exist'
       write(LIS_logunit,*) '[ERR] program stopping ...'
       call LIS_endrun
    endif
#endif
  end subroutine read_albedoclimo
end module LIS_albedoMod
