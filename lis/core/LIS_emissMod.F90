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
module LIS_emissMod
!BOP
!
! !MODULE: LIS_emissMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read emissivity
!  data.
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the
!  emissivity climatology data and allows the users to
!  specify the frequency of climatology (in months).
!  The climatological data is temporally interpolated
!  between months to the current simulation date.
!
! !REVISION HISTORY:
!
!  08 Aug 2005: Sujay Kumar; Initial implementation
!
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_histDataMod
  use LIS_fileIOMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE

!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LIS_emiss_setup    !allocates memory for required structures
  public :: LIS_read_emiss     !reads the emiss data
  public :: LIS_diagnoseemiss      !diagnoses emiss data for history output
  public :: LIS_emiss_finalize !cleanup allocated structures
  public :: LIS_emiss_reset    !resets alarms and data structures

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LIS_emiss !data structure containing emiss fraction data.

!EOP
  type, public :: emiss_type_dec
     character(len=LIS_CONST_PATH_LEN) :: emissfile
     character*50  :: emissIntervalType
     real          :: emissInterval
     logical       :: firstInstance
     real, allocatable :: emiss(:)
     real, allocatable :: emissp1(:)
     real, allocatable :: emissp2(:)
!     real*8        :: alarmTime
     real, allocatable    :: rlat1(:)
     real, allocatable    :: rlon1(:)
     integer, allocatable :: n111(:), n121(:)
     integer, allocatable :: n211(:), n221(:)
     real, allocatable    :: w111(:), w121(:)
     real, allocatable    :: w211(:), w221(:)
  end type emiss_type_dec

  type(emiss_type_dec), allocatable :: LIS_emiss(:)

contains

!BOP
!
! !ROUTINE: LIS_emiss_setup
! \label{LIS_emiss_setup}
!
! !INTERFACE:
  subroutine LIS_emiss_setup
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
! !DESCRIPTION:
!
! Allocates memory for data structures for reading
! the emiss fraction datasets
!
!  The routines invoked are:
!  \begin{description}
!   \item[emissivitysetup](\ref{emissivitysetup}) \newline
!    calls the registry to invoke the emiss setup methods.
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
    logical :: file_exists

    TRACE_ENTER("emiss_setup")
    ndoms = 0
    do n=1,LIS_rc%nnest
       if(LIS_rc%useemissmap(n).ne."none") then
          ndoms = ndoms+1
       endif
    enddo

    if(ndoms.gt.0) then
       allocate(LIS_emiss(LIS_rc%nnest))

       do n=1,LIS_rc%nnest
          LIS_emiss(n)%firstInstance = .true.
       enddo

       do n=1,LIS_rc%nnest
          allocate(LIS_emiss(n)%emiss(LIS_rc%ntiles(n)))
          allocate(LIS_emiss(n)%emissp1(LIS_rc%ntiles(n)))
          allocate(LIS_emiss(n)%emissp2(LIS_rc%ntiles(n)))
          LIS_emiss(n)%emissp1 = 0.0
          LIS_emiss(n)%emissp2 = 0.0
          LIS_emiss(n)%emiss = 0.0

          if(LIS_rc%useemissmap(n).eq."LDT") then

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

             inquire(file=LIS_rc%paramfile(n), exist=file_exists)
             if(file_exists) then

                ios = nf90_open(path=trim(LIS_rc%paramfile(n)),&
                     mode=NF90_NOWRITE,ncid=nid)
                call LIS_verify(ios,'Error in nf90_open in read_emissclimo')

                ios = nf90_get_att(nid, NF90_GLOBAL, &
                     'EMISSIVITY_DATA_INTERVAL', &
                     LIS_emiss(n)%emissIntervalType)
                call LIS_verify(ios,'Error in nf90_get_att in read_emissclimo')
                ios = nf90_close(nid)
                call LIS_verify(ios,'Error in nf90_close in read_emissclimo')
             else
                write(LIS_logunit,*) LIS_rc%paramfile(n), ' does not exist'
                write(LIS_logunit,*) 'program stopping ...'
                call LIS_endrun
             endif
#endif
             if(LIS_emiss(n)%emissIntervalType.eq."monthly") then
                LIS_emiss(n)%emissInterval = 2592000
             endif

!The intervaltype and interval is set in the plugin, now
!register the alarm.
             call LIS_registerAlarm("LIS emissivity read alarm",LIS_rc%ts, &
                  LIS_emiss(n)%emissInterval,&
                  intervalType=LIS_emiss(n)%emissIntervalType)

             call LIS_computeTemporalWeights(LIS_rc,&
                  LIS_emiss(n)%emissIntervalType, &
                  t1,t2,wt1,wt2)

             allocate(value1(LIS_rc%ntiles(n)))
             allocate(value2(LIS_rc%ntiles(n)))

             value1 = LIS_rc%udef
             value2 = LIS_rc%udef

             call read_emissclimo(n,t1,value1)
             call read_emissclimo(n,t2,value2)

             do i=1,LIS_rc%ntiles(n)
                LIS_emiss(n)%emissp1(i) = value1(i)
                LIS_emiss(n)%emissp2(i) = value2(i)
             enddo

          else
             call emissivitysetup(trim(LIS_rc%useemissmap(n))//char(0),n)
             call reademissivity(trim(LIS_rc%useemissmap(n))//char(0),&
               n,wt1,wt2,LIS_emiss(n)%emissp1,LIS_emiss(n)%emissp2)
          endif

          do i=1,LIS_rc%ntiles(n)
             LIS_emiss(n)%emiss(i) = (wt1*LIS_emiss(n)%emissp1(i))+&
                  (wt2*LIS_emiss(n)%emissp2(i))
          enddo
       enddo
    endif
    TRACE_EXIT("emiss_setup")

  end subroutine LIS_emiss_setup

!BOP
!
! !ROUTINE: LIS_read_emiss
! \label{LIS_read_emiss}
!
! !INTERFACE:
  subroutine LIS_read_emiss(n)
! !USES:

    implicit none
! !ARGUMENTS:
    integer, intent(in) :: n
!
! !DESCRIPTION:
!
!  Reads the emiss fraction climalotogy and temporally interpolates
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
!   \item[reademissivity](\ref{reademissivity}) \newline
!    invokes the generic method in the registry to read the
!    emiss climatology data
!  \end{description}
!
!EOP

    logical :: emissAlarmCheck
    integer :: i
    integer         :: t1, t2
    real            :: wt1, wt2
    real, allocatable   :: value1(:) ! temporary value holder for t1
    real, allocatable   :: value2(:) ! temporary value holder for t2

    TRACE_ENTER("emiss_read")
    if(LIS_rc%useemissmap(n).ne."none") then
       if(LIS_rc%useemissmap(n).eq."LDT") then
          emissAlarmCheck = LIS_isAlarmRinging(LIS_rc,&
               "LIS emissivity read alarm",&
               LIS_emiss(n)%emissIntervalType)

          call LIS_computeTemporalWeights(LIS_rc,&
               LIS_emiss(n)%emissIntervalType, &
               t1,t2,wt1,wt2)

          if(emissAlarmCheck) then
             allocate(value1(LIS_rc%ntiles(n)))
             allocate(value2(LIS_rc%ntiles(n)))

             call read_emissclimo(n,t1,value1)
             call read_emissclimo(n,t2,value2)

             do i=1,LIS_rc%ntiles(n)
                LIS_emiss(n)%emissp1(i) = value1(i)
                LIS_emiss(n)%emissp2(i) = value2(i)
             enddo
          endif
       else
          call reademissivity(trim(LIS_rc%useemissmap(n))//char(0),&
               n, wt1, wt2, LIS_emiss(n)%emissp1, LIS_emiss(n)%emissp2)
       endif

       do i=1,LIS_rc%ntiles(n)
          LIS_emiss(n)%emiss(i) = (wt1*LIS_emiss(n)%emissp1(i))+&
               (wt2*LIS_emiss(n)%emissp2(i))
       end do
    endif
    TRACE_EXIT("emiss_read")
  end subroutine LIS_read_emiss

!BOP
!
! !ROUTINE: LIS_emiss_finalize
! \label{LIS_emiss_finalize}
!
! !INTERFACE:
  subroutine LIS_emiss_finalize
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
       if(LIS_rc%useemissmap(n).ne."none") ndoms = ndoms+1
    enddo

    if(ndoms.gt.0) then
       do n=1,LIS_rc%nnest
          deallocate(LIS_emiss(n)%emiss)
          deallocate(LIS_emiss(n)%emissp1)
          deallocate(LIS_emiss(n)%emissp2)
       enddo
       deallocate(LIS_emiss)
    endif
  end subroutine LIS_emiss_finalize

!BOP
!
! !ROUTINE: LIS_diagnoseemiss
! \label{LIS_diagnoseemiss}
!
! !INTERFACE:
  subroutine LIS_diagnoseemiss(n)
! !USES:

! !ARGUMENTS:
    implicit none
    integer, intent(in)   :: n

! !DESCRIPTION:
!  This routine maps the emiss data to the LIS history writer.
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
    real, allocatable    :: temp(:)

    TRACE_ENTER("emiss_diag")
    allocate(temp(LIS_rc%ntiles(n)))
    if(LIS_rc%useemissmap(n).ne."none") then
       temp = LIS_rc%udef
       do t=1,LIS_rc%ntiles(n)
          if(LIS_domain(n)%tile(t)%index.ne.-1) then
             temp(t) = LIS_emiss(n)%emiss(t)
          endif
          call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_EMISSFORC,vlevel=1,&
               value=temp(t),unit="-", &
               direction="-")
       enddo
    endif
    deallocate(temp)
    TRACE_EXIT("emiss_diag")

  end subroutine LIS_diagnoseemiss

!BOP
!
! !ROUTINE: LIS_emiss_reset
! \label{LIS_emiss_reset}
!
! !INTERFACE:
  subroutine LIS_emiss_reset
! !USES:

! !DESCRIPTION:
!
! Resets the data structures for reading
! the emiss fraction datasets
!
!  The routines invoked are:
!  \begin{description}
!   \item[emissivitysetup](\ref{emissivitysetup}) \newline
!    calls the registry to invoke the emiss data reading methods.
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

    TRACE_ENTER("emiss_reset")
    ndoms = 0
    do n=1,LIS_rc%nnest
       if(LIS_rc%useemissmap(n).ne."none") then
          ndoms = ndoms+1
       endif
    enddo

    if(ndoms.gt.0) then

       do n=1,LIS_rc%nnest

          LIS_emiss(n)%emissp1 = 0.0
          LIS_emiss(n)%emissp2 = 0.0
          LIS_emiss(n)%emiss = 0.0

          if(LIS_rc%useemissmap(n).ne."LDT") then
             call emissivitysetup(trim(LIS_rc%useemissmap(n))//char(0),n)
          else
!             LIS_emiss(n)%emissIntervalType = "monthly"
          endif

          !Read the data for the first time
          call LIS_computeTemporalWeights(LIS_rc,&
               LIS_emiss(n)%emissIntervalType, &
               t1,t2,wt1,wt2)

          allocate(value1(LIS_rc%ntiles(n)))
          allocate(value2(LIS_rc%ntiles(n)))

          if(LIS_rc%useemissmap(n).eq."LDT") then
             call read_emissclimo(n,t1,value1)
             call read_emissclimo(n,t2,value2)
          else
             ! EMK Fixed calls
             !call reademissivity(trim(LIS_rc%useemissmap(n))//char(0),&
             !     n,t1,value1)
             !call reademissivity(trim(LIS_rc%useemissmap(n))//char(0),&
             !     n,t2,value2)
             call reademissivity(trim(LIS_rc%useemissmap(n))//char(0), &
                  n, wt1, wt2, value1, value2)
          endif

          do i=1,LIS_rc%ntiles(n)
             LIS_emiss(n)%emissp1(i) = value1(i)
             LIS_emiss(n)%emissp2(i) = value2(i)
          enddo
          deallocate(value1)
          deallocate(value2)

          do i=1,LIS_rc%ntiles(n)
             LIS_emiss(n)%emiss(i) = (wt1*LIS_emiss(n)%emissp1(i))+&
                  (wt2*LIS_emiss(n)%emissp2(i))
          enddo
       end do
    endif
    TRACE_EXIT("emiss_reset")
  end subroutine LIS_emiss_reset


!BOP
!
! !ROUTINE: read_emissclimo
!  \label{read_emissclimo}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
  subroutine read_emissclimo(n,time,array)
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

    implicit none
! !ARGUMENTS:
    integer, intent(in)         :: n
    integer, intent(in)         :: time
    real, intent(inout)         :: array(LIS_rc%ntiles(n))
! !DESCRIPTION:
!  This subroutine reads the emiss data climatology
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \end{description}
!
!EOP

    integer :: ios1
    integer :: ios,nid,emissid,ncId, nrId,mid
    integer :: nc,nr,t,months
    integer :: mo
    real, allocatable :: emiss(:,:,:)
    real, allocatable :: localemiss(:,:)
    logical :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    mo = time
    allocate(localemiss(LIS_rc%lnc(n),LIS_rc%lnr(n)))

    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then

       write(LIS_logunit,*)'Reading emissivity map for month ', mo
       ios = nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in read_emissclimo')

       ios = nf90_inq_dimid(nid,"east_west",ncId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_emissclimo')

       ios = nf90_inq_dimid(nid,"north_south",nrId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_emissclimo')

       ios = nf90_inq_dimid(nid,"month",mId)
       call LIS_verify(ios,'Error in nf90_inq_dimid in read_emissclimo')

       ios = nf90_inquire_dimension(nid,ncId, len=nc)
       call LIS_verify(ios, &
            'Error in nf90_inquire_dimension in read_emissclimo')

       ios = nf90_inquire_dimension(nid,nrId, len=nr)
       call LIS_verify(ios, &
            'Error in nf90_inquire_dimension in read_emissclimo')

       ios = nf90_inquire_dimension(nid,mId, len=months)
       call LIS_verify(ios, &
            'Error in nf90_inquire_dimension in read_emissclimo')

       allocate(emiss(LIS_rc%gnc(n),LIS_rc%gnr(n),months))

       ios = nf90_inq_varid(nid,'EMISS',emissid)
       call LIS_verify(ios,'EMISS field not found in the LIS param file')

       ios = nf90_get_var(nid,emissid,emiss)
       call LIS_verify(ios,'Error in nf90_get_var in read_emissclimo')

       ios = nf90_close(nid)
       call LIS_verify(ios,'Error in nf90_close in read_emissclimo')

       localemiss(:,:) = &
            emiss(LIS_ews_halo_ind(n,LIS_localPet+1):&
            LIS_ewe_halo_ind(n,LIS_localPet+1), &
            LIS_nss_halo_ind(n,LIS_localPet+1): &
            LIS_nse_halo_ind(n,LIS_localPet+1),mo)

       deallocate(emiss)

       do t=1,LIS_rc%ntiles(n)
          array(t) = localemiss(LIS_domain(n)%tile(t)%col,&
               LIS_domain(n)%tile(t)%row)
       enddo
    else
       write(LIS_logunit,*) 'emiss map: ',LIS_rc%paramfile(n), &
            ' does not exist'
       write(LIS_logunit,*) 'program stopping ...'
       call LIS_endrun
    endif
    deallocate(localemiss)
#endif
  end subroutine read_emissclimo

end module LIS_emissMod
