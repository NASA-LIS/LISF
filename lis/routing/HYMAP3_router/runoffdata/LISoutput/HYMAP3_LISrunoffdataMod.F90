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

module HYMAP3_LISrunoffdataMod
!BOP
!
! !MODULE: LISrunoffdataMod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  8 Jan 2016: Sujay Kumar, initial implementation
! 17 Mar 2016: Augusto Getirana, Save in memory input file name and
! surface runoff and baseflow variables - this will reduce the number of
! times input files are read
!
! !USES:
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: HYMAP3_LISrunoffdata_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------

  public :: HYMAP3_LISrunoffdata_struc

  type, public :: HYMAP3_LISrunoffdatadec
     real                 :: outInterval
     character(len=LIS_CONST_PATH_LEN)        :: odir
     character(len=LIS_CONST_PATH_LEN)        :: domfile
     integer, allocatable :: n11(:)
     !ag - 17Mar2016
     logical             :: domaincheck
     integer             :: nc,nr
     character(len=LIS_CONST_PATH_LEN)       :: previous_filename
     real                :: datares
     real, allocatable   :: qs(:,:),qsb(:,:),evap(:,:)
  end type HYMAP3_LISrunoffdatadec

  type(HYMAP3_LISrunoffdatadec), allocatable :: &
       HYMAP3_LISrunoffdata_struc(:)

contains

!BOP
!
! !ROUTINE: HYMAP3_LISrunoffdata_init
! \label{HYMAP3_LISrunoffdata_init}
!
  subroutine HYMAP3_LISrunoffdata_init

    !USES:
    use LIS_coreMod
    use LIS_logMod
    use LIS_timeMgrMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

    implicit none

    integer              :: n
    integer              :: status
    integer              :: ftn
    character*100        :: lis_map_proj
    character*10         :: time
    integer              :: nc,nr
    integer              :: ncId, nrId
    real                 :: lat1,lat2
    real                 :: lon1,lon2
    real                 :: dx,dy
    real                 :: gridDesc(50)

    external :: neighbor_interp_input
    external :: upscaleByAveraging_input

    allocate(HYMAP3_LISrunoffdata_struc(LIS_rc%nnest))

    call ESMF_ConfigFindLabel(LIS_config,&
         "LIS runoff output directory:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config, &
            HYMAP3_LISrunoffdata_struc(n)%odir,rc=status)
       call LIS_verify(status,&
            "LIS runoff output directory: not defined")
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "LIS runoff output interval:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
       call LIS_verify(status,&
            "LIS runoff output interval: not defined")
       call LIS_parseTimeString(time, &
            HYMAP3_LISrunoffdata_struc(n)%outInterval)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "LIS runoff output domain file:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config, &
            HYMAP3_LISrunoffdata_struc(n)%domfile,rc=status)
       call LIS_verify(status,&
            "LIS runoff output domain file: not defined")
    enddo

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    do n=1,LIS_rc%nnest
       call LIS_verify(nf90_open( &
            path=HYMAP3_LISrunoffdata_struc(n)%domfile,&
            mode=NF90_NOWRITE,ncid=ftn),&
            'Error opening file '// &
            trim(HYMAP3_LISrunoffdata_struc(n)%domfile))

       call LIS_verify(nf90_get_att(ftn,NF90_GLOBAL,'MAP_PROJECTION',&
            lis_map_proj),&
            'Error in nf90_get_att: MAP_PROJECTION')

       if(lis_map_proj.eq."EQUIDISTANT CYLINDRICAL") then

          call LIS_verify(nf90_inq_dimid(ftn,'east_west',ncId),&
               'Error in nf90_inq_dimid: east_west')

          call LIS_verify(nf90_inq_dimid(ftn,'north_south',nrId),&
               'Error in nf90_inq_dimid: north_south')

          call LIS_verify(nf90_inquire_dimension(ftn,ncId, &
               len=HYMAP3_LISrunoffdata_struc(n)%nc),&
               'Error in nf90_inquire_dimension: ncId')

          call LIS_verify(nf90_inquire_dimension(ftn,nrId, &
               len=HYMAP3_LISrunoffdata_struc(n)%nr),&
               'Error in nf90_inquire_dimension: nrId')

          HYMAP3_LISrunoffdata_struc(n)%domainCheck = .true.

          if(HYMAP3_LISrunoffdata_struc(n)%nc.ne.LIS_rc%lnc(n).or.&
               HYMAP3_LISrunoffdata_struc(n)%nr.ne.LIS_rc%lnr(n)) then

             HYMAP3_LISrunoffdata_struc(n)%domainCheck = .false.
             call LIS_verify(nf90_get_att(ftn,NF90_GLOBAL, &
                  'SOUTH_WEST_CORNER_LAT',&
                  lat1),&
                  'Error in nf90_get_att: SOUTH_WEST_CORNER_LAT')

             call LIS_verify(nf90_get_att(ftn,NF90_GLOBAL, &
                  'SOUTH_WEST_CORNER_LON',&
                  lon1),&
                  'Error in nf90_get_att: SOUTH_WEST_CORNER_LON')

             call LIS_verify(nf90_get_att(ftn,NF90_GLOBAL,'DX',&
                  dx),&
                  'Error in nf90_get_att: DX')

             call LIS_verify(nf90_get_att(ftn,NF90_GLOBAL,'DY',&
                  dy),&
                  'Error in nf90_get_att: DY')

             gridDesc = 0

             lat2 = (HYMAP3_LISrunoffdata_struc(n)%nr-1)*dx + lat1
             lon2 = (HYMAP3_LISrunoffdata_struc(n)%nc-1)*dy + lon1

             gridDesc(1) = 0
             gridDesc(2) = HYMAP3_LISrunoffdata_struc(n)%nc
             gridDesc(3) = HYMAP3_LISrunoffdata_struc(n)%nr
             gridDesc(4) = lat1
             gridDesc(5) = lon1
             gridDesc(7) = lat2
             gridDesc(8) = lon2
             gridDesc(6) = 128
             gridDesc(9) = dx
             gridDesc(10) = dy
             gridDesc(20) = 64

             HYMAP3_LISrunoffdata_struc(n)%datares = min(dx,dy)

             if(LIS_isAtAfinerResolution(n, &
                  HYMAP3_LISrunoffdata_struc(n)%datares)) then
                allocate(HYMAP3_LISrunoffdata_struc(n)%n11( &
                     LIS_rc%lnc(n)*LIS_rc%lnr(n)))
                call neighbor_interp_input(n,gridDesc,&
                     HYMAP3_LISrunoffdata_struc(n)%n11)
             else
                nc = HYMAP3_LISrunoffdata_struc(n)%nc
                nr = HYMAP3_LISrunoffdata_struc(n)%nr
                allocate(HYMAP3_LISrunoffdata_struc(n)%n11(nc*nr))
                call upscaleByAveraging_input(gridDesc,&
                     LIS_rc%gridDesc(n,:),&
                     nc*nr,&
                     LIS_rc%lnc(n)*LIS_rc%lnr(n),&
                     HYMAP3_LISrunoffdata_struc(n)%n11)
             endif
          endif
       else
          write(LIS_logunit,*) &
               '[ERR] currently only LIS data in lat/lon projection' // &
               ' is supported'
          call LIS_endrun()
       endif
    enddo
#endif

    !ag - 17Mar2016
    do n=1, LIS_rc%nnest
       HYMAP3_LISrunoffdata_struc(n)%previous_filename='none'
       allocate(HYMAP3_LISrunoffdata_struc(n)%qs( &
            LIS_rc%lnc(n),LIS_rc%lnr(n)))
       allocate(HYMAP3_LISrunoffdata_struc(n)%qsb( &
            LIS_rc%lnc(n),LIS_rc%lnr(n)))
       allocate(HYMAP3_LISrunoffdata_struc(n)%evap( &
            LIS_rc%lnc(n),LIS_rc%lnr(n)))
       HYMAP3_LISrunoffdata_struc(n)%qs=LIS_rc%udef
       HYMAP3_LISrunoffdata_struc(n)%qsb=LIS_rc%udef
       HYMAP3_LISrunoffdata_struc(n)%evap=LIS_rc%udef
    enddo

  end subroutine HYMAP3_LISrunoffdata_init
end module HYMAP3_LISrunoffdataMod
