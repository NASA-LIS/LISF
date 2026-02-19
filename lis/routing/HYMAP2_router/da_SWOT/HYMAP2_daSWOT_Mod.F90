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
module HYMAP2_daSWOT_Mod
!BOP
!
! !MODULE: HYMAP2_daSWOT_Mod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
! 15 Apr 24: Yeosang Yoon; Initial specification;
!                          copied from HYMAP2_daWL_Mod
!
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_dataAssimMod
  use LIS_logMod

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: HYMAP2_daSWOT_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: HYMAP2_daSWOT_struc
!EOP

 type, public :: daSWOT_dec

     integer                :: nbins
     integer                :: ntimes
     integer                :: scal
     integer                :: nsites
     integer                :: useLocalUpd
     integer                :: localUpdDX
     real, allocatable      :: sites(:,:)
     real, allocatable      :: localWeight(:,:,:)

  end type daSWOT_dec

  type(daSWOT_dec), allocatable :: HYMAP2_daSWOT_struc(:)

contains
!BOP
!
! !ROUTINE: HYMAP2_daSWOT_init
! \label{HYMAP2_daSWOT_init}
!
! !INTERFACE:
  subroutine HYMAP2_daSWOT_init(k)
! !USES:
    use HYMAP2_routingMod
    use HYMAP2_initMod
    use LIS_constantsMod, only: LIS_CONST_PATH_LEN
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
!
! !DESCRIPTION:
!
!  This routine initializes the DA related data structures in HYMAP2
!  and the runtime DA configuration options
!
!EOP

    implicit none

    integer                :: k
    integer                :: n
    integer                :: nid
    integer                :: status
    integer                :: ncId, nrId, siteId
    integer                :: nc,nr,drainid,sid
    logical                :: local_upd_flag
    character(len=LIS_CONST_PATH_LEN) :: localWeightMap
    logical                :: file_exists

    if(.not.allocated(HYMAP2_daSWOT_struc)) then
       allocate(HYMAP2_daSWOT_struc(LIS_rc%nnest))

       call ESMF_ConfigFindLabel(LIS_config,&
            "HYMAP2 use localization update in DA:",rc=status)
       do n=1, LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config,&
               HYMAP2_daSWOT_struc(n)%useLocalUpd,rc=status)
          call LIS_verify(status,&
               "HYMAP2 use localization update in DA: not defined")
       enddo

       local_upd_flag = .false.
       do n=1,LIS_rc%nnest
          if(HYMAP2_daSWOT_struc(n)%useLocalUpd.eq.1) then
             local_upd_flag = .true.
          endif
       enddo

       if(local_upd_flag) then
          call ESMF_ConfigFindLabel(LIS_config,&
               "HYMAP2 localization update window size:",rc=status)
          do n=1, LIS_rc%nnest
             call ESMF_ConfigGetAttribute(LIS_config,&
                  HYMAP2_daSWOT_struc(n)%localupdDX,rc=status)
             call LIS_verify(status,&
                  "HYMAP2 localization update window size: not defined")
          enddo

          call ESMF_ConfigFindLabel(LIS_config,&
               "HYMAP2 localization weight map:",rc=status)
          do n=1, LIS_rc%nnest
             call ESMF_ConfigGetAttribute(LIS_config,&
                  localWeightMap,rc=status)
             call LIS_verify(status,&
                  "HYMAP2 localization weight map: not defined")
#if (defined USE_NETCDF3 || defined USE_NETCDF4)

             inquire(file=localWeightMap, exist=file_exists)
             if(file_exists) then

                write(LIS_logunit,*)'[INFO] Reading localization map from '&
                     //trim(localWeightMap)
                status = nf90_open(path=localWeightMap,&
                     mode=NF90_NOWRITE,ncid=nid)
                call LIS_verify(status, &
                     'Error in nf90_open in HYMAP2_daSWOT_Mod')

                status = nf90_inq_dimid(nid,"east_west",ncId)
                call LIS_verify(status, &
                     'Error in nf90_inq_dimid in HYMAP2_daSWOT_Mod')

                status = nf90_inq_dimid(nid,"north_south",nrId)
                call LIS_verify(status, &
                     'Error in nf90_inq_dimid in HYMAP2_daSWOT_Mod')

                status = nf90_inq_dimid(nid,"count",siteId)
                call LIS_verify(status, &
                     'Error in nf90_inq_dimid in HYMAP2_daSWOT_Mod')

                status = nf90_inquire_dimension(nid,ncId, len=nc)
                call LIS_verify(status, &
                     'Error in nf90_inquire_dimension in HYMAP2_daSWOT_Mod')

                status = nf90_inquire_dimension(nid,nrId, len=nr)
                call LIS_verify(status, &
                     'Error in nf90_inquire_dimension in HYMAP2_daSWOT_Mod')

                status = nf90_inquire_dimension(nid,siteId, &
                     len=HYMAP2_daSWOT_struc(n)%nsites)
                call LIS_verify(status, &
                     'Error in nf90_inquire_dimension in HYMAP2_daSWOT_Mod')

                status = nf90_inq_varid(nid,'weight',drainid)
                call LIS_verify(status, &
                     'distance field not found in the localization weight file')

                status = nf90_inq_varid(nid,'sites',sid)
                call LIS_verify(status, &
                     'sites field not found in the localization weight file')
                allocate(HYMAP2_daSWOT_struc(n)%localWeight(&
                     LIS_rc%gnc(n), LIS_rc%gnr(n), &
                     HYMAP2_daSWOT_struc(n)%nsites))

                allocate(HYMAP2_daSWOT_struc(n)%sites( &
                     LIS_rc%lnc(n), LIS_rc%lnr(n)))

                status = nf90_get_var(nid,drainid,&
                     HYMAP2_daSWOT_struc(n)%localWeight)
                call LIS_verify(status, &
                     'Error in nf90_get_var in HYMAP2_daSWOT_Mod')

                status = nf90_get_var(nid,sid,&
                     HYMAP2_daSWOT_struc(n)%sites, &
                     start=(/LIS_ews_halo_ind(n,LIS_localPet+1),&
                     LIS_nss_halo_ind(n,LIS_localPet+1)/),&
                     count = (/LIS_ewe_halo_ind(n,LIS_localPet+1) - &
                     LIS_ews_halo_ind(n,LIS_localPet+1)+1, &
                     LIS_nse_halo_ind(n,LIS_localPet+1) - &
                     LIS_nss_halo_ind(n,LIS_localPet+1)+1/))

                call LIS_verify(status, &
                     'Error in nf90_get_var in HYMAP2_daSWOT_Mod')

                status = nf90_close(nid)
                call LIS_verify(status, &
                     'Error in nf90_close in HYMAP2_daSWOT_Mod')

             else
                write(LIS_logunit,*) '[ERR] localization map: ', &
                     trim(localWeightMap), ' does not exist'
                write(LIS_logunit,*) '[ERR] program stopping ...'
                call LIS_endrun
             endif
#endif
          enddo
       endif
    endif
  end subroutine HYMAP2_daSWOT_init
end module HYMAP2_daSWOT_Mod
