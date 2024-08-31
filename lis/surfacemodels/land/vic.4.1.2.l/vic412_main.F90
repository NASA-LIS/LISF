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
!BOP
!
! !ROUTINE: vic412_main
! \label{vic412_main}
!
! !REVISION HISTORY:
! 02 Aug 2011; James Geiger, Initial implementation of VIC 4.1.1 into LIS.
! 03 Apr 2012; Shugong Wang, add support to restart file with MPI-safe tile id 
! 03 Sep 2012; Shugong Wang, add support to output VIC PET variables
! !INTERFACE:
subroutine vic412_main(n)
! !USES:
  use LIS_coreMod,        only : LIS_rc, LIS_domain
  use LIS_timeMgrMod,     only : LIS_isAlarmRinging 
  use LIS_histDataMod
  use vic412_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)  :: n
!
! !DESCRIPTION:
!  This is the entry point for calling the VIC 4.1.1 LSM physics.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP

   integer :: t, ts, nest
   integer :: iret
   logical :: alarmCheck
   integer :: vt_scheme ! added by Shugong Wang to support reading restart file
   integer :: vegclass  ! added by Shugong Wang to calculate MPI-safe tile id 
   character*3 :: fnest

   write(fnest,'(i3.3)') n
   alarmCheck = LIS_isAlarmRinging(LIS_rc,"VIC412 model alarm "//trim(fnest))
   vt_scheme = vic412_struc(n)%veg_tiling_scheme ! added by Shugong Wang 04/30/2012
   if ( alarmCheck ) then 
      ts = vic412_struc(n)%tscount
      nest = n
      do t = 1, LIS_rc%npatch(n,LIS_rc%lsm_index)
         ! retrieve landuse type  Shugong Wang 05/07/2012
         if (vt_scheme == 1) then ! LIS-based tiling
             vegclass = LIS_domain(n)%tile(t)%vegt
         else                ! VIC-based tiling
             vegclass = -1
         endif
         call vic412_run(t, ts, vt_scheme,  & ! vt_scheme is added by Shugong on 04/30/12
                         vegclass,          & ! vegclass is added by Shugong on 05/07/12
                         nest,              &
                         LIS_MOC_EVAP,      &
                         LIS_MOC_QS,        &
                         LIS_MOC_QSB,       &
                         LIS_MOC_CANOPINT,  &
                         LIS_MOC_SMLIQFRAC, &
                         LIS_MOC_SWNET,     &
                         LIS_MOC_ECANOP,    &
                         LIS_MOC_TVEG,      &
                         LIS_MOC_ESOIL,     &
                         LIS_MOC_ARESIST,   &
                         LIS_MOC_AVGSURFT,  &
                         LIS_MOC_ALBEDO,    &
                         LIS_MOC_LWNET,     &
                         LIS_MOC_SUBSNOW,   &
                         LIS_MOC_RADT,      &
                         LIS_MOC_QLE,       &
                         LIS_MOC_QH,        &
                         LIS_MOC_QG,        &
                         LIS_MOC_QF,        &
                         LIS_MOC_ACOND,         &
                         LIS_MOC_BARESOILT,     &
                         LIS_MOC_DELCOLDCONT,   &
                         LIS_MOC_DELINTERCEPT,  &
                         LIS_MOC_DELSOILMOIST,  &
                         LIS_MOC_DELSURFSTOR,   &
                         LIS_MOC_QFZ,           &
                         LIS_MOC_QSM,           &
                         LIS_MOC_QV,            &
                         LIS_MOC_RAINF,         &
                         LIS_MOC_ROOTMOIST,     &
                         LIS_MOC_SMFROZFRAC,    &
                         LIS_MOC_SNFRALBEDO,    &
                         LIS_MOC_SNOWCOVER,     &
                         LIS_MOC_SNOWDEPTH,     &
                         LIS_MOC_SNOWF,         &
                         LIS_MOC_SNOWT,         &
                         LIS_MOC_SOILMOIST,     &
                         LIS_MOC_SWDOWNFORC,    &
                         LIS_MOC_SWE,           &
                         LIS_MOC_TAIRFORC,      &
                         LIS_MOC_TOTALPRECIP,   &
                         LIS_MOC_VIC_PET_SATSOIL, &
                         LIS_MOC_VIC_PET_H2OSURF, &
                         LIS_MOC_VIC_PET_SHORT,   &
                         LIS_MOC_VIC_PET_TALL,    &
                         LIS_MOC_VIC_PET_NATVEG,  &
                         LIS_MOC_VIC_PET_VEGNOCR, &
                         LIS_MOC_SOILTEMP,        &
                         LIS_MOC_WRSI,            &
                         LIS_MOC_SWI,             &
                         LIS_MOC_WATERTABLED)

      enddo
      vic412_struc(n)%tscount = vic412_struc(n)%tscount + 1
   endif
end subroutine vic412_main

