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
! !ROUTINE: geowrsi2_dynsetup
! \label{geowrsi2_dynsetup}
! 
! !INTERFACE:
subroutine geowrsi2_dynsetup(n)

! !USES:
  use LIS_coreMod,    only : LIS_rc, LIS_masterproc
  use LIS_historyMod, only : LIS_readvar_restart
  use LIS_logMod,     only : LIS_logunit, LIS_endrun, LIS_verify
  use geowrsi2_lsmMod
  use geowrsi2_module
  use geowrsi2_physics_module, only : offsetTStepByNumTSteps, stepnFromTStep, &
           gTStepsBeforeSeasonStarts, gTStepsAfterSeasonEnds, gTimeStepsPerYear, &
           DifferenceOf2TSteps

  implicit none
! !ARGUMENTS: 
  integer, intent(in)  :: n
!
! !DESCRIPTION:
!  This routine sets up and ingests time-dependent (dynamic) variables 
!   for GeoWRSI model, including yearly SOS/SOSa files. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP

! Other arguments
  character*20  :: wformat

  integer       :: i, t, rc
  type(geowrsi2dec), pointer :: geowrsi2Pt
  logical, allocatable, dimension(:), target :: mask_default
  logical,     pointer, dimension(:)         :: mask

  integer*4 :: SOS
  integer*4 :: nTsteps
  integer   :: year_temp, ts_temp, offset_temp

! _______________________________________________

   if( LIS_rc%npatch(n,LIS_rc%lsm_index) == 0 ) then
!      write(LIS_logunit,*) " MSG: THIS PROCESSING ELEMENT,",LIS_localPet,&
!         " HAS ZERO-LAND TILES WITHIN DOMAIN FOR NEST, ",n,".  SKIPPING ..."
      return
   endif


!- Read-in SOS-files :
   if( geowrsi2_lsmRunMode == "WRSI" .and. LIS_rc%nensem(n) == 1  ) then

     if( geowrsi2_struc(n)%eos_alarm ) then

        wformat="bil"   ! Hard-coded to bil for now
        write(LIS_logunit,*) "... READING IN SOS/SOSa (.BIL) FILES "
!        write(*,*) "... READING IN SOS/SOSa (.BIL) FILES "

      ! Mask Default (Set to true):
        allocate(mask_default(LIS_rc%npatch(n,LIS_rc%lsm_index)))
        mask_default = .true.
        mask => mask_default

     !- Reinitialize SOS fields:
        geowrsi2_struc(n)%wrsi%SOS  = LIS_rc%udef     ! 0 dekad init
        geowrsi2_struc(n)%wrsi%SOSa = 0

        year_temp = geowrsi2_struc(n)%lastSOScalcOfSeasonYr 

        call geowrsi2_readSOS_Bilfile( n, &
                     LIS_rc%npatch(n,LIS_rc%lsm_index),    &
                     year_temp, geowrsi2_struc(n)%wrsi%SOS,&
                     geowrsi2_struc(n)%wrsi%SOSa,          &
                     geowrsi2_struc(n)%wrsi%SOS_offset )


     !- Update Final Year and Final TStep based on SOS+LGP timesteps:
        nTsteps = 0.
        do t = 1, LIS_rc%npatch(n,LIS_rc%lsm_index)
          if( mask(t) .eqv. .true. ) then
            geowrsi2Pt => geowrsi2_struc(n)%wrsi(t)

            geowrsi2Pt%FinalYear = geowrsi2Pt%InitialYear
            if( geowrsi2Pt%FinalTStep < geowrsi2Pt%InitialTStep ) then
              geowrsi2Pt%FinalYear = geowrsi2Pt%FinalYear + 1
            endif
            nTsteps = stepnFromTStep( &
                      int(geowrsi2Pt%FinalYear,4), geowrsi2Pt%FinalTStep, &
                      int(geowrsi2Pt%InitialYear, 4), geowrsi2Pt%InitialTStep, &
                      (gTStepsBeforeSeasonStarts+gTStepsAfterSeasonEnds) )

            if( (geowrsi2Pt%SOS == int(LIS_rc%udef))       .or. &
                (geowrsi2Pt%LGP_TIMESTEPS == int(LIS_rc%udef))  &
              ) cycle

            geowrsi2Pt%FinalYear = geowrsi2Pt%FinalYearFromDataFile

            geowrsi2Pt%FinalYear  = geowrsi2Pt%InitialYear
            geowrsi2Pt%FinalTStep = int(LIS_rc%udef)     ! Reinitialize

            if( geowrsi2Pt%SOS < geowrsi2Pt%InitialTStep ) then
               geowrsi2Pt%FinalYear = geowrsi2Pt%FinalYear + 1
            endif

            geowrsi2Pt%FinalTStep = geowrsi2Pt%SOS

            call offsetTStepByNumTSteps(       &
                     geowrsi2Pt%FinalTStep,    &
                     geowrsi2Pt%LGP_TIMESTEPS, &
                     geowrsi2Pt%FinalYear)

            if( (geowrsi2Pt%LGP_TIMESTEPS+gTStepsBeforeSeasonStarts) >= nTsteps ) then
              geowrsi2Pt%FinalTStep = geowrsi2Pt%FinalTStep - &
                   ((geowrsi2Pt%LGP_TIMESTEPS+gTStepsBeforeSeasonStarts+1)-nTsteps)
              do while( geowrsi2Pt%FinalTStep <= 0 )
                geowrsi2Pt%FinalTStep = geowrsi2Pt%FinalTStep + gTimeStepsPerYear
                geowrsi2Pt%FinalYear = geowrsi2Pt%FinalYear - 1
              enddo
            endif

          endif
        end do
        deallocate(mask_default)

        geowrsi2_struc(n)%eos_alarm = .false.

      endif
   end if

end subroutine geowrsi2_dynsetup
