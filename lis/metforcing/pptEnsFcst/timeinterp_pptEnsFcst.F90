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
! !ROUTINE: timeinterp_pptEnsFcst
! \label{timeinterp_pptEnsFcst}
!
! !REVISION HISTORY: 
! 27Sep2016 : K. Arsenault; implemented into LIS-7
!
! !INTERFACE:
subroutine timeinterp_pptEnsFcst(n, findex)

! !USES:
  use ESMF
  use LIS_coreMod,         only : LIS_rc,LIS_domain
  use LIS_constantsMod,    only : LIS_CONST_SOLAR
  use LIS_logMod,          only : LIS_logunit, LIS_verify, LIS_endrun
  use LIS_timeMgrMod,      only : LIS_time2date
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod,   only : LIS_FORC_Base_State
  use pptEnsFcst_forcingMod, only : pptensfcst_struc
  use pptEnsFcst_VariablesMod, only : forcopts

! !ARGUMENTS: 
  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: findex

! !DESCRIPTION: 
!
!  Temporally interpolates the forcing data to the current model 
!  timestep. If needed, downward shortwave radiation is interpolated using a
!  zenith-angled based approach, depending on input forcing. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[findex]
!   index of the forcing
!  \end{description}
! 
!EOP

! Time-interp weights:
  real    :: wt1, wt2
! Zterp terms:
  real    :: gmt1, gmt2
  integer :: bdoy,byr,bmo,bda,bhr,bmn
  real*8  :: btime
  integer :: zdoy
  real    :: zw1,zw2,czb,cze,czm

  integer :: t, m, k
  integer :: index1
  integer :: mfactor, tid
  integer :: status

! ESMF Fields and Pointer Arrays:
  type(ESMF_Field)   :: prcpField, cpcpField
  real,pointer       :: prcp(:), cpcp(:)
! __________________________________________________________

! Get Meteorological Field  - ESMF State Get:
  if( forcopts%read_rainf ) then
    call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Rainf%varname(1),prcpField,&
         rc=status)
    call LIS_verify(status, 'Error: Enable Rainf in the forcing variables list')
    call ESMF_FieldGet(prcpField,localDE=0,farrayPtr=prcp,rc=status)
    call LIS_verify(status)
    prcp  = 0.0
  endif

  if( forcopts%read_cpcp ) then
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_CRainf%varname(1),cpcpField,&
          rc=status)
     call LIS_verify(status, 'Error: Enable CRainf in the forcing variables list')
     call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
     call LIS_verify(status)
     cpcp = 0.0
  endif

#if 0
!- Time interpolate or apply average rate:
   wt1 = (pptensfcst_struc%metforc_time2 - LIS_rc%time) / &
         (pptensfcst_struc%metforc_time2 - pptensfcst_struc%metforc_time1)
   wt2 = 1.0 - wt1
#endif 

  ! Ensemble member factor for this ensemble forcing dataset:
  mfactor = LIS_rc%nensem(n)/pptensfcst_struc%max_ens_members

!- Time Averaged Longwave, Block Interpolation

  do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
     do m=1,pptensfcst_struc%max_ens_members
        do k=1,mfactor
           tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
           index1 = LIS_domain(n)%tile(tid)%index

          if( forcopts%read_rainf ) then 
            if((pptensfcst_struc%metdata1(forcopts%index_rainf,m,index1).ne.LIS_rc%udef).and.&
               (pptensfcst_struc%metdata2(forcopts%index_rainf,m,index1).ne.LIS_rc%udef)) then
              prcp(tid) = pptensfcst_struc%metdata2(forcopts%index_rainf,m,index1)

              if( prcp(tid) < 0.0 ) then ! mm/s
                prcp(tid) = 0.0
              endif
              if( prcp(tid) > 1.0 ) then ! mm/s
                prcp(tid) = 0.025
              endif

            endif
          endif
 
          if( forcopts%read_cpcp ) then 
            if((pptensfcst_struc%metdata1(forcopts%index_cpcp,m,index1).ne.LIS_rc%udef).and.&
               (pptensfcst_struc%metdata2(forcopts%index_cpcp,m,index1).ne.LIS_rc%udef)) then
              cpcp(tid) = pptensfcst_struc%metdata2(forcopts%index_cpcp,m,index1)

              if( cpcp(tid) < 0.0 ) then
                cpcp(tid) = 0.0
              endif

            endif     
          endif

       end do
     enddo
   enddo   ! End tile loop

end subroutine timeinterp_pptEnsFcst
