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
#include "LDT_NetCDF_inc.h"

module LDT_SimpleWeight_TdscaleMod
!BOP
!
! !ROUTINE: LDT_SimpleWeight_TdscaleMod
!  \label{LDT_SimpleWeight_TdscaleMod}
!
! !REVISION HISTORY:
!  15  Dec 2014: KR Arsenault; Implemented Simple Weight as Module
!
! !INTERFACE:

! USE:

  implicit none 
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: LDT_SimpleWeight_Tdscale

Contains

!
!BOP
! !ROUTINE: LDT_SimpleWeight_Tdscale
! \label{LDT_SimpleWeight_Tdscale}
!
! !INTERFACE: 
  subroutine LDT_SimpleWeight_Tdscale(n)

! USE:
  use ESMF
  use LDT_FORC_AttributesMod
  use LDT_metforcingMod, only : LDT_forc, LDT_FORC_State, LDT_FORC_Base_State
  use LDT_coreMod,       only : LDT_rc, LDT_domain
  use LDT_logMod,        only : LDT_verify
  use LDT_metforcDscaleMod, only : LDT_metDscale, LDT_FORC_Dscale_State


   implicit none
! ARGUMENTS:     
   integer, intent(in) :: n
!
! !DESCRIPTION:
!  This subroutine temporally disaggregates coarser forcing datasets
!  (e.g., daily) using finer time scale forcing datasets (e.g., 6-hrly)
!  by estimating weights using the finer-scale data.
!
!  Important note:  Works best with accumulated variables, like 
!     precipitation fields.
!   
!  04 Dec 2014:  KR Arsenault: Added new downscaling routine
!
!EOP
    integer                    :: i,t,m,fobjcount
    character*100, allocatable :: forcobjs(:)
    type(ESMF_Field)           :: base1Field, base2Field
    type(ESMF_Field)           :: sumDscale1Field
    real,          pointer     :: forcdata_base1(:), forcdata_base2(:)
    real,          pointer     :: forcdata_sum_dscale1(:)
    integer                    :: status, status1, status2
    real                       :: factor1, factor2
! __________________________________________________________________

    call ESMF_StateGet(LDT_FORC_State(n),itemCount=fobjcount,rc=status)
    call LDT_verify(status,'ESMF_StateGet failed for objcount in LDT_SimpleWeight_TdscaleMod')

    allocate(forcobjs(fobjcount))

    call ESMF_StateGet(LDT_FORC_State(n),itemNameList=forcobjs,rc=status)
    call LDT_verify(status,'ESMF_StateGet failed for forcobjs in LDT_SimpleWeight_TdscaleMod')

    do i=1,fobjcount

       call ESMF_StateGet(LDT_FORC_Base_State(n,1),trim(forcobjs(i)),base1Field,&
            rc=status1)
       call ESMF_StateGet(LDT_FORC_Base_State(n,2),trim(forcobjs(i)),base2Field,&
            rc=status1)

       call ESMF_StateGet(LDT_FORC_Dscale_State(n),"Sum_"//trim(forcobjs(i))//"_1",&
            sumDscale1Field,rc=status1)

       if( status1.eq.0 ) then
          call ESMF_FieldGet(base1Field,localDE=0,farrayPtr=forcdata_base1, &
               rc=status2)
          call LDT_verify(status2,'ESMF_FieldGet failed in LDT_SimpleWeight_TdscaleMod')

          call ESMF_FieldGet(base2Field,localDE=0,farrayPtr=forcdata_base2, &
               rc=status2)
          call LDT_verify(status2,'ESMF_FieldGet failed in LDT_SimpleWeight_TdscaleMod')

          call ESMF_FieldGet(sumdscale1Field,localDE=0,farrayPtr=forcdata_sum_dscale1, &
               rc=status2)
          call LDT_verify(status2,'ESMF_FieldGet failed in LDT_SimpleWeight_TdscaleMod')

          factor1 = 1.0
          factor2 = 1.0
          if( index(forcobjs(i),"Rainfall") > 0. ) then
            factor1 = LDT_metDscale%fine_ts
            factor2 = LDT_metDscale%coarse_ts
          endif

        ! Search for present forcing fields:
          do t=1,LDT_rc%ntiles(n)
             if( forcdata_base1(t) .ne. LDT_rc%udef .and. &
                 forcdata_base2(t) .ne. LDT_rc%udef .and. &
                 forcdata_sum_dscale1(t).gt.0. ) then

             ! Calculate temporally downscaled coarser forcing dataset:
               forcdata_base2(t) = forcdata_base2(t)*factor2 &
                       *((forcdata_base1(t)*factor1)/forcdata_sum_dscale1(t))

             ! Convert coarser forcing dataset to a rate:
               forcdata_base2(t) = forcdata_base2(t) / factor1

             endif
          enddo

       endif
    enddo     ! End forcing type loop

    deallocate(forcobjs)

  end subroutine LDT_SimpleWeight_Tdscale

end module LDT_SimpleWeight_TdscaleMod
