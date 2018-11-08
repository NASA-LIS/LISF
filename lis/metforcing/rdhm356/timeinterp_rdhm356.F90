!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: timeinterp_rdhm356
! \label{timeinterp_rdhm356}
! 
! !REVISION HISTORY:
!  25May2006: Kristi Arsenault; Initial implementation
!  10 Oct 2006:   Sujay Kumar: Switched to using ESMF_State for storing
!                              forcing data. 
!  03 May 2010: Soni Yatheendradas; Precip and Temper. input grids now 
!               can have different extents/directories and different from 
!               the run-domain extent, as per the new input grids posted
!               onto the DMIP2 website for Sierra Nevada
!  19 Dec 2013; Shugong Wang; RDHM356
!
! !INTERFACE:
  subroutine timeinterp_rdhm356(n, findex)
! !USES:
    use ESMF
    use LIS_coreMod, only: LIS_rc, LIS_domain
    use LIS_FORC_AttributesMod
    use LIS_metforcingMod, only : LIS_FORC_Base_State, LIS_forc
    use LIS_logMod,         only : LIS_verify, LIS_endrun 
    use rdhm356_forcingMod, only : rdhm356_struc_precip, &
                                   rdhm356_struc_temper, &
                                   const_wind 

    implicit none

! !ARGUMENTS: 
    integer, intent(in) :: n
    integer, intent(in) :: findex

!
! !DESCRIPTION:
!  Temporally interpolates the forcing data to the current model 
!  timestep. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!EOP
    integer :: index1
    integer :: c,wt1,wt2

    integer          :: status
    type(ESMF_Field) :: pcpField, tmpField, uField, vField
    real, pointer    :: pcp(:), tmp(:), uwind(:), vwind(:)
    character*32     :: fname
    real,  allocatable :: ratio(:)

    allocate(ratio(LIS_rc%ntiles(n)))

    wt1=(rdhm356_struc_precip(n)%rdhm356time2 - LIS_rc%time)/                            &
         (rdhm356_struc_precip(n)%rdhm356time2 - rdhm356_struc_precip(n)%rdhm356time1)
    wt2=1.0-wt1

    call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Tair%varname(1),tmpField,&
         rc=status)
    call LIS_verify(status, 'Error: Enable Tair in the forcing variables list')

    call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Rainf%varname(1),pcpField,&
         rc=status)
    call LIS_verify(status, 'Error: Enable Rainf in the forcing variables list')
  
    call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_E%varname(1),uField,&
         rc=status)
    call LIS_verify(status, 'Error: Enable Wind_E in the forcing variables list')
  
    call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_N%varname(1),vField,&
         rc=status)
    call LIS_verify(status, 'Error: Enable Wind_N in the forcing variables list')

    call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
    call LIS_verify(status)
    
    call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
    call LIS_verify(status)
    
    call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
    call LIS_verify(status)
      
    call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
    call LIS_verify(status)
 
!    write(fname, *), LIS_rc%ts
!    open(unit=1001, file='time_interp', status='unknown')
!    do c = 1,LIS_rc%ntiles(n)
!       index1 = LIS_domain(n)%tile(c)%index
!       if ( LIS_forc(n,findex)%metdata1(1,index1) .ne. LIS_rc%udef .and. &
!            LIS_forc(n,findex)%metdata2(1,index1) .ne. LIS_rc%udef ) then
!           pcp(c) = wt1*LIS_forc(n,findex)%metdata1(1,index1) + wt2*LIS_forc(n,findex)%metdata2(1,index1) 
!       endif
!       if ( LIS_forc(n,findex)%metdata1(2,index1) .ne. LIS_rc%udef .and. &
!            LIS_forc(n,findex)%metdata2(2,index1) .ne. LIS_rc%udef ) then
!          tmp(c)  = wt1*LIS_forc(n,findex)%metdata1(2,index1) + wt2*LIS_forc(n,findex)%metdata2(2,index1)
!       endif
!
!       write(1001, '(2I8,6F15.6)'), index1, findex, wt1, wt2,              &
!                                    LIS_forc(n,findex)%metdata1(1,index1), & 
!                                    LIS_forc(n,findex)%metdata2(1,index1), &
!                                    LIS_forc(n,findex)%metdata1(2,index1), &
!                                    LIS_forc(n,findex)%metdata2(2,index1)
!       ! set uwind to be the constant wind speed
!       uwind(c) = const_wind(n)  
!       ! set vwind to be zero
!       vwind(c) = 0.0 
!    enddo
!    close(1001)
    
    do c = 1,LIS_rc%ntiles(n)
       index1 = LIS_domain(n)%tile(c)%index
       pcp(c) = rdhm356_struc_precip(n)%metdata2(index1) 
       tmp(c) = rdhm356_struc_temper(n)%metdata2(index1)

       ! set uwind to be the constant wind speed
       uwind(c) = const_wind(n)  
       ! set vwind to be zero
       vwind(c) = 0.0 
    enddo
    close(1001)
    deallocate(ratio)

    !call LIS_endrun()
  end subroutine timeinterp_rdhm356
