!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
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
    use LDT_coreMod, only: LDT_rc, LDT_domain
    use LDT_FORC_AttributesMod
    use LDT_metforcingMod, only : LDT_FORC_Base_State, LDT_forc
    use LDT_logMod,         only : LDT_verify, LDT_endrun 
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


    wt1=(rdhm356_struc_precip(n)%rdhm356time2 - LDT_rc%time)/                            &
         (rdhm356_struc_precip(n)%rdhm356time2 - rdhm356_struc_precip(n)%rdhm356time1)
    wt2=1.0-wt1

    call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_Tair%varname(1),tmpField,&
         rc=status)
    call LDT_verify(status, 'Error: Enable Tair in the forcing variables list')

    call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_Rainf%varname(1),pcpField,&
         rc=status)
    call LDT_verify(status, 'Error: Enable Rainf in the forcing variables list')
  
    call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_Wind_E%varname(1),uField,&
         rc=status)
    call LDT_verify(status, 'Error: Enable Wind_E in the forcing variables list')
  
    call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_Wind_N%varname(1),vField,&
         rc=status)
    call LDT_verify(status, 'Error: Enable Wind_N in the forcing variables list')

    call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
    call LDT_verify(status)
    
    call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
    call LDT_verify(status)
    
    call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
    call LDT_verify(status)
      
    call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
    call LDT_verify(status)
 
!    write(fname, *), LDT_rc%ts
!    open(unit=1001, file='time_interp', status='unknown')
!    do c = 1,LDT_rc%ntiles(n)
!       index1 = LDT_domain(n)%tile(c)%index
!       if ( LDT_forc(n,findex)%metdata1(1,index1) .ne. LDT_rc%udef .and. &
!            LDT_forc(n,findex)%metdata2(1,index1) .ne. LDT_rc%udef ) then
!           pcp(c) = wt1*LDT_forc(n,findex)%metdata1(1,index1) + wt2*LDT_forc(n,findex)%metdata2(1,index1) 
!       endif
!       if ( LDT_forc(n,findex)%metdata1(2,index1) .ne. LDT_rc%udef .and. &
!            LDT_forc(n,findex)%metdata2(2,index1) .ne. LDT_rc%udef ) then
!          tmp(c)  = wt1*LDT_forc(n,findex)%metdata1(2,index1) + wt2*LDT_forc(n,findex)%metdata2(2,index1)
!       endif
!
!       write(1001, '(2I8,6F15.6)'), index1, findex, wt1, wt2,              &
!                                    LDT_forc(n,findex)%metdata1(1,index1), & 
!                                    LDT_forc(n,findex)%metdata2(1,index1), &
!                                    LDT_forc(n,findex)%metdata1(2,index1), &
!                                    LDT_forc(n,findex)%metdata2(2,index1)
!       ! set uwind to be the constant wind speed
!       uwind(c) = const_wind(n)  
!       ! set vwind to be zero
!       vwind(c) = 0.0 
!    enddo
!    close(1001)
    
    do c = 1,LDT_rc%ntiles(n)
       index1 = LDT_domain(n)%tile(c)%index
       pcp(c) = LDT_forc(n,findex)%metdata2(1,index1) 
       tmp(c) = LDT_forc(n,findex)%metdata2(2,index1)

       ! set uwind to be the constant wind speed
       uwind(c) = const_wind(n)  
       ! set vwind to be zero
       vwind(c) = 0.0 
    enddo
    close(1001)

  end subroutine timeinterp_rdhm356
