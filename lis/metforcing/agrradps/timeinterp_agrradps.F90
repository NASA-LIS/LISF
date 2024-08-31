!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: timeinterp_agrradps
! \label{timeinterp_agrradps}
! 
! !REVISION HISTORY: 
!  20 Jan 2006:   Yudong Tian: Initial Implementation
!  10 Oct 2006:   Sujay Kumar: Switched to using ESMF_State for storing
!                              forcing data. 
! !INTERFACE:
subroutine timeinterp_agrradps(n,findex)
! !USES:
    use ESMF
    use LIS_coreMod,         only : LIS_rc, LIS_domain
    use LIS_FORC_AttributesMod 
    use LIS_metforcingMod,   only : LIS_FORC_Base_State, LIS_forc
    use LIS_logMod,          only : LIS_verify
    use agrradps_forcingMod, only : agrradps_struc

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n
    integer, intent(in) :: findex

! !DESCRIPTION:
!  Temporally interpolates the forcing data to the current model 
!  timestep. 
!  All variables except precipitation is linearly interpolated. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!
!EOP
    integer :: t,index
    integer :: f
    integer :: tid, tid1, tid2
    real    :: wt1,wt2
    
    integer          :: status
    type(ESMF_Field) :: swdField
    type(ESMF_Field) :: lwdField
    real, pointer    :: swd(:),lwd(:)
!EOP
  
  wt1 = (agrradps_struc(n)%agrtime2-LIS_rc%time) / & 
        (agrradps_struc(n)%agrtime2-agrradps_struc(n)%agrtime1)
  wt2 = 1.0 - wt1
  

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),    &
                     trim(LIS_FORC_SWdown%varname(1)), &
                     swdField,                         &
                     rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),    &
                     trim(LIS_FORC_LWdown%varname(1)), &
                     lwdField,                         &
                     rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
  call LIS_verify(status)

  do t = 1,LIS_rc%ntiles(n)
     index = LIS_domain(n)%tile(t)%index
     if((agrradps_struc(n)%metdata1(1,index).ne.LIS_rc%udef) .and. &
          (agrradps_struc(n)%metdata2(1,index).ne.LIS_rc%udef)) then  
        swd(t) = &
             wt1 * agrradps_struc(n)%metdata1(1,index) +  & 
             wt2 *agrradps_struc(n)%metdata2(1,index)
     endif
  enddo
  
  do t = 1,LIS_rc%ntiles(n)
     index = LIS_domain(n)%tile(t)%index
     if((agrradps_struc(n)%metdata1(2,index).ne.LIS_rc%udef) .and. &
          (agrradps_struc(n)%metdata2(2,index).ne.LIS_rc%udef)) then  
        lwd(t) = &
             wt1 * agrradps_struc(n)%metdata1(2,index) +  & 
             wt2 * agrradps_struc(n)%metdata2(2,index)
     endif
  enddo

end subroutine timeinterp_agrradps
