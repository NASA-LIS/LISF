!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: timeinterp_agrrad
! \label{timeinterp_agrrad}
! 
!
! !REVISION HISTORY:
! 08 Dec 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine timeinterp_agrrad(n, findex)
! !USES:
  use ESMF
  use LIS_coreMod,        only : LIS_rc,LIS_domain
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod,  only : LIS_FORC_Base_State, LIS_forc
  use LIS_logMod,         only : LIS_verify
  use agrrad_forcingMod,  only : agrrad_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
  
! !DESCRIPTION:
! 
! Temporally interpolates the AGRMET radiation forcing data to the 
! model timestep.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[findex]
!   index of the forcing source
!  \end{description}
! 
!EOP
  
  real    :: wt1,wt2
  integer :: t
  integer :: index1

  integer          :: status
  type(ESMF_Field) :: swdField
  type(ESMF_Field) :: lwdField
  real, pointer    :: swd(:),lwd(:)

  
  wt1 = (agrrad_struc(n)%agrtime2-LIS_rc%time) / & 
        (agrrad_struc(n)%agrtime2-agrrad_struc(n)%agrtime1)
  wt2 = 1.0 - wt1
  

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),    &
                     trim(LIS_FORC_SWdown%varname(1)), &
                     swdField,                         &
                     rc=status)
  call LIS_verify(status, 'Error: Enable SWdown in forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),    &
                     trim(LIS_FORC_LWdown%varname(1)), &
                     lwdField,                         &
                     rc=status)
  call LIS_verify(status, 'Error: enables LWdown in forcing variables list')

  call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
  call LIS_verify(status)

  do t = 1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     if ( (agrrad_struc(n)%metdata1(1,index1).ne.LIS_rc%udef) .and. &
          (agrrad_struc(n)%metdata2(1,index1).ne.LIS_rc%udef) ) then  
        swd(t) = &
             wt1 * agrrad_struc(n)%metdata1(1,index1) + & 
             wt2 * agrrad_struc(n)%metdata2(1,index1)
     endif
  enddo

  do t = 1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     if ( (agrrad_struc(n)%metdata1(2,index1).ne.LIS_rc%udef) .and. &
          (agrrad_struc(n)%metdata2(2,index1).ne.LIS_rc%udef) ) then  
        lwd(t) = &
             wt1 * agrrad_struc(n)%metdata1(2,index1) + & 
             wt2 * agrrad_struc(n)%metdata2(2,index1)
     endif
  enddo

end subroutine timeinterp_agrrad
