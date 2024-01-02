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
! !ROUTINE: timeinterp_petusgs
! \label{timeinterp_petusgs}
! 
! !REVISION HISTORY: 
!  20 Jan 2006:  Yudong Tian: Initial Implementation
!  10 Oct 2006:  Sujay Kumar: Switched to using ESMF_State for storing
!                              forcing data. 
!  10 Mar 2012:  K. Arsenault: Applied to USGS PET dataset
!  25 Oct 2013:  K. Arsenault: added PET USGS to LIS7
!
! !INTERFACE:
subroutine timeinterp_petusgs( n, findex )

! !USES:
  use ESMF
  use LIS_coreMod,        only : LIS_rc,LIS_domain
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod,  only : LIS_FORC_Base_State, LIS_forc
  use LIS_logMod,         only : LIS_verify
  use petusgs_forcingMod, only : petusgs_struc
  use LIS_forecastMod, only : LIS_get_iteration_index

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Temporally interpolates the forcing data to the current model 
!   timestep. 
!
!  Also, daily PET values are scaled up (i.e. multiplied) by a factor 
!   of 100 to preserve the precision, so here the read-in values are
!   divided by 100 to obtain actual PET summary.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!
!EOP

   integer          :: t, index1
   integer          :: mfactor, m, k, kk
   integer          :: status
   type(ESMF_Field) :: petfield
!   type(ESMF_Field) :: rETfield
   real, pointer    :: pet(:)
!   real, pointer    :: rET(:)

! ----------------------------------------

    call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_PET%varname(1)),petfield,&
          rc=status)
    call LIS_verify(status, 'Error: enable PET in forcing variables list')

!    call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_RefET%varname(1)),rETfield,&
!          rc=status)
!    call LIS_verify(status, 'Error: enable RefET in forcing variables list')

    call ESMF_FieldGet(petfield,localDE=0,farrayptr=pet,rc=status)
    call LIS_verify(status)

!    call ESMF_FieldGet(retField,localDE=0,farrayptr=rET,rc=status)
!    call LIS_verify(status)

!------------------------------------------------------------------------

   mfactor = LIS_rc%nensem(n)/petusgs_struc(n)%nIter

   do k=1,LIS_rc%ntiles(n)/mfactor
     do m=1,mfactor
       t = m + (k-1)*mfactor
       index1 = LIS_domain(n)%tile(t)%index
       kk = LIS_get_iteration_index(n, k, index1, mfactor)

       if( petusgs_struc(n)%metdata2(kk,1,index1) .ne. -1.0 ) &

       !- Convert scaled daily PET to actual values (precision: 0.01):
          pet(t) = petusgs_struc(n)%metdata2(kk,1,index1) / 100.     

       !- Convert mm/day to rate with seconds:
          pet(t) = pet(t) / 86400.     ! mm/s

!       if (petusgs_struc(n)%metdata2(kk,1,index1) .ne. -1.0) &
!        rET(t) = petusgs_struc(n)%metdata2(kk,1,index1)
      enddo
    enddo

end subroutine timeinterp_petusgs
  
