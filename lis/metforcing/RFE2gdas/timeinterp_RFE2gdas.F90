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
! !ROUTINE: timeinterp_RFE2gdas
! \label{timeinterp_RFE2gdas}
!
! !REVISION HISTORY:
!
! 2 June 2010: Soni Yatheendradas; Initial LIS version for FEWSNET
! 8 Mar 2012: Updated for LIS7
!
! !INTERFACE:
subroutine timeinterp_RFE2gdas(n,findex)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_domain
  use LIS_metforcingMod, only : LIS_FORC_Base_State, LIS_forc
  use LIS_FORC_AttributesMod
  use LIS_logMod, only : LIS_verify,LIS_logunit
  use RFE2gdas_forcingMod, only : RFE2gdas_struc
  use LIS_forecastMod, only : LIS_get_iteration_index

  implicit none

! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Temporally interpolates the forcing data to the current model
!  timestep (evenly distributed over sub-daily timesteps for now). 
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    index of the supplemental forcing source
!  \item[suppdata1]
!    previous forcing values
!  \item[suppdata2]
!    next forcing values
!  \end{description}
!
!EOP

  integer :: t,m,index1,mfactor,kk,k
  
  integer          :: status
  type(ESMF_Field) :: pcpField, cpcpField
  real, pointer    :: pcp(:), cpcp(:)
  real,  allocatable :: ratio(:)


  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),&
       trim(LIS_FORC_Rainf%varname(1)),pcpField,&
       rc=status)
  call LIS_verify(status, 'Error: enable Rainf in forcing variables list')
  
  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LIS_verify(status)

#if 0
!- Applying a convective rainfall amount to observed precip field: 
   if( LIS_FORC_CRainf%selectOpt == 1 ) then
     allocate(ratio(LIS_rc%ntiles(n)))

     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),&
          trim(LIS_FORC_CRainf%varname(1)),cpcpField,&
          rc=status)
     call LIS_verify(status, 'Error: enable CRainf in forcing variables list')
  
     call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
     call LIS_verify(status)
  !------------------------------------------------------------------------
  ! Compute ratio between convective model precip and total model precip
  ! so that it can be applied to the observed global precip
  !------------------------------------------------------------------------
     do t = 1,LIS_rc%ntiles(n)
        if( cpcp(t) > 0. ) then 
        endif
        if( pcp(t) .ne. 0.0 .and.  &
            pcp(t) .ne. LIS_rc%udef .and.  &
           cpcp(t) .ne. LIS_rc%udef) then
          ratio(t) = cpcp(t) / pcp(t)
          if( cpcp(t) > 0. ) then
          endif
          if (ratio(t) .gt. 1.0) ratio(t) = 1.0
          if (ratio(t) .lt. 0.0) ratio(t) = 0.0
        else
          ratio(t) = 0.0
        endif
     enddo
   endif
#endif

!- Convert precipitation sum to rate:
  mfactor = LIS_rc%nensem(n)/RFE2gdas_struc(n)%nIter
  do k=1,LIS_rc%ntiles(n)/mfactor
     do m=1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        kk = LIS_get_iteration_index(n, k, index1, mfactor)
        if( RFE2gdas_struc(n)%metdata2(kk,1,index1) .ge. 0.0 ) then
           pcp(t) = RFE2gdas_struc(n)%metdata2(kk,1,index1)
        else
           pcp(t) = LIS_rc%udef
!         cpcp(t) = LIS_rc%udef
        endif
     enddo
  enddo

#if 0
   if( LIS_FORC_CRainf%selectOpt == 1 ) then
     deallocate(ratio)
   endif
#endif

end subroutine timeinterp_RFE2gdas
