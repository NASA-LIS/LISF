!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp36_qcsoilmLAIveg
! \label{noahmp36_qcsoilmLAIveg}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 1 Aug 2016: Mahdi Navari; Modified for NoahMP36 
! 18 Jun 2021: Michel Bechtold: SM and LAI updating with S1 backscatter w/ WCM
!
! !INTERFACE:
subroutine noahmp36_qcsoilmLAI(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use noahmp36_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the soilmoisture related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)       :: sm1Field
  type(ESMF_Field)       :: laiField
  integer                :: t
  integer                :: status
  real, pointer          :: soilm1(:)
  real, pointer          :: lai(:)
  real                   :: smmax1!,smmax2,smmax3,smmax4
  real                   :: smmin1!,smmin2,smmin3,smmin4
  real                   :: laimax
  real                   :: laimin
  integer                :: gid
  real                   :: laitmp

  logical                :: update_flag(LIS_rc%ngrid(n))
  real                   :: perc_violation(LIS_rc%ngrid(n))
  real                   :: laimean(LIS_rc%ngrid(n))
  integer                :: nlaimean(LIS_rc%ngrid(n))
  integer                :: N_ens
  real                   :: state_tmp(LIS_rc%nensem(n)),state_mean

  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet for Soil Moisture Layer 1 failed in noahmp36_qcsoilmLAI")
 
  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet for Soil Moisture Layer 1 failed in noahmp36_qcsoilmLAI")

  call ESMF_AttributeGet(sm1Field,"Max Value",smmax1,rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: Max Value failed in noahmp36_qcsoilmLAI")

  call ESMF_AttributeGet(sm1Field,"Min Value",smmin1,rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: Min Value failed in noahmp36_qcsoilmLAI")

  call ESMF_StateGet(LSM_State,"LAI",laiField,rc=status)
  call LIS_verify(status,&
           "ESMF_StateGet for LAI failed in noahmp36_qcsoilmLAI")
  call ESMF_FieldGet(laiField,localDE=0,farrayPtr=lai,rc=status)
  call LIS_verify(status,&
           "ESMF_FieldGet for LAI failed in noahmp36_qcsoilmLAI")

  call ESMF_AttributeGet(laiField,"Max Value",laimax,rc=status)
  call LIS_verify(status,&
           "ESMF_AttributeGet for LAI Max Value failed in noahmp36_qcsoilmLAI")
  call ESMF_AttributeGet(laiField,"Min Value",laimin,rc=status)
  call LIS_verify(status,&
           "ESMF_AttributeGet for LAI Min Value failed in noahmp36_qcsoilmLAI")



  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     if(soilm1(t).gt.smmax1) soilm1(t) = smmax1
     if(soilm1(t).lt.smmin1) soilm1(t) = smmin1
  enddo

  update_flag    = .true.
  perc_violation = 0.0
  laimean       = 0.0
  nlaimean      = 0

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     laitmp =  lai(t)

     if(laitmp.lt.laimin.or.laitmp.gt.laimax) then
        update_flag(gid) = .false.
        perc_violation(gid) = perc_violation(gid) +1
     endif

  enddo

  do gid=1,LIS_rc%ngrid(n)
     perc_violation(gid) = perc_violation(gid)/LIS_rc%nensem(n)
  enddo

! For ensembles that are unphysical, compute the
! ensemble average after excluding them. This
! is done only if the majority of the ensemble
! members are good (>60%)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)
     if(.not.update_flag(gid)) then
        if(perc_violation(gid).lt.0.8) then
           if((lai(t).gt.laimin).and.&
                (lai(t).lt.laimax)) then 
              laimean(gid) = laimean(gid) + &
                   lai(t) 
              nlaimean(gid) = nlaimean(gid) + 1
           endif
        endif
     endif
  enddo
  
  do gid=1,LIS_rc%ngrid(n)
     if(nlaimean(gid).gt.0) then
        laimean(gid) = laimean(gid)/nlaimean(gid)
     endif
  enddo


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     laitmp =  lai(t)

! If the update is unphysical, simply set to the average of
! the good ensemble members. If all else fails, do not
! update.

     if(update_flag(gid)) then
        lai(t) = laitmp
     elseif(perc_violation(gid).lt.0.8) then
        if(laitmp.lt.laimin.or.laitmp.gt.laimax) then
           lai(t) = laimean(gid)
        else
           lai(t) = lai(t) 
        endif
     endif
  enddo

#if 0 
  N_ens = LIS_rc%nensem(n)
  do t=1,N_ens
     state_tmp(t) = lai(t)
  enddo
  state_mean =sum(state_tmp)/N_ens

  write(113,fmt='(i4.4,i2.2,i2.2,i2.2,i2.2,i2.2,21F8.3)') &
       LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, &
       LIS_rc%mn, LIS_rc%ss, &
       state_mean, &
       state_tmp
#endif
end subroutine noahmp36_qcsoilmLAI

