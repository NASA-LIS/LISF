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
! !ROUTINE: noahmp36_qcalbedo
! \label{noahmp36_qcalbedo}
!
! !REVISION HISTORY:
! 22 Dec 2017: Sujay Kumar; Initial Specification
! !INTERFACE:
subroutine noahmp36_qcalbedo(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod
  use noahmp36_lsmMod
  use LIS_logMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  QC's the related state prognostic variable objects for
!  veg data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)       :: albdField, albiField
  integer                :: t
  integer                :: status
  real, pointer          :: albd(:),albi(:)

  real                   :: albdmax,albimax
  real                   :: albdmin,albimin

  integer                :: gid
  real                   :: albdtmp,albitmp

  logical                :: update_flag(LIS_rc%ngrid(n))
  real                   :: perc_violation(LIS_rc%ngrid(n))

  real                   :: albdmean(LIS_rc%ngrid(n))
  integer                :: nalbdmean(LIS_rc%ngrid(n))

  real                   :: albimean(LIS_rc%ngrid(n))
  integer                :: nalbimean(LIS_rc%ngrid(n))
 
  integer                :: N_ens
  real                   :: state_tmp(LIS_rc%nensem(n)),state_mean

  call ESMF_StateGet(LSM_State,"ALBD",albdField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"ALBI",albiField,rc=status)
  call LIS_verify(status)
 
  call ESMF_FieldGet(albdField,localDE=0,farrayPtr=albd,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(albiField,localDE=0,farrayPtr=albi,rc=status)
  call LIS_verify(status)

  call ESMF_AttributeGet(albdField,"Max Value",albdmax,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(albdField,"Min Value",albdmin,rc=status)
  call LIS_verify(status)

  call ESMF_AttributeGet(albiField,"Max Value",albimax,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(albiField,"Min Value",albimin,rc=status)
  call LIS_verify(status)

  update_flag    = .true.
  perc_violation = 0.0
  albdmean       = 0.0
  nalbdmean      = 0

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     albdtmp =  albd(t)

     if(albdtmp.lt.albdmin.or.albdtmp.gt.albdmax) then
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
           if((albd(t).gt.albdmin).and.&
                (albd(t).lt.albdmax)) then 
              albdmean(gid) = albdmean(gid) + &
                   albd(t) 
              nalbdmean(gid) = nalbdmean(gid) + 1
           endif
        endif
     endif
  enddo
  
  do gid=1,LIS_rc%ngrid(n)
     if(nalbdmean(gid).gt.0) then
        albdmean(gid) = albdmean(gid)/nalbdmean(gid)
     endif
  enddo


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     albdtmp =  albd(t)

! If the update is unphysical, simply set to the average of
! the good ensemble members. If all else fails, do not
! update.

     if(update_flag(gid)) then
        albd(t) = albdtmp
     elseif(perc_violation(gid).lt.0.8) then
        if(albdtmp.lt.albdmin.or.albdtmp.gt.albdmax) then
           albd(t) = albdmean(gid)
        else
           albd(t) = albd(t) 
        endif
     endif
  enddo


  update_flag    = .true.
  perc_violation = 0.0
  albdmean       = 0.0
  nalbdmean      = 0

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     albitmp =  albi(t)

     if(albitmp.lt.albimin.or.albitmp.gt.albimax) then
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
           if((albi(t).gt.albimin).and.&
                (albi(t).lt.albimax)) then 
              albimean(gid) = albimean(gid) + &
                   albi(t) 
              nalbimean(gid) = nalbimean(gid) + 1
           endif
        endif
     endif
  enddo
  
  do gid=1,LIS_rc%ngrid(n)
     if(nalbimean(gid).gt.0) then
        albimean(gid) = albimean(gid)/nalbimean(gid)
     endif
  enddo


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     albitmp =  albi(t)

! If the update is unphysical, simply set to the average of
! the good ensemble members. If all else fails, do not
! update.

     if(update_flag(gid)) then
        albi(t) = albitmp
     elseif(perc_violation(gid).lt.0.8) then
        if(albitmp.lt.albimin.or.albitmp.gt.albimax) then
           albi(t) = albimean(gid)
        else
           albi(t) = albi(t) 
        endif
     endif
  enddo

end subroutine noahmp36_qcalbedo

