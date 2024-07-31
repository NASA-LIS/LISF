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
! !ROUTINE: noah36_qcsnodep
! \label{noah36_qcsnodep}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!  02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
!  21 Jul 2011: James Geiger; Modified for Noah 3.2
!  30 Jan 2015: Yuqiong Liu; added additional QC
! !INTERFACE:
subroutine noah36_qcsnodep(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod
  use noah36_lsmMod
  use LIS_logMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  QC's the related state prognostic variable objects for
!  SNODEP data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)       :: sweField
  type(ESMF_Field)       :: snodField
  integer                :: t,gid
  integer                :: status
  real, pointer          :: swe(:)
  real, pointer          :: snod(:)

  real                   :: swemax,snodmax
  real                   :: swemin,snodmin
  real                   :: sndens
  logical                :: update_flag(LIS_rc%ngrid(n))
  real                   :: perc_violation(LIS_rc%ngrid(n))

  real                   :: snodmean(LIS_rc%ngrid(n))
  integer                :: nsnodmean(LIS_rc%ngrid(n))
 

  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Snowdepth",snodField,rc=status)
  call LIS_verify(status)
 
  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snodField,localDE=0,farrayPtr=snod,rc=status)
  call LIS_verify(status)

  call ESMF_AttributeGet(sweField,"Max Value",swemax,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(sweField,"Min Value",swemin,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(snodField,"Max Value",snodmax,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(snodField,"Min Value",snodmin,rc=status)
  call LIS_verify(status)

  update_flag    = .true. 
  perc_violation = 0.0
  snodmean       = 0.0 
  nsnodmean      = 0
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row) 

     if((snod(t).lt.snodmin)) then 
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
        if(perc_violation(gid).lt.0.4) then
           if(snod(t).ge.0) then 
              snodmean(gid) = snodmean(gid) + &
                   snod(t)
              nsnodmean(gid) = nsnodmean(gid) + 1
           endif
        endif
     endif
  enddo

  do gid=1,LIS_rc%ngrid(n)
     if(nsnodmean(gid).gt.0) then 
        snodmean(gid) = snodmean(gid)/nsnodmean(gid)
     endif
  enddo

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row) 
          
!Use the model's snow density from the previous timestep
     sndens = 0.0
     if(noah36_struc(n)%noah(t)%snowh.gt.0) then 
        sndens = noah36_struc(n)%noah(t)%sneqv/noah36_struc(n)%noah(t)%snowh
     endif
        
! If the update is unphysical, simply set to the average of 
! the good ensemble members. If all else fails, do not 
! update. 

     if(update_flag(gid)) then 
        snod(t) = snod(t)
     elseif(perc_violation(gid).lt.0.3) then 
        if(snod(t).lt.snodmin) then 
           snod(t) = snodmean(gid)
        else
           snod(t) =noah36_struc(n)%noah(t)%snowh
        endif
     endif

     if(snod(t).gt.snodmax) then 
        snod(t) = snodmax
     endif

     swe(t) = snod(t)*sndens

  enddo
end subroutine noah36_qcsnodep

