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
! !ROUTINE: noah36_updatesnodep
! \label{noah36_updatesnodep}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!  02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
!  21 Jul 2011: James Geiger; Modified for Noah 3.2
! 01 May 2014: Yuqiong Liu; modified for better QC
!
! !INTERFACE:
subroutine noah36_updatesnodep(n, LSM_State, LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use noah36_lsmMod
  use LIS_logMod,   only : LIS_logunit, LIS_verify

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
!
! !DESCRIPTION:
!
!  Updates the related state prognostic variable objects for
!  SNODEP data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \item[LSM\_Incr\_State] ESMF State container for LSM state increments \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)       :: sweField, sweIncrField
  type(ESMF_Field)       :: snodField, snodIncrField

  integer                :: t
  integer                :: status
  integer                :: gid
  real, pointer          :: swe(:), sweincr(:)
  real, pointer          :: snod(:), snodincr(:)
  real                   :: swetmp, snodtmp,sndens
  logical                :: update_flag(LIS_rc%ngrid(n))
  real                   :: perc_violation(LIS_rc%ngrid(n))

  real                   :: snodmean(LIS_rc%ngrid(n))
  integer                :: nsnodmean(LIS_rc%ngrid(n))

#if 0 
  integer :: svk_col,svk_row,ii,jj
  real    :: svk_statebf(LIS_rc%lnc(n),LIS_rc%lnr(n))
#endif

  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Snowdepth",snodField,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LSM_Incr_State,"SWE",sweIncrField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_Incr_State,"Snowdepth",snodIncrField,rc=status)
  call LIS_verify(status)

 
  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snodField,localDE=0,farrayPtr=snod,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweIncrField,localDE=0,farrayPtr=sweincr,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snodIncrField,localDE=0,farrayPtr=snodincr,rc=status)
  call LIS_verify(status)

  update_flag    = .true. 
  perc_violation = 0.0
  snodmean       = 0.0 
  nsnodmean      = 0
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row) 

     swetmp = swe(t) + sweincr(t)
     snodtmp = snod(t) + snodincr(t)
     
!     if((swetmp.lt.0.or.snodtmp.lt.0).or.&
!         (swetmp.gt.snodtmp)) then 

     if((snodtmp.lt.0)) then 
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
           if(snod(t)+snodincr(t).ge.0) then 
              snodmean(gid) = snodmean(gid) + &
                   snod(t)+snodincr(t)
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
          
     snodtmp = snod(t) + snodincr(t)

!Use the model's snow density from the previous timestep
     sndens = 0.0
     if(noah36_struc(n)%noah(t)%snowh.gt.0) then 
        sndens = noah36_struc(n)%noah(t)%sneqv/noah36_struc(n)%noah(t)%snowh
     endif
        
! If the update is unphysical, simply set to the average of 
! the good ensemble members. If all else fails, do not 
! update. 

     if(update_flag(gid)) then 
        snod(t) = snodtmp
     elseif(perc_violation(gid).lt.0.3) then 
        if(snodtmp.lt.0) then 
           snod(t) = snodmean(gid)
        else
           snod(t) = snod(t) + snodincr(t)
        endif
     endif
     swe(t) = snod(t)*sndens

  enddo
  
#if 0 
  svk_statebf = 0.0
  
  do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     svk_col = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col
     svk_row = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row
     
     svk_statebf(svk_col,svk_row) =  svk_statebf(svk_col,svk_row) + &
          snod(t)
  enddo

  do jj=1,LIS_rc%lnr(n)
     do ii=1,LIS_rc%lnc(n)
        if(svk_statebf(ii,jj).gt.0) then 
           svk_statebf(ii,jj) = svk_statebf(ii,jj)/LIS_rc%nensem(n)
        endif
     enddo
  enddo

  open(100,file='stateupd.bin',form='unformatted')
  write(100) svk_statebf
  close(100)
  stop
#endif

end subroutine noah36_updatesnodep

