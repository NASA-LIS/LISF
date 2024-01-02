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
! !ROUTINE: clsmf25_qc_snowobs
! \label{clsmf25_qc_snowobs}
!
! !REVISION HISTORY:
! 25Feb2008: Sujay Kumar: Initial Specification
!
! !INTERFACE:
subroutine clsmf25_qc_snowobs(n,OBS_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod,  only : LIS_verify
  use LIS_constantsMod, only : LIS_CONST_TKFRZ
  use LIS_vegDataMod, only : LIS_lai,LIS_gfrac
  use clsmf25_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)      :: n
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION:
!
!  This subroutine performs any model-based QC of the observation 
!  prior to data assimilation. Here the snow observations
!  are flagged when LSM indicates that (1) rain is falling (2)
!  ground is fully or partially covered with snow. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF state container for observations \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)         :: obs_snow_field

  real, pointer            :: snowobs(:)
  integer                  :: t
  integer                  :: gid
  integer                  :: status

  
  call ESMF_StateGet(OBS_State,"Observation01",obs_snow_field,&
       rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet failed in clsmf25_qc_soilmobs")
  call ESMF_FieldGet(obs_snow_field,localDE=0,farrayPtr=snowobs,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet failed in clsmf25_qc_soilmobs")
  
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid  = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
!     if(LIS_rc%mo.eq.5.or.LIS_rc%mo.eq.6.or.LIS_rc%mo.eq.7.or.&
!          LIS_rc%mo.eq.8) then 
!        snowobs(gid) = LIS_rc%udef
!     endif

     if(snowobs(gid).ne.LIS_rc%udef) then 
        if(LIS_lai(n)%tlai(t).gt.2) then 
           snowobs(gid) = LIS_rc%udef        
!           print*, 'LAI too high ',gid,LIS_lai(n)%tlai(t)
        elseif(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt.le.3) then !forest types
           snowobs(gid) = LIS_rc%udef
!           print*, 'veg too high ',gid,LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt
           !assume that snow will not form at 5 deg. celcius or higher ground temp. 
        elseif(clsmf25_struc(n)%cat_diagn(t)%tsurf.ge.278.15) then 
           snowobs(gid) = LIS_rc%udef
!           print*, 'ground too warm ',gid,clsmf25_struc(n)%cat_diagn(n)%tsurf
        elseif(clsmf25_struc(n)%cat_diagn(t)%tp(1)+273.16.ge.278.15) then 
           snowobs(gid) = LIS_rc%udef
!           print*, 'tp too warm ',gid,clsmf25_struc(n)%cat_diagn(n)%tp(1)+273.16
        endif
     endif
  enddo

end subroutine clsmf25_qc_snowobs

