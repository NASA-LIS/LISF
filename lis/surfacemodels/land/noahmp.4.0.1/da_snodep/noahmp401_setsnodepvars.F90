!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp401_setsnodepvars
! \label{noahmp401_setsnodepvars}
!
! !REVISION HISTORY:
! 15 Aug 2017: Sujay Kumar; Initial Specification
! 03 Oct 2018: Yeosang Yoon; Modified for NoahMP 3.6
!
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
! 21 Jul 2011: James Geiger; Modified for Noah 3.2
! 03 Oct 2018: Yeosang Yoon; Modified for NoahMP 3.6
! 14 Dec 2018: Yeosang Yoon; Modified for NoahMP 4.0.1
!
! !INTERFACE:
subroutine noahmp401_setsnodepvars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_surface
  use LIS_snowMod, only : LIS_snow_struc
  use LIS_logMod, only : LIS_logunit, LIS_verify
  use noahmp401_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
! 
! !DESCRIPTION:
! 
!  This routine assigns the snow progognostic variables to noah's
!  model space. The state vector consists of total SWE and snow depth. 
!  This routine also updates other model prognostics (snice, snliq,
!  snow thickness, snow temperature) based on the update. 
! 
!EOP
  type(ESMF_Field)       :: sweField
  type(ESMF_Field)       :: snodField
  real, pointer          :: swe(:)
  real, pointer          :: snod(:)
  real                   :: dsneqv,dsnowh
  integer                :: t
  integer                :: status
  
  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Snowdepth",snodField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snodField,localDE=0,farrayPtr=snod,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     dsneqv = swe(t)*1000.0 - noahmp401_struc(n)%noahmp401(t)%sneqv !in mm
     dsnowh = snod(t) - noahmp401_struc(n)%noahmp401(t)%snowh  !in m

!NOTE: if snow becomes unphysical after applying the deltas, then 
! the update is rejected. This may cause inconsistent ensemble 
! updates. TBD

     call noahmp401_snodep_update(n, t, dsneqv, dsnowh)

     if(noahmp401_struc(n)%noahmp401(t)%sneqv.lt.0.or.&
          noahmp401_struc(n)%noahmp401(t)%snowh.lt.0) then 
        print*, dsneqv, dsnowh
        print*, swe(t), snod(t)
        stop
     endif
  enddo
end subroutine noahmp401_setsnodepvars


