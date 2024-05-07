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
! !ROUTINE: noahmp401_qcsnodep
! \label{noahmp401_qcsnow}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
! 21 Jul 2011: James Geiger; Modified for Noah 3.2
! 30 Jan 2015: Yuqiong Liu; added additional QC
! 03 Oct 2018: Yeosang Yoon; Modified for NoahMP 3.6
! 14 Dec 2018: Yeosang Yoon; Modified for NoahMP 4.0.1
! 09 Jan 2020: Yeosang Yoon; update QC
!
! !INTERFACE:
subroutine noahmp401_qcsnodep(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod
  use noahmp401_lsmMod
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
  integer                :: t, gid
  integer                :: status
  real, pointer          :: swe(:)
  real, pointer          :: snod(:)

  real                   :: swemax,snodmax
  real                   :: swemin,snodmin

  real                   :: sndens
  logical                :: update_flag(LIS_rc%ngrid(n))
 
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
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     if((snod(t).lt.snodmin) .or. swe(t).lt.swemin) then
        update_flag(gid) = .false.
     endif

  enddo

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

!Use the model's snow density from the previous timestep
     sndens = 0.0
     if(noahmp401_struc(n)%noahmp401(t)%snowh.gt.0) then
       sndens = noahmp401_struc(n)%noahmp401(t)%sneqv/noahmp401_struc(n)%noahmp401(t)%snowh
     endif

!If the update is unphysical, do not update.
     if(update_flag(gid)) then
        snod(t) = snod(t)
        swe(t)  = snod(t)*sndens
     else ! do not update
        snod(t) = noahmp401_struc(n)%noahmp401(t)%snowh
        swe(t)  = noahmp401_struc(n)%noahmp401(t)%sneqv
     end if

     if(swe(t).gt.swemax) then
        swe(t) = swemax
     endif
     if(snod(t).gt.snodmax) then
        snod(t) = snodmax
     endif

  end do

end subroutine noahmp401_qcsnodep

