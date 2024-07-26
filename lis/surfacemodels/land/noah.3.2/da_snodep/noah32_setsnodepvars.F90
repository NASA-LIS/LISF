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
! !ROUTINE: noah32_setsnodepvars
! \label{noah32_setsnodepvars}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!  02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
!  21 Jul 2011: James Geiger; Modified for Noah 3.2
!
! !INTERFACE:
subroutine noah32_setsnodepvars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_domain
  use noah32_lsmMod
  use LIS_snowMod, only : LIS_snow_struc
  use LIS_logMod, only : LIS_logunit, LIS_verify

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
! 
! !DESCRIPTION:
! 
!  This routine assigns the snow progognostic variables to noah's
!  model space. 
! 
!EOP
  type(ESMF_State)       :: LSM_State
  type(ESMF_Field)       :: sweField
  type(ESMF_Field)       :: snodField
  real, allocatable          :: swe(:)
  real, allocatable          :: snod(:)
  real                   :: npts(LIS_rc%ngrid(n))
  integer                :: t,c,r,gid
  integer                :: status 
!  real                   :: temp(LIS_rc%lnc(n), LIS_rc%lnr(n))

  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Snowdepth",snodField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snodField,localDE=0,farrayPtr=snod,rc=status)
  call LIS_verify(status)

  LIS_snow_struc(n)%sneqv = 0.0
  LIS_snow_struc(n)%snowdepth = 0.0
  npts = 0 

  do t=1,LIS_rc%ntiles(n)
     noah32_struc(n)%noah(t)%sneqv = swe(t)
     noah32_struc(n)%noah(t)%snowh = snod(t)

     c = LIS_domain(n)%tile(t)%col
     r = LIS_domain(n)%tile(t)%row
     gid = LIS_domain(n)%gindex(c,r)     
   
     LIS_snow_struc(n)%sneqv(gid)   = LIS_snow_struc(n)%sneqv(gid)+swe(t)
     LIS_snow_struc(n)%snowdepth(t) = LIS_snow_struc(n)%snowdepth(gid)+snod(t)
     npts(gid) = npts(gid)+1     
  enddo

!  temp = -9999.0
  do t=1,LIS_rc%ngrid(n)
!     if(LIS_snow_struc(n)%sneqv(t).gt.50) print*, 'here',t,npts(t),&
!          LIS_snow_struc(n)%snowdepth(t),LIS_snow_struc(n)%sneqv(t)
!     temp(LIS_domain(n)%grid(t)%col, LIS_domain(n)%grid(t)%row) = LIS_snow_struc(n)%sneqv(t)
     LIS_snow_struc(n)%sneqv(t)     = LIS_snow_struc(n)%sneqv(t)/npts(t)
     LIS_snow_struc(n)%snowdepth(t) = LIS_snow_struc(n)%snowdepth(t)/npts(t)
  enddo
  
!  open(100,file='varfield.bin',form='unformatted')
!  write(100) temp
!  close(100)
!  stop
end subroutine noah32_setsnodepvars

