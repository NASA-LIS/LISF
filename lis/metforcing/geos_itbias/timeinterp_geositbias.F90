!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: timeinterp_geositbias
! \label{timeinterp_geositbias}
!
! !REVISION HISTORY:
! 02 Oct 2025: Fadji Maina, initial code (based on geos-it)
!
! !INTERFACE:
subroutine timeinterp_geositbias(n,findex)

! !USES:
  use ESMF
  use LIS_FORC_AttributesMod
  use LIS_coreMod,       only : LIS_rc,LIS_domain,LIS_localPet
  use LIS_metforcingMod, only : LIS_forc,LIS_FORC_Base_State
  use LIS_constantsMod,  only : LIS_CONST_SOLAR
  use LIS_timeMgrMod,    only : LIS_tick
  use LIS_logMod,        only : LIS_verify,LIS_endrun
  use geositbias_forcingMod,   only : geositbias_struc
  use LIS_forecastMod,   only : LIS_get_iteration_index
  use LIS_ran2_gasdev

  implicit none

! !ARGUMENTS:
  integer, intent(in):: n
  integer, intent(in):: findex
!
! !DESCRIPTION:
!  Temporally interpolates the forcing data to the current model
!  timestep.  Because all variables are a 1-hourly time average,
!  no interpolation in time is actually performed.  The identical
!  data is used for all timesteps within the hourly input file.
!
!  The routines invoked are:
!  \begin{description}
!   \item[LIS\_time2date](\ref{LIS_time2date}) \newline
!    converts the time to a date format
!   \item[zterp](\ref{zterp}) \newline
!    zenith-angle based interpolation
!  \end{description}
!EOP
  integer :: t,k,kk
  integer :: index1
  integer :: bdoy,byr,bmo
  integer :: bda,bhr,bmn
  integer :: bss
  real*8  :: btime
  real    :: gmt1,gmt2
  real    :: bts
  integer          :: status
  integer          :: mfactor,m
  type(ESMF_Field) :: tairField,qairField,psField,lwgabField
  real,pointer     :: tair(:),qair(:),ps(:),lwgab(:)

  btime = geositbias_struc(n)%geositbiastime1
  byr = LIS_rc%yr
  bmo = LIS_rc%mo
  bda = LIS_rc%da
  bhr = LIS_rc%hr
  bmn = 30
  bss = 0
  if (LIS_rc%mn.lt.30) then
     bts = -(60*60)
  else
     bts = 0
  endif
  call LIS_tick(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn,bss,bts)

  btime = geositbias_struc(n)%geositbiastime2
  byr = LIS_rc%yr             !next hour
  bmo = LIS_rc%mo
  bda = LIS_rc%da
  bhr = LIS_rc%hr
  bmn = 30
  bss = 0
  if (LIS_rc%mn.lt.30) then
     bts = 0
  else
     bts = 60*60
  endif
  call LIS_tick(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn,bss,bts)

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),                &
       LIS_FORC_Tair%varname(1),tairField,                         &
       rc=status)
  call LIS_verify(status,                                          &
       'Error: Enable Tair in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),                &
       LIS_FORC_Qair%varname(1),QairField,                         &
       rc=status)
  call LIS_verify(status,                                          &
       'Error: Enable Qair in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),                &
       LIS_FORC_LWdown%varname(1),lwgabField,                      &
       rc=status)
  call LIS_verify(status,                                          &
       'Error: Enable LWdown in the forcing variables list')


  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),                &
       LIS_FORC_Psurf%varname(1),psField,                          &
       rc=status)
  call LIS_verify(status,                                          &
       'Error: Enable Psurf in the forcing variables list')


  call ESMF_FieldGet(tairField,localDE=0,farrayPtr=tair,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(qairField,localDE=0,farrayPtr=qair,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(lwgabField,localDE=0,farrayPtr=lwgab,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(psField,localDE=0,farrayPtr=ps,rc=status)
  call LIS_verify(status)

  mfactor = LIS_rc%nensem(n)/geositbias_struc(n)%nIter

  do k = 1,LIS_rc%ntiles(n)/mfactor
     do m = 1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        kk = LIS_get_iteration_index(n,k,index1,mfactor)
        tair(t) = geositbias_struc(n)%metdata1(kk,1,index1)
        qair(t) = geositbias_struc(n)%metdata1(kk,2,index1)
        lwgab(t) = geositbias_struc(n)%metdata1(kk,4,index1)
        ps(t) = geositbias_struc(n)%metdata1(kk,7,index1)
     enddo
  enddo
end subroutine timeinterp_geositbias

