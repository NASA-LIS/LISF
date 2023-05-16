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
! !ROUTINE: timeinterp_geosit
! \label{timeinterp_geosit}
!
! !REVISION HISTORY:
! 18 Mar 2015: James Geiger, initial code (based on merra-land)
! 20 Apr 2023: David Mocko,  initial code (based on merra2)
!
! !INTERFACE:
      subroutine timeinterp_geosit(n,findex)

! !USES:
      use ESMF
      use LIS_FORC_AttributesMod
      use LIS_coreMod,       only : LIS_rc,LIS_domain,LIS_localPet
      use LIS_metforcingMod, only : LIS_forc,LIS_FORC_Base_State
      use LIS_constantsMod,  only : LIS_CONST_SOLAR
      use LIS_timeMgrMod,    only : LIS_tick
      use LIS_logMod,        only : LIS_logunit,LIS_verify,LIS_endrun
      use geosit_forcingMod,   only : geosit_struc
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
      integer :: t,zdoy,k,kk
      integer :: index1
      integer :: bdoy,byr,bmo
      integer :: bda,bhr,bmn
      integer :: bss
      real*8  :: btime
      real    :: wt1,wt2,czb,cze,czm,gmt1,gmt2
      real    :: zw1,zw2,bts
      integer          :: status
      integer          :: mfactor,m
      type(ESMF_Field) :: tairField,qairField,swgdnField,lwgabField
      type(ESMF_Field) :: uwindField,vwindField,psField
      type(ESMF_Field) :: prectotField,precconField,precsnoField
      type(ESMF_Field) :: swlandField,pardrField,pardfField,hlmlField
      real,pointer     :: tair(:),qair(:),swgdn(:),lwgab(:)
      real,pointer     :: uwind(:),vwind(:),ps(:)
      real,pointer     :: prectot(:),preccon(:),precsno(:)
      real,pointer     :: swland(:),pardr(:),pardf(:),hlml(:)

      btime = geosit_struc(n)%geosittime1
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

      btime = geosit_struc(n)%geosittime2
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
                         LIS_FORC_Tair%varname(1),tairField,           &
                         rc=status)
      call LIS_verify(status,                                          &
              'Error: Enable Tair in the forcing variables list')

      call ESMF_StateGet(LIS_FORC_Base_State(n,findex),                &
                         LIS_FORC_Qair%varname(1),QairField,           &
                         rc=status)
      call LIS_verify(status,                                          &
              'Error: Enable Qair in the forcing variables list')

      call ESMF_StateGet(LIS_FORC_Base_State(n,findex),                &
                         LIS_FORC_SWdown%varname(1),swgdnField,        &
                         rc=status)
      call LIS_verify(status,                                          &
              'Error: Enable SWdown in the forcing variables list')

      call ESMF_StateGet(LIS_FORC_Base_State(n,findex),                &
                         LIS_FORC_LWdown%varname(1),lwgabField,        &
                         rc=status)
      call LIS_verify(status,                                          &
              'Error: Enable LWdown in the forcing variables list')

      call ESMF_StateGet(LIS_FORC_Base_State(n,findex),                &
                         LIS_FORC_Wind_E%varname(1),uwindField,        &
                         rc=status)
      call LIS_verify(status,                                          &
              'Error: Enable Wind_E in the forcing variables list')

      call ESMF_StateGet(LIS_FORC_Base_State(n,findex),                &
                         LIS_FORC_Wind_N%varname(1),vwindField,        &
                         rc=status)
      call LIS_verify(status,                                          &
              'Error: Enable Wind_N in the forcing variables list')

      call ESMF_StateGet(LIS_FORC_Base_State(n,findex),                &
                         LIS_FORC_Psurf%varname(1),psField,            &
                         rc=status)
      call LIS_verify(status,                                          &
              'Error: Enable Psurf in the forcing variables list')

      call ESMF_StateGet(LIS_FORC_Base_State(n,findex),                &
                         LIS_FORC_Rainf%varname(1),prectotField,       &
                         rc=status)
      call LIS_verify(status,                                          &
              'Error: Enable Rainf in the forcing variables list')

      if (LIS_FORC_CRainf%selectOpt.eq.1) then
         call ESMF_StateGet(LIS_FORC_Base_State(n,findex),             &
                            LIS_FORC_CRainf%varname(1),precconField,   &
                            rc=status)
         call LIS_verify(status,                                       &
                 'Error: Enable CRainf in the forcing variables list')
      endif

      if (LIS_FORC_Snowf%selectOpt.eq.1) then
         call ESMF_StateGet(LIS_FORC_Base_State(n,findex),             &
                            LIS_FORC_Snowf%varname(1),precsnoField,    &
                            rc=status)
         call LIS_verify(status,                                       &
                 'Error: Enable Snowf in the forcing variables list')
      endif

      if (LIS_FORC_Swnet%selectOpt.eq.1) then
         call ESMF_StateGet(LIS_FORC_Base_State(n,findex),             &
                            LIS_FORC_Swnet%varname(1),swlandField,&
                            rc=status)
         call LIS_verify(status,                                       &
                 'Error: Enable SWnet in the forcing variables list')
      endif

      if (LIS_FORC_Pardr%selectOpt.eq.1) then
         call ESMF_StateGet(LIS_FORC_Base_State(n,findex),             &
                            LIS_FORC_Pardr%varname(1),pardrField,      &
                            rc=status)
         call LIS_verify(status,                                       &
                 'Error: Enable PARDR in the forcing variables list')
      endif

      if (LIS_FORC_Pardf%selectOpt.eq.1) then
         call ESMF_StateGet(LIS_FORC_Base_State(n,findex),             &
                            LIS_FORC_Pardf%varname(1),pardfField,      &
                            rc=status)
         call LIS_verify(status,                                       &
                 'Error: Enable PARDF in the forcing variables list')
      endif

      if (LIS_FORC_Forc_Hgt%selectOpt.eq.1) then
         call ESMF_StateGet(LIS_FORC_Base_State(n,findex),             &
                            LIS_FORC_Forc_Hgt%varname(1),hlmlField,    &
                            rc=status)
         call LIS_verify(status,                                       &
                 'Error: Enable Forc_Hgt in the forcing variables list')
      endif

      call ESMF_FieldGet(tairField,localDE=0,farrayPtr=tair,rc=status)
      call LIS_verify(status)

      call ESMF_FieldGet(qairField,localDE=0,farrayPtr=qair,rc=status)
      call LIS_verify(status)

      call ESMF_FieldGet(swgdnField,localDE=0,farrayPtr=swgdn,rc=status)
      call LIS_verify(status)

      call ESMF_FieldGet(lwgabField,localDE=0,farrayPtr=lwgab,rc=status)
      call LIS_verify(status)

      call ESMF_FieldGet(uwindField,localDE=0,farrayPtr=uwind,rc=status)
      call LIS_verify(status)

      call ESMF_FieldGet(vwindField,localDE=0,farrayPtr=vwind,rc=status)
      call LIS_verify(status)

      call ESMF_FieldGet(psField,localDE=0,farrayPtr=ps,rc=status)
      call LIS_verify(status)

      call ESMF_FieldGet(prectotField,localDE=0,farrayPtr=prectot,rc=status)
      call LIS_verify(status)

      if (LIS_FORC_CRainf%selectOpt.eq.1) then
         call ESMF_FieldGet(precconField,localDE=0,farrayPtr=preccon,rc=status)
         call LIS_verify(status)
      endif

      if (LIS_FORC_Snowf%selectOpt.eq.1) then
         call ESMF_FieldGet(precsnoField,localDE=0,farrayPtr=precsno,rc=status)
         call LIS_verify(status)
      endif

      if (LIS_FORC_Swnet%selectOpt.eq.1) then
         call ESMF_FieldGet(swlandField,localDE=0,farrayPtr=swland,rc=status)
         call LIS_verify(status)
      endif

      if (LIS_FORC_Pardr%selectOpt.eq.1) then
         call ESMF_FieldGet(pardrField,localDE=0,farrayPtr=pardr,rc=status)
         call LIS_verify(status)
      endif

      if (LIS_FORC_Pardf%selectOpt.eq.1) then
         call ESMF_FieldGet(pardfField,localDE=0,farrayPtr=pardf,rc=status)
         call LIS_verify(status)
      endif

      if (LIS_FORC_Forc_Hgt%selectOpt.eq.1) then
         call ESMF_FieldGet(hlmlField,localDE=0,farrayPtr=hlml,rc=status)
         call LIS_verify(status)
      endif

      mfactor = LIS_rc%nensem(n)/geosit_struc(n)%nIter

      do k = 1,LIS_rc%ntiles(n)/mfactor
         do m = 1,mfactor
            t = m + (k-1)*mfactor
            index1 = LIS_domain(n)%tile(t)%index
            kk = LIS_get_iteration_index(n,k,index1,mfactor)

            tair(t) = geosit_struc(n)%metdata1(kk,1,index1)
            qair(t) = geosit_struc(n)%metdata1(kk,2,index1)
            swgdn(t) = geosit_struc(n)%metdata1(kk,3,index1)
            if (swgdn(t).gt.LIS_CONST_SOLAR) then
               swgdn(t) = LIS_CONST_SOLAR
            endif
            lwgab(t) = geosit_struc(n)%metdata1(kk,4,index1)
            uwind(t) = geosit_struc(n)%metdata1(kk,5,index1)
            vwind(t) = geosit_struc(n)%metdata1(kk,6,index1)
            ps(t) = geosit_struc(n)%metdata1(kk,7,index1)
            prectot(t) = geosit_struc(n)%metdata1(kk,8,index1)
            if (prectot(t).lt.0.0) then
               prectot(t) = 0.0
            endif
            if (LIS_FORC_CRainf%selectOpt.eq.1) then
               preccon(t) = geosit_struc(n)%metdata1(kk,9,index1)
            endif
            if (LIS_FORC_Snowf%selectOpt.eq.1) then
               precsno(t) = geosit_struc(n)%metdata1(kk,10,index1)
            endif
            if (LIS_FORC_Swnet%selectOpt.eq.1) then
               swland(t) = geosit_struc(n)%metdata1(kk,11,index1)
            endif
            if (LIS_FORC_Pardr%selectOpt.eq.1) then
               pardr(t) = geosit_struc(n)%metdata1(kk,12,index1)
            endif
            if (LIS_FORC_Pardf%selectOpt.eq.1) then
               pardf(t) = geosit_struc(n)%metdata1(kk,13,index1)
            endif
            if (LIS_FORC_Forc_Hgt%selectOpt.eq.1) then
               hlml(t) = geosit_struc(n)%metdata1(kk,14,index1)
            endif
         enddo
      enddo

      end subroutine timeinterp_geosit

