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
! !ROUTINE: timeinterp_nldas20
! \label{timeinterp_nldas20}
!
! !REVISION HISTORY:
! 11 Jul 2024: David Mocko, Initial Specification
!                           (derived from timeinterp_nldas2.F90)
!
! !INTERFACE:
subroutine timeinterp_nldas20(n,findex)
! !USES:
  use ESMF
  use LIS_FORC_AttributesMod
  use LIS_coreMod,        only : LIS_rc,LIS_domain
  use LIS_constantsMod,   only : LIS_CONST_SOLAR
  use LIS_metforcingMod,  only : LIS_forc,LIS_FORC_Base_State
  use LIS_timeMgrMod,     only : LIS_tick,LIS_time2date
  use LIS_logMod,         only : LIS_logunit,LIS_verify,LIS_endrun
  use nldas20_forcingMod, only : nldas20_struc
  use LIS_forecastMod,    only : LIS_get_iteration_index

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Temporally interpolates the forcing data to the current model
!  timestep.  Downward shortwave radiation is interpolated using
!  a zenith-angled based approach.  Precipitation fields are not
!  temporally interpolated, as it is backward-accumulated in the
!  forcing files.  All other variables are linearly interpolated
!  between the hourly blocks.
!
!  The routines invoked are:
!  \begin{description}
!   \item[LIS\_time2date](\ref{LIS_time2date}) \newline
!    converts the time to a date format
!   \item[LIS\_tick](\ref{LIS_tick}) \newline
!    advances or retracts time by the specified amount
!   \item[zterp](\ref{zterp}) \newline
!    zenith-angle based interpolation
!  \end{description}
!
!EOP
  integer :: zdoy
  real    :: zw1,zw2
  real    :: czm,cze,czb
  real    :: wt1,wt2,swt1,swt2
  real    :: gmt1,gmt2,tempbts
  integer :: t,index1
  integer :: bdoy,byr,bmo,bda,bhr,bmn
  real*8  :: btime,newtime1,newtime2
  real    :: tempgmt1,tempgmt2
  integer :: tempbdoy,tempbyr,tempbmo,tempbda,tempbhr,tempbmn
  integer :: tempbss

  integer          :: status
  type(ESMF_Field) :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field) :: psurfField,pcpField,cpcpField,fhgtField,acondField
  type(ESMF_Field) :: PETField,CAPEField
  real, pointer    :: tmp(:),q2(:),uwind(:),vwind(:)
  real, pointer    :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)
  real, pointer    :: fheight(:),acond(:),pet(:),cape(:)
  logical          :: forcing_z,forcing_ch,forcing_pet,forcing_cape
  integer          :: mfactor,m,k,kk

!________________________________________

  btime = nldas20_struc(n)%nldas20time1
  call LIS_time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)

  tempbdoy = bdoy
  tempgmt1 = gmt1
  tempbyr = byr
  tempbmo = bmo
  tempbda = bda
  tempbhr = bhr
  if (tempbhr.eq.24) tempbhr = 0
  tempbmn = bmn
  tempbss = 0
  tempbts = 0
  call LIS_tick(newtime1,tempbdoy,tempgmt1,tempbyr,tempbmo,tempbda, &
       tempbhr,tempbmn,tempbss,tempbts)

  btime = nldas20_struc(n)%nldas20time2
  call LIS_time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)
  tempbdoy = bdoy
  tempgmt2 = gmt2
  tempbyr = byr
  tempbmo = bmo
  tempbda = bda
  tempbhr = bhr
  if (tempbhr.eq.24) tempbhr = 0
  tempbmn = bmn
  tempbss = 0
  tempbts = 0
  call LIS_tick(newtime2,tempbdoy,tempgmt2,tempbyr,tempbmo,tempbda, &
       tempbhr,tempbmn,tempbss,tempbts)

!=== Interpolate Data in time
  wt1 = (nldas20_struc(n)%nldas20time2 - LIS_rc%time) /            &
       (nldas20_struc(n)%nldas20time2 - nldas20_struc(n)%nldas20time1)
  wt2 = 1.0 - wt1
  swt1 = (newtime2 - LIS_rc%time) / (newtime2 - newtime1)
  swt2 = 1.0 - swt1

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),                &
       LIS_FORC_Tair%varname(1),tmpField,rc=status)
  call LIS_verify(status,                                          &
       "[ERR] Enable Tair in the forcing variables list")

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),                &
       LIS_FORC_Qair%varname(1),q2Field,rc=status)
  call LIS_verify(status,                                          &
       "[ERR] Enable Qair in the forcing variables list")

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),                &
       LIS_FORC_SWdown%varname(1),swdField,rc=status)
  call LIS_verify(status,                                          &
       "[ERR] Enable SWdown in the forcing variables list")

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),                &
       LIS_FORC_LWdown%varname(1),lwdField,rc=status)
  call LIS_verify(status,                                          &
       "[ERR] Enable LWdown in the forcing variables list")

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),                &
       LIS_FORC_Wind_E%varname(1),uField,rc=status)
  call LIS_verify(status,                                          &
       "[ERR] Enable Wind_E in the forcing variables list")

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),                &
       LIS_FORC_Wind_N%varname(1),vField,rc=status)
  call LIS_verify(status,                                          &
       "[ERR] Enable Wind_N in the forcing variables list")

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),                &
       LIS_FORC_Psurf%varname(1),psurfField,rc=status)
  call LIS_verify(status,                                          &
       "[ERR] Enable Psurf in the forcing variables list")

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),                &
       LIS_FORC_Rainf%varname(1),pcpField,rc=status)
  call LIS_verify(status,                                          &
       "[ERR] Enable Rainf in the forcing variables list")

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),                &
       LIS_FORC_CRainf%varname(1),cpcpField,rc=status)
  call LIS_verify(status,                                          &
       "[ERR] Enable CRainf in the forcing variables list")

  forcing_z = .false.
  forcing_ch = .false.
  if (nldas20_struc(n)%model_level_data.eq.1) then
     if (LIS_FORC_Forc_Hgt%selectOpt.eq.1) then
        call ESMF_StateGet(LIS_FORC_Base_State(n,findex),          &
             LIS_FORC_Forc_Hgt%varname(1),fhgtField,rc=status)
        call LIS_verify(status,                                    &
             "[ERR] Enable Forc_Hgt in the forcing variables list")
        forcing_z = .true.
     else
        write(LIS_logunit,*)                                       &
             "[ERR] Trying to use the NLDAS-2 FORB forcing data"
        write(LIS_logunit,*)                                       &
             "[ERR] without turning on Forc_Hgt variable in the"
        write(LIS_logunit,*)                                       &
             "[ERR] Forcing variable list file.  Please turn the"
        write(LIS_logunit,*)                                       &
             "[ERR] Forc_Hgt variable on when using FORB data."
        write(LIS_logunit,*) "[ERR] Program stopping ..."
        call LIS_endrun()
     endif

     if (LIS_FORC_Ch%selectOpt.eq.1) then
        call ESMF_StateGet(LIS_FORC_Base_State(n,findex),          &
             LIS_FORC_Ch%varname(1),acondField,rc=status)
        call LIS_verify(status,                                    &
             "[ERR] Enable Ch in the forcing variables list")
        forcing_ch = .true.
     else
        write(LIS_logunit,*)                                       &
             "[ERR] Trying to use the NLDAS-2 FORB forcing data"
        write(LIS_logunit,*)                                       &
             "[ERR] without turning on Ch variable in the"
        write(LIS_logunit,*)                                       &
             "[ERR] Forcing variable list file.  Please turn the"
        write(LIS_logunit,*)                                       &
             "[ERR] Ch variable on when using FORB data."
        write(LIS_logunit,*) "[ERR] Program stopping ..."
        call LIS_endrun()
     endif
  endif

  if (LIS_FORC_PET%selectOpt.eq.1) then
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),             &
          LIS_FORC_PET%varname(1),PETField,rc=status)
     call LIS_verify(status,                                       &
          "[ERR] Enable PET in the forcing variables list")
     forcing_pet = .true.
  else
     forcing_pet = .false.
  endif

  if (LIS_FORC_CAPE%selectOpt.eq.1) then
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),             &
          LIS_FORC_CAPE%varname(1),CAPEField,rc=status)
     call LIS_verify(status,                                       &
          "[ERR] Enable CAPE in the forcing variables list")
     forcing_cape = .true.
  else
     forcing_cape = .false.
  endif

  call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
  call LIS_verify(status)

! Loop over number of forcing ensembles:
  mfactor = LIS_rc%nensem(n) / nldas20_struc(n)%nIter

  do k = 1,(LIS_rc%ntiles(n)/mfactor)
     do m = 1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        zdoy = LIS_rc%doy
        ! Compute and apply zenith angle weights
        call zterp(1,LIS_domain(n)%grid(index1)%lat,               &
             LIS_domain(n)%grid(index1)%lon,                       &
             gmt1,gmt2,LIS_rc%gmt,zdoy,zw1,zw2,                    &
             czb,cze,czm,LIS_rc)

        kk = LIS_get_iteration_index(n,k,index1,mfactor)

        if ((nldas20_struc(n)%metdata1(kk,3,index1).ne.LIS_rc%udef).and. &
             (nldas20_struc(n)%metdata2(kk,3,index1).ne.LIS_rc%udef)) then
           swd(t) = (nldas20_struc(n)%metdata1(kk,3,index1)*zw1) + &
                (nldas20_struc(n)%metdata2(kk,3,index1)*zw2)
           ! In cases of small cos(zenith) angles, use linear weighting
           ! to avoid overly large weights.
           if ((swd(t).gt.nldas20_struc(n)%metdata1(kk,3,index1)).and. &
                (swd(t).gt.nldas20_struc(n)%metdata2(kk,3,index1)).and. &
                ((czb.lt.0.1).or.(cze.lt.0.1))) then
              swd(t) = (nldas20_struc(n)%metdata1(kk,3,index1)*swt1) + &
                   (nldas20_struc(n)%metdata2(kk,3,index1)*swt2)
           endif
        endif

        if (swd(t).gt.LIS_CONST_SOLAR) then
           write(unit=LIS_logunit,fmt=*)                           &
                "[WARN] sw radiation too high!!"
           write(unit=LIS_logunit,fmt=*)                           &
                "[WARN] it is: ",swd(t)
           write(unit=LIS_logunit,fmt=*)                           &
                "[WARN] data1 =",nldas20_struc(n)%metdata1(kk,3,index1)
           write(unit=LIS_logunit,fmt=*)                           &
                "[WARN] data2 =",nldas20_struc(n)%metdata2(kk,3,index1)
           write(unit=LIS_logunit,fmt=*)                           &
                "[WARN] zw1 =",zw1,"zw2 =",zw2
           write(unit=LIS_logunit,fmt=*)                           &
                "[WARN] swt1 =",swt1,"swt2 =",swt2
        endif
     enddo
  enddo                     ! End for SWdown

! Do block precipitation interpolation
  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LIS_verify(status)

  do k = 1,(LIS_rc%ntiles(n)/mfactor)
     do m = 1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        kk = LIS_get_iteration_index(n,k,index1,mfactor)
        if (nldas20_struc(n)%metdata2(kk,8,index1).ne.LIS_rc%udef) then
           pcp(t) = nldas20_struc(n)%metdata2(kk,8,index1)
           pcp(t) = pcp(t) / (60.0*60.0)
        endif
     enddo
  enddo

  call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
  call LIS_verify(status)

! Input is actually convective precip fraction; Calc actual CPCP below
  do k = 1,(LIS_rc%ntiles(n)/mfactor)
     do m = 1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        kk = LIS_get_iteration_index(n,k,index1,mfactor)
        if (nldas20_struc(n)%metdata2(kk,9,index1).ne.LIS_rc%udef) then
           cpcp(t) = nldas20_struc(n)%metdata2(kk,9,index1)
           if (nldas20_struc(n)%model_pcp_data.gt.0) then
              cpcp(t) = cpcp(t) / (60.0*60.0)
           else
              cpcp(t) = cpcp(t) * pcp(t)
           endif
        endif
     enddo
  enddo

! Linearly interpolate everything else
  call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
  call LIS_verify(status)

  do k = 1,(LIS_rc%ntiles(n)/mfactor)
     do m = 1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        kk = LIS_get_iteration_index(n,k,index1,mfactor)
        if ((nldas20_struc(n)%metdata1(kk,1,index1).ne.LIS_rc%udef).and. &
             (nldas20_struc(n)%metdata2(kk,1,index1).ne.LIS_rc%udef)) then
           tmp(t) = (nldas20_struc(n)%metdata1(kk,1,index1)*wt1) + &
                (nldas20_struc(n)%metdata2(kk,1,index1)*wt2)
        endif
     enddo
  enddo

  call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
  call LIS_verify(status)

  do k = 1,(LIS_rc%ntiles(n)/mfactor)
     do m = 1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        kk = LIS_get_iteration_index(n,k,index1,mfactor)
        if ((nldas20_struc(n)%metdata1(kk,2,index1).ne.LIS_rc%udef).and. &
             (nldas20_struc(n)%metdata2(kk,2,index1).ne.LIS_rc%udef)) then
           q2(t) = (nldas20_struc(n)%metdata1(kk,2,index1)*wt1) + &
                (nldas20_struc(n)%metdata2(kk,2,index1)*wt2)
        endif
     enddo
  enddo

  call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
  call LIS_verify(status)

  do k = 1,(LIS_rc%ntiles(n)/mfactor)
     do m = 1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        kk = LIS_get_iteration_index(n,k,index1,mfactor)
        if ((nldas20_struc(n)%metdata1(kk,4,index1).ne.LIS_rc%udef).and. &
             (nldas20_struc(n)%metdata2(kk,4,index1).ne.LIS_rc%udef)) then
           lwd(t) = (nldas20_struc(n)%metdata1(kk,4,index1)*wt1) + &
                (nldas20_struc(n)%metdata2(kk,4,index1)*wt2)
        endif
     enddo
  enddo

  call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
  call LIS_verify(status)

  do k = 1,(LIS_rc%ntiles(n)/mfactor)
     do m = 1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        kk = LIS_get_iteration_index(n,k,index1,mfactor)
        if ((nldas20_struc(n)%metdata1(kk,5,index1).ne.LIS_rc%udef).and. &
             (nldas20_struc(n)%metdata2(kk,5,index1).ne.LIS_rc%udef)) then
           uwind(t) = (nldas20_struc(n)%metdata1(kk,5,index1)*wt1) + &
                (nldas20_struc(n)%metdata2(kk,5,index1)*wt2)
        endif
     enddo
  enddo

  call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
  call LIS_verify(status)

  do k = 1,(LIS_rc%ntiles(n)/mfactor)
     do m = 1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        kk = LIS_get_iteration_index(n,k,index1,mfactor)
        if ((nldas20_struc(n)%metdata1(kk,6,index1).ne.LIS_rc%udef).and. &
             (nldas20_struc(n)%metdata2(kk,6,index1).ne.LIS_rc%udef)) then
           vwind(t) = (nldas20_struc(n)%metdata1(kk,6,index1)*wt1) + &
                (nldas20_struc(n)%metdata2(kk,6,index1)*wt2)
        endif
     enddo
  enddo

  call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
  call LIS_verify(status)

  do k = 1,(LIS_rc%ntiles(n)/mfactor)
     do m = 1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        kk = LIS_get_iteration_index(n,k,index1,mfactor)
        if ((nldas20_struc(n)%metdata1(kk,7,index1).ne.LIS_rc%udef).and. &
             (nldas20_struc(n)%metdata2(kk,7,index1).ne.LIS_rc%udef)) then
           psurf(t) = (nldas20_struc(n)%metdata1(kk,7,index1)*wt1) + &
                (nldas20_struc(n)%metdata2(kk,7,index1)*wt2)
        endif
     enddo
  enddo

  if (forcing_PET) then
     call ESMF_FieldGet(PETField,localDE=0,farrayPtr=pet,rc=status)
     call LIS_verify(status)

     do k = 1,LIS_rc%ntiles(n)/mfactor
        do m = 1,mfactor
           t = m + (k-1)*mfactor
           index1 = LIS_domain(n)%tile(t)%index
           kk = LIS_get_iteration_index(n,k,index1,mfactor)
           if ((nldas20_struc(n)%metdata1(kk,10,index1).ne.LIS_rc%udef).and. &
                (nldas20_struc(n)%metdata2(kk,10,index1).ne.LIS_rc%udef)) then
              pet(t) = (nldas20_struc(n)%metdata1(kk,10,index1)*wt1) + &
                   (nldas20_struc(n)%metdata2(kk,10,index1)*wt2)
              ! Convert NLDAS-2 PET from kg/m^2 to kg/m^2/sec - dmm
              pet(t) = pet(t)/(60.0*60.0)
           endif
        enddo
     enddo
  endif

  if (forcing_CAPE) then
     call ESMF_FieldGet(CAPEField,localDE=0,farrayPtr=cape,rc=status)
     call LIS_verify(status)

     do k = 1,LIS_rc%ntiles(n)/mfactor
        do m = 1,mfactor
           t = m + (k-1)*mfactor
           index1 = LIS_domain(n)%tile(t)%index
           kk = LIS_get_iteration_index(n,k,index1,mfactor)
           if ((nldas20_struc(n)%metdata1(kk,11,index1).ne.LIS_rc%udef) .and. &
                (nldas20_struc(n)%metdata2(kk,11,index1).ne.LIS_rc%udef)) then
              cape(t) = (nldas20_struc(n)%metdata1(kk,11,index1)*wt1) + &
                   (nldas20_struc(n)%metdata2(kk,11,index1)*wt2)
           endif
        enddo
     enddo
  endif

  if (forcing_z) then
     call ESMF_FieldGet(fhgtField,localDE=0,farrayPtr=fheight,rc=status)
     call LIS_verify(status)

     do k = 1,LIS_rc%ntiles(n)/mfactor
        do m = 1,mfactor
           t = m + (k-1)*mfactor
           index1 = LIS_domain(n)%tile(t)%index
           kk = LIS_get_iteration_index(n,k,index1,mfactor)
           if ((nldas20_struc(n)%metdata1(kk,12,index1).ne.LIS_rc%udef).and. &
                (nldas20_struc(n)%metdata2(kk,12,index1).ne.LIS_rc%udef)) then
              fheight(t) = (nldas20_struc(n)%metdata1(kk,12,index1)*wt1) + &
                   (nldas20_struc(n)%metdata2(kk,12,index1)*wt2)
           endif
        enddo
     enddo
  endif

  if (forcing_ch) then
     call ESMF_FieldGet(acondField,localDE=0,farrayPtr=acond,rc=status)
     call LIS_verify(status)

     do k = 1,LIS_rc%ntiles(n)/mfactor
        do m = 1,mfactor
           t = m + (k-1)*mfactor
           index1 = LIS_domain(n)%tile(t)%index
           kk = LIS_get_iteration_index(n,k,index1,mfactor)
           if ((nldas20_struc(n)%metdata1(kk,13,index1).ne.LIS_rc%udef).and. &
                (nldas20_struc(n)%metdata2(kk,13,index1).ne.LIS_rc%udef)) then
              acond(t) = (nldas20_struc(n)%metdata1(kk,13,index1)*wt1) + &
                   (nldas20_struc(n)%metdata2(kk,13,index1)*wt2)
           endif
        enddo
     enddo
  endif

end subroutine timeinterp_nldas20

