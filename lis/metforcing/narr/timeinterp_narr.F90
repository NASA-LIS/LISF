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
! !ROUTINE: timeinterp_narr
! \label{timeinterp_narr}
!
! !REVISION HISTORY:
!  30 APR 2009: Sujay Kumar, Initial Specification
!
! !INTERFACE:
subroutine timeinterp_narr(n, findex)
! !USES:
  use ESMF
  use LIS_coreMod,    only : LIS_rc, LIS_domain
  use LIS_logMod,     only : LIS_verify
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod, only : LIS_FORC_Base_State, LIS_forc
  use narr_forcingMod,  only : narr_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex

!
! !DESCRIPTION:
!  This routine temporally interpolates the forcing variables from 
!  NARR data to the LIS timestep. 
!
!EOP

  type(ESMF_Field)   :: tmpField, q2Field, lprsField,prsField,swdownField
  type(ESMF_Field)   :: uwindField, vwindField, lwdownField, rainfField,crainfField
 
  real,   pointer    :: tmp(:), q2(:), lprs(:), prs(:),uwind(:),vwind(:)
  real,   pointer    :: swdown(:),lwdown(:),rainf(:), crainf(:)
  integer            :: t,k
  integer            :: index1
  integer            :: status
  integer            :: ivar
  real               :: wt1, wt2
  real               :: thislevelval, nextlevelval  
  wt1 = (narr_struc(n)%narrtime2-LIS_rc%time)/&
       (narr_struc(n)%narrtime2-narr_struc(n)%narrtime1)
  wt2 = 1-wt1

  ivar = 0

!TMP
  do k=1,narr_struc(n)%nlevels+1  ! '+1' because of surface level
     ivar = ivar+1
     if (k .eq. 1) then !surface level and no averaging
        call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Tair%varname(k),&
             tmpField,rc=status)
        call LIS_verify(status, 'Error: Enable Tair in the forcing variables list')
        call ESMF_FieldGet(tmpField, localDE=0,farrayPtr=tmp, rc=status)
        call LIS_verify(status, 'Error: tmp fieldget in timeinterp_narr')

        do t=1,LIS_rc%ntiles(n)
           index1 = LIS_domain(n)%tile(t)%index
           if((narr_struc(n)%metdata1(ivar,index1).ne.LIS_rc%udef)&
                .and.&
                (narr_struc(n)%metdata2(ivar,index1).ne.LIS_rc%udef)) then 
              tmp(t) = narr_struc(n)%metdata1(ivar,index1)*wt1+&
                   narr_struc(n)%metdata2(ivar,index1)*wt2
           endif
        enddo
     else if (k .lt. narr_struc(n)%nlevels+1) then
        call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Tair%varname(k),tmpField,&
             rc=status)
        call LIS_verify(status, 'Error: Enable Tair in the forcing variables list')
        call ESMF_FieldGet(tmpField, localDE=0,farrayPtr=tmp, rc=status)
        call LIS_verify(status, 'Error: tmp fieldget in timeinterp_narr')
        do t=1,LIS_rc%ntiles(n)
           index1 = LIS_domain(n)%tile(t)%index
           if(((narr_struc(n)%metdata1(ivar,index1).ne.LIS_rc%udef).and.&
                (narr_struc(n)%metdata2(ivar,index1).ne.LIS_rc%udef)) .and. &
                ((narr_struc(n)%metdata1(ivar+1,index1).ne.LIS_rc%udef).and.&
                (narr_struc(n)%metdata2(ivar+1,index1).ne.LIS_rc%udef)))  then 
              thislevelval = narr_struc(n)%metdata1(ivar,index1)*wt1+&
                   narr_struc(n)%metdata2(ivar,index1)*wt2
              nextlevelval = narr_struc(n)%metdata1(ivar+1,index1)*wt1+&
                   narr_struc(n)%metdata2(ivar+1,index1)*wt2
              tmp(t) = (thislevelval+nextlevelval)/2
           endif
        enddo
     else ! k .eq. n+1
           !index = uppermost pressure level
           !ignore; just advance ivar by one
     endif
  enddo

!Q2
  do k=1,narr_struc(n)%nlevels+1  ! '+1' because of surface level
     ivar = ivar+1
     if (k .eq. 1) then !surface level and no averaging
        call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Qair%varname(k),q2Field,&
             rc=status)
        call LIS_verify(status, 'Error: Enable Qair in the forcing variables list')
        call ESMF_FieldGet(q2Field, localDE=0,farrayPtr=q2, rc=status)
        call LIS_verify(status, 'Error: q2 fieldget in timeinterp_narr')
        do t=1,LIS_rc%ntiles(n)
           index1 = LIS_domain(n)%tile(t)%index
           if((narr_struc(n)%metdata1(ivar,index1).ne.LIS_rc%udef).and.&
                (narr_struc(n)%metdata2(ivar,index1).ne.LIS_rc%udef)) then 
              q2(t) = narr_struc(n)%metdata1(ivar,index1)*wt1+&
                   narr_struc(n)%metdata2(ivar,index1)*wt2
           endif
        enddo
     else if (k .lt. narr_struc(n)%nlevels+1) then
        call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Qair%varname(k),q2Field,&
             rc=status)
        call LIS_verify(status, 'Error: Enable Qair in the forcing variables list')
        call ESMF_FieldGet(q2Field, localDE=0,farrayPtr=q2, rc=status)
        call LIS_verify(status, 'Error: q2 fieldget in timeinterp_narr')
        do t=1,LIS_rc%ntiles(n)
           index1 = LIS_domain(n)%tile(t)%index
           if(((narr_struc(n)%metdata1(ivar,index1).ne.LIS_rc%udef).and.&
                (narr_struc(n)%metdata2(ivar,index1).ne.LIS_rc%udef)) .and. &
                ((narr_struc(n)%metdata1(ivar+1,index1).ne.LIS_rc%udef).and. &
                (narr_struc(n)%metdata2(ivar+1,index1).ne.LIS_rc%udef)))  then 
              thislevelval = narr_struc(n)%metdata1(ivar,index1)*wt1+&
                   narr_struc(n)%metdata2(ivar,index1)*wt2
              nextlevelval = narr_struc(n)%metdata1(ivar+1,index1)*wt1+&
                   narr_struc(n)%metdata2(ivar+1,index1)*wt2
              q2(t) = (thislevelval+nextlevelval)/2
           endif
        enddo
     else ! k .eq. n+1
           !index = uppermost pressure level
           !ignore; just advance ivar by one
     endif
  end do
!PRS and LPRS
!recall level pressure is stored in k = 2::30 of corresp. LIS_FORC, ie, no surf field
  do k=1,narr_struc(n)%nlevels+1
     ivar = ivar+1
     if (k .eq. 1) then
        !put surface pressure into (1) of Psurf
        !do nothing for Level pressure

        call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Psurf%varname(k),prsField,&
             rc=status)
        call LIS_verify(status, 'Error: Enable Psurf in the forcing variables list')
        call ESMF_FieldGet(prsField, localDE=0,farrayPtr=prs, rc=status)
        call LIS_verify(status, 'Error: prs fieldget in timeinterp_narr')
        
        do t=1,LIS_rc%ntiles(n)
           if((narr_struc(n)%metdata1(ivar,index1).ne.LIS_rc%udef).and.&
                (narr_struc(n)%metdata2(ivar,index1).ne.LIS_rc%udef)) then 
              prs(t) = narr_struc(n)%metdata1(ivar,index1)*wt1+&
                   narr_struc(n)%metdata2(ivar,index1)*wt2
           endif
        enddo
     else if (k .lt. narr_struc(n)%nlevels+1) then
         !put avg pressure (of k and k+1) into Psurf(k)
         !put pressure into LPressure(k)

        call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_LPressure%varname(k),lprsField,&
             rc=status)
        call LIS_verify(status, 'Error: Enable LPressure in the forcing variables list')
        call ESMF_FieldGet(lprsField, localDE=0,farrayPtr=lprs, rc=status)
        call LIS_verify(status, 'Error: lprs fieldget in timeinterp_narr')
        
        call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Psurf%varname(k),prsField,&
             rc=status)
        call LIS_verify(status, 'Error: Enable Psurf in the forcing variables list')
        call ESMF_FieldGet(prsField, localDE=0,farrayPtr=prs, rc=status)
        call LIS_verify(status, 'Error: prs fieldget in timeinterp_narr')
        do t=1,LIS_rc%ntiles(n)
           if(((narr_struc(n)%metdata1(ivar,index1).ne.LIS_rc%udef).and. &
                (narr_struc(n)%metdata2(ivar,index1).ne.LIS_rc%udef)) .and. &
                ((narr_struc(n)%metdata1(ivar+1,index1).ne.LIS_rc%udef).and. &
                (narr_struc(n)%metdata2(ivar+1,index1).ne.LIS_rc%udef))) then 
              thislevelval = narr_struc(n)%metdata1(ivar,index1)*wt1+&
                   narr_struc(n)%metdata2(ivar,index1)*wt2
              nextlevelval = narr_struc(n)%metdata1(ivar,index1)*wt1+&
                   narr_struc(n)%metdata2(ivar+1,index1)*wt2
              lprs(t) = thislevelval
              prs(t) = (thislevelval+nextlevelval)/2
           endif
        enddo
     else  !k = narr_struc(n)%nlevels+1
           !put pressure into LPressure(k)
        call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_LPressure%varname(k),lprsField,&
             rc=status)
        call LIS_verify(status, 'Error: Enable LPressure in the forcing variables list')
        call ESMF_FieldGet(lprsField, localDE=0,farrayPtr=lprs, rc=status)
        call LIS_verify(status, 'Error: lprs fieldget in timeinterp_narr')
        
        do t=1,LIS_rc%ntiles(n)
           index1 = LIS_domain(n)%tile(t)%index
           if((narr_struc(n)%metdata1(ivar,index1).ne.LIS_rc%udef).and.&
                (narr_struc(n)%metdata2(ivar,index1).ne.LIS_rc%udef)) then 
              lprs(t) = narr_struc(n)%metdata1(ivar,index1)*wt1+&
                   narr_struc(n)%metdata2(ivar,index1)*wt2
           endif
        enddo
     endif
  enddo
  

! replace this with zterp interpolation?
  ivar = ivar+1
  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_SWdown%varname(1),swdownField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable SWdown in the forcing variables list')
  
  call ESMF_FieldGet(swdownField, localDE=0,farrayPtr=swdown, rc=status)
  call LIS_verify(status, 'Error: swdown fieldget in timeinterp_narr')

  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     if((narr_struc(n)%metdata1(ivar,index1).ne.LIS_rc%udef).and.&
          (narr_struc(n)%metdata2(ivar,index1).ne.LIS_rc%udef)) then 
        swdown(t) = narr_struc(n)%metdata1(ivar,index1)*wt1+&
             narr_struc(n)%metdata2(ivar,index1)*wt2
     endif
  enddo
  
  ivar = ivar+1
  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Lwdown%varname(1),lwdownField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Lwdown in the forcing variables list')
  
  call ESMF_FieldGet(lwdownField, localDE=0,farrayPtr=lwdown, rc=status)
  call LIS_verify(status, 'Error: lwdown fieldget in timeinterp_narr')
  
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     if((narr_struc(n)%metdata1(ivar,index1).ne.LIS_rc%udef).and.&
          (narr_struc(n)%metdata2(ivar,index1).ne.LIS_rc%udef)) then 
        lwdown(t) = narr_struc(n)%metdata1(ivar,index1)*wt1+&
             narr_struc(n)%metdata2(ivar,index1)*wt2
     endif
  enddo
  
  ivar = ivar+1
  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_E%varname(1),uwindField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Wind_E in the forcing variables list')
  
  call ESMF_FieldGet(uwindField, localDE=0,farrayPtr=uwind, rc=status)
  call LIS_verify(status, 'Error: wind_e fieldget in timeinterp_narr')
  
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     if((narr_struc(n)%metdata1(ivar,index1).ne.LIS_rc%udef).and.&
          (narr_struc(n)%metdata2(ivar,index1).ne.LIS_rc%udef)) then 
        uwind(t) = narr_struc(n)%metdata1(ivar,index1)*wt1+&
             narr_struc(n)%metdata2(ivar,index1)*wt2
     endif
  enddo

  ivar = ivar+1
  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_N%varname(1),vwindField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Wind_N in the forcing variables list')
  
  call ESMF_FieldGet(vwindField, localDE=0,farrayPtr=vwind, rc=status)
  call LIS_verify(status, 'Error: wind_n fieldget in timeinterp_narr')

  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     if((narr_struc(n)%metdata1(ivar,index1).ne.LIS_rc%udef).and.&
          (narr_struc(n)%metdata2(ivar,index1).ne.LIS_rc%udef)) then 
        vwind(t) = narr_struc(n)%metdata1(ivar,index1)*wt1+&
             narr_struc(n)%metdata2(ivar,index1)*wt2
     endif
  enddo

  ivar = ivar+1
  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Rainf%varname(1),rainfField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Rainf in the forcing variables list')
  
  call ESMF_FieldGet(rainfField, localDE=0,farrayPtr=rainf, rc=status)
  call LIS_verify(status, 'Error: Rainf fieldget in timeinterp_narr')
  
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     if(narr_struc(n)%metdata2(ivar,index1).ne.LIS_rc%udef) then 
        rainf(t) = narr_struc(n)%metdata2(ivar,index1)
     endif
  enddo

  ivar = ivar+1
  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_CRainf%varname(1),crainfField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable CRainf in the forcing variables list')
  
  call ESMF_FieldGet(crainfField, localDE=0,farrayPtr=crainf, rc=status)
  call LIS_verify(status, 'Error: CRainf fieldget in timeinterp_narr')

  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     if(narr_struc(n)%metdata2(ivar,index1).ne.LIS_rc%udef) then 
        crainf(t) = narr_struc(n)%metdata2(ivar,index1)
     endif
  enddo
  
end subroutine timeinterp_narr
