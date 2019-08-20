!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: timeinterp_gdas3d
! \label{timeinterp_gdas3d}
!
! !REVISION HISTORY:
!
! !INTERFACE:
subroutine timeinterp_gdas3d(n, findex)
! !USES:
  use ESMF
  use LIS_coreMod,    only : LIS_rc, LIS_domain
  use LIS_logMod,     only : LIS_verify
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod, only : LIS_FORC_Base_State, LIS_forc
  use gdas3d_forcingMod,  only : gdas3d_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!
!EOP

  type(ESMF_Field)   :: tmpField, q2Field, lprsField,prsField,o3Field
  type(ESMF_Field)   :: uwindField, vwindField
 
  real,   pointer    :: tmp(:), q2(:), lprs(:), prs(:),o3(:),uwind(:),vwind(:)
  integer            :: t,k
  integer            :: index1
  integer            :: status
  integer            :: ivar
  real               :: wt1, wt2
  
  wt1 = (gdas3d_struc(n)%gdastime2-LIS_rc%time)/&
       (gdas3d_struc(n)%gdastime2-gdas3d_struc(n)%gdastime1)
  wt2 = 1-wt1

  ivar = 0
! The GDAS layers are counted from TOA to surface. We need to reverse 
! the indices



!   kwh 4/8/09: added '+1'
  do k=1,gdas3d_struc(n)%nlayer +1
     ivar = ivar+1

     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),&
          LIS_FORC_LPressure%varname(k),lprsField,&
          rc=status)
     call LIS_verify(status, 'Error: Enable LPressure in the forcing variables list')
     
     call ESMF_FieldGet(lprsField, localDE=0,farrayPtr=lprs, rc=status)
     call LIS_verify(status, 'Error: fieldget in timeinterp_gdas3d')

     
     do t=1,LIS_rc%ntiles(n)
        index1 = LIS_domain(n)%tile(t)%index
        if((gdas3d_struc(n)%metdata1(ivar,index1).ne.LIS_rc%udef).and.&
             (gdas3d_struc(n)%metdata2(ivar,index1).ne.LIS_rc%udef)) then 
           lprs(t) = gdas3d_struc(n)%metdata1(ivar,index1)*wt1+&
                gdas3d_struc(n)%metdata2(ivar,index1)*wt2
        endif
     enddo
  enddo

  do k=1,gdas3d_struc(n)%nlayer
     ivar = ivar+1

     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Psurf%varname(k),prsField,&
          rc=status)
     call LIS_verify(status, 'Error: Enable Psurf in the forcing variables list')
     
     call ESMF_FieldGet(prsField, localDE=0,farrayPtr=prs, rc=status)
     call LIS_verify(status, 'Error: fieldget in timeinterp_gdas3d')

     do t=1,LIS_rc%ntiles(n)
        index1 = LIS_domain(n)%tile(t)%index
        if((gdas3d_struc(n)%metdata1(ivar,index1).ne.LIS_rc%udef).and.&
             (gdas3d_struc(n)%metdata2(ivar,index1).ne.LIS_rc%udef)) then 
           prs(t) = (gdas3d_struc(n)%metdata1(ivar,index1)*wt1+&
                gdas3d_struc(n)%metdata2(ivar,index1)*wt2)*100.0
        endif
     enddo
  enddo

  do k=1,gdas3d_struc(n)%nlayer
     ivar = ivar+1
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Tair%varname(k),tmpField,&
          rc=status)
     call LIS_verify(status, 'Error: Enable Tair in the forcing variables list')
     
     call ESMF_FieldGet(tmpField, localDE=0,farrayPtr=tmp, rc=status)
     call LIS_verify(status, 'Error: fieldget in timeinterp_gdas3d')
     
     do t=1,LIS_rc%ntiles(n)
        index1 = LIS_domain(n)%tile(t)%index
        if((gdas3d_struc(n)%metdata1(ivar,index1).ne.LIS_rc%udef).and.&
             (gdas3d_struc(n)%metdata2(ivar,index1).ne.LIS_rc%udef)) then 
           tmp(t) = gdas3d_struc(n)%metdata1(ivar,index1)*wt1+&
                gdas3d_struc(n)%metdata2(ivar,index1)*wt2
        endif
     enddo
  enddo

  do k=1,gdas3d_struc(n)%nlayer
     ivar = ivar+1

     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Qair%varname(k),q2Field,&
          rc=status)
     call LIS_verify(status, 'Error: Enable Qair in the forcing variables list')
     
     call ESMF_FieldGet(q2Field, localDE=0,farrayPtr=q2, rc=status)
     call LIS_verify(status, 'Error: fieldget in timeinterp_gdas3d')

     do t=1,LIS_rc%ntiles(n)
        index1 = LIS_domain(n)%tile(t)%index
        if((gdas3d_struc(n)%metdata1(ivar,index1).ne.LIS_rc%udef).and.&
             (gdas3d_struc(n)%metdata2(ivar,index1).ne.LIS_rc%udef)) then 
           q2(t) = gdas3d_struc(n)%metdata1(ivar,index1)*wt1+&
                gdas3d_struc(n)%metdata2(ivar,index1)*wt2
        endif
     enddo
  end do

  do k=1,gdas3d_struc(n)%nlayer
     ivar = ivar+1
     
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_O3%varname(k),o3Field,&
          rc=status)
     call LIS_verify(status, 'Error: Enable O3 in the forcing variables list')
     
     call ESMF_FieldGet(o3Field, localDE=0,farrayPtr=o3, rc=status)
     call LIS_verify(status, 'Error: fieldget in timeinterp_gdas3d')
          
     do t=1,LIS_rc%ntiles(n)
        index1 = LIS_domain(n)%tile(t)%index
        if((gdas3d_struc(n)%metdata1(ivar,index1).ne.LIS_rc%udef).and.&
             (gdas3d_struc(n)%metdata2(ivar,index1).ne.LIS_rc%udef)) then 
           o3(t) = gdas3d_struc(n)%metdata1(ivar,index1)*wt1+&
                gdas3d_struc(n)%metdata2(ivar,index1)*wt2
        endif
     enddo
  enddo

  ivar = ivar+1
  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_N%varname(1),vwindField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Wind_N in the forcing variables list')
  
  call ESMF_FieldGet(vwindField, localDE=0,farrayPtr=vwind, rc=status)
  call LIS_verify(status, 'Error: fieldget in timeinterp_gdas3d')
  
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     if((gdas3d_struc(n)%metdata1(ivar,index1).ne.LIS_rc%udef).and.&
          (gdas3d_struc(n)%metdata2(ivar,index1).ne.LIS_rc%udef)) then 
        vwind(t) = gdas3d_struc(n)%metdata1(ivar,index1)*wt1+&
             gdas3d_struc(n)%metdata2(ivar,index1)*wt2
     endif
  enddo

  ivar = ivar+1
  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_E%varname(1),uwindField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Wind_E in the forcing variables list')
  
  call ESMF_FieldGet(uwindField, localDE=0,farrayPtr=uwind, rc=status)
  call LIS_verify(status, 'Error: fieldget in timeinterp_gdas3d')
  
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     if((gdas3d_struc(n)%metdata1(ivar,index1).ne.LIS_rc%udef).and.&
          (gdas3d_struc(n)%metdata2(ivar,index1).ne.LIS_rc%udef)) then 
        uwind(t) = gdas3d_struc(n)%metdata1(ivar,index1)*wt1+&
             gdas3d_struc(n)%metdata2(ivar,index1)*wt2
     endif
  enddo
end subroutine timeinterp_gdas3d
