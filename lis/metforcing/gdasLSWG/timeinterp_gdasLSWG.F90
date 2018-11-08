!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: timeinterp_gdasLSWG
! \label{timeinterp_gdasLSWG}
!
! !REVISION HISTORY:
!  20 Oct 2009: Sujay Kumar, Initial Specification
!
! !INTERFACE:
subroutine timeinterp_gdasLSWG(n, findex)
! !USES:
  use ESMF
  use LIS_coreMod,    only : LIS_rc, LIS_domain
  use LIS_logMod,     only : LIS_verify
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod, only : LIS_FORC_Base_State
  use gdasLSWG_forcingMod,  only : gdasLSWG_struc
#if (defined RTMS)
  USE units_conversion
  USE Type_Kinds, ONLY: fp
#endif
  
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
  real                :: metdata1(LIS_rc%met_nf(findex), LIS_rc%ngrid(n))
  real                :: metdata2(LIS_rc%met_nf(findex), LIS_rc%ngrid(n))  
!
! !DESCRIPTION:
!
!   This routine performs the temporal interpolation of the GDAS data to the
!   LIS model timestep
!EOP

  type(ESMF_Field)   :: tmpField, q2Field, lprsField,prsField
 
  real,   pointer    :: tmp(:), q2(:), lprs(:), prs(:)
  integer            :: t,k
  integer            :: index1
  integer            :: status
  integer            :: ivar
  integer            :: tid,tid1,tid2
  real               :: wt1, wt2
  real               :: thislevelval, nextlevelval  
  real, parameter    :: qsmall=1.0e-6

  wt1 = (gdasLSWG_struc(n)%time2-LIS_rc%time)/&
       (gdasLSWG_struc(n)%time2-gdasLSWG_struc(n)%time1)
  wt2 = 1-wt1

!TMP
  ivar = 0
  do k=2,26  ! '+1' because of surface level
     ivar = ivar+1
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Tair%varname(k),tmpField,&
          rc=status)
     call LIS_verify(status, 'Error: Enable Tair in the forcing variables list')
     call ESMF_FieldGet(tmpField, localDE=0,farrayPtr=tmp, rc=status)
     call LIS_verify(status, 'Error: tmp fieldget in timeinterp_gdasLSWG')
     do t=1,LIS_rc%ntiles(n)
        index1 = LIS_domain(n)%tile(t)%index
        if(((metdata1(ivar,index1).ne.LIS_rc%udef).and.&
             (metdata2(ivar,index1).ne.LIS_rc%udef)) .and. &
             ((metdata1(ivar+1,index1).ne.LIS_rc%udef).and.&
             (metdata2(ivar+1,index1).ne.LIS_rc%udef)))  then 
           thislevelval = metdata1(ivar,index1)*wt1+&
                metdata2(ivar,index1)*wt2
           nextlevelval = metdata1(ivar+1,index1)*wt1+&
                metdata2(ivar+1,index1)*wt2
           tmp(t) = (thislevelval+nextlevelval)/2
        endif
     enddo
  enddo

!PRS and LPRS
!recall level pressure is stored in k = 2::27 of corresp. LIS_FORC, ie, no surf field
  ivar = 52
  do k=2,26+1
     ivar = ivar+1
     
     if (k .lt. 26+1) then
        !put avg pressure (of k and k+1) into Psurf(k)
        !put pressure into LPressure(k)
        
        call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_LPressure%varname(k),lprsField,&
             rc=status)
        call LIS_verify(status, 'Error: Enable LPressure in the forcing variables list')
        call ESMF_FieldGet(lprsField, localDE=0,farrayPtr=lprs, rc=status)
        call LIS_verify(status, 'Error: lprs fieldget in timeinterp_gdasLSWG')
        
        call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Psurf%varname(k),prsField,&
             rc=status)
        call LIS_verify(status, 'Error: Enable Psurf in the forcing variables list')
        call ESMF_FieldGet(prsField, localDE=0,farrayPtr=prs, rc=status)
        call LIS_verify(status, 'Error: prs fieldget in timeinterp_gdasLSWG')
        do t=1,LIS_rc%ntiles(n)
           if(((metdata1(ivar,index1).ne.LIS_rc%udef).and. &
                (metdata2(ivar,index1).ne.LIS_rc%udef)) .and. &
                ((metdata1(ivar+1,index1).ne.LIS_rc%udef).and. &
                (metdata2(ivar+1,index1).ne.LIS_rc%udef))) then 
              thislevelval = metdata1(ivar,index1)*wt1+&
                   metdata2(ivar,index1)*wt2
              nextlevelval = metdata1(ivar,index1)*wt1+&
                   metdata2(ivar+1,index1)*wt2
              lprs(t) = thislevelval
              prs(t) = (thislevelval+nextlevelval)/2
           endif
        enddo
     else  !k = gdasLSWG_struc(n)%nlevels+1
           !put pressure into LPressure(k)

        call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_LPressure%varname(k),lprsField,&
             rc=status)
        call LIS_verify(status, 'Error: Enable LPressure in the forcing variables list')
        call ESMF_FieldGet(lprsField, localDE=0,farrayPtr=lprs, rc=status)
        call LIS_verify(status, 'Error: lprs fieldget in timeinterp_gdasLSWG')
        
        do t=1,LIS_rc%ntiles(n)
           index1 = LIS_domain(n)%tile(t)%index
           if((metdata1(ivar,index1).ne.LIS_rc%udef).and.&
                (metdata2(ivar,index1).ne.LIS_rc%udef)) then 
              lprs(t) = metdata1(ivar,index1)*wt1+&
                   metdata2(ivar,index1)*wt2
           endif
        enddo
     endif
  enddo

!Q2
  ivar = 26  !order reversed because rel humidity conversion first requires temp and press to compute
  do k=2,26  
     ivar = ivar+1
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Qair%varname(k),q2Field,&
          rc=status)
     call LIS_verify(status, 'Error: Enable Qair in the forcing variables list')
     call ESMF_FieldGet(q2Field, localDE=0,farrayPtr=q2, rc=status)
     call LIS_verify(status, 'Error: q2 fieldget in timeinterp_gdasLSWG')

     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Tair%varname(k),tmpField,&
          rc=status)
     call LIS_verify(status, 'Error: Enable Tair in the forcing variables list')
     call ESMF_FieldGet(tmpField, localDE=0,farrayPtr=tmp, rc=status)
     call LIS_verify(status, 'Error: tmp fieldget in timeinterp_gdasLSWG')

     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Psurf%varname(k),prsField,&
          rc=status)
     call LIS_verify(status, 'Error: Enable Psurf in the forcing variables list')
     call ESMF_FieldGet(prsField, localDE=0,farrayPtr=prs, rc=status)
     call LIS_verify(status, 'Error: prs fieldget in timeinterp_gdasLSWG')

     do t=1,LIS_rc%ntiles(n)
        index1 = LIS_domain(n)%tile(t)%index
        if(((metdata1(ivar,index1).ne.LIS_rc%udef).and.&
             (metdata2(ivar,index1).ne.LIS_rc%udef)) .and. &
             ((metdata1(ivar+1,index1).ne.LIS_rc%udef).and. &
             (metdata2(ivar+1,index1).ne.LIS_rc%udef)))  then 
           thislevelval = metdata1(ivar,index1)*wt1+&
                metdata2(ivar,index1)*wt2
           nextlevelval = metdata1(ivar+1,index1)*wt1+&
                metdata2(ivar+1,index1)*wt2
           q2(t) = (thislevelval+nextlevelval)/2
           
           !Q2--change Relative Humidity (RH) into LIS units, specific humidity (kg/kg)
           !use CRTM Profile Utility conversion to go from RH(%) to MR(g/kg)
           ! and then MR(g/kg) to SH(g/kg) and then mult by (1kg/1000g) to get to LIS unit
           !conversion assumes hPa and not Pa; hence prs division by 100
           q2(t) = max(q2(t),qsmall)
#if(defined RTMS)
           q2(t) = RH_to_MR( prs(t)/100.0_fp,tmp(t)*1.0_fp,q2(t)*1.0_fp)
           q2(t) = MR_to_SA( q2(t)*1.0_fp)/1000.0
#endif
        endif
     enddo
  enddo
end subroutine timeinterp_gdasLSWG
