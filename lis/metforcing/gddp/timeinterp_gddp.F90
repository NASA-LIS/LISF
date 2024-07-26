!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: timeinterp_gddp
! \label{timeinterp_gddp}
!
! !REVISION HISTORY:
!
! 03 Feb 2022: Sujay Kumar, Initial specification
!
! !INTERFACE:
subroutine timeinterp_gddp(n,findex)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_constantsMod
  use LIS_metforcingMod
  use LIS_FORC_AttributesMod
  use LIS_timeMgrMod
  use LIS_logMod
  use gddp_forcingMod
  use LIS_forecastMod
 
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Temporally interpolates the forcing data to the current model 
!  timestep. Downward shortwave radiation is interpolated using a
!  zenith-angled based approach. Precipitation and longwave radiation
!  are not temporally interpolated, and the previous 3 hourly value
!  is used. All other variables are linearly interpolated between 
!  the 3 hourly blocks. 
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
!EOP
  integer :: zdoy
  real    :: zw1, zw2
  real    :: czm, cze, czb
  real    :: wt1, wt2,swt1,swt2
  real    :: gmt1, gmt2, tempbts
  integer :: t,index1
  integer :: bdoy,byr,bmo,bda,bhr,bmn
  real*8  :: btime,newtime1,newtime2
  real    :: tempgmt1,tempgmt2
  integer :: tempbdoy,tempbyr,tempbmo,tempbda,tempbhr,tempbmn
  integer :: tempbss
  integer            :: status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField,fhgtField,acondField
  type(ESMF_Field)   :: PETField,CAPEField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)
  real,pointer       :: fheight(:),acond(:),pet(:),cape(:)
  logical            :: forcing_z, forcing_ch, forcing_pet, forcing_cape
  integer            :: mfactor, m, k, kk
! ________________________________________

  btime=gddp_struc(n)%gddptime1
  call LIS_time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)
  
  tempbdoy=bdoy
  tempgmt1=gmt1
  tempbyr=byr
  tempbmo=bmo
  tempbda=bda
  tempbhr=bhr
  if (tempbhr.eq.24) tempbhr=0
  tempbmn=bmn
  tempbss=0
  tempbts=0
  call LIS_tick(newtime1,tempbdoy,tempgmt1,& 
       tempbyr,tempbmo,tempbda,tempbhr,tempbmn, & 
       tempbss,tempbts)
  
  btime=gddp_struc(n)%gddptime2
  call LIS_time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)
  tempbdoy=bdoy
  tempgmt2=gmt2
  tempbyr=byr
  tempbmo=bmo
  tempbda=bda
  tempbhr=bhr
  if (tempbhr.eq.24) tempbhr=0
  tempbmn=bmn
  tempbss=0
  tempbts=0
  call LIS_tick(newtime2,tempbdoy,tempgmt2,&
       tempbyr,tempbmo,tempbda,tempbhr,tempbmn,&
       tempbss,tempbts)
  
!=== Interpolate Data in time      
  wt1=(gddp_struc(n)%gddptime2-LIS_rc%time)/ & 
       (gddp_struc(n)%gddptime2-gddp_struc(n)%gddptime1)


  wt2=1.0-wt1
  swt1=(newtime2-LIS_rc%time)/(newtime2-newtime1)
  swt2=1.0-swt1


  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Tair%varname(1),tmpField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Tair in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Qair%varname(1),q2Field,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Qair in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_SWdown%varname(1),swdField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable SWdown in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_LWdown%varname(1),lwdField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable LWdown in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_E%varname(1),uField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Wind_E in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_N%varname(1),vField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Wind_N in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Psurf%varname(1),psurfField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Psurf in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Rainf%varname(1),pcpField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Rainf in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_CRainf%varname(1),cpcpField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable CRainf in the forcing variables list')

  if(LIS_FORC_Forc_Hgt%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Forc_Hgt%varname(1),&
          fhgtField, rc=status)
     call LIS_verify(status, 'Error: Enable Forc_Hgt in the forcing variables list')
     forcing_z = .true.
  else
     forcing_z = .false.
  endif

  if(LIS_FORC_Ch%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Ch%varname(1),acondField,&
          rc=status)
     call LIS_verify(status, 'Error: Enable Ch in the forcing variables list')
     forcing_ch = .true.
  else
     forcing_ch = .false.
  endif

  if(LIS_FORC_PET%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_PET%varname(1),PETField,&
          rc=status)
     call LIS_verify(status, 'Error: Enable PET in the forcing variables list')
     forcing_pet = .true.
  else
     forcing_pet = .false.
  endif

  if(LIS_FORC_CAPE%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_CAPE%varname(1),CAPEField,&
          rc=status)
     call LIS_verify(status, 'Error: Enable CAPE in the forcing variables list')
     forcing_cape = .true.
  else
     forcing_cape = .false.
  endif

  call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
  call LIS_verify(status)

! Loop over number of forcing ensembles:
  mfactor = LIS_rc%nensem(n)

  do k=1,LIS_rc%ntiles(n)/mfactor
     do m=1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index


        zdoy=LIS_rc%doy
        
           ! compute and apply zenith angle weights
        call zterp(1,LIS_domain(n)%grid(index1)%lat,&
             LIS_domain(n)%grid(index1)%lon,&
             gmt1,gmt2,LIS_rc%gmt,zdoy,zw1,zw2,czb,cze,czm,LIS_rc)
        
        if(gddp_struc(n)%metdata1(m,3,index1).ne.LIS_rc%udef.and.&
             gddp_struc(n)%metdata2(m,3,index1).ne.LIS_rc%udef) then 
           swd(t) = gddp_struc(n)%metdata1(m,3,index1)*zw1+&
                gddp_struc(n)%metdata2(m,3,index1)*zw2
          ! In cases of small cos(zenith) angles, use linear weighting
          !  to avoid overly large weights

           if((swd(t).gt.gddp_struc(n)%metdata1(m,3,index1).and. & 
                swd(t).gt.gddp_struc(n)%metdata2(m,3,index1)).and. & 
                (czb.lt.0.1.or.cze.lt.0.1))then
              swd(t) = gddp_struc(n)%metdata1(m,3,index1)*swt1+ & 
                   gddp_struc(n)%metdata2(m,3,index1)*swt2
           endif
        endif
        
        if(swd(t).lt.0.0) then
           write(unit=LIS_logunit,fmt=*)'[ERR] sw radiation is unphysical'
           write(unit=LIS_logunit,fmt=*)'[ERR] it is', LIS_localPet, t, swd(t)
           write(unit=LIS_logunit,fmt=*)'[ERR] data1=',gddp_struc(n)%metdata1(m,3,index1)
           write(unit=LIS_logunit,fmt=*)'[ERR] data2=',gddp_struc(n)%metdata2(m,3,index1)
           write(unit=LIS_logunit,fmt=*)'[ERR] wts=',wt1,wt2
           call LIS_endrun()

        endif
     enddo  ! End for SWdown
  enddo
  
  ! do block precipitation interpolation
  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LIS_verify(status)

  do k=1,LIS_rc%ntiles(n)/mfactor
     do m=1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        if(gddp_struc(n)%metdata2(m,7,index1).ne.LIS_rc%udef) then 
           pcp(t)=gddp_struc(n)%metdata2(m,7,index1)
           if(pcp(t).lt.0) then
              pcp(t) = 0.0
           endif
        endif
     enddo
  enddo

  !linearly interpolate everything else

  call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
  call LIS_verify(status)
  
  do k=1,LIS_rc%ntiles(n)/mfactor
     do m=1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        if((gddp_struc(n)%metdata1(m,1,index1).ne.LIS_rc%udef).and.&
             (gddp_struc(n)%metdata2(m,1,index1).ne.LIS_rc%udef)) then 
           tmp(t) = gddp_struc(n)%metdata1(m,1,index1)*wt1+ & 
                gddp_struc(n)%metdata2(m,1,index1)*wt2
        endif
     enddo
  enddo
  
  call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
  call LIS_verify(status)

  do k=1,LIS_rc%ntiles(n)/mfactor
     do m=1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        if((gddp_struc(n)%metdata1(m,2,index1).ne.LIS_rc%udef).and.&
             (gddp_struc(n)%metdata2(m,2,index1).ne.LIS_rc%udef)) then 
           q2(t) =gddp_struc(n)%metdata1(m,2,index1)*wt1+ & 
                gddp_struc(n)%metdata2(m,2,index1)*wt2
        endif
     enddo
  enddo

  call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
  call LIS_verify(status)

  do k=1,LIS_rc%ntiles(n)/mfactor
     do m=1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        if((gddp_struc(n)%metdata1(m,4,index1).ne.LIS_rc%udef).and.&
             (gddp_struc(n)%metdata2(m,4,index1).ne.LIS_rc%udef)) then 
           lwd(t) =gddp_struc(n)%metdata1(m,4,index1)*wt1+ & 
                gddp_struc(n)%metdata2(m,4,index1)*wt2
        endif
     enddo
  enddo

  call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
  call LIS_verify(status)

  do k=1,LIS_rc%ntiles(n)/mfactor
     do m=1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index

        if((gddp_struc(n)%metdata1(m,5,index1).ne.LIS_rc%udef).and.&
             (gddp_struc(n)%metdata2(m,5,index1).ne.LIS_rc%udef)) then 
           uwind(t) = gddp_struc(n)%metdata1(m,5,index1)*wt1+ & 
                gddp_struc(n)%metdata2(m,5,index1)*wt2
        endif
     enddo
  enddo
  
  call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
  call LIS_verify(status)

  do k=1,LIS_rc%ntiles(n)/mfactor
     do m=1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        
        vwind(t) = 0.0           

     enddo
  enddo
  call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
  call LIS_verify(status)

  do k=1,LIS_rc%ntiles(n)/mfactor
     do m=1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        
        if((gddp_struc(n)%metdata1(m,6,index1).ne.LIS_rc%udef).and.&
             (gddp_struc(n)%metdata2(m,6,index1).ne.LIS_rc%udef)) then 
           psurf(t) =gddp_struc(n)%metdata1(m,6,index1)*wt1+ & 
                gddp_struc(n)%metdata2(m,6,index1)*wt2
        endif
        
     enddo
  enddo
  
end subroutine timeinterp_gddp
 
