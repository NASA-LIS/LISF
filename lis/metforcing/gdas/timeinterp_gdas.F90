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
! !ROUTINE: timeinterp_gdas
! \label{timeinterp_gdas}
!
! !REVISION HISTORY:
!  30 Nov 2000: Jon Radakovich; Initial code based on geteta.f
!  17 Apr 2001: Jon Gottschalck; A few changes to allow model init.  
!  13 Aug 2001: Urszula Jambor; Introduced missing data replacement.     
!   5 Nov 2001: Urszula Jambor; Reset tiny negative SW values to zero. 
!  10 Oct 2006: Sujay Kumar: Switched to using ESMF_State for storing
!               forcing data. 
!  25 Jan 2012: Sujay Kumar; added the tweaks to enable ensemble forcings
! !INTERFACE:
subroutine timeinterp_gdas(n, findex)
! !USES:
  use ESMF
  use LIS_coreMod,        only  : LIS_rc, LIS_domain
  use LIS_FORC_AttributesMod 
  use LIS_timeMgrMod,     only  : LIS_time2date, LIS_tick
  use LIS_constantsMod,   only  : LIS_CONST_SOLAR
  use LIS_metforcingMod, only  : LIS_forc, LIS_FORC_Base_State
  use LIS_logMod,         only  : LIS_logunit, LIS_verify, LIS_endrun
  use gdas_forcingMod,    only  : gdas_struc

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
!  If the ensemble of forcing option is specified, then the code 
!  picks out the ensemble member (in each subgrid tile) to assign the
!  forcing. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
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
  integer :: t,zdoy
  integer :: index1
  integer :: bdoy,byr,bmo
  integer :: bda,bhr,bmn
  real*8  :: btime
  real    :: wt1,wt2,czb,cze,czm,gmt1,gmt2
  real    :: zw1,zw2
  integer           :: status
  type(ESMF_Field)  :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)  :: psurfField,pcpField,cpcpField
  real,pointer      :: tmp(:),q2(:),uwind(:),vwind(:)
  real,pointer      :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)

  btime=gdas_struc(n)%gdastime1
  call LIS_time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)
  btime=gdas_struc(n)%gdastime2
  call LIS_time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)
  
  wt1 = (gdas_struc(n)%gdastime2-LIS_rc%time) / & 
       (gdas_struc(n)%gdastime2-gdas_struc(n)%gdastime1)
  wt2 = 1.0 - wt1

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_Tair%varname(1)),tmpField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Tair in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_Qair%varname(1)),q2Field,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Qair in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_SWdown%varname(1)),swdField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable SWdown in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_LWdown%varname(1)),lwdField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable LWdown in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_Wind_E%varname(1)),uField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Wind_E in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_Wind_N%varname(1)),vField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Wind_N in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_Psurf%varname(1)),psurfField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Psurf in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_Rainf%varname(1)),pcpField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Rainf in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_CRainf%varname(1)),cpcpField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable CRainf in the forcing variables list')

  
  call ESMF_FieldGet(swdField,localDE=0, farrayPtr=swd,rc=status)
  call LIS_verify(status)
  
  call ESMF_FieldGet(lwdField,localDE=0, farrayPtr=lwd,rc=status)
  call LIS_verify(status)
        
  call ESMF_FieldGet(pcpField,localDE=0, farrayPtr=pcp,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(cpcpField,localDE=0, farrayPtr=cpcp,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(tmpField,localDE=0, farrayPtr=tmp,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(q2Field,localDE=0, farrayPtr=q2,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(uField,localDE=0, farrayPtr=uwind,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(vField,localDE=0, farrayPtr=vwind,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(psurfField,localDE=0, farrayPtr=psurf,rc=status)
  call LIS_verify(status)

  pcp  = 0.0
  cpcp = 0.0

!only the particular ensemble member is initialized. 
  do t = 1, LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     zdoy = LIS_rc%doy
     call zterp( 0, LIS_domain(n)%grid(index1)%lat,&
          LIS_domain(n)%grid(index1)%lon, gmt1, gmt2, & 
          LIS_rc%gmt,zdoy,zw1,zw2,czb,cze,czm,LIS_rc)
     
     swd(t) = zw1 * gdas_struc(n)%metdata2(3,index1)
     if (swd(t) < 0) then
        write(LIS_logunit,*) '[WARN] SW radiation is negative!!'
        write(LIS_logunit,*) '[WARN] sw=', swd(t), '... negative'
        write(LIS_logunit,*) '[WARN] gdas2=', gdas_struc(n)%metdata2(3,index1)
        call LIS_endrun
     end if
     
     if (swd(t).gt.LIS_CONST_SOLAR) then
        swd(t)=gdas_struc(n)%metdata2(3,index1)
     endif
  end do

  do t = 1, LIS_rc%ntiles(n)     
     index1 = LIS_domain(n)%tile(t)%index
     lwd(t) = gdas_struc(n)%metdata2(4,index1)     

     ! pcp comes from the averaged GDAS forcing data.
     ! In some cases, it is derived from the hrf06 and hrf03 values
     ! ( 2*hrf06 - hrf03 ).
     !
     ! One would expect the 2*hrf06 precip to be greater than or equal to
     ! the hrf03 precip, but in some cases it is slightly smaller.
     ! When the 2*hrf06 precip is slightly smaller than the hrf03, one gets
     ! a small negative precip rate.
     !
     ! The hrf06 file is a 6 hour forecast from hr.  So it contains the events
     ! that also occurred during the 3 hour forecast from hr, namely, hrf03.
     ! Thus hrf06 contains the precip from hrf03 plus any additional precip that
     ! occurred in the following 3 hours.
     !
     ! Should (2*hrf06 - hrf03) yield negative values for precip, then we
     ! assume that there was no additional precip for the last 3 hours of
     ! that 6-hour forecast and that the run of the 6-hour forecast produced
     ! slightly different results for the first 3 hours (as compared against
     ! the hrf03 forecast).
     !
     ! When LIS is using values from an hrf06 file, then LIS has already
     ! run over the period that corresponds to the first 3 hours of that
     ! 6-hour forecast.  LIS needs data corresponding to the last 3 hours.
     ! Thus always reset negative precip rates to 0.
     pcp(t) = gdas_struc(n)%metdata2(8,index1)    
     if ( pcp(t) < 0.0 ) then
        pcp(t) = 0.0
     endif

     cpcp(t) = gdas_struc(n)%metdata2(9,index1)     
     tmp(t) = wt1 * gdas_struc(n)%metdata1(1,index1) & 
          + wt2 *gdas_struc(n)%metdata2(1,index1)
     q2(t) = wt1 * gdas_struc(n)%metdata1(2,index1) & 
          + wt2 *gdas_struc(n)%metdata2(2,index1)
     uwind(t) = wt1 * gdas_struc(n)%metdata1(5,index1) & 
          + wt2 *gdas_struc(n)%metdata2(5,index1)     
     vwind(t) = wt1 * gdas_struc(n)%metdata1(6,index1) & 
          + wt2 *gdas_struc(n)%metdata2(6,index1)     
     psurf(t) = wt1 * gdas_struc(n)%metdata1(7,index1) & 
          + wt2 *gdas_struc(n)%metdata2(7,index1)
  enddo

end subroutine timeinterp_gdas
