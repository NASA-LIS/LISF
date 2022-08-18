!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: timeinterp_gdasbc
! \label{timeinterp_gdasbc}
!
! !REVISION HISTORY:
!
! 02 Feb 2004: Sujay Kumar; Initial Specification
! 24 Aug 2007: Chuck Alonge; Modified for use with NLDAS-2 data
! 14 Mar 2014: David Mocko: Added CAPE and PET forcing from NLDAS-2
!
! !INTERFACE:
subroutine timeinterp_gdasbc(n,findex)
! !USES:
  use ESMF
  use LIS_coreMod,    only : LIS_rc,LIS_domain
  use LIS_constantsMod, only  : LIS_CONST_SOLAR
  use LIS_metforcingMod, only : LIS_forc, LIS_FORC_Base_State
  use LIS_FORC_AttributesMod
  use LIS_timeMgrMod, only : LIS_tick, LIS_time2date
  use LIS_logMod,     only : LIS_logunit, LIS_verify
  use gdasbc_forcingMod, only : gdasbc_struc
  use LIS_forecastMod, only : LIS_get_iteration_index
 
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
  type(ESMF_Field)   :: pcpField
  real,pointer       :: pcp(:)
  integer            :: mfactor, m, k, kk
! ________________________________________

  btime=gdasbc_struc(n)%gdasbctime1
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
  
  btime=gdasbc_struc(n)%gdasbctime2
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
  wt1=(gdasbc_struc(n)%gdasbctime2-LIS_rc%time)/ & 
       (gdasbc_struc(n)%gdasbctime2-gdasbc_struc(n)%gdasbctime1)
  wt2=1.0-wt1
  swt1=(newtime2-LIS_rc%time)/(newtime2-newtime1)
  swt2=1.0-swt1

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Rainf%varname(1),pcpField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Rainf in the forcing variables list')


! Loop over number of forcing ensembles:
  mfactor = LIS_rc%nensem(n)/gdasbc_struc(n)%nIter

  ! do block precipitation interpolation
  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LIS_verify(status)

  do k=1,LIS_rc%ntiles(n)/mfactor
     do m=1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        kk = LIS_get_iteration_index(n, k, index1, mfactor)
        if(gdasbc_struc(n)%metdata2(kk,1,index1).ne.LIS_rc%udef) then 
          pcp(t)=gdasbc_struc(n)%metdata2(kk,1,index1)
        endif
     end do
  enddo

end subroutine timeinterp_gdasbc
 
