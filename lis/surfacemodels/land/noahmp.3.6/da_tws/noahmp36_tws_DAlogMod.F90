!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! 16Feb12  Ben Zaitchik; Initial Specification
!

module noahmp36_tws_DAlogMod
  
  use LIS_constantsMod,  only : LIS_CONST_RHOFW
  use ESMF
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------
  public :: noahmp36_tws_DAlog
!-----------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------
  public :: NOAHMPpred_struc
!EOP

  type, public ::NOAHMPpred_dec
     
     real,allocatable ::clmnwater(:,:)
     
  end type NOAHMPpred_dec
  
  type (NOAHMPpred_dec),allocatable :: NOAHMPpred_struc(:)
  
contains 
  
  subroutine noahmp36_tws_DAlog(n)
    
    ! USES:
    use LIS_coreMod, only : LIS_rc,LIS_surface
    use LIS_timeMgrMod
    use noahmp36_lsmMod
    use LIS_logMod, only : LIS_logunit, LIS_verify
    !      use smootherDA_runMod, only : smootherDA_increments_mode
      
    ! ARGUMENTS:  
    integer, intent(in)    :: n 
      
    ! DESCRIPTION:
    ! Calculates total column water storage three times per month, to
    ! approximate the GRACE return frequency
    
    integer                  :: i,m,t,gid,d          
    integer                  :: yr,mo,da,hr,mn,ss
    integer                  :: yr1, mo1, da1
    integer                  :: yr2, mo2, da2
    integer                  :: yr3, mo3, da3
    integer                  :: tw_tmp1, tw_tmp2
    type(ESMF_Time)          :: tTime1,tTime2,tTime3
    type(ESMF_TimeInterval)  :: tw1, tw2
    integer                  :: status

    if(LIS_rc%DAincrMode(n).eq.0) then
       if(LIS_rc%twInterval.eq.2592000.0) then 
          if((LIS_rc%da.eq.1).and.(LIS_rc%hr.eq.12).and.(LIS_rc%mn.eq.0)) then
             if(.not.allocated(NOAHMPpred_struc)) then 
                allocate(NOAHMPpred_struc(LIS_rc%nnest))
                allocate(NOAHMPpred_struc(n)%clmnwater(3,&
                     LIS_rc%npatch(n,LIS_rc%lsm_index)))
             endif
             NOAHMPpred_struc(n)%clmnwater = 0.0
          end if
          
          if(((LIS_rc%da.eq.5).and.(LIS_rc%hr.eq.4).and.(LIS_rc%mn.eq.0)).or. &
               ((LIS_rc%da.eq.15).and.(LIS_rc%hr.eq.12).and.(LIS_rc%mn.eq.0)).or. &
               ((LIS_rc%da.eq.25).and.(LIS_rc%hr.eq.18).and.(LIS_rc%mn.eq.0))) then
             
             d = (LIS_rc%da+5)/10
             write(LIS_logunit,*)'[INFO] logging obspred data for GRACE-DA'

             NOAHMPpred_struc(n)%clmnwater(d,:) = 0.0

             do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
                NOAHMPpred_struc(n)%clmnwater(d,t)= &
                     NOAHMPpred_struc(n)%clmnwater(d,t) + &
                     NOAHMP36_struc(n)%noahmp36(t)%sneqv   +         &
                     (NOAHMP36_struc(n)%noahmp36(t)%canliq  +         &
                     NOAHMP36_struc(n)%noahmp36(t)%canice) +         &
                     (NOAHMP36_struc(n)%noahmp36(t)%smc(1)  *         &
                     NOAHMP36_struc(n)%sldpth(1)*LIS_CONST_RHOFW)        +         &
                     (NOAHMP36_struc(n)%noahmp36(t)%smc(2)  *         &
                     NOAHMP36_struc(n)%sldpth(2)*LIS_CONST_RHOFW)        +         &
                     (NOAHMP36_struc(n)%noahmp36(t)%smc(3)  *         &
                     NOAHMP36_struc(n)%sldpth(3)*LIS_CONST_RHOFW)        +         &
                     (NOAHMP36_struc(n)%noahmp36(t)%smc(4)  *         &
                     NOAHMP36_struc(n)%sldpth(4)*LIS_CONST_RHOFW)        +         &
                     NOAHMP36_struc(n)%noahmp36(t)%wa                   
             enddo
          endif
       
       else
          call ESMF_TimeGet(LIS_twMidTime, yy = yr, &
               mm = mo, &
               dd = da, &
               h  = hr, &
               m  = mn,& 
               s  = ss, & 
               calendar = LIS_calendar, & 
               rc = status)

          if((LIS_rc%mo.eq.mo).and.(LIS_rc%da.eq.da) &
               .and.(LIS_rc%hr.eq.12).and.(LIS_rc%mn.eq.0)) then
             if(.not.allocated(NOAHMPpred_struc)) then 
                allocate(NOAHMPpred_struc(LIS_rc%nnest))
                allocate(NOAHMPpred_struc(n)%clmnwater(3,&
                     LIS_rc%npatch(n,LIS_rc%lsm_index)))
             endif
             NOAHMPpred_struc(n)%clmnwater = 0.0
          end if
          
          tw_tmp1 = nint(LIS_rc%twInterval/3.0)
          tw_tmp2 = tw_tmp1/2
          call ESMF_TimeIntervalSet(tw1,s=tw_tmp1,rc=status)
          call ESMF_TimeIntervalSet(tw2,s=tw_tmp2,rc=status)
          
          tTime1 = LIS_twMidTime + tw2
          tTime2 = tTime1 + tw1
          tTime3 = tTime2 + tw1
          
          call ESMF_TimeGet(tTime1,yy=yr1,mm=mo1,dd=da1,calendar=LIS_calendar,&
               rc=status)
          call ESMF_TimeGet(tTime2,yy=yr2,mm=mo2,dd=da2,calendar=LIS_calendar,&
               rc=status)
          call ESMF_TimeGet(tTime3,yy=yr3,mm=mo3,dd=da3,calendar=LIS_calendar,&
               rc=status)
          
          if(&
               ((LIS_rc%yr.eq.yr1).and.(LIS_rc%mo.eq.mo1).and.(LIS_rc%da.eq.da1)&
               .and.(LIS_rc%hr.eq.12).and.(LIS_rc%mn.eq.0)).or. &
               ((LIS_rc%yr.eq.yr2).and.(LIS_rc%mo.eq.mo2).and.(LIS_rc%da.eq.da2)&
               .and.(LIS_rc%hr.eq.12).and.(LIS_rc%mn.eq.0)).or. &
               ((LIS_rc%yr.eq.yr3).and.(LIS_rc%mo.eq.mo3).and.(LIS_rc%da.eq.da3)&
               .and.(LIS_rc%hr.eq.12).and.(LIS_rc%mn.eq.0)) & 
               ) then 
               
             d = -1
             if((LIS_rc%yr.eq.yr1).and.(LIS_rc%mo.eq.mo1).and.(LIS_rc%da.eq.da1)&
                  .and.(LIS_rc%hr.eq.12).and.(LIS_rc%mn.eq.0)) then 
                d = 1
             elseif((LIS_rc%yr.eq.yr2).and.(LIS_rc%mo.eq.mo2).and.(LIS_rc%da.eq.da2)&
                  .and.(LIS_rc%hr.eq.12).and.(LIS_rc%mn.eq.0)) then 
                d = 2
             elseif((LIS_rc%yr.eq.yr3).and.(LIS_rc%mo.eq.mo3).and.(LIS_rc%da.eq.da3)&
                  .and.(LIS_rc%hr.eq.12).and.(LIS_rc%mn.eq.0)) then 
                d = 3
             endif
             write(LIS_logunit,*)'[INFO] logging obspred data for GRACE-DA'
             NOAHMPpred_struc(n)%clmnwater(d,:) = 0.0

             do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
                NOAHMPpred_struc(n)%clmnwater(d,t)= &
                     NOAHMPpred_struc(n)%clmnwater(d,t) + &
                     NOAHMP36_struc(n)%noahmp36(t)%sneqv   +         &
                     (NOAHMP36_struc(n)%noahmp36(t)%canliq  +         &
                     NOAHMP36_struc(n)%noahmp36(t)%canice) +         &
                     (NOAHMP36_struc(n)%noahmp36(t)%smc(1)  *         &
                     NOAHMP36_struc(n)%sldpth(1)*LIS_CONST_RHOFW)        +         &
                     (NOAHMP36_struc(n)%noahmp36(t)%smc(2)  *         &
                     NOAHMP36_struc(n)%sldpth(2)*LIS_CONST_RHOFW)        +         &
                     (NOAHMP36_struc(n)%noahmp36(t)%smc(3)  *         &
                     NOAHMP36_struc(n)%sldpth(3)*LIS_CONST_RHOFW)        +         &
                     (NOAHMP36_struc(n)%noahmp36(t)%smc(4)  *         &
                     NOAHMP36_struc(n)%sldpth(4)*LIS_CONST_RHOFW)        +         &
                     NOAHMP36_struc(n)%noahmp36(t)%wa                   

             enddo
          
          endif
          
       endif
    endif

  end subroutine noahmp36_tws_DAlog
  
end module noahmp36_tws_DAlogMod
