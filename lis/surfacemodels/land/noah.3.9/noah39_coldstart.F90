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
! !ROUTINE: noah39_coldstart
! \label{noah39_coldstart}
!
! !REVISION HISTORY:
!  28 Apr 2002: Sujay Kumar; Initial version
!   8 May 2009: Sujay Kumar; additions for Noah3.1
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
!   9 Sep 2011: David Mocko, changes for Noah3.3 in LIS6.1
!  30 Oct 2014: David Mocko, added Noah-3.6 into LIS-7
! 
! !INTERFACE:
subroutine noah39_coldstart(mtype)
! !USES:
   use LIS_coreMod,      only : LIS_rc
   use LIS_logMod,       only : LIS_logunit
   use LIS_timeMgrMod,   only : LIS_date2time
   use noah39_lsmMod

!
! !DESCRIPTION:
!  
!  This routine initializes the Noah-3.9 state variables with some 
!  predefined values uniformly for the entire domain.  These initial
!  values will be overwritten by the values read from the supplied 
!  Noah-3.9 model restart file. 
! 
!EOP
   implicit none
   integer :: mtype
   integer :: t,l,n,i
   real, allocatable :: smcmax(:,:)
   real, allocatable :: psisat(:,:)
   real, allocatable :: bexpp(:,:)
   real, allocatable :: dksat(:,:)
   real, allocatable :: quartz(:,:)
   real, allocatable :: dwsat(:,:)
   real, allocatable :: rsmin(:,:), rgl(:,:), hs(:,:), z0(:,:), lai(:,:)
   real, allocatable :: cfactr(:,:), cmcmax(:,:),sbeta(:,:),rsmax(:,:)
   real, allocatable :: topt(:,:),refdk(:,:),fxexp(:,:),refkdt(:,:)
   real, allocatable :: czil(:,:),csoil(:,:),frzk(:,:),snup(:,:)
   real, allocatable :: smcref(:,:),smcdry(:,:),smcwlt(:,:),f1(:,:)
   real, allocatable :: slope(:,:),emiss(:,:)
   real        :: cb_gridDesc(6)
   integer     :: c,r

   do n=1,LIS_rc%nnest
      if (trim(LIS_rc%startcode).eq."coldstart") then
         write(LIS_logunit,*)                                          &
               '[INFO] noah39_coldstart -- cold-starting Noah-3.9'
         do t=1,LIS_rc%npatch(n,mtype)
            noah39_struc(n)%noah(t)%sigma_sm=0.04 ! needed for MCMC uncertainty analysis
            noah39_struc(n)%noah(t)%t1=noah39_struc(n)%initskintemp
            noah39_struc(n)%noah(t)%cmc=noah39_struc(n)%initcanopywater
            noah39_struc(n)%noah(t)%snowh=noah39_struc(n)%initsnowdepth
            noah39_struc(n)%noah(t)%sneqv=noah39_struc(n)%initsnowequiv
            noah39_struc(n)%noah(t)%ch=1.E-4
            noah39_struc(n)%noah(t)%cm=1.E-4
            do l=1,noah39_struc(n)%nslay
               noah39_struc(n)%noah(t)%stc(l)=noah39_struc(n)%inittemp(l)
               noah39_struc(n)%noah(t)%smc(l)=noah39_struc(n)%initsm(l)
               noah39_struc(n)%noah(t)%sh2o(l)=noah39_struc(n)%initsmliq(l)
            enddo
            noah39_struc(n)%noah(t)%emiss = 0.96
            noah39_struc(n)%noah(t)%snotime1 = 0.0
         enddo
      endif
      LIS_rc%yr=LIS_rc%syr
      LIS_rc%mo=LIS_rc%smo 
      LIS_rc%da=LIS_rc%sda
      LIS_rc%hr=LIS_rc%shr
      LIS_rc%mn=LIS_rc%smn
      LIS_rc%ss=LIS_rc%sss
      
      call LIS_date2time(LIS_rc%time,LIS_rc%doy,LIS_rc%gmt,LIS_rc%yr,&
           LIS_rc%mo,LIS_rc%da,LIS_rc%hr,LIS_rc%mn,LIS_rc%ss) 
      write(LIS_logunit,*) '[INFO] noah39_coldstart -- ',              &
                         'Using the specified start time ',LIS_rc%time
   enddo
end subroutine noah39_coldstart
