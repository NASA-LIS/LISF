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
! !ROUTINE: noah32_coldstart
! \label{noah32_coldstart}
!
! !REVISION HISTORY:
!  28 Apr 2002: Sujay Kumar; Initial version
!   8 May 2009: Sujay Kumar; additions for Noah3.1
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
! 
! !INTERFACE:
subroutine noah32_coldstart(mtype)
! !USES:
   use LIS_coreMod,      only : LIS_rc
   use LIS_logMod,       only : LIS_logunit
   use LIS_timeMgrMod,   only : LIS_date2time
   use noah32_lsmMod

!
! !DESCRIPTION:
!  
!  This routine initializes the Noah3.2 state variables with some 
!  predefined values uniformly for the entire domain.  These initial
!  values will be overwritten by the values read from the supplied 
!  Noah3.2 model restart file. 
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
               'MSG: noah32_coldstart -- cold-starting Noah3.2'
         do t=1,LIS_rc%npatch(n,mtype)
            noah32_struc(n)%noah(t)%t1=noah32_struc(n)%initskintemp
            noah32_struc(n)%noah(t)%cmc=noah32_struc(n)%initcanopywater
            noah32_struc(n)%noah(t)%snowh=noah32_struc(n)%initsnowdepth
            noah32_struc(n)%noah(t)%sneqv=noah32_struc(n)%initsnowequiv
            noah32_struc(n)%noah(t)%ch=1.E-4
            noah32_struc(n)%noah(t)%cm=1.E-4
            do l=1,noah32_struc(n)%nslay
               noah32_struc(n)%noah(t)%stc(l)=noah32_struc(n)%inittemp(l)
               noah32_struc(n)%noah(t)%smc(l)=noah32_struc(n)%initsm(l)
               noah32_struc(n)%noah(t)%sh2o(l)=noah32_struc(n)%initsmliq(l)
            enddo
            noah32_struc(n)%noah(t)%emiss = 0.96
            noah32_struc(n)%noah(t)%snotime1 = 0.0
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
      write(LIS_logunit,*) 'MSG: noah32_coldstart -- ',              &
                         'Using the specified start time ',LIS_rc%time
   enddo
end subroutine noah32_coldstart
