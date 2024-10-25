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
! !ROUTINE: noah271_coldstart
! \label{noah271_coldstart}
!
! !REVISION HISTORY:
!  28 Apr 2002: Sujay Kumar; Initial version
!  27 Oct 2010: David Mocko, changes for Noah2.7.1 in LIS6.1
! 
! !INTERFACE:
subroutine noah271_coldstart()
! !USES:
   use LIS_coreMod
   use LIS_logMod,       only : LIS_logunit
   use LIS_timeMgrMod,   only : LIS_date2time
   use noah271_lsmMod
   use LIS_fileIOMod,      only : LIS_readData
!
! !DESCRIPTION:
!  
!  This routine initializes the Noah2.7.1 state variables with some 
!  predefined values uniformly for the entire domain.  These initial
!  values will be overwritten by the values read from the supplied 
!  Noah2.7.1 model restart file. 
! 
!EOP
   implicit none
   integer :: t,l,n

   integer :: i
   real, allocatable :: smcmax(:,:)
   real, allocatable :: psisat(:,:)
   real, allocatable :: bexpp(:,:)
   real, allocatable :: dksat(:,:)
   real, allocatable :: quartz(:,:)
   real, allocatable :: dwsat(:,:)
   real, allocatable :: smc1(:,:), smc2(:,:), smc3(:,:), smc4(:,:)
   real, allocatable :: sh2o1(:,:), sh2o2(:,:), sh2o3(:,:), sh2o4(:,:)
   real, allocatable :: rsmin(:,:), rgl(:,:), hs(:,:), z0(:,:), lai(:,:)
   real, allocatable :: cfactr(:,:), cmcmax(:,:),sbeta(:,:),rsmax(:,:)
   real, allocatable :: topt(:,:),refdk(:,:),fxexp(:,:),refkdt(:,:)
   real, allocatable :: czil(:,:),csoil(:,:),frzk(:,:),snup(:,:)
   real        :: cb_gridDesc(6)
   integer     :: c,r

   do n=1,LIS_rc%nnest
      if (trim(LIS_rc%startcode).eq."coldstart") then
         write(LIS_logunit,*)                                          &
               'MSG: noah271_coldstart -- cold-starting Noah2.7.1'
         do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
            noah271_struc(n)%noah(t)%t1=noah271_struc(n)%initskintemp
            noah271_struc(n)%noah(t)%cmc=noah271_struc(n)%initcanopywater
            noah271_struc(n)%noah(t)%snowh=noah271_struc(n)%initsnowdepth
            noah271_struc(n)%noah(t)%sneqv=noah271_struc(n)%initsnowequiv
!            noah271_struc(n)%noah(t)%ch=0.0150022404
!            noah271_struc(n)%noah(t)%cm=0.0205970779
            noah271_struc(n)%noah(t)%ch=1.E-4
            noah271_struc(n)%noah(t)%cm=1.E-4
            do l=1,noah271_struc(n)%nslay
               noah271_struc(n)%noah(t)%stc(l)=noah271_struc(n)%inittemp(l)
               noah271_struc(n)%noah(t)%smc(l)=noah271_struc(n)%initsm(l)
               noah271_struc(n)%noah(t)%sh2o(l)=noah271_struc(n)%initsmliq(l)
            enddo
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
      write(LIS_logunit,*) 'MSG: noah271_coldstart -- ',             &
                         'Using the specified start time ',LIS_rc%time

!begin hack for ga simulation
#if 0
      write(LIS_logunit,*) 'reading calibrated parameters...'
      do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
         noah271_struc(n)%noah(i)%smcmax = 2.01301366091e-01  
         noah271_struc(n)%noah(i)%psisat = exp(-5.50178766251e+00) 
         noah271_struc(n)%noah(i)%dksat  = exp(-1.53239078522e+01) 
         noah271_struc(n)%noah(i)%bexp   = exp(2.67149925232e+00)  
      enddo				 
#endif
! for mcmc run
#if 0 
      write(LIS_logunit,*) 'reading calibrated parameters...'
      open(111, file='smcmax_mcmc.bin', form="unformatted")        
      allocate(smcmax(LIS_rc%npatch(n,LIS_rc%lsm_index)))
      read(111) smcmax
      close(111)
      do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
         noah271_struc(n)%noah(i)%smcmax = smcmax(i)
      enddo
      deallocate(smcmax)

      open(111, file='psisat_mcmc.bin', form="unformatted")        
      allocate(psisat(LIS_rc%npatch(n,LIS_rc%lsm_index)))
      read(111) psisat
      close(111)
      do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
         noah271_struc(n)%noah(i)%psisat = psisat(i)
      enddo
      deallocate(psisat)

      open(111, file='dksat_mcmc.bin', form="unformatted")        
      allocate(dksat(LIS_rc%npatch(n,LIS_rc%lsm_index)))
      read(111) dksat
      close(111)
      do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
         noah271_struc(n)%noah(i)%dksat = dksat(i)
      enddo
      deallocate(dksat)

      open(111, file='bexp_mcmc.bin', form="unformatted")        
      allocate(bexpp(LIS_rc%npatch(n,LIS_rc%lsm_index)))
      read(111) bexpp
      close(111)
      do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
         noah271_struc(n)%noah(i)%bexp = bexpp(i)

      enddo
      deallocate(bexpp)
#endif
#if 0 
      open(111, file='quartz_mcmc.bin', form="unformatted")        
      allocate(quartz(LIS_rc%npatch(n,LIS_rc%lsm_index)))
      read(111) quartz
      close(111)
      do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
         noah271_struc(n)%noah(i)%quartz = quartz(i)
      enddo
      deallocate(quartz)
#endif
#if 0 
      cb_gridDesc(1) = 30.5
      cb_gridDesc(2) = -124.50
      cb_gridDesc(3) = 50.5
      cb_gridDesc(4) = -75.5
      cb_gridDesc(5) = 1.0
      cb_gridDesc(6) = 1.0 

      open(111, file='smcmax_calib.1gd4r', access='direct',status='old', &
           form="unformatted", recl=4)        
      allocate(smcmax(LIS_rc%lnc(n),LIS_rc%lnr(n)))    
      call LIS_readData(n, 111,  cb_gridDesc,smcmax)
!      read(111,rec=1) smcmax
      close(111)
      do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
         noah271_struc(n)%noah(i)%smcmax = &
              smcmax(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
      enddo
      deallocate(smcmax)
     allocate(psisat(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='psisat_calib.1gd4r',access='direct',status='old', &
           form="unformatted", recl=4)       
     call LIS_readData(n, 111,  cb_gridDesc,psisat)
!     read(111,rec=1) psisat
     close(111)
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%psisat = &
             psisat(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(psisat)

     allocate(dksat(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='dksat_calib.1gd4r',access='direct',status='old', &
           form="unformatted", recl=4)  
     call LIS_readData(n, 111,  cb_gridDesc,dksat)
!     read(111,rec=1) dksat
     close(111)
     
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%dksat = &
             dksat(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(dksat)

     allocate(dwsat(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='dwsat_calib.1gd4r',access='direct',status='old', &
           form="unformatted", recl=4) 
     call LIS_readData(n, 111,  cb_gridDesc,dwsat)
!     read(111,rec=1) dwsat
     close(111)
     
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%dwsat = &
             dwsat(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(dwsat)

     allocate(bexpp(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='bexp_calib.1gd4r',access='direct',status='old', &
           form="unformatted", recl=4) 
     call LIS_readData(n, 111,  cb_gridDesc,bexpp)
!     read(111,rec=1) bexpp
     close(111)
     
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%bexp = &
             bexpp(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(bexpp)

     allocate(quartz(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='quartz_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,quartz)
!     read(111,rec=1) quartz
     close(111)
     
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%quartz = &
             quartz(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(quartz)


     allocate(smc1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='smc1_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,smc1)
!     read(111,rec=1) smc1
     close(111)
     
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%smc(1) = &
             smc1(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(smc1)

     allocate(smc2(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='smc2_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,smc2)
!     read(111,rec=1) smc2
     close(111)
     
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%smc(2) = &
             smc2(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(smc2)

     allocate(smc3(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='smc3_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,smc3)
!     read(111,rec=1) smc3
     close(111)
     
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%smc(3) = &
             smc3(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(smc3)

     allocate(smc4(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='smc4_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,smc4)
!     read(111,rec=1) smc4
     close(111)
     
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%smc(4) = &
             smc4(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(smc4)

     allocate(sh2o1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='sh2o1_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,sh2o1)
!     read(111,rec=1) sh2o1
     close(111)
     
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%sh2o(1) = &
             sh2o1(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(sh2o1)

     allocate(sh2o2(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='sh2o2_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,sh2o2)
!     read(111,rec=1) sh2o2
     close(111)
     
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%sh2o(2) = &
             sh2o2(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(sh2o2)

     allocate(sh2o3(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='sh2o3_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,sh2o3)
!     read(111,rec=1) sh2o3
     close(111)
     
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%sh2o(3) = &
             sh2o3(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(sh2o3)

     allocate(sh2o4(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='sh2o4_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,sh2o4)
!     read(111,rec=1) sh2o4
     close(111)
     
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%sh2o(4) = &
             sh2o4(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(sh2o4)

     allocate(rsmin(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='rsmin_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,rsmin)
!     read(111,rec=1) rsmin
     close(111)
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%rsmin = &
             rsmin(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(rsmin)


     allocate(rgl(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='rgl_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,rgl)
!     read(111,rec=1) rgl
     close(111)
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%rgl = &
             rgl(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(rgl)

     allocate(hs(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='hs_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,hs)
!     read(111,rec=1) hs
     close(111)
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%hs = &
             hs(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(hs)

     allocate(z0(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='z0_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,z0)
!     read(111,rec=1) z0
     close(111)
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%z0 = &
             z0(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(z0)

     allocate(lai(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='lai_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,lai)
!     read(111,rec=1) lai
     close(111)
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%lai = &
             lai(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(lai)

     allocate(cfactr(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='cfactr_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,cfactr)
!     read(111,rec=1) cfactr
     close(111)
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%cfactr = &
             cfactr(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(cfactr)

     allocate(cmcmax(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='cmcmax_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,cmcmax)
!     read(111,rec=1) cmcmax
     close(111)
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%cmcmax = &
             cmcmax(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(cmcmax)

     allocate(sbeta(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='sbeta_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,sbeta)
!     read(111,rec=1) sbeta
     close(111)
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%sbeta = &
             sbeta(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(sbeta)

     allocate(rsmax(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='rsmax_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,rsmax)
!     read(111,rec=1) rsmax
     close(111)
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%rsmax = &
             rsmax(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(rsmax)

     allocate(topt(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='topt_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,topt)
!     read(111,rec=1) topt
     close(111)
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%topt = &
             topt(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(topt)

     allocate(refdk(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='refdk_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,refdk)
!     read(111,rec=1) refdk
     close(111)
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%refdk = &
             refdk(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(refdk)

     allocate(fxexp(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='fxexp_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,fxexp)
!     read(111,rec=1) fxexp
     close(111)
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%fxexp = &
             fxexp(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(fxexp)

     allocate(refkdt(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='refkdt_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,refkdt)
!     read(111,rec=1) refkdt
     close(111)
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%refkdt = &
             refkdt(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(refkdt)

     allocate(czil(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='czil_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,czil)
!     read(111,rec=1) czil
     close(111)
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%czil = &
             czil(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
        noah271_struc(n)%noah(i)%czil = 0.1
     enddo
     deallocate(czil)
     
     allocate(csoil(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='csoil_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,csoil)
!     read(111,rec=1) csoil
     close(111)
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%csoil = &
             csoil(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(csoil)

     allocate(frzk(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='frzk_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,frzk)
!     read(111,rec=1) frzk
     close(111)
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%frzk = &
             frzk(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(frzk)

     allocate(snup(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='snup_calib.1gd4r',access='direct',status='old',form='unformatted',recl=4)
     call LIS_readData(n, 111,  cb_gridDesc,snup)
!     read(111,rec=1) snup
     close(111)
     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        noah271_struc(n)%noah(i)%snup = &
             snup(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
     enddo
     deallocate(snup)
#endif
!end hack

   enddo
end subroutine noah271_coldstart
