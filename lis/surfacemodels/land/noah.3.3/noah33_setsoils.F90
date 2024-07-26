!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!#include "LIS_misc.h"
!BOP
!
! !ROUTINE: noah33_setsoils
!  \label{noah33_setsoils}
!
! !REVISION HISTORY:
!  28 Apr 2002: Kristi Arsenault;  Added Noah2.5 LSM, Initial Code
!  13 Oct 2003: Sujay Kumar; Domain independent modifications
!  15 Sep 2005: Matthew Garcia; update for ARMS project (Windows capability)
!  04 Oct 2005: Matthew Garcia; additions for Noah Unified
!   9 Jul 2008: Chris Franks/AFWA; Save wilting point for calculation of relative soil moisture
!   8 May 2009: Sujay Kumar; additions for Noah3.1
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
!   9 Sep 2011: David Mocko, changes for Noah3.3 in LIS6.1
!  20 Jan 2011: David Mocko, added slope type
!
! !INTERFACE:
subroutine noah33_setsoils(mtype)
! !USES:
  use LIS_coreMod,   only : LIS_rc, LIS_domain, LIS_surface
  use LIS_soilsMod,  only : LIS_soils
  use LIS_logMod,    only : LIS_logunit, LIS_endrun, &
                            LIS_getNextUnitNumber, LIS_releaseUnitNumber
  use LIS_fileIOMod, only : LIS_read_param
  use noah33_lsmMod     

!
! !DESCRIPTION:
!  This subroutine sets the soil parameters in Noah3.3.  Noah3.3 uses a
!  lookup table based on the input soil texture to derive these parameters. 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[LIS\_mapSoilType](\ref{LIS_mapSoilType}) \newline
!     method to derive the soil texture type from the sand, silt, and 
!     clay fractions.
!   \end{description}
!EOP      
  implicit none
  integer :: mtype
  integer :: gid
  integer :: jj,t,i,k,c,r,n
!  real :: maxdt,mindz,tiledt,psis,ks,b
!  integer :: newdt,idiv

  real, allocatable :: smcmax(:,:)
  real, allocatable :: psisat(:,:)
  real, allocatable :: dwsat(:,:)
  real, allocatable :: dksat(:,:)

  integer           :: slcats
  integer           :: iindex
  character*4       :: sltype
  real, allocatable     :: bb(:)
  real, allocatable     :: drysmc(:)
  real, allocatable     :: f11(:)
  real, allocatable     :: maxsmc(:)
  real, allocatable     :: refsmc(:)
  real, allocatable     :: satpsi(:)
  real, allocatable     :: satdk(:)
  real, allocatable     :: satdw(:)
  real, allocatable     :: wltsmc(:)
  real, allocatable     :: qtz(:)
  real, allocatable     :: fxexp(:)
  real, allocatable     :: sbeta(:)
  
  real              :: frzfact
  real, allocatable :: slope_data(:)
  real ::  SBETA_DATA,FXEXP_DATA,CSOIL_DATA,SALP_DATA,REFDK_DATA,      &
           REFKDT_DATA,FRZK_DATA,ZBOT_DATA,SMLOW_DATA,SMHIGH_DATA,     &
           CZIL_DATA,LVCOEF_DATA
  integer :: NUM_SLOPE
  real, allocatable :: placeslopetype(:,:)
  real              :: zsoil
  real, allocatable :: dsoil(:,:)
  real, allocatable :: droot(:,:)
  integer       :: ftn
  logical       :: file_exists 

  do n=1,LIS_rc%nnest

     write(unit=LIS_logunit,fmt=*)                                     &
                            '[INFO] reading soil files in Noah LSM'

     if(LIS_rc%usetexturemap(n).eq."none") then
        do i=1,LIS_rc%npatch(n,mtype)
           call LIS_mapSoilType(noah33_struc(n)%soilscheme,&
                LIS_surface(n,mtype)%tile(i)%sand, &
                LIS_surface(n,mtype)%tile(i)%clay, &
                LIS_surface(n,mtype)%tile(i)%silt,&
                LIS_surface(n,mtype)%tile(i)%soilt)
        enddo
     endif

     if (noah33_struc(n)%fixedsoiltype.ne.0) then
        write(LIS_logunit,*) '[INFO] Fixing Noah3.3 soil index to type: ',    &
                              noah33_struc(n)%fixedsoiltype
     endif
     do i=1,LIS_rc%npatch(n,mtype)
        if (noah33_struc(n)%fixedsoiltype.eq.0) then
           noah33_struc(n)%noah(i)%soiltype = &
                LIS_surface(n,mtype)%tile(i)%soilt
        else
           noah33_struc(n)%noah(i)%soiltype = noah33_struc(n)%fixedsoiltype
        endif
        !     if(noah33_struc(n)%noah(i)%soiltype.eq.-9999) noah33_struc(n)%noah(i)%soiltype = 14     
        if (noah33_struc(n)%noah(i)%soiltype .eq. -9999.or.&
             noah33_struc(n)%noah(i)%soiltype.eq.0) then
           if ( noah33_struc(n)%soilscheme == 2 ) then ! statsgo
              noah33_struc(n)%noah(i)%soiltype = 12 ! clay
           else ! zobler
              noah33_struc(n)%noah(i)%soiltype = 6 ! clay loam
           endif
        endif
        if ( noah33_struc(n)%noah(i)%soiltype .eq. 14 .and. &
             noah33_struc(n)%noah(i)%vegt .ne. LIS_rc%waterclass ) then
           noah33_struc(n)%noah(i)%soiltype = 7
        endif
     enddo
!Hardcoded for STATSGO classes only for soil types found to be water
! at a land point

     if (noah33_struc(n)%fixedslopetype.ne.0) then
        write(LIS_logunit,*) '[INFO] Fixing Noah3.3 slope index to type: ',   &
                              noah33_struc(n)%fixedslopetype
     else
        write(LIS_logunit,*) '[INFO] Reading Noah3.3 slope index from maps'
        allocate(placeslopetype(LIS_rc%lnc(n),LIS_rc%lnr(n)))
        call LIS_read_param(n,"SLOPETYPE",placeslopetype)
     endif

     if(noah33_struc(n)%usedsoilmap.ne.0) then 
        allocate(dsoil(LIS_rc%lnc(n),LIS_rc%lnr(n)))

        call LIS_read_param(n,"SOILDEPTH",dsoil)
        
        do t=1,LIS_rc%npatch(n,mtype)
           allocate(noah33_struc(n)%noah(t)%lyrthk(noah33_struc(n)%nslay))
           zsoil = 0.0
           c = LIS_surface(n,mtype)%tile(t)%col
           r = LIS_surface(n,mtype)%tile(t)%row

           do k=1,noah33_struc(n)%nslay
              if(k.ne.noah33_struc(n)%nslay) then 
                 noah33_struc(n)%noah(t)%lyrthk(k) = noah33_struc(n)%lyrthk(k)
                 zsoil = zsoil +  noah33_struc(n)%lyrthk(k)
              else
                 if((dsoil(c,r)-zsoil).gt.0) then
                    noah33_struc(n)%noah(t)%lyrthk(k) = &
                         dsoil(c,r) - zsoil
                 else
                    write(LIS_logunit,*) '[ERR] The total soil depth is shorter'
                    write(LIS_logunit,*) '[ERR] than the layer structure. Please'
                    write(LIS_logunit,*) '[ERR] setup a thinner layer structure'
                    call LIS_endrun()
                 endif
              endif
           enddo
        enddo
        deallocate(dsoil)
     endif

     if(noah33_struc(n)%usedrootmap.ne.0) then 
        allocate(droot(LIS_rc%lnc(n),LIS_rc%lnr(n)))

        call LIS_read_param(n,"ROOTDEPTH",droot)
        
        do t=1,LIS_rc%npatch(n,mtype)

           c = LIS_surface(n,mtype)%tile(t)%col
           r = LIS_surface(n,mtype)%tile(t)%row
           
           zsoil = 0.0
           noah33_struc(n)%noah(t)%nroot = 1
           do k=1,noah33_struc(n)%nslay
              zsoil = zsoil + noah33_struc(n)%noah(t)%lyrthk(k)
              if(zsoil.ge.droot(c,r)) then 
                 exit; 
              else
                 noah33_struc(n)%noah(t)%nroot = &
                      noah33_struc(n)%noah(t)%nroot + 1
              endif
           enddo
        enddo
        deallocate(droot)
     endif

     do i=1,LIS_rc%npatch(n,mtype)
        if (noah33_struc(n)%fixedslopetype.eq.0) then
           if (placeslopetype(LIS_surface(n,mtype)%tile(i)%col,               &
                LIS_surface(n,mtype)%tile(i)%row).gt.0) then
              noah33_struc(n)%noah(i)%slopetype =                      &
                       int(placeslopetype(LIS_surface(n,mtype)%tile(i)%col,   &
                                          LIS_surface(n,mtype)%tile(i)%row))
           else
              noah33_struc(n)%noah(i)%slopetype = 3
              write(LIS_logunit,*) '[INFO] Noah3.3 slope type -- ',           &
                                   'not defined for point ',i
              write(LIS_logunit,*) '[INFO] Noah3.3 slope type -- ',           &
                                   'set to default value of 3'
           endif
        else
           noah33_struc(n)%noah(i)%slopetype = noah33_struc(n)%fixedslopetype
        endif
     enddo
     if (noah33_struc(n)%fixedslopetype.eq.0) then
        deallocate(placeslopetype)
     endif

!-----------------------------------------------------------------------
! Read in the Noah3.3 Soil Parameter File
!-----------------------------------------------------------------------
     inquire(file=noah33_struc(n)%sfile,exist=file_exists)
     if(.not.file_exists) then
        write(LIS_logunit,*) '[ERR] Noah3.3 soil parameter table file, ',  &
             trim(noah33_struc(n)%sfile),', does not exist. '
        write(LIS_logunit,*) '[ERR] Program stopping ...'
       call LIS_endrun()
     endif

     ftn = LIS_getNextUnitNumber()       
     write(LIS_logunit,*) '[INFO] Reading Noah3.3 soil parameter file: ',     &
                                             trim(noah33_struc(n)%sfile)
     open(unit=ftn,file=noah33_struc(n)%sfile,status='old',             &
                                                    access='sequential')
     read(ftn,*) 
!     read(ftn,*) (noah33_struc(n)%lyrthk(i), i=1,noah33_struc(n)%nslay)
!     read(ftn,*)
     read(ftn,fmt='(A4)') SLTYPE
     read(ftn,*) SLCATS,IINDEX

     allocate(bb(slcats))
     allocate(drysmc(slcats))
     allocate(f11(slcats))
     allocate(maxsmc(slcats))
     allocate(refsmc(slcats))
     allocate(satpsi(slcats))
     allocate(satdk(slcats))
     allocate(satdw(slcats))
     allocate(wltsmc(slcats))
     allocate(qtz(slcats))
     allocate(fxexp(slcats))
     allocate(sbeta(slcats))

     do i=1,SLCATS
        read(ftn,*) IINDEX, bb(i), drysmc(i), f11(i), maxsmc(i), &
             refsmc(i), satpsi(i), satdk(i), satdw(i), wltsmc(i), qtz(i),&
             fxexp(i),sbeta(i)
     enddo    
     close(ftn)
!-----------------------------------------------------------------------
! Assign SOIL Parameters to each tile based on the
! type of Zobler soil class present in that tile.
!-----------------------------------------------------------------------
     do i = 1,LIS_rc%npatch(n,mtype)
        noah33_struc(n)%noah(i)%smcwlt = &
             wltsmc(noah33_struc(n)%noah(i)%soiltype)
        noah33_struc(n)%noah(i)%smcdry = &
             drysmc(noah33_struc(n)%noah(i)%soiltype)
        noah33_struc(n)%noah(i)%smcref = &
             refsmc(noah33_struc(n)%noah(i)%soiltype)
        noah33_struc(n)%noah(i)%dwsat = &
             satdw(noah33_struc(n)%noah(i)%soiltype)
     end do     
! Set tile soil porosity
     if (LIS_rc%useporositymap(n).eq."none") then ! default, from look-up table
        do i = 1,LIS_rc%npatch(n,mtype)
           noah33_struc(n)%noah(i)%smcmax = &
                maxsmc(noah33_struc(n)%noah(i)%soiltype)
        end do
     else
        do i = 1,LIS_rc%npatch(n,mtype)
           noah33_struc(n)%noah(i)%smcmax = &
                LIS_soils(n)%porosity(LIS_surface(n,mtype)%tile(i)%col,&
                LIS_surface(n,mtype)%tile(i)%row,1)
        end do
     end if
     
! Set tile soil Psi-sat
     if (LIS_rc%usepsisatmap(n).eq."none") then ! default, from look-up table
        do i = 1,LIS_rc%npatch(n,mtype)
           noah33_struc(n)%noah(i)%psisat = &
                satpsi(noah33_struc(n)%noah(i)%soiltype)
        end do
     else
        do i = 1,LIS_rc%npatch(n,mtype)
           noah33_struc(n)%noah(i)%psisat = &
                LIS_soils(n)%psisat(LIS_surface(n,mtype)%tile(i)%col,&
                LIS_surface(n,mtype)%tile(i)%row)
        end do
     end if

! Set tile soil K-sat
     if (LIS_rc%useksatmap(n).eq."none") then ! default, from look-up table
        do i = 1,LIS_rc%npatch(n,mtype)
           noah33_struc(n)%noah(i)%dksat = &
                satdk(noah33_struc(n)%noah(i)%soiltype)
        end do
     else
        do i = 1,LIS_rc%npatch(n,mtype)
           noah33_struc(n)%noah(i)%dksat = &
                LIS_soils(n)%ksat(LIS_surface(n,mtype)%tile(i)%col,&
                LIS_surface(n,mtype)%tile(i)%row)
        end do
     end if
     
! Set tile soil b-parameter
     if (LIS_rc%usebexpmap(n).eq."none") then ! default, from look-up table
        do i = 1,LIS_rc%npatch(n,mtype)
           noah33_struc(n)%noah(i)%bexp = &
                bb(noah33_struc(n)%noah(i)%soiltype)
        end do
     else
        do i = 1,LIS_rc%npatch(n,mtype)
           noah33_struc(n)%noah(i)%bexp = &
                LIS_soils(n)%bexp(LIS_surface(n,mtype)%tile(i)%col,&
                LIS_surface(n,mtype)%tile(i)%row)
        end do
     end if
     
! Set tile soil Qz content
     if (LIS_rc%usequartzmap(n).eq."none") then ! default
        do i = 1,LIS_rc%npatch(n,mtype)
           noah33_struc(n)%noah(i)%quartz = &
                qtz(noah33_struc(n)%noah(i)%soiltype)
        end do
     else 
        do i = 1,LIS_rc%npatch(n,mtype)
           noah33_struc(n)%noah(i)%quartz = &
                LIS_soils(n)%quartz(LIS_surface(n,mtype)%tile(i)%col,&
                LIS_surface(n,mtype)%tile(i)%row)
        end do
     end if
     do i = 1,LIS_rc%npatch(n,mtype)
        noah33_struc(n)%noah(i)%f1 = &
             f11(noah33_struc(n)%noah(i)%soiltype)
        noah33_struc(n)%noah(i)%fxexp = &
             fxexp(noah33_struc(n)%noah(i)%soiltype)
        noah33_struc(n)%noah(i)%sbeta = &
             sbeta(noah33_struc(n)%noah(i)%soiltype)
     enddo
     deallocate(bb)
     deallocate(drysmc)
     deallocate(f11)
     deallocate(maxsmc)
     deallocate(refsmc)
     deallocate(satpsi)
     deallocate(satdk)
     deallocate(satdw)
     deallocate(wltsmc)
     deallocate(qtz)
     deallocate(fxexp)
     deallocate(sbeta)

     ! These are parameters from the GENPARM.TBL
     inquire(file=noah33_struc(n)%gfile,exist=file_exists)
     if(.not.file_exists) then
        write(LIS_logunit,*) '[ERR] Noah3.3 general parameter table file, ',  &
             trim(noah33_struc(n)%gfile),', does not exist. '
        write(LIS_logunit,*) '[ERR] Program stopping ...'
       call LIS_endrun()
     endif

     write(LIS_logunit,*) '[INFO] Reading Noah3.3 general parameter file: ',  &
                           trim(noah33_struc(n)%gfile)
     open(unit=ftn,file=noah33_struc(n)%gfile,status='old',             &
                       access='sequential')
     read(ftn,*)
     read(ftn,*)
     read(ftn,*) NUM_SLOPE
     allocate(slope_data(NUM_SLOPE))
     do jj=1,NUM_SLOPE
        read(ftn,*) slope_data(jj)
     enddo
!     read(ftn,*)
!     read(ftn,*)SBETA_DATA
!     read(ftn,*)
!     read(ftn,*)FXEXP_DATA
     read(ftn,*)
     read(ftn,*)CSOIL_DATA
     read(ftn,*)
     read(ftn,*)SALP_DATA
     read(ftn,*)
     read(ftn,*)REFDK_DATA
     read(ftn,*)
     read(ftn,*)REFKDT_DATA
     read(ftn,*)
     read(ftn,*)FRZK_DATA
     read(ftn,*)
     read(ftn,*)ZBOT_DATA
!     read(ftn,*)
!     read(ftn,*)CZIL_DATA
     read(ftn,*)
     read(ftn,*)SMLOW_DATA
     read(ftn,*)
     read(ftn,*)SMHIGH_DATA
     read(ftn,*)
     read(ftn,*)LVCOEF_DATA
     close(ftn)
     call LIS_releaseUnitNumber(ftn)
        
     do i = 1,LIS_rc%npatch(n,mtype)
        if (LIS_rc%useslopemap(n).eq."none") then
           noah33_struc(n)%noah(i)%slope =                             &
                slope_data(noah33_struc(n)%noah(i)%slopetype)
        else
           gid = LIS_surface(n,mtype)%tile(i)%tile_id
           noah33_struc(n)%noah(i)%slope = LIS_domain(n)%tile(gid)%slope
        endif
!        noah33_struc(n)%noah(i)%sbeta = SBETA_DATA
!        noah33_struc(n)%noah(i)%fxexp = FXEXP_DATA
        noah33_struc(n)%noah(i)%csoil = CSOIL_DATA
        noah33_struc(n)%noah(i)%salp = SALP_DATA
        noah33_struc(n)%noah(i)%refdk = REFDK_DATA
        noah33_struc(n)%noah(i)%refkdt = REFKDT_DATA
        noah33_struc(n)%noah(i)%frzk = FRZK_DATA
        noah33_struc(n)%noah(i)%zbot = ZBOT_DATA
!        noah33_struc(n)%noah(i)%czil = CZIL_DATA
        noah33_struc(n)%noah(i)%lvcoef = LVCOEF_DATA
     enddo

! overwrite default values with PTF 
     if(noah33_struc(n)%useptf.eq.1) then !use cosby et al.(1984) PTFs
        do i = 1,LIS_rc%npatch(n,mtype)
           
           noah33_struc(n)%noah(i)%smcmax = 0.489 - &
                0.00126*LIS_surface(n,mtype)%tile(i)%sand*100.0

           noah33_struc(n)%noah(i)%psisat = (10.0*(10.0**(1.88-0.0131*&
                LIS_surface(n,mtype)%tile(i)%sand*100.0)))/1000.0

           noah33_struc(n)%noah(i)%dksat =  (0.0070556*(10.0**(-0.884+0.0153*&
                LIS_surface(n,mtype)%tile(i)%sand*100.0)))/1000.0

           noah33_struc(n)%noah(i)%bexp = 2.91 + (0.159*&
                LIS_surface(n,mtype)%tile(i)%clay*100.0)

        end do
     endif
     
     deallocate(slope_data)

!begin hack
#if 0 
!     allocate(smcmax(LIS_rc%lnc(n),LIS_rc%lnr(n)))
!     open(111,file='smcmax_calib.bin',form='unformatted')
!     read(111) smcmax
!     close(111)
     
!     do i=1,LIS_rc%npatch(n,mtype)
!        noah33_struc(n)%noah(i)%smcmax = &
!             smcmax(LIS_surface(n,mtype)%tile(i)%col, LIS_surface(n,mtype)%tile(i)%row)
!     enddo
!     deallocate(smcmax)

     allocate(z0(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='z0_calib.bin',form='unformatted')
     read(111) z0
     close(111)
     
     do i=1,LIS_rc%npatch(n,mtype)
        noah33_struc(n)%noah(i)%z0 = &
             z0(LIS_surface(n,mtype)%tile(i)%col, LIS_surface(n,mtype)%tile(i)%row)
     enddo
     write(LIS_logunit,*) 'z0 ',z0
     deallocate(z0)

     allocate(rsmin(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='rsmin_calib.bin',form='unformatted')
     read(111) rsmin
     close(111)
     
     do i=1,LIS_rc%npatch(n,mtype)
        noah33_struc(n)%noah(i)%rsmin = &
             rsmin(LIS_surface(n,mtype)%tile(i)%col, LIS_surface(n,mtype)%tile(i)%row)
     enddo
     write(LIS_logunit,*) 'rsmin ',rsmin
     deallocate(rsmin)

     allocate(emiss(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='emiss_calib.bin',form='unformatted')
     read(111) emiss
     close(111)
     
     do i=1,LIS_rc%npatch(n,mtype)
        noah33_struc(n)%noah(i)%emiss = &
             emiss(LIS_surface(n,mtype)%tile(i)%col, LIS_surface(n,mtype)%tile(i)%row)
     enddo
     write(LIS_logunit,*) 'emiss ',emiss
     deallocate(emiss)

     allocate(czil(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='czil_calib.bin',form='unformatted')
     read(111) czil
     close(111)
     
     do i=1,LIS_rc%npatch(n,mtype)
        noah33_struc(n)%noah(i)%czil = &
             czil(LIS_surface(n,mtype)%tile(i)%col, LIS_surface(n,mtype)%tile(i)%row)
     enddo
     write(LIS_logunit,*) 'czil ',czil
     deallocate(czil)

     allocate(psisat(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='psisat_calib.bin',form='unformatted')
     read(111) psisat
     close(111)
     
     do i=1,LIS_rc%npatch(n,mtype)
        noah33_struc(n)%noah(i)%psisat = &
             psisat(LIS_surface(n,mtype)%tile(i)%col, LIS_surface(n,mtype)%tile(i)%row)
     enddo
     write(LIS_logunit,*) 'psisat ',psisat
     deallocate(psisat)

     allocate(dksat(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='dksat_calib.bin',form='unformatted')
     read(111) dksat
     close(111)
     
     do i=1,LIS_rc%npatch(n,mtype)
        noah33_struc(n)%noah(i)%dksat = &
             dksat(LIS_surface(n,mtype)%tile(i)%col, LIS_surface(n,mtype)%tile(i)%row)
     enddo
     write(LIS_logunit,*) 'dksat ',dksat
     deallocate(dksat)

     allocate(dwsat(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='dwsat_calib.bin',form='unformatted')
     read(111) dwsat
     close(111)
     
     do i=1,LIS_rc%npatch(n,mtype)
        noah33_struc(n)%noah(i)%dwsat = &
             dwsat(LIS_surface(n,mtype)%tile(i)%col, LIS_surface(n,mtype)%tile(i)%row)
     enddo
     write(LIS_logunit,*) 'dwsat ',dwsat
     deallocate(dwsat)

     allocate(quartz(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='quartz_calib.bin',form='unformatted')
     read(111) quartz
     close(111)
     
     do i=1,LIS_rc%npatch(n,mtype)
        noah33_struc(n)%noah(i)%quartz = &
             quartz(LIS_surface(n,mtype)%tile(i)%col, LIS_surface(n,mtype)%tile(i)%row)
     enddo
     write(LIS_logunit,*) 'quartz ',quartz
     deallocate(quartz)

     allocate(sbeta(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='sbeta_calib.bin',form='unformatted')
     read(111) sbeta
     close(111)
     
     do i=1,LIS_rc%npatch(n,mtype)
        noah33_struc(n)%noah(i)%sbeta = &
             sbeta(LIS_surface(n,mtype)%tile(i)%col, LIS_surface(n,mtype)%tile(i)%row)
     enddo
     write(LIS_logunit,*) 'sbeta ',sbeta
     deallocate(sbeta)

     allocate(fxexp(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     open(111,file='fxexp_calib.bin',form='unformatted')
     read(111) fxexp
     close(111)
     
     do i=1,LIS_rc%npatch(n,mtype)
        noah33_struc(n)%noah(i)%fxexp = &
             fxexp(LIS_surface(n,mtype)%tile(i)%col, LIS_surface(n,mtype)%tile(i)%row)

        noah33_struc(n)%noah(i)%kdt = noah33_struc(n)%refkdt* noah33_struc(n)%noah(i)%dksat/&
             noah33_struc(n)%refdk

        frzfact = &
             (noah33_struc(n)%noah(i)%smcmax/noah33_struc(n)%noah(i)%smcref)*(0.412/0.468)

        noah33_struc(n)%noah(i)%frzx = noah33_struc(n)%frzk*frzfact 

     enddo
     write(LIS_logunit,*) 'fxexp ',fxexp
     deallocate(fxexp)

#endif
!end hack

!  call soils_finalize()

!-----------------------------------------------------------------------
! Check numerical stability at saturation on basis of timestep and 
! soil parameters.  From formula derived by Yihua Wu, NASA/GSFC HSB
!-----------------------------------------------------------------------
!This routine as of now is "unsafe" in a parallel mode since each processor
!can change the timestep. There needs to be some gather/scatter...
#if 0 
     maxdt = LIS_rc%ts
     mindz = noahattrib%lyrthk(1) 
     do i = 2,noah33_struc(n)%nslay ! soil layer loop
        mindz = min(mindz,noahattrib%lyrthk(i))
     end do ! i
     do i = 1,LIS_rc%npatch(n,mtype)          ! tile loop
! NOTE: depends on *positive* Psi-sat values, as in Noah3.3 LSM
        psis = noah33_struc(n)%noah(i)%psisat
        ks = noah33_struc(n)%noah(i)%dksat
        b = noah33_struc(n)%noah(i)%bexp
        tiledt = 2 * mindz**2 / (b * psis * ks)
        maxdt = min(maxdt,tiledt) 
     end do ! i
     if (LIS_rc%ts.gt.maxdt) then
        write(LIS_logunit,*) '[WARN] noah33_setsoils -- timestep may',   &
                             '[WARN] be too large for numerical stability'
        write(LIS_logunit,*) '[WARN] noah33_setsoils -- soil parameters',&
                             '[WARN] suggest dt .le.',maxdt
        newdt = 3600
        idiv = 1
        do while (newdt.gt.maxdt)
           if (mod(3600,idiv).eq.0) then
              newdt = 3600 / idiv
           endif
           idiv = idiv + 1
           if (idiv.gt.60) exit
        enddo
        if (newdt.gt.maxdt) then
           write(LIS_logunit,*) ' noah33_setsoils -- new timestep',&
                                'is less than 1 minute -- stopping'
           call LIS_endrun
        else
           write(LIS_logunit,*) ' noah33_setsoils -- setting',     &
                                'timestep to ',newdt,' sec'
           LIS_rc%ts = newdt
        end if
     else
        write(LIS_logunit,*) ' noah33_setsoils -- timestep OK',    &
                      'for soil hydraulic stability (',LIS_rc%ts,' sec)'
     end if
#endif
  enddo
  return

end subroutine noah33_setsoils

!BOP
!
! !ROUTINE: noah33_resetsoils
!  \label{noah33_resetsoils}
!
! !REVISION HISTORY:
!  28 Apr 2002: Kristi Arsenault;  Added Noah2.5 LSM, Initial Code
!  13 Oct 2003: Sujay Kumar; Domain independent modifications
!  15 Sep 2005: Matthew Garcia; update for ARMS project (Windows capability)
!  04 Oct 2005: Matthew Garcia; additions for Noah Unified
!   9 Jul 2008: Chris Franks/AFWA; Save wilting point for calculation of relative soil moisture
!   8 May 2009: Sujay Kumar; additions for Noah3.1
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
!   9 Sep 2011: David Mocko, changes for Noah3.3 in LIS6.1
!  20 Jan 2011: David Mocko, added slope type
!
! !INTERFACE:
subroutine noah33_resetsoils(mtype)
! !USES:
  use LIS_logMod,    only : LIS_logunit
!EOP

  implicit none
  integer :: mtype

  write(LIS_logunit,*)                        &
       'Noah3.3 resetting soils'
  

!!!! !USES:
!!!  use LIS_coreMod,   only : LIS_rc, LIS_domain, LIS_surface
!!!  use LIS_soilsMod,  only : LIS_soils
!!!  use LIS_logMod,    only : LIS_logunit, LIS_endrun, &
!!!                            LIS_getNextUnitNumber, LIS_releaseUnitNumber
!!!  use LIS_fileIOMod, only : LIS_read_param
!!!  use noah33_lsmMod     
!!!
!!!!
!!!! !DESCRIPTION:
!!!!  This subroutine resets the soil parameters in Noah3.3.  Noah3.3 uses a
!!!!  lookup table based on the input soil texture to derive these parameters. 
!!!! 
!!!!  The routines invoked are: 
!!!!  \begin{description}
!!!!   \item[LIS\_mapSoilType](\ref{LIS_mapSoilType}) \newline
!!!!     method to derive the soil texture type from the sand, silt, and 
!!!!     clay fractions.
!!!!   \end{description}
!!!!EOP      
!!!  implicit none
!!!  integer :: mtype
!!!  integer :: gid
!!!  integer :: jj,t,i,k,c,r,n
!!!!  real :: maxdt,mindz,tiledt,psis,ks,b
!!!!  integer :: newdt,idiv
!!!
!!!  real, allocatable :: smcmax(:,:)
!!!  real, allocatable :: psisat(:,:)
!!!  real, allocatable :: dwsat(:,:)
!!!  real, allocatable :: dksat(:,:)
!!!
!!!  integer           :: slcats
!!!  integer           :: iindex
!!!  character*4       :: sltype
!!!  real, allocatable     :: bb(:)
!!!  real, allocatable     :: drysmc(:)
!!!  real, allocatable     :: f11(:)
!!!  real, allocatable     :: maxsmc(:)
!!!  real, allocatable     :: refsmc(:)
!!!  real, allocatable     :: satpsi(:)
!!!  real, allocatable     :: satdk(:)
!!!  real, allocatable     :: satdw(:)
!!!  real, allocatable     :: wltsmc(:)
!!!  real, allocatable     :: qtz(:)
!!!  
!!!  real              :: frzfact
!!!  real, allocatable :: slope_data(:)
!!!  real ::  SBETA_DATA,FXEXP_DATA,CSOIL_DATA,SALP_DATA,REFDK_DATA,      &
!!!           REFKDT_DATA,FRZK_DATA,ZBOT_DATA,SMLOW_DATA,SMHIGH_DATA,     &
!!!           CZIL_DATA,LVCOEF_DATA
!!!  integer :: NUM_SLOPE
!!!  real, allocatable :: placeslopetype(:,:)
!!!  real              :: zsoil
!!!  real, allocatable :: dsoil(:,:)
!!!  real, allocatable :: droot(:,:)
!!!  integer       :: ftn
!!!
!!!  do n=1,LIS_rc%nnest
!!!
!!!     write(unit=LIS_logunit,fmt=*)                                     &
!!!                            'noah33_resetsoils -- reading soil files'
!!!
!!!     if(LIS_rc%usetexturemap(n).eq."none") then
!!!        do i=1,LIS_rc%npatch(n,mtype)
!!!           call LIS_mapSoilType(noah33_struc(n)%soilscheme,&
!!!                LIS_surface(n,mtype)%tile(i)%sand, &
!!!                LIS_surface(n,mtype)%tile(i)%clay, &
!!!                LIS_surface(n,mtype)%tile(i)%silt,&
!!!                LIS_surface(n,mtype)%tile(i)%soilt)
!!!        enddo
!!!     endif
!!!
!!!     if (noah33_struc(n)%fixedsoiltype.ne.0) then
!!!        write(LIS_logunit,*) 'Fixing Noah3.3 soil index to type: ',    &
!!!                              noah33_struc(n)%fixedsoiltype
!!!     endif
!!!     do i=1,LIS_rc%npatch(n,mtype)
!!!        if (noah33_struc(n)%fixedsoiltype.eq.0) then
!!!           noah33_struc(n)%noah(i)%soiltype = &
!!!                LIS_surface(n,mtype)%tile(i)%soilt
!!!        else
!!!           noah33_struc(n)%noah(i)%soiltype = noah33_struc(n)%fixedsoiltype
!!!        endif
!!!        !     if(noah33_struc(n)%noah(i)%soiltype.eq.-9999) noah33_struc(n)%noah(i)%soiltype = 14     
!!!        if (noah33_struc(n)%noah(i)%soiltype .eq. -9999.or.&
!!!             noah33_struc(n)%noah(i)%soiltype.eq.0) then
!!!           if ( noah33_struc(n)%soilscheme == 2 ) then ! statsgo
!!!              noah33_struc(n)%noah(i)%soiltype = 12 ! clay
!!!           else ! zobler
!!!              noah33_struc(n)%noah(i)%soiltype = 6 ! clay loam
!!!           endif
!!!        endif
!!!        if ( noah33_struc(n)%noah(i)%soiltype .eq. 14 .and. &
!!!             noah33_struc(n)%noah(i)%vegt .ne. LIS_rc%waterclass ) then
!!!           noah33_struc(n)%noah(i)%soiltype = 7
!!!        endif
!!!     enddo
!!!!Hardcoded for STATSGO classes only for soil types found to be water
!!!! at a land point
!!!
!!!     if (noah33_struc(n)%fixedslopetype.ne.0) then
!!!        write(LIS_logunit,*) 'Fixing Noah3.3 slope index to type: ',   &
!!!                              noah33_struc(n)%fixedslopetype
!!!     else
!!!        write(LIS_logunit,*) 'Reading Noah3.3 slope index from maps'
!!!        allocate(placeslopetype(LIS_rc%lnc(n),LIS_rc%lnr(n)))
!!!        call LIS_read_param(n,"SLOPETYPE",placeslopetype)
!!!     endif
!!!
!!!     if(noah33_struc(n)%usedsoilmap.ne.0) then 
!!!        allocate(dsoil(LIS_rc%lnc(n),LIS_rc%lnr(n)))
!!!
!!!        call LIS_read_param(n,"SOILDEPTH",dsoil)
!!!        
!!!        do t=1,LIS_rc%npatch(n,mtype)
!!!           allocate(noah33_struc(n)%noah(t)%lyrthk(noah33_struc(n)%nslay))
!!!           zsoil = 0.0
!!!           c = LIS_surface(n,mtype)%tile(t)%col
!!!           r = LIS_surface(n,mtype)%tile(t)%row
!!!
!!!           do k=1,noah33_struc(n)%nslay
!!!              if(k.ne.noah33_struc(n)%nslay) then 
!!!                 noah33_struc(n)%noah(t)%lyrthk(k) = noah33_struc(n)%lyrthk(k)
!!!                 zsoil = zsoil +  noah33_struc(n)%lyrthk(k)
!!!              else
!!!                 if((dsoil(c,r)-zsoil).gt.0) then
!!!                    noah33_struc(n)%noah(t)%lyrthk(k) = &
!!!                         dsoil(c,r) - zsoil
!!!                 else
!!!                    write(LIS_logunit,*) 'The total soil depth is shorter'
!!!                    write(LIS_logunit,*) 'than the layer structure. Please'
!!!                    write(LIS_logunit,*) 'resetup a thinner layer structure'
!!!                    call LIS_endrun()
!!!                 endif
!!!              endif
!!!           enddo
!!!        enddo
!!!        deallocate(dsoil)
!!!     endif
!!!
!!!     if(noah33_struc(n)%usedrootmap.ne.0) then 
!!!        allocate(droot(LIS_rc%lnc(n),LIS_rc%lnr(n)))
!!!
!!!        call LIS_read_param(n,"ROOTDEPTH",droot)
!!!        
!!!        do t=1,LIS_rc%npatch(n,mtype)
!!!
!!!           c = LIS_surface(n,mtype)%tile(t)%col
!!!           r = LIS_surface(n,mtype)%tile(t)%row
!!!           
!!!           zsoil = 0.0
!!!           noah33_struc(n)%noah(t)%nroot = 1
!!!           do k=1,noah33_struc(n)%nslay
!!!              zsoil = zsoil + noah33_struc(n)%noah(t)%lyrthk(k)
!!!              if(zsoil.ge.droot(c,r)) then 
!!!                 exit; 
!!!              else
!!!                 noah33_struc(n)%noah(t)%nroot = &
!!!                      noah33_struc(n)%noah(t)%nroot + 1
!!!              endif
!!!           enddo
!!!        enddo
!!!        deallocate(droot)
!!!     endif
!!!
!!!     do i=1,LIS_rc%npatch(n,mtype)
!!!        if (noah33_struc(n)%fixedslopetype.eq.0) then
!!!           if (placeslopetype(LIS_surface(n,mtype)%tile(i)%col,               &
!!!                LIS_surface(n,mtype)%tile(i)%row).gt.0) then
!!!              noah33_struc(n)%noah(i)%slopetype =                      &
!!!                       int(placeslopetype(LIS_surface(n,mtype)%tile(i)%col,   &
!!!                                          LIS_surface(n,mtype)%tile(i)%row))
!!!           else
!!!              noah33_struc(n)%noah(i)%slopetype = 3
!!!              write(LIS_logunit,*) 'Noah3.3 slope type -- ',           &
!!!                                   'not defined for point ',i
!!!              write(LIS_logunit,*) 'Noah3.3 slope type -- ',           &
!!!                                   'reset to default value of 3'
!!!           endif
!!!        else
!!!           noah33_struc(n)%noah(i)%slopetype = noah33_struc(n)%fixedslopetype
!!!        endif
!!!     enddo
!!!     if (noah33_struc(n)%fixedslopetype.eq.0) then
!!!        deallocate(placeslopetype)
!!!     endif
!!!
!!!!-----------------------------------------------------------------------
!!!! Read in the Noah3.3 Soil Parameter File
!!!!-----------------------------------------------------------------------
!!!     ftn = LIS_getNextUnitNumber()       
!!!     write(LIS_logunit,*) 'Reading Noah3.3 soil parameter file: ',     &
!!!                                             trim(noah33_struc(n)%sfile)
!!!     open(unit=ftn,file=noah33_struc(n)%sfile,status='old',             &
!!!                                                    access='sequential')
!!!     read(ftn,*) 
!!!!     read(ftn,*) (noah33_struc(n)%lyrthk(i), i=1,noah33_struc(n)%nslay)
!!!!     read(ftn,*)
!!!     read(ftn,fmt='(A4)') SLTYPE
!!!     read(ftn,*) SLCATS,IINDEX
!!!
!!!     allocate(bb(slcats))
!!!     allocate(drysmc(slcats))
!!!     allocate(f11(slcats))
!!!     allocate(maxsmc(slcats))
!!!     allocate(refsmc(slcats))
!!!     allocate(satpsi(slcats))
!!!     allocate(satdk(slcats))
!!!     allocate(satdw(slcats))
!!!     allocate(wltsmc(slcats))
!!!     allocate(qtz(slcats))
!!!
!!!     do i=1,SLCATS
!!!        read(ftn,*) IINDEX, bb(i), drysmc(i), f11(i), maxsmc(i), &
!!!             refsmc(i), satpsi(i), satdk(i), satdw(i), wltsmc(i), qtz(i)
!!!     enddo    
!!!     close(ftn)
!!!     
!!!!-----------------------------------------------------------------------
!!!! Assign SOIL Parameters to each tile based on the
!!!! type of Zobler soil class present in that tile.
!!!!-----------------------------------------------------------------------
!!!     do i = 1,LIS_rc%npatch(n,mtype)
!!!        noah33_struc(n)%noah(i)%smcwlt = &
!!!             wltsmc(noah33_struc(n)%noah(i)%soiltype)
!!!        noah33_struc(n)%noah(i)%smcdry = &
!!!             drysmc(noah33_struc(n)%noah(i)%soiltype)
!!!        noah33_struc(n)%noah(i)%smcref = &
!!!             refsmc(noah33_struc(n)%noah(i)%soiltype)
!!!        noah33_struc(n)%noah(i)%dwsat = &
!!!             satdw(noah33_struc(n)%noah(i)%soiltype)
!!!     end do     
!!!! Reset tile soil porosity
!!!     if (LIS_rc%useporositymap(n).eq."none") then ! default, from look-up table
!!!        do i = 1,LIS_rc%npatch(n,mtype)
!!!           noah33_struc(n)%noah(i)%smcmax = &
!!!                maxsmc(noah33_struc(n)%noah(i)%soiltype)
!!!        end do
!!!     else
!!!        do i = 1,LIS_rc%npatch(n,mtype)
!!!           noah33_struc(n)%noah(i)%smcmax = &
!!!                LIS_soils(n)%porosity(LIS_surface(n,mtype)%tile(i)%col,&
!!!                LIS_surface(n,mtype)%tile(i)%row,1)
!!!        end do
!!!     end if
!!!     
!!!! Reset tile soil Psi-sat
!!!     if (LIS_rc%usepsisatmap(n).eq."none") then ! default, from look-up table
!!!        do i = 1,LIS_rc%npatch(n,mtype)
!!!           noah33_struc(n)%noah(i)%psisat = &
!!!                satpsi(noah33_struc(n)%noah(i)%soiltype)
!!!        end do
!!!     else
!!!        do i = 1,LIS_rc%npatch(n,mtype)
!!!           noah33_struc(n)%noah(i)%psisat = &
!!!                LIS_soils(n)%psisat(LIS_surface(n,mtype)%tile(i)%col,&
!!!                LIS_surface(n,mtype)%tile(i)%row)
!!!        end do
!!!     end if
!!!
!!!! Reset tile soil K-sat
!!!     if (LIS_rc%useksatmap(n).eq."none") then ! default, from look-up table
!!!        do i = 1,LIS_rc%npatch(n,mtype)
!!!           noah33_struc(n)%noah(i)%dksat = &
!!!                satdk(noah33_struc(n)%noah(i)%soiltype)
!!!        end do
!!!     else
!!!        do i = 1,LIS_rc%npatch(n,mtype)
!!!           noah33_struc(n)%noah(i)%dksat = &
!!!                LIS_soils(n)%ksat(LIS_surface(n,mtype)%tile(i)%col,&
!!!                LIS_surface(n,mtype)%tile(i)%row)
!!!        end do
!!!     end if
!!!     
!!!! Reset tile soil b-parameter
!!!     if (LIS_rc%usebexpmap(n).eq."none") then ! default, from look-up table
!!!        do i = 1,LIS_rc%npatch(n,mtype)
!!!           noah33_struc(n)%noah(i)%bexp = &
!!!                bb(noah33_struc(n)%noah(i)%soiltype)
!!!        end do
!!!     else
!!!        do i = 1,LIS_rc%npatch(n,mtype)
!!!           noah33_struc(n)%noah(i)%bexp = &
!!!                LIS_soils(n)%bexp(LIS_surface(n,mtype)%tile(i)%col,&
!!!                LIS_surface(n,mtype)%tile(i)%row)
!!!        end do
!!!     end if
!!!     
!!!! Reset tile soil Qz content
!!!     if (LIS_rc%usequartzmap(n).eq."none") then ! default
!!!        do i = 1,LIS_rc%npatch(n,mtype)
!!!           noah33_struc(n)%noah(i)%quartz = &
!!!                qtz(noah33_struc(n)%noah(i)%soiltype)
!!!        end do
!!!     else 
!!!        do i = 1,LIS_rc%npatch(n,mtype)
!!!           noah33_struc(n)%noah(i)%quartz = &
!!!                LIS_soils(n)%quartz(LIS_surface(n,mtype)%tile(i)%col,&
!!!                LIS_surface(n,mtype)%tile(i)%row)
!!!        end do
!!!     end if
!!!     do i = 1,LIS_rc%npatch(n,mtype)
!!!        noah33_struc(n)%noah(i)%f1 = &
!!!             f11(noah33_struc(n)%noah(i)%soiltype)
!!!     enddo
!!!     
!!!     deallocate(bb)
!!!     deallocate(drysmc)
!!!     deallocate(f11)
!!!     deallocate(maxsmc)
!!!     deallocate(refsmc)
!!!     deallocate(satpsi)
!!!     deallocate(satdk)
!!!     deallocate(satdw)
!!!     deallocate(wltsmc)
!!!     deallocate(qtz)
!!!
!!!! These are parameters from the GENPARM.TBL
!!!     write(LIS_logunit,*) 'Reading Noah3.3 general parameter file: ',  &
!!!                           trim(noah33_struc(n)%gfile)
!!!     open(unit=ftn,file=noah33_struc(n)%gfile,status='old',             &
!!!                       access='sequential')
!!!     read(ftn,*)
!!!     read(ftn,*)
!!!     read(ftn,*) NUM_SLOPE
!!!     allocate(slope_data(NUM_SLOPE))
!!!     do jj=1,NUM_SLOPE
!!!        read(ftn,*) slope_data(jj)
!!!     enddo
!!!     read(ftn,*)
!!!     read(ftn,*)SBETA_DATA
!!!     read(ftn,*)
!!!     read(ftn,*)FXEXP_DATA
!!!     read(ftn,*)
!!!     read(ftn,*)CSOIL_DATA
!!!     read(ftn,*)
!!!     read(ftn,*)SALP_DATA
!!!     read(ftn,*)
!!!     read(ftn,*)REFDK_DATA
!!!     read(ftn,*)
!!!     read(ftn,*)REFKDT_DATA
!!!     read(ftn,*)
!!!     read(ftn,*)FRZK_DATA
!!!     read(ftn,*)
!!!     read(ftn,*)ZBOT_DATA
!!!     read(ftn,*)
!!!     read(ftn,*)CZIL_DATA
!!!     read(ftn,*)
!!!     read(ftn,*)SMLOW_DATA
!!!     read(ftn,*)
!!!     read(ftn,*)SMHIGH_DATA
!!!     read(ftn,*)
!!!     read(ftn,*)LVCOEF_DATA
!!!     close(ftn)
!!!     call LIS_releaseUnitNumber(ftn)
!!!        
!!!     do i = 1,LIS_rc%npatch(n,mtype)
!!!        if (LIS_rc%useslopemap(n).eq."none") then
!!!           noah33_struc(n)%noah(i)%slope =                             &
!!!                slope_data(noah33_struc(n)%noah(i)%slopetype)
!!!        else
!!!           gid = LIS_surface(n,mtype)%tile(i)%tile_id
!!!           noah33_struc(n)%noah(i)%slope = LIS_domain(n)%tile(gid)%slope
!!!        endif
!!!        noah33_struc(n)%noah(i)%sbeta = SBETA_DATA
!!!        noah33_struc(n)%noah(i)%fxexp = FXEXP_DATA
!!!        noah33_struc(n)%noah(i)%csoil = CSOIL_DATA
!!!        noah33_struc(n)%noah(i)%salp = SALP_DATA
!!!        noah33_struc(n)%noah(i)%refdk = REFDK_DATA
!!!        noah33_struc(n)%noah(i)%refkdt = REFKDT_DATA
!!!        noah33_struc(n)%noah(i)%frzk = FRZK_DATA
!!!        noah33_struc(n)%noah(i)%zbot = ZBOT_DATA
!!!        noah33_struc(n)%noah(i)%czil = CZIL_DATA
!!!        noah33_struc(n)%noah(i)%lvcoef = LVCOEF_DATA
!!!     enddo
!!!
!!!! overwrite default values with PTF 
!!!     if(noah33_struc(n)%useptf.eq.1) then !use cosby et al.(1984) PTFs
!!!        do i = 1,LIS_rc%npatch(n,mtype)
!!!           
!!!           noah33_struc(n)%noah(i)%smcmax = 0.489 - &
!!!                0.00126*LIS_surface(n,mtype)%tile(i)%sand*100.0
!!!
!!!           noah33_struc(n)%noah(i)%psisat = (10.0*(10.0**(1.88-0.0131*&
!!!                LIS_surface(n,mtype)%tile(i)%sand*100.0)))/1000.0
!!!
!!!           noah33_struc(n)%noah(i)%dksat =  (0.0070556*(10.0**(-0.884+0.0153*&
!!!                LIS_surface(n,mtype)%tile(i)%sand*100.0)))/1000.0
!!!
!!!           noah33_struc(n)%noah(i)%bexp = 2.91 + (0.159*&
!!!                LIS_surface(n,mtype)%tile(i)%clay*100.0)
!!!
!!!        end do
!!!     endif
!!!     
!!!     deallocate(slope_data)
!!!
!!!!begin hack
!!!#if 0 
!!!!     allocate(smcmax(LIS_rc%lnc(n),LIS_rc%lnr(n)))
!!!!     open(111,file='smcmax_calib.bin',form='unformatted')
!!!!     read(111) smcmax
!!!!     close(111)
!!!     
!!!!     do i=1,LIS_rc%npatch(n,mtype)
!!!!        noah33_struc(n)%noah(i)%smcmax = &
!!!!             smcmax(LIS_surface(n,mtype)%tile(i)%col, LIS_surface(n,mtype)%tile(i)%row)
!!!!     enddo
!!!!     deallocate(smcmax)
!!!
!!!     allocate(z0(LIS_rc%lnc(n),LIS_rc%lnr(n)))
!!!     open(111,file='z0_calib.bin',form='unformatted')
!!!     read(111) z0
!!!     close(111)
!!!     
!!!     do i=1,LIS_rc%npatch(n,mtype)
!!!        noah33_struc(n)%noah(i)%z0 = &
!!!             z0(LIS_surface(n,mtype)%tile(i)%col, LIS_surface(n,mtype)%tile(i)%row)
!!!     enddo
!!!     write(LIS_logunit,*) 'z0 ',z0
!!!     deallocate(z0)
!!!
!!!     allocate(rsmin(LIS_rc%lnc(n),LIS_rc%lnr(n)))
!!!     open(111,file='rsmin_calib.bin',form='unformatted')
!!!     read(111) rsmin
!!!     close(111)
!!!     
!!!     do i=1,LIS_rc%npatch(n,mtype)
!!!        noah33_struc(n)%noah(i)%rsmin = &
!!!             rsmin(LIS_surface(n,mtype)%tile(i)%col, LIS_surface(n,mtype)%tile(i)%row)
!!!     enddo
!!!     write(LIS_logunit,*) 'rsmin ',rsmin
!!!     deallocate(rsmin)
!!!
!!!     allocate(emiss(LIS_rc%lnc(n),LIS_rc%lnr(n)))
!!!     open(111,file='emiss_calib.bin',form='unformatted')
!!!     read(111) emiss
!!!     close(111)
!!!     
!!!     do i=1,LIS_rc%npatch(n,mtype)
!!!        noah33_struc(n)%noah(i)%emiss = &
!!!             emiss(LIS_surface(n,mtype)%tile(i)%col, LIS_surface(n,mtype)%tile(i)%row)
!!!     enddo
!!!     write(LIS_logunit,*) 'emiss ',emiss
!!!     deallocate(emiss)
!!!
!!!     allocate(czil(LIS_rc%lnc(n),LIS_rc%lnr(n)))
!!!     open(111,file='czil_calib.bin',form='unformatted')
!!!     read(111) czil
!!!     close(111)
!!!     
!!!     do i=1,LIS_rc%npatch(n,mtype)
!!!        noah33_struc(n)%noah(i)%czil = &
!!!             czil(LIS_surface(n,mtype)%tile(i)%col, LIS_surface(n,mtype)%tile(i)%row)
!!!     enddo
!!!     write(LIS_logunit,*) 'czil ',czil
!!!     deallocate(czil)
!!!
!!!     allocate(psisat(LIS_rc%lnc(n),LIS_rc%lnr(n)))
!!!     open(111,file='psisat_calib.bin',form='unformatted')
!!!     read(111) psisat
!!!     close(111)
!!!     
!!!     do i=1,LIS_rc%npatch(n,mtype)
!!!        noah33_struc(n)%noah(i)%psisat = &
!!!             psisat(LIS_surface(n,mtype)%tile(i)%col, LIS_surface(n,mtype)%tile(i)%row)
!!!     enddo
!!!     write(LIS_logunit,*) 'psisat ',psisat
!!!     deallocate(psisat)
!!!
!!!     allocate(dksat(LIS_rc%lnc(n),LIS_rc%lnr(n)))
!!!     open(111,file='dksat_calib.bin',form='unformatted')
!!!     read(111) dksat
!!!     close(111)
!!!     
!!!     do i=1,LIS_rc%npatch(n,mtype)
!!!        noah33_struc(n)%noah(i)%dksat = &
!!!             dksat(LIS_surface(n,mtype)%tile(i)%col, LIS_surface(n,mtype)%tile(i)%row)
!!!     enddo
!!!     write(LIS_logunit,*) 'dksat ',dksat
!!!     deallocate(dksat)
!!!
!!!     allocate(dwsat(LIS_rc%lnc(n),LIS_rc%lnr(n)))
!!!     open(111,file='dwsat_calib.bin',form='unformatted')
!!!     read(111) dwsat
!!!     close(111)
!!!     
!!!     do i=1,LIS_rc%npatch(n,mtype)
!!!        noah33_struc(n)%noah(i)%dwsat = &
!!!             dwsat(LIS_surface(n,mtype)%tile(i)%col, LIS_surface(n,mtype)%tile(i)%row)
!!!     enddo
!!!     write(LIS_logunit,*) 'dwsat ',dwsat
!!!     deallocate(dwsat)
!!!
!!!     allocate(quartz(LIS_rc%lnc(n),LIS_rc%lnr(n)))
!!!     open(111,file='quartz_calib.bin',form='unformatted')
!!!     read(111) quartz
!!!     close(111)
!!!     
!!!     do i=1,LIS_rc%npatch(n,mtype)
!!!        noah33_struc(n)%noah(i)%quartz = &
!!!             quartz(LIS_surface(n,mtype)%tile(i)%col, LIS_surface(n,mtype)%tile(i)%row)
!!!     enddo
!!!     write(LIS_logunit,*) 'quartz ',quartz
!!!     deallocate(quartz)
!!!
!!!     allocate(sbeta(LIS_rc%lnc(n),LIS_rc%lnr(n)))
!!!     open(111,file='sbeta_calib.bin',form='unformatted')
!!!     read(111) sbeta
!!!     close(111)
!!!     
!!!     do i=1,LIS_rc%npatch(n,mtype)
!!!        noah33_struc(n)%noah(i)%sbeta = &
!!!             sbeta(LIS_surface(n,mtype)%tile(i)%col, LIS_surface(n,mtype)%tile(i)%row)
!!!     enddo
!!!     write(LIS_logunit,*) 'sbeta ',sbeta
!!!     deallocate(sbeta)
!!!
!!!     allocate(fxexp(LIS_rc%lnc(n),LIS_rc%lnr(n)))
!!!     open(111,file='fxexp_calib.bin',form='unformatted')
!!!     read(111) fxexp
!!!     close(111)
!!!     
!!!     do i=1,LIS_rc%npatch(n,mtype)
!!!        noah33_struc(n)%noah(i)%fxexp = &
!!!             fxexp(LIS_surface(n,mtype)%tile(i)%col, LIS_surface(n,mtype)%tile(i)%row)
!!!
!!!        noah33_struc(n)%noah(i)%kdt = noah33_struc(n)%refkdt* noah33_struc(n)%noah(i)%dksat/&
!!!             noah33_struc(n)%refdk
!!!
!!!        frzfact = &
!!!             (noah33_struc(n)%noah(i)%smcmax/noah33_struc(n)%noah(i)%smcref)*(0.412/0.468)
!!!
!!!        noah33_struc(n)%noah(i)%frzx = noah33_struc(n)%frzk*frzfact 
!!!
!!!     enddo
!!!     write(LIS_logunit,*) 'fxexp ',fxexp
!!!     deallocate(fxexp)
!!!
!!!#endif
!!!!end hack
!!!
!!!!  call soils_finalize()
!!!
!!!!-----------------------------------------------------------------------
!!!! Check numerical stability at saturation on basis of timestep and 
!!!! soil parameters.  From formula derived by Yihua Wu, NASA/GSFC HSB
!!!!-----------------------------------------------------------------------
!!!!This routine as of now is "unsafe" in a parallel mode since each processor
!!!!can change the timestep. There needs to be some gather/scatter...
!!!#if 0 
!!!     maxdt = LIS_rc%ts
!!!     mindz = noahattrib%lyrthk(1) 
!!!     do i = 2,noah33_struc(n)%nslay ! soil layer loop
!!!        mindz = min(mindz,noahattrib%lyrthk(i))
!!!     end do ! i
!!!     do i = 1,LIS_rc%npatch(n,mtype)          ! tile loop
!!!! NOTE: depends on *positive* Psi-sat values, as in Noah3.3 LSM
!!!        psis = noah33_struc(n)%noah(i)%psisat
!!!        ks = noah33_struc(n)%noah(i)%dksat
!!!        b = noah33_struc(n)%noah(i)%bexp
!!!        tiledt = 2 * mindz**2 / (b * psis * ks)
!!!        maxdt = min(maxdt,tiledt) 
!!!     end do ! i
!!!     if (LIS_rc%ts.gt.maxdt) then
!!!        write(LIS_logunit,*) 'WRN: noah33_resetsoils -- timestep may',   &
!!!                             'be too large for numerical stability'
!!!        write(LIS_logunit,*) 'WRN: noah33_resetsoils -- soil parameters',&
!!!                             'suggest dt .le.',maxdt
!!!        newdt = 3600
!!!        idiv = 1
!!!        do while (newdt.gt.maxdt)
!!!           if (mod(3600,idiv).eq.0) then
!!!              newdt = 3600 / idiv
!!!           endif
!!!           idiv = idiv + 1
!!!           if (idiv.gt.60) exit
!!!        enddo
!!!        if (newdt.gt.maxdt) then
!!!           write(LIS_logunit,*) ' noah33_resetsoils -- new timestep',&
!!!                                'is less than 1 minute -- stopping'
!!!           call LIS_endrun
!!!        else
!!!           write(LIS_logunit,*) ' noah33_resetsoils -- setting',     &
!!!                                'timestep to ',newdt,' sec'
!!!           LIS_rc%ts = newdt
!!!        end if
!!!     else
!!!        write(LIS_logunit,*) ' noah33_resetsoils -- timestep OK',    &
!!!                      'for soil hydraulic stability (',LIS_rc%ts,' sec)'
!!!     end if
!!!#endif
!!!  enddo
!!!  return
!!!
end subroutine noah33_resetsoils
