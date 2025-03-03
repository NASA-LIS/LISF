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
! !ROUTINE: noah36_setsoils
!  \label{noah36_setsoils}
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
!  30 Oct 2014: David Mocko, added Noah-3.6 into LIS-7
!
! !INTERFACE:
subroutine noah36_setsoils(mtype)
! !USES:
  use LIS_coreMod,   only : LIS_rc, LIS_domain, LIS_surface
  use LIS_soilsMod,  only : LIS_soils
  use LIS_logMod,    only : LIS_logunit, LIS_endrun, &
                            LIS_getNextUnitNumber, LIS_releaseUnitNumber
  use LIS_fileIOMod, only : LIS_read_param
  use noah36_lsmMod     

!
! !DESCRIPTION:
!  This subroutine sets the soil parameters in Noah-3.6.  Noah-3.6 uses a
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
!  real, allocatable     :: fxexp(:)
!  real, allocatable     :: sbeta(:)
  
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

  do n=1,LIS_rc%nnest

     write(unit=LIS_logunit,fmt=*)                                     &
                            '[INFO] noah36_setsoils -- reading soil files'

     if(LIS_rc%usetexturemap(n).eq."none") then
        do i=1,LIS_rc%npatch(n,mtype)
           call LIS_mapSoilType(noah36_struc(n)%soilscheme,&
                LIS_surface(n,mtype)%tile(i)%sand, &
                LIS_surface(n,mtype)%tile(i)%clay, &
                LIS_surface(n,mtype)%tile(i)%silt,&
                LIS_surface(n,mtype)%tile(i)%soilt)
        enddo
     endif

     if (noah36_struc(n)%fixedsoiltype.ne.0) then
        write(LIS_logunit,*) '[INFO] Fixing Noah-3.6 soil index to type: ',   &
                              noah36_struc(n)%fixedsoiltype
     endif
     do i=1,LIS_rc%npatch(n,mtype)
        if (noah36_struc(n)%fixedsoiltype.eq.0) then
           noah36_struc(n)%noah(i)%soiltype = &
                LIS_surface(n,mtype)%tile(i)%soilt
        else
           noah36_struc(n)%noah(i)%soiltype = noah36_struc(n)%fixedsoiltype
        endif
        !     if(noah36_struc(n)%noah(i)%soiltype.eq.-9999) noah36_struc(n)%noah(i)%soiltype = 14     
        if (noah36_struc(n)%noah(i)%soiltype .eq. -9999.or.&
             noah36_struc(n)%noah(i)%soiltype.eq.0) then
           if ( noah36_struc(n)%soilscheme == 2 ) then ! statsgo
              noah36_struc(n)%noah(i)%soiltype = 12 ! clay
           else ! zobler
              noah36_struc(n)%noah(i)%soiltype = 6 ! clay loam
           endif
        endif
        if ( noah36_struc(n)%noah(i)%soiltype .eq. 14 .and. &
             noah36_struc(n)%noah(i)%vegt .ne. LIS_rc%waterclass ) then
           noah36_struc(n)%noah(i)%soiltype = 7
        endif
     enddo
!Hardcoded for STATSGO classes only for soil types found to be water
! at a land point

     if (noah36_struc(n)%fixedslopetype.ne.0) then
        write(LIS_logunit,*) '[INFO] Fixing Noah-3.6 slope index to type: ',  &
                              noah36_struc(n)%fixedslopetype
     else
        write(LIS_logunit,*) '[INFO] Reading Noah-3.6 slope index from maps'
        allocate(placeslopetype(LIS_rc%lnc(n),LIS_rc%lnr(n)))
        call LIS_read_param(n,"SLOPETYPE",placeslopetype)
     endif

     if(noah36_struc(n)%usedsoilmap.ne.0) then 
        allocate(dsoil(LIS_rc%lnc(n),LIS_rc%lnr(n)))

        call LIS_read_param(n,"SOILDEPTH",dsoil)
        
        do t=1,LIS_rc%npatch(n,mtype)
           allocate(noah36_struc(n)%noah(t)%lyrthk(noah36_struc(n)%nslay))
           zsoil = 0.0
           c = LIS_surface(n,mtype)%tile(t)%col
           r = LIS_surface(n,mtype)%tile(t)%row

           do k=1,noah36_struc(n)%nslay
              if(k.ne.noah36_struc(n)%nslay) then 
                 noah36_struc(n)%noah(t)%lyrthk(k) = noah36_struc(n)%lyrthk(k)
                 zsoil = zsoil +  noah36_struc(n)%lyrthk(k)
              else
                 if((dsoil(c,r)-zsoil).gt.0) then
                    noah36_struc(n)%noah(t)%lyrthk(k) = &
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

     if(noah36_struc(n)%usedrootmap.ne.0) then 
        allocate(droot(LIS_rc%lnc(n),LIS_rc%lnr(n)))

        call LIS_read_param(n,"ROOTDEPTH",droot)
        
        do t=1,LIS_rc%npatch(n,mtype)

           c = LIS_surface(n,mtype)%tile(t)%col
           r = LIS_surface(n,mtype)%tile(t)%row
           
           zsoil = 0.0
           noah36_struc(n)%noah(t)%nroot = 1
           do k=1,noah36_struc(n)%nslay
              zsoil = zsoil + noah36_struc(n)%noah(t)%lyrthk(k)
              if(zsoil.ge.droot(c,r)) then 
                 exit; 
              else
                 noah36_struc(n)%noah(t)%nroot = &
                      noah36_struc(n)%noah(t)%nroot + 1
              endif
           enddo
        enddo
        deallocate(droot)
     endif

     do i=1,LIS_rc%npatch(n,mtype)
        if (noah36_struc(n)%fixedslopetype.eq.0) then
           if (placeslopetype(LIS_surface(n,mtype)%tile(i)%col,               &
                LIS_surface(n,mtype)%tile(i)%row).gt.0) then
              noah36_struc(n)%noah(i)%slopetype =                      &
                       int(placeslopetype(LIS_surface(n,mtype)%tile(i)%col,   &
                                          LIS_surface(n,mtype)%tile(i)%row))
           else
              noah36_struc(n)%noah(i)%slopetype = 3
              write(LIS_logunit,*) '[WARN] Noah-3.6 slope type -- ',          &
                                   '[WARN] not defined for point ',i
              write(LIS_logunit,*) '[WARN] Noah-3.6 slope type -- ',          &
                                   '[WARN] set to default value of 3'
           endif
        else
           noah36_struc(n)%noah(i)%slopetype = noah36_struc(n)%fixedslopetype
        endif
     enddo
     if (noah36_struc(n)%fixedslopetype.eq.0) then
        deallocate(placeslopetype)
     endif

!-----------------------------------------------------------------------
! Read in the Noah-3.6 Soil Parameter File
!-----------------------------------------------------------------------
     ftn = LIS_getNextUnitNumber()       
     write(LIS_logunit,*) '[INFO] Reading Noah-3.6 soil parameter file: ',    &
                                             trim(noah36_struc(n)%sfile)
     open(unit=ftn,file=noah36_struc(n)%sfile,status='old',            &
                                                    access='sequential')
     read(ftn,*) 
!     read(ftn,*) (noah36_struc(n)%lyrthk(i), i=1,noah36_struc(n)%nslay)
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
!     allocate(fxexp(slcats))
!     allocate(sbeta(slcats))

     do i=1,SLCATS
        read(ftn,*) IINDEX, bb(i), drysmc(i), f11(i), maxsmc(i), &
             refsmc(i), satpsi(i), satdk(i), satdw(i), wltsmc(i), qtz(i)
!             refsmc(i), satpsi(i), satdk(i), satdw(i), wltsmc(i), qtz(i),&
!             fxexp(i),sbeta(i)
     enddo    
     close(ftn)
!-----------------------------------------------------------------------
! Assign SOIL Parameters to each tile based on the
! type of Zobler soil class present in that tile.
!-----------------------------------------------------------------------
     do i = 1,LIS_rc%npatch(n,mtype)
        noah36_struc(n)%noah(i)%smcwlt = &
             wltsmc(noah36_struc(n)%noah(i)%soiltype)
        noah36_struc(n)%noah(i)%smcdry = &
             drysmc(noah36_struc(n)%noah(i)%soiltype)
        noah36_struc(n)%noah(i)%smcref = &
             refsmc(noah36_struc(n)%noah(i)%soiltype)
        noah36_struc(n)%noah(i)%dwsat = &
             satdw(noah36_struc(n)%noah(i)%soiltype)
     end do     
! Set tile soil porosity
     if (LIS_rc%useporositymap(n).eq."none") then ! default, from look-up table
        do i = 1,LIS_rc%npatch(n,mtype)
           noah36_struc(n)%noah(i)%smcmax = &
                maxsmc(noah36_struc(n)%noah(i)%soiltype)
        end do
     else
        do i = 1,LIS_rc%npatch(n,mtype)
           noah36_struc(n)%noah(i)%smcmax = &
                LIS_soils(n)%porosity(LIS_surface(n,mtype)%tile(i)%col,&
                LIS_surface(n,mtype)%tile(i)%row,1)
        end do
     end if
     
! Set tile soil Psi-sat
     if (LIS_rc%usepsisatmap(n).eq."none") then ! default, from look-up table
        do i = 1,LIS_rc%npatch(n,mtype)
           noah36_struc(n)%noah(i)%psisat = &
                satpsi(noah36_struc(n)%noah(i)%soiltype)
        end do
     else
        do i = 1,LIS_rc%npatch(n,mtype)
           noah36_struc(n)%noah(i)%psisat = &
                LIS_soils(n)%psisat(LIS_surface(n,mtype)%tile(i)%col,&
                LIS_surface(n,mtype)%tile(i)%row)
        end do
     end if

! Set tile soil K-sat
     if (LIS_rc%useksatmap(n).eq."none") then ! default, from look-up table
        do i = 1,LIS_rc%npatch(n,mtype)
           noah36_struc(n)%noah(i)%dksat = &
                satdk(noah36_struc(n)%noah(i)%soiltype)
        end do
     else
        do i = 1,LIS_rc%npatch(n,mtype)
           noah36_struc(n)%noah(i)%dksat = &
                LIS_soils(n)%ksat(LIS_surface(n,mtype)%tile(i)%col,&
                LIS_surface(n,mtype)%tile(i)%row)
        end do
     end if
     
! Set tile soil b-parameter
     if (LIS_rc%usebexpmap(n).eq."none") then ! default, from look-up table
        do i = 1,LIS_rc%npatch(n,mtype)
           noah36_struc(n)%noah(i)%bexp = &
                bb(noah36_struc(n)%noah(i)%soiltype)
        end do
     else
        do i = 1,LIS_rc%npatch(n,mtype)
           noah36_struc(n)%noah(i)%bexp = &
                LIS_soils(n)%bexp(LIS_surface(n,mtype)%tile(i)%col,&
                LIS_surface(n,mtype)%tile(i)%row)
        end do
     end if
     
! Set tile soil Qz content
     if (LIS_rc%usequartzmap(n).eq."none") then ! default
        do i = 1,LIS_rc%npatch(n,mtype)
           noah36_struc(n)%noah(i)%quartz = &
                qtz(noah36_struc(n)%noah(i)%soiltype)
        end do
     else 
        do i = 1,LIS_rc%npatch(n,mtype)
           noah36_struc(n)%noah(i)%quartz = &
                LIS_soils(n)%quartz(LIS_surface(n,mtype)%tile(i)%col,&
                LIS_surface(n,mtype)%tile(i)%row)
        end do
     end if
     do i = 1,LIS_rc%npatch(n,mtype)
        noah36_struc(n)%noah(i)%f1 = &
             f11(noah36_struc(n)%noah(i)%soiltype)
!        noah36_struc(n)%noah(i)%fxexp = &
!             fxexp(noah36_struc(n)%noah(i)%soiltype)
!        noah36_struc(n)%noah(i)%sbeta = &
!             sbeta(noah36_struc(n)%noah(i)%soiltype)
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
!     deallocate(fxexp)
!     deallocate(sbeta)

! These are parameters from the GENPARM.TBL
     write(LIS_logunit,*) '[INFO] Reading Noah-3.6 general parameter file: ',  &
                           trim(noah36_struc(n)%gfile)
     open(unit=ftn,file=noah36_struc(n)%gfile,status='old',             &
                       access='sequential')
     read(ftn,*)
     read(ftn,*)
     read(ftn,*) NUM_SLOPE
     allocate(slope_data(NUM_SLOPE))
     do jj=1,NUM_SLOPE
        read(ftn,*) slope_data(jj)
     enddo
     read(ftn,*)
     read(ftn,*)SBETA_DATA
     read(ftn,*)
     read(ftn,*)FXEXP_DATA
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
     read(ftn,*)
     read(ftn,*)CZIL_DATA
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
           noah36_struc(n)%noah(i)%slope =                             &
                slope_data(noah36_struc(n)%noah(i)%slopetype)
        else
           gid = LIS_surface(n,mtype)%tile(i)%tile_id
           noah36_struc(n)%noah(i)%slope = LIS_domain(n)%tile(gid)%slope
        endif
        noah36_struc(n)%noah(i)%sbeta = SBETA_DATA
        noah36_struc(n)%noah(i)%fxexp = FXEXP_DATA
        noah36_struc(n)%noah(i)%csoil = CSOIL_DATA
        noah36_struc(n)%noah(i)%salp = SALP_DATA
        noah36_struc(n)%noah(i)%refdk = REFDK_DATA
        noah36_struc(n)%noah(i)%refkdt = REFKDT_DATA
        noah36_struc(n)%noah(i)%frzk = FRZK_DATA
        noah36_struc(n)%noah(i)%zbot = ZBOT_DATA
        noah36_struc(n)%noah(i)%czil = CZIL_DATA
        noah36_struc(n)%noah(i)%lvcoef = LVCOEF_DATA
     enddo

! overwrite default values with PTF 
     if(noah36_struc(n)%useptf.eq.1) then !use cosby et al.(1984) PTFs
        do i = 1,LIS_rc%npatch(n,mtype)
           
           noah36_struc(n)%noah(i)%smcmax = 0.489 - &
                0.00126*LIS_surface(n,mtype)%tile(i)%sand*100.0

           noah36_struc(n)%noah(i)%psisat = (10.0*(10.0**(1.88-0.0131*&
                LIS_surface(n,mtype)%tile(i)%sand*100.0)))/1000.0

           noah36_struc(n)%noah(i)%dksat =  (0.0070556*(10.0**(-0.884+0.0153*&
                LIS_surface(n,mtype)%tile(i)%sand*100.0)))/1000.0

           noah36_struc(n)%noah(i)%bexp = 2.91 + (0.159*&
                LIS_surface(n,mtype)%tile(i)%clay*100.0)

        end do
     endif
     
     deallocate(slope_data)

  enddo
  return

end subroutine noah36_setsoils
