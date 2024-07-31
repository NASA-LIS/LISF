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
!
! !ROUTINE: cable_setvegparms
! \label{cable_setvegparms}
!
! !REVISION HISTORY:
!  25 Jul 2011: David Mocko, CABLE LSM implementation in LISv6.+
!  13 Sep 2011: Claire Carouge (ccc), CABLE LSM improvements
!  23 May 2013: David Mocko, latest CABLE v1.4b version for LIS6.2
!
! !INTERFACE:
subroutine cable_setvegparms
! !USES:
  use LIS_coreMod,        only : LIS_rc, LIS_surface
  use LIS_logMod,         only : LIS_logunit, LIS_endrun
  use LIS_vegDataMod,         only : LIS_lai
  use cable_dimensions,   only : ms,ncp,ncs
  use cable_lsmMod
!
! !DESCRIPTION:
!  This subroutine retrieves CABLE vegetation parameters.
!  The current implementation uses a table-based lookup based
!  on vegetation classes to initialize the following parameters.
!
!EOP
  implicit none
  
  character(len=70), dimension(:), allocatable :: veg_desc
  character(len=80) :: comments
  character(len=10) :: vegtypetmp
  character(len=25) :: vegnametmp
  integer :: n,t,tid,a,jveg,nvegt
  real :: notused,tdepth
  
  real, dimension(:),   allocatable :: vegin_canst1
  real, dimension(:),   allocatable :: vegin_dleaf
  real, dimension(:),   allocatable :: vegin_vcmax
  real, dimension(:),   allocatable :: vegin_ejmax
  real, dimension(:),   allocatable :: vegin_hc
  real, dimension(:),   allocatable :: vegin_xfang
  real, dimension(:),   allocatable :: vegin_rp20
  real, dimension(:),   allocatable :: vegin_rpcoef
  real, dimension(:),   allocatable :: vegin_rs20
  real, dimension(:),   allocatable :: vegin_shelrb
  real, dimension(:),   allocatable :: vegin_frac4
  real, dimension(:),   allocatable :: vegin_wai
  real, dimension(:),   allocatable :: vegin_vegcf
  real, dimension(:),   allocatable :: vegin_extkn
  real, dimension(:),   allocatable :: vegin_tminvj
  real, dimension(:),   allocatable :: vegin_tmaxvj
  real, dimension(:),   allocatable :: vegin_vbeta
  real, dimension(:),   allocatable :: vegin_rootbeta
  real, dimension(:,:), allocatable :: vegin_froot
  real, dimension(:,:), allocatable :: vegin_cplant
  real, dimension(:,:), allocatable :: vegin_csoil
  real, dimension(:,:), allocatable :: vegin_ratecp
  real, dimension(:,:), allocatable :: vegin_ratecs
  character*100                 :: lcmap

  do n=1,LIS_rc%nnest
     
     write(LIS_logunit,*)                                          &
          'MSG: cable_setvegparms -- reading vegetation files'
     
!-----------------------------------------------------------------------
! Set CABLE vegetation type at tile from the LIS domain
!-----------------------------------------------------------------------

     lcmap = LIS_rc%lcscheme

     if (cable_struc(n)%fixedvegtype.ne.0) then
        write(LIS_logunit,*) 'Fixing CABLE vegetation to type: ',  &
             cable_struc(n)%fixedvegtype
     endif
     
     do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        if (cable_struc(n)%fixedvegtype.eq.0) then
           cable_struc(n)%cable(t)%vegtype =                       &
                LIS_surface(n,LIS_rc%lsm_index)%tile(t)%vegt
! If using the modified IGBP vegetation map (which has 20 types),
! re-map types 18-20 to types within the standard 17 IGBP types. - dmm
! This kluge is temporary until the standard 17-type IGBP map is put
! into LIS.  For now, do a (modified) reverse procedure as shown here:
!     ftp://ftp.emc.ncep.noaa.gov/mmb/gcp/ldas/noahlsm/README
           if (lcmap.eq."MODIS") then
              if (cable_struc(n)%cable(t)%vegtype.eq.18)           &
                   cable_struc(n)%cable(t)%vegtype = 7
              if (cable_struc(n)%cable(t)%vegtype.ge.19)           &
                   cable_struc(n)%cable(t)%vegtype = 16
           endif
        else
           cable_struc(n)%cable(t)%vegtype =                       &
                cable_struc(n)%fixedvegtype
        endif
! map to MODIS classes if ECOCLIMAP is used. 
        if(lcmap.eq."ECOCLIMAP2") then 
           if(cable_struc(n)%cable(t)%vegtype.eq.1) & 
                cable_struc(n)%cable(t)%vegtype = 16
           if(cable_struc(n)%cable(t)%vegtype.eq.2) & 
                cable_struc(n)%cable(t)%vegtype = 16
           if(cable_struc(n)%cable(t)%vegtype.eq.3) & 
                cable_struc(n)%cable(t)%vegtype = 15
           if(cable_struc(n)%cable(t)%vegtype.eq.4) & 
                cable_struc(n)%cable(t)%vegtype = 4
           if(cable_struc(n)%cable(t)%vegtype.eq.5) & 
                cable_struc(n)%cable(t)%vegtype = 1
           if(cable_struc(n)%cable(t)%vegtype.eq.6) & 
                cable_struc(n)%cable(t)%vegtype = 2
           if(cable_struc(n)%cable(t)%vegtype.eq.7) & 
                cable_struc(n)%cable(t)%vegtype = 12
           if(cable_struc(n)%cable(t)%vegtype.eq.8) & 
                cable_struc(n)%cable(t)%vegtype = 12
           if(cable_struc(n)%cable(t)%vegtype.eq.9) & 
                cable_struc(n)%cable(t)%vegtype = 14
           if(cable_struc(n)%cable(t)%vegtype.eq.10) & 
                cable_struc(n)%cable(t)%vegtype = 10
           if(cable_struc(n)%cable(t)%vegtype.eq.11) & 
                cable_struc(n)%cable(t)%vegtype = 8
           if(cable_struc(n)%cable(t)%vegtype.eq.12) & 
                cable_struc(n)%cable(t)%vegtype = 14
           lcmap = "MODIS"
        endif
     enddo
     
!-----------------------------------------------------------------------
! Read in the CABLE Vegetation Parameter File
!-----------------------------------------------------------------------
     write(LIS_logunit,*)                                          &
          'Reading CABLE vegetation parameter file: ',  &
          trim(cable_struc(n)%vfile)
     open(unit=11,file=cable_struc(n)%vfile,status='old',          &
          access='sequential')
     
     if (lcmap.eq."MODIS") then
! Assume using IGBP vegetation types
        read(11,*) comments
        read(11,*) nvegt
     elseif (lcmap.eq."UMD") then
! Assume using UMD vegetation types?? - dmm
        read(11,*)
        read(11,*)
        read(11,*) nvegt    ! read # vegetation types
        read(11,*)
        read(11,*)
! Still need to check this code against CASA (13-type?) parameter file - dmm
        comments = 'UMD'
     else
        write(LIS_logunit,*) 'Not currently a valid vegetation',   &
             ' type map for CABLE -- stopping'
        call LIS_endrun
     endif
     write(LIS_logunit,*) 'CABLE Landuse type ',trim(comments),    &
          ' - found ',nvegt,' categories'
     ! Allocate memory for type-specific vegetation parameters:
     allocate(veg_desc(nvegt))
     allocate(vegin_canst1(nvegt),vegin_dleaf(nvegt))
     allocate(vegin_vcmax(nvegt),vegin_ejmax(nvegt))
     allocate(vegin_hc(nvegt),vegin_xfang(nvegt))
     allocate(vegin_rp20(nvegt),vegin_rpcoef(nvegt))
     allocate(vegin_rs20(nvegt),vegin_shelrb(nvegt))
     allocate(vegin_frac4(nvegt),vegin_wai(nvegt))
     allocate(vegin_vegcf(nvegt),vegin_extkn(nvegt))
     allocate(vegin_tminvj(nvegt),vegin_tmaxvj(nvegt))
     allocate(vegin_vbeta(nvegt),vegin_rootbeta(nvegt))
     allocate(vegin_froot(ms,nvegt),vegin_cplant(ncp,nvegt))
     allocate(vegin_csoil(ncs,nvegt),vegin_ratecp(ncp,nvegt))
     allocate(vegin_ratecs(ncs,nvegt))

     if (lcmap.eq."MODIS") then
        ! Added to read new format (BP dec 2007)
        do a = 1,nvegt
           read(11,*) jveg, vegtypetmp, vegnametmp
           if (jveg.gt.nvegt) then
              write(LIS_logunit,*) 'vegetation type index out of', &
                   'range in parameter file -- stopping'
              call LIS_endrun
           endif
           veg_desc(jveg) = vegnametmp
           read(11,*) vegin_hc(jveg),vegin_xfang(jveg),notused,    &
                vegin_dleaf(jveg)
           read(11,*) ! rholeaf not used
           read(11,*) ! tauleaf not used
           read(11,*) notused, notused, notused, notused
           ! rhosoil not used
           read(11,*) notused, vegin_wai(jveg),                    &
                vegin_canst1(jveg), vegin_shelrb(jveg),      &
                vegin_vegcf(jveg), vegin_extkn(jveg)
           read(11,*) vegin_vcmax(jveg), vegin_rp20(jveg),         &
                vegin_rpcoef(jveg), vegin_rs20(jveg)
           read(11,*) vegin_tminvj(jveg), vegin_tmaxvj(jveg),      &
                vegin_vbeta(jveg), vegin_rootbeta(jveg)
           read(11,*) vegin_cplant(1:ncp,jveg),                    &
                vegin_csoil(1:ncs,jveg)
           ! rates not currently set to vary with veg type
           read(11,*) vegin_ratecp(1:ncp,jveg),                    &
                vegin_ratecs(1:ncs,jveg)
        enddo
        close(11)
     elseif (lcmap.eq."AVHRR") then
        do a = 1,nvegt
           ! Read description of each veg type
           read(11,'(8X,A70)') veg_desc(a)
        enddo
        read(11,*)
        read(11,*)
        read(11,*) vegin_canst1
        read(11,*)
        read(11,*) vegin_dleaf
        read(11,*) vegin_vcmax
        read(11,*) vegin_hc
        read(11,*) vegin_xfang
        read(11,*) vegin_rp20
        read(11,*) vegin_rpcoef
        read(11,*) vegin_rs20
        read(11,*) vegin_shelrb
        read(11,*) vegin_frac4
        read(11,*) vegin_wai
        read(11,*) vegin_vegcf
        read(11,*) vegin_extkn
        read(11,*) vegin_tminvj
        read(11,*) vegin_tmaxvj
        read(11,*) vegin_vbeta
        read(11,*)
        read(11,*) vegin_rootbeta
        read(11,*) vegin_cplant(1,:)
        read(11,*) vegin_cplant(2,:)
        read(11,*) vegin_cplant(3,:)
        read(11,*) vegin_csoil(1,:)
        read(11,*) vegin_csoil(2,:)
        read(11,*)
        read(11,*) vegin_ratecp(:,1)
        ! Set ratecp to be the same for all veg types:
        vegin_ratecp(1,:)=vegin_ratecp(1,1)
        vegin_ratecp(2,:)=vegin_ratecp(2,1)
        vegin_ratecp(3,:)=vegin_ratecp(3,1)
        read(11,*)
        read(11,*) vegin_ratecs(:,1)
        vegin_ratecs(1,:)=vegin_ratecs(1,1)
        vegin_ratecs(2,:)=vegin_ratecs(2,1)
        close(11)
     endif

     cable_struc(n)%nvegt = nvegt

     !-----------------------------------------------------------------------
     ! Assign vegetation parameters to each tile based on
     ! the type of vegetation class present in that tile.
     !-----------------------------------------------------------------------
     do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        ! Parameters that are not spatially dependent
        ! Fix in-canopy turbulence scheme globally:
        cable_struc(n)%cable(t)%meth = 1
        ! Parameters that vary spatially by vegetation type
        cable_struc(n)%cable(t)%canst1 =                           &
             vegin_canst1(cable_struc(n)%cable(t)%vegtype)
        cable_struc(n)%cable(t)%dleaf =                            &
             vegin_dleaf(cable_struc(n)%cable(t)%vegtype)
        cable_struc(n)%cable(t)%vcmax =                            &
             vegin_vcmax(cable_struc(n)%cable(t)%vegtype)
        cable_struc(n)%cable(t)%ejmax = 2.0 *                      &
             vegin_vcmax(cable_struc(n)%cable(t)%vegtype)
        cable_struc(n)%cable(t)%hc =                               &
             vegin_hc(cable_struc(n)%cable(t)%vegtype)
        cable_struc(n)%cable(t)%xfang =                            &
             vegin_xfang(cable_struc(n)%cable(t)%vegtype)
        cable_struc(n)%cable(t)%vbeta =                            &
             vegin_vbeta(cable_struc(n)%cable(t)%vegtype)
        cable_struc(n)%cable(t)%rp20 =                             &
             vegin_rp20(cable_struc(n)%cable(t)%vegtype)
        cable_struc(n)%cable(t)%rpcoef =                           &
             vegin_rpcoef(cable_struc(n)%cable(t)%vegtype)
        cable_struc(n)%cable(t)%rs20 =                             &
             vegin_rs20(cable_struc(n)%cable(t)%vegtype)
        cable_struc(n)%cable(t)%shelrb =                           &
             vegin_shelrb(cable_struc(n)%cable(t)%vegtype)
        cable_struc(n)%cable(t)%wai =                              &
             vegin_wai(cable_struc(n)%cable(t)%vegtype)
        cable_struc(n)%cable(t)%vegcf =                            &
             vegin_vegcf(cable_struc(n)%cable(t)%vegtype)
        cable_struc(n)%cable(t)%extkn =                            &
             vegin_extkn(cable_struc(n)%cable(t)%vegtype)
        cable_struc(n)%cable(t)%tminvj =                           &
             vegin_tminvj(cable_struc(n)%cable(t)%vegtype)
        cable_struc(n)%cable(t)%tmaxvj =                           &
             vegin_tmaxvj(cable_struc(n)%cable(t)%vegtype)
        do a = 1,ncp
           cable_struc(n)%cable(t)%cplant_init(a) =                &
                vegin_cplant(a,cable_struc(n)%cable(t)%vegtype)
           cable_struc(n)%cable(t)%ratecp(a) =                     &
                vegin_ratecp(a,cable_struc(n)%cable(t)%vegtype)
        enddo
        do a = 1,ncs
           cable_struc(n)%cable(t)%csoil_init(a) =                 &
                vegin_csoil(a,cable_struc(n)%cable(t)%vegtype)
           cable_struc(n)%cable(t)%ratecs(a) =                     &
                vegin_ratecs(a,cable_struc(n)%cable(t)%vegtype)
        enddo
        tdepth = 0.0
        do a = 1,ms
           tdepth = tdepth + (cable_struc(n)%cable(t)%zse(a)*100.0)
           cable_struc(n)%cable(t)%froot(a) = min(1.0,1.0 -        &
                vegin_rootbeta(cable_struc(n)%cable(t)%vegtype)**tdepth)
        enddo
        do a = ms,2,-1
           cable_struc(n)%cable(t)%froot(a) =                      &
                cable_struc(n)%cable(t)%froot(a) - &
                cable_struc(n)%cable(t)%froot(a-1)
        enddo
        ! This kluge assigns the value of frac4 as a function of vegetation type.
        ! frac4 is read in as a f(vegetation type) for the 13-type CASA/UMD map.
        ! For the 17-type IGBP vegetation data text file, frac4 is not available
        ! and is instead normally read in from CCAM grid.  For now, just assign
        ! frac4 to zero, or to one for savannah/grassland types, similar to how
        ! these values are set in the 13-type CASA/UMD map. - dmm
        if (lcmap.eq."MODIS") then
           cable_struc(n)%cable(t)%frac4 = 0.0
           if ((cable_struc(n)%cable(t)%vegtype.ge.8).and.         &
                (cable_struc(n)%cable(t)%vegtype.le.10)) then
              cable_struc(n)%cable(t)%frac4 = 1.0
           endif
        endif
     enddo

     write(LIS_logunit,*) 'Reading CABLE LAI from LAI maps'
     if (LIS_rc%uselaimap(n).ne."none") then
        do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
           cable_struc(n)%cable(t)%lai = LIS_lai(n)%tlai(tid)
           ! KLUGE: CABLE doesn't work with LAI=0. Problem in cable_canopy.F90 with gbhu. ccc
	        if (abs(cable_struc(n)%cable(t)%lai) < 0.009) then
              cable_struc(n)%cable(t)%lai = 0.009
           endif
        enddo
     endif

     deallocate(veg_desc)
     deallocate(vegin_canst1,vegin_dleaf)
     deallocate(vegin_vcmax,vegin_ejmax)
     deallocate(vegin_hc,vegin_xfang)
     deallocate(vegin_rp20,vegin_rpcoef)
     deallocate(vegin_rs20,vegin_shelrb)
     deallocate(vegin_frac4,vegin_wai)
     deallocate(vegin_vegcf,vegin_extkn)
     deallocate(vegin_tminvj,vegin_tmaxvj)
     deallocate(vegin_vbeta,vegin_rootbeta)
     deallocate(vegin_froot,vegin_cplant)
     deallocate(vegin_csoil,vegin_ratecp)
     deallocate(vegin_ratecs)
  enddo

end subroutine cable_setvegparms
