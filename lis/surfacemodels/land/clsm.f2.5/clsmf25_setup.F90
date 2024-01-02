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
!
! !ROUTINE: clsmf25_setup
! \label{clsmf25_setup}
!
! !REVISION HISTORY:
! 16 Dec 2005: Sujay Kumar; Initial Code
! 23 Nov 2012: David Mocko, Updates for Catchment Fortuna-2.5
! 
! !INTERFACE:
subroutine clsmf25_setup()
! !USES:
  use LIS_coreMod,  only : LIS_rc,LIS_domain, LIS_surface
  use LIS_logMod,   only : LIS_logunit, LIS_endrun
  use LIS_soilsMod, only : LIS_soils
  use clsmf25_lsmMod
  use LIS_pluginIndices

! !DESCRIPTION: 
! 
!  This routine is the entry point to set up the parameters
!  required for the Catchment Fortuna-2.5 LSM.
!  These include the soils, greenness, albedo, etc.
!  
! The routines invoked are: 
! \begin{description}
! \item[clsmf25\_read\_land\_parameter](\ref{clsmf25_read_land_parameter}) \newline
!   reads specified 2d land parameter
! \item[clsmf25\_read\_land3d\_parameter](\ref{clsmf25_read_land3d_parameter}) \newline
!   reads specified 3d land parameter
! \end{description}
!EOP
  implicit none
  integer           :: n,k,c,r,kk
  real, allocatable :: pvalue(:,:)
  real, allocatable :: pvalue3d(:,:,:)
  

  write(LIS_logunit,*) '[INFO] Reading in CLSM F2.5 Parameters ... '

  do n = 1,LIS_rc%nnest

     ! Map LIS landcover to CLSM
     allocate(pvalue(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     pvalue = -9999.0

     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        if (LIS_rc%lcscheme.eq."UMD") then
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.0) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 10
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.1) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 3
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.2) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 1
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.3) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 3
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.4) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 2
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.5) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 2
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.6) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 2
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.7) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 4
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.8) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 5
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.9) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 5
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.10) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 4
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.11) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 4
     ! Types 7/8 no longer exist in Catchment; Map to type 5 instead. - dmm
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.12) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 5
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.13) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 5
           pvalue(LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col, &
                LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row) = &
                clsmf25_struc(n)%cat_param(k)%vegcls

        elseif(LIS_rc%lcscheme.eq."MOSAIC") then
           clsmf25_struc(n)%cat_param(k)%vegcls = &
                LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt

        elseif(LIS_rc%lcscheme.eq."IGBPNCEP") then  !yliu17: add support for MODIS land cover
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.1) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 3
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.2) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 1
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.3) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 3
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.4) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 2
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.5) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 2
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.6) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 5 !2
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.7) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 5 !4
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.8) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 4 !5
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.9) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 4 !5
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.10) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 4
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.11) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 4 !4
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.12) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 4 !5
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.13) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 5
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.eq.14) &
                clsmf25_struc(n)%cat_param(k)%vegcls = 4
           if (LIS_surface(n,LIS_rc%lsm_index)%tile(k)%vegt.ge.15) &  !yliu: assign all other MODIS land cover
                clsmf25_struc(n)%cat_param(k)%vegcls = 5              ! types to 5 for now
           pvalue(LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col, &
                LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row) = &
                clsmf25_struc(n)%cat_param(k)%vegcls
        else
            write(LIS_logunit,*) '[ERR] Land cover data ', LIS_rc%lcscheme, &
                                 '  is currently not supported with CLSM F2.5.'
            write(LIS_logunit,*) ' LIS run stopping ...'
            call LIS_endrun
        endif
     enddo
     clsmf25_struc(n)%catopen = 0

!B parameter
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"BEXP",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%bee = & 
             pvalue(c,r)
     enddo

!psisat
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"PSISAT",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%psis = & 
             pvalue(c,r)
    enddo
!porosity
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"POROSITY",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%poros = & 
             pvalue(c,r)
    enddo

!porosity
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"KSAT",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%cond = & 
             pvalue(c,r)
    enddo
!wpwet
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"WPWET",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%wpwet = & 
             pvalue(c,r)
    enddo
!dpth
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"BEDROCKDEPTH",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%dpth = & 
             pvalue(c,r)
    enddo

!ATAU
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"ATAUCLSM",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%atau = & 
             pvalue(c,r)
    enddo

!BTAU
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"BTAUCLSM",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%btau = & 
             pvalue(c,r)
    enddo

!GNU
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"GNUCLSM",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%gnu = & 
             pvalue(c,r)
    enddo
!ARS1
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"ARS1CLSM",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%ars1 = & 
             pvalue(c,r)
    enddo
!ARS2
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"ARS2CLSM",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%ars2 = & 
             pvalue(c,r)
    enddo
!ARS3
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"ARS3CLSM",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%ars3 = & 
             pvalue(c,r)
    enddo

!ARA1
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"ARA1CLSM",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%ara1 = & 
             pvalue(c,r)
    enddo

!ARA2
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"ARA2CLSM",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%ara2 = & 
             pvalue(c,r)
    enddo

!ARA3
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"ARA3CLSM",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%ara3 = & 
             pvalue(c,r)
    enddo
!ARA4
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"ARA4CLSM",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%ara4 = & 
             pvalue(c,r)
    enddo
!ARW1
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"ARW1CLSM",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%arw1 = & 
             pvalue(c,r)
    enddo

!ARW2
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"ARW2CLSM",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%arw2 = & 
             pvalue(c,r)
    enddo

!ARW3
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"ARW3CLSM",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%arw3 = & 
             pvalue(c,r)
    enddo

!ARW4
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"ARW4CLSM",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%arw4 = & 
             pvalue(c,r)
    enddo

!BF1
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"BF1CLSM",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%bf1 = & 
             pvalue(c,r)
    enddo
!BF2
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"BF2CLSM",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%bf2 = & 
             pvalue(c,r)
    enddo
!BF3
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"BF3CLSM",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%bf3 = & 
             pvalue(c,r)
    enddo
!TSA1
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"TSA1CLSM",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%tsa1 = & 
             pvalue(c,r)
    enddo
!TSA2
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"TSA2CLSM",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%tsa2 = & 
             pvalue(c,r)
    enddo
!TSB1
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"TSB1CLSM",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%tsb1 = & 
             pvalue(c,r)
    enddo
!TSB2
     pvalue = -9999.0
     call clsmf25_read_land_parameter(n,"TSB2CLSM",pvalue)
     do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
        r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
        clsmf25_struc(n)%cat_param(k)%tsb2 = & 
             pvalue(c,r)
    enddo
    deallocate(pvalue)
   
    call clsmf25_compute_land_parameters(n)

    if(clsmf25_struc(n)%usemodisalbflag.eq.1) then 
       allocate(pvalue3d(LIS_rc%lnc(n),LIS_rc%lnr(n),12))
       
       pvalue3d = -9999.0       
       call clsmf25_read_land3d_parameter(n,"ALBNIRDIFF",pvalue3d)
       do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
          r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
          do kk=1,12
             clsmf25_struc(n)%modis_alb_param(k,kk)%albnf = &
                  pvalue3d(c,r,kk)
          enddo
       enddo

       pvalue3d = -9999.0       
       call clsmf25_read_land3d_parameter(n,"ALBNIRDIR",pvalue3d)
       do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
          r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
          do kk=1,12
             clsmf25_struc(n)%modis_alb_param(k,kk)%albnr = &
                  pvalue3d(c,r,kk)
          enddo
       enddo

       pvalue3d = -9999.0       
       call clsmf25_read_land3d_parameter(n,"ALBVISDIR",pvalue3d)
       do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
          r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
          do kk=1,12
             clsmf25_struc(n)%modis_alb_param(k,kk)%albvr = &
                  pvalue3d(c,r,kk)
          enddo
       enddo
       pvalue3d = -9999.0       
       call clsmf25_read_land3d_parameter(n,"ALBVISDIFF",pvalue3d)
       do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          c = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%col
          r = LIS_surface(n,LIS_rc%lsm_index)%tile(k)%row
          do kk=1,12
             clsmf25_struc(n)%modis_alb_param(k,kk)%albvf = &
                  pvalue3d(c,r,kk)
          enddo

       enddo
       deallocate(pvalue3d)
    endif

    write(LIS_logunit,*) "[INFO] Done in Catchment Fortuna-2.5 setup"
 enddo

end subroutine clsmf25_setup


!BOP
!
! !ROUTINE: clsmf25_read_land_parameter
!  \label{clsmf25_read_land_parameter}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
  subroutine clsmf25_read_land_parameter(n,vname,pvalue)
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    use LIS_coreMod,        only : LIS_rc, LIS_localPet,&
         LIS_ews_ind, LIS_ewe_ind,&
         LIS_nss_ind, LIS_nse_ind, LIS_ews_halo_ind,LIS_ewe_halo_ind, &
         LIS_nss_halo_ind, LIS_nse_halo_ind
    use LIS_logMod,         only : LIS_logunit, LIS_getNextUnitNumber, &
         LIS_releaseUnitNumber, LIS_endrun, LIS_verify
    use LIS_fileIOMod
    
    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n
    character(len=*)    :: vname
    real                :: pvalue(LIS_rc%lnc(n),LIS_rc%lnr(n))
  
! !DESCRIPTION:
!  This subroutine reads the CLSM model parameter
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[locallc]
!    landlc for the region of interest
!   \end{description}
!
!EOP      

    integer          :: ios1
    integer          :: ios,nid,ntypesId, varid,ncId, nrId
    integer          :: nc,nr
    integer          :: c,r,t
    real             :: glb_pvalue(LIS_rc%gnc(n),LIS_rc%gnr(n))
    logical          :: file_exists
! _______________________________________________________

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    pvalue = -9999.0
    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then 
       
       write(LIS_logunit,*)'[INFO] Reading map of '//trim(vname)
       ios = nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in clsmf25_read_land_parameter')
       
       ios = nf90_inq_varid(nid,trim(vname),varid)
       call LIS_verify(ios,&
            trim(vname)//' field not found in the LIS param file')
       
       ios = nf90_get_var(nid,varid,glb_pvalue)
       call LIS_verify(ios,&
            'Error in nf90_get_var in clsmf25_read_land_parameter')
       
       pvalue(:,:) = glb_pvalue(&
            LIS_ews_halo_ind(n,LIS_localPet+1):&         
            LIS_ewe_halo_ind(n,LIS_localPet+1), &
            LIS_nss_halo_ind(n,LIS_localPet+1): &
            LIS_nse_halo_ind(n,LIS_localPet+1))
       
       ios = nf90_close(nid)
       call LIS_verify(ios,'Error in nf90_close in clsmf25_read_land_parameter')
       
    else
       write(LIS_logunit,*) "[ERR] ",trim(LIS_rc%paramfile(n)), &
            ' does not exist.'
       write(LIS_logunit,*) '[ERR] program stopping ...'
       call LIS_endrun
    endif
#endif

  end subroutine clsmf25_read_land_parameter

!BOP
!
! !ROUTINE: clsmf25_read_land3d_parameter
!  \label{clsmf25_read_land3d_parameter}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
  subroutine clsmf25_read_land3d_parameter(n,vname,pvalue)
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    use LIS_coreMod,        only : LIS_rc, LIS_localPet,&
         LIS_ews_ind, LIS_ewe_ind,&
         LIS_nss_ind, LIS_nse_ind, LIS_ews_halo_ind,LIS_ewe_halo_ind, &
         LIS_nss_halo_ind, LIS_nse_halo_ind
    use LIS_logMod,         only : LIS_logunit, LIS_getNextUnitNumber, &
         LIS_releaseUnitNumber, LIS_endrun, LIS_verify
    use LIS_fileIOMod
    
    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n
    character(len=*)    :: vname
    real                :: pvalue(LIS_rc%lnc(n),LIS_rc%lnr(n),12)
  
! !DESCRIPTION:
!  This subroutine reads the CLSM model parameter
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[locallc]
!    landlc for the region of interest
!   \end{description}
!
!EOP      

    integer          :: ios1
    integer          :: ios,nid,ntypesId, varid,ncId, nrId
    integer          :: nc,nr
    integer          :: c,r,t
    real             :: glb_pvalue(LIS_rc%gnc(n),LIS_rc%gnr(n),12)
    logical          :: file_exists

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    pvalue = -9999.0
    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then 
       
       write(LIS_logunit,*)'[INFO] Reading map of '//trim(vname)
       ios = nf90_open(path=LIS_rc%paramfile(n),&
            mode=NF90_NOWRITE,ncid=nid)
       call LIS_verify(ios,'Error in nf90_open in clsmf25_read_land3d_parameter')
       
       ios = nf90_inq_varid(nid,trim(vname),varid)
       call LIS_verify(ios,&
            trim(vname)//' field not found in the LIS param file')
       
       ios = nf90_get_var(nid,varid,glb_pvalue)
       call LIS_verify(ios,&
            'Error in nf90_get_var in clsmf25_read_land3d_parameter')
       
       pvalue(:,:,:) = glb_pvalue(&
            LIS_ews_halo_ind(n,LIS_localPet+1):&         
            LIS_ewe_halo_ind(n,LIS_localPet+1), &
            LIS_nss_halo_ind(n,LIS_localPet+1): &
            LIS_nse_halo_ind(n,LIS_localPet+1),:)
       
       ios = nf90_close(nid)
       call LIS_verify(ios,'Error in nf90_close in clsmf25_read_land3d_parameter')
       
    else
       write(LIS_logunit,*) "[ERR] ",trim(LIS_rc%paramfile(n)), &
            " does not exist."
       write(LIS_logunit,*) "[ERR] LIS run stopping ..."
       call LIS_endrun
    endif
#endif

  end subroutine clsmf25_read_land3d_parameter
