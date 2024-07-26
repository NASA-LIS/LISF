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
! !ROUTINE: NoahMP36_setup
! \label{NoahMP36_setup}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   9/4/14: Shugong Wang; initial implementation for LIS 7 and NoahMP36
!   2/7/18: Soni Yatheendradas; code added for OPTUE to work 
!
! !INTERFACE:
subroutine NoahMP36_setup()
! !USES:
    use LIS_logMod,    only: LIS_logunit, LIS_verify, LIS_endrun
    use LIS_fileIOMod, only: LIS_read_param
    use LIS_coreMod,   only: LIS_rc, LIS_surface
    !use NOAHMP_VEG_PARAMETERS_36, only: read_mp_veg_parameters
    !use MODULE_SF_NOAHMPLSM_36, only: read_mp_veg_parameters
    use MODULE_SF_NOAHMPLSM_36, only: read_mp_veg_parameters, &
           SLCATS, LUCATS, CSOIL_DATA, BB, SATDK, SATDW, &
           SATPSI, QTZ, MAXSMC, REFSMC, WLTSMC, &
           CZIL_DATA, FRZK_DATA, REFDK_DATA, REFKDT_DATA, SLOPE_DATA, &
           TOPT_DATA, RGLTBL, RSMAX_DATA, RSTBL, HSTBL, NROTBL, &
           CH2OP, DLEAF, Z0MVT, HVT, HVB, RC, RHOL, RHOS, TAUL, TAUS, &
           XL, CWPVT, C3PSN, KC25, AKC, KO25, AKO, AVCMX, AQE, &         
           LTOVRC,  DILEFC,  DILEFW,  RMF25,  SLA,  FRAGR,  TMIN, &
           VCMX25,  TDLEF,  BP, MP, QE25, RMS25, RMR25, ARM, &
           FOLNMX, WDPOOL, WRRAT, MRP, DRYSMC ! SY: adding 
           ! calibratable parameters + DVEG to this list for OPTUE to work by 
           ! updating further below the corresponding values in NoahMP36_module 
           ! after call to subroutine SOIL_VEG_GEN_PARM_36   
           ! SY: Not used by NoahMP3.6 from REDPRM:  F11, DRYSMC            
           ! SY: Not used by NoahMP3.6 from read_mp_veg_parameters : AQE and 
           !     SLAREA (since subroutine BVOCFLUX not used), DEN
           ! SY: Also added DRYSMC used as constraint by OPTUE
 

    use NoahMP36_lsmMod

!
! !DESCRIPTION:
!
!  This routine is the entry point to set up the parameters
!  required for NoahMP36.  These include: 
!    vegetype     - land cover type index [-]
!    soiltype     - soil type index [-]
!    slopetype    - slope type for Noah baseflow [-]
!    tbot         - deep-layer soil temperature [K]
!    pblh         - planetary boundary layer height [m]
! 
!  The routines invoked are:
!  \begin{description}
!  \item[LIS\_read\_param](\ref{LIS_read_param}) \newline
!    retrieves LIS parameter data from NetCDF file
!  \item[NOAHMP36\_read\_MULTILEVEL\_param](\ref{NOAHMP36_read_MULTILEVEL_param}) \newline
!    retrieves MULTILEVEL spatial parameter from NetCDF file
!  \end{description}
!EOP
    implicit none
    integer           :: mtype
    integer           :: t, k, n
    integer           :: col, row
    real, allocatable :: placeholder(:,:)
    
    mtype = LIS_rc%lsm_index
    
    do n=1, LIS_rc%nnest
        ! allocate memory for place holder for #n nest
        allocate(placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n)))
        
        !------------------------------------!
        ! reading spatial spatial parameters !
        !------------------------------------!
        ! vegetype takes value from the LIS built-in parameter vegt
        !TODO: convert vegetation data source into vegetation types
        if(LIS_rc%uselcmap(n) .ne. 'none') then
            write(LIS_logunit,*) "NoahMP36: retrieve parameter VEGETYPE from LIS"
            do t=1, LIS_rc%npatch(n, mtype)
                NOAHMP36_struc(n)%noahmp36(t)%vegetype= LIS_surface(n, mtype)%tile(t)%vegt
            enddo
        else 
            ! read: vegetype
            write(LIS_logunit,*) "NoahMP36: reading parameter VEGETYPE from ", trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NOAHMP36_struc(n)%LDT_ncvar_vegetype), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NOAHMP36_struc(n)%noahmp36(t)%vegetype = placeholder(col, row)
            enddo 
        endif
        ! read: soiltype
        !write(LIS_logunit,*) "NoahMP36: reading parameter SOILTYPE from ", trim(LIS_rc%paramfile(n))
        !call LIS_read_param(n, trim(NOAHMP36_struc(n)%LDT_ncvar_soiltype), placeholder)
        !do t = 1, LIS_rc%npatch(n, mtype)
        !    col = LIS_surface(n, mtype)%tile(t)%col
        !    row = LIS_surface(n, mtype)%tile(t)%row
        !    NOAHMP36_struc(n)%noahmp36(t)%soiltype = placeholder(col, row)
        !enddo 

        ! soiltype takes value from the LIS built-in parameter soilt
        !TODO: convert soil texture into soil types according to scheme
        if(LIS_rc%usetexturemap(n) .ne. 'none') then
            write(LIS_logunit,*) "NoahMP36: retrieve parameter SOILTYPE from LIS"
            do t=1, LIS_rc%npatch(n, mtype)
                NOAHMP36_struc(n)%noahmp36(t)%soiltype= LIS_surface(n, mtype)%tile(t)%soilt
            enddo
        else 
            ! read: soiltype
            write(LIS_logunit,*) "NoahMP36: reading parameter SOILTYPE from ", trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(NOAHMP36_struc(n)%LDT_ncvar_soiltype), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NOAHMP36_struc(n)%noahmp36(t)%soiltype = placeholder(col, row)
            enddo 
        endif

        ! read: slopetype
        write(LIS_logunit,*) "NoahMP36: reading parameter SLOPETYPE from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(NOAHMP36_struc(n)%LDT_ncvar_slopetype), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            NOAHMP36_struc(n)%noahmp36(t)%slopetype = placeholder(col, row)
        enddo 

        ! read: tbot
        write(LIS_logunit,*) "NoahMP36: reading parameter TBOT from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(NOAHMP36_struc(n)%LDT_ncvar_tbot), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            NOAHMP36_struc(n)%noahmp36(t)%tbot = placeholder(col, row)
        enddo 

        ! read: pblh
        write(LIS_logunit,*) "NoahMP36: reading parameter PBLH from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(NOAHMP36_struc(n)%LDT_ncvar_pblh), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            NOAHMP36_struc(n)%noahmp36(t)%pblh = placeholder(col, row)
        enddo 

        !----------------------------------------------!
        ! MULTILEVEL reading spatial spatial parameters !
        !----------------------------------------------!
        ! read: shdfac_monthly
        write(LIS_logunit,*) "NoahMP36: reading parameter SHDFAC_MONTHLY from ", trim(LIS_rc%paramfile(n))
        do k = 1, 12
            call NOAHMP36_read_MULTILEVEL_param(n, NOAHMP36_struc(n)%LDT_ncvar_shdfac_monthly, k, placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                NOAHMP36_struc(n)%noahmp36(t)%shdfac_monthly(k) = placeholder(col, row)
            enddo 
        enddo 

        ! read: smceq for (opt_run=5)  Miguez-Macho & Fan groundwater with equilibrium water table
        if(NOAHMP36_struc(n)%run_opt .eq. 5) then
            write(LIS_logunit,*) "NoahMP36: reading parameter SMCEQ from ", trim(LIS_rc%paramfile(n))
            do k = 1, NOAHMP36_struc(n)%nsoil
                call NOAHMP36_read_MULTILEVEL_param(n, NOAHMP36_struc(n)%LDT_ncvar_smceq, k, placeholder)
                do t = 1, LIS_rc%npatch(n, mtype)
                    col = LIS_surface(n, mtype)%tile(t)%col
                    row = LIS_surface(n, mtype)%tile(t)%row
                    NOAHMP36_struc(n)%noahmp36(t)%smceq(k) = placeholder(col, row)
                enddo 
            enddo 
        endif

        deallocate(placeholder)
        call SOIL_VEG_GEN_PARM_36(NOAHMP36_struc(n)%landuse_tbl_name,   & 
                                  NOAHMP36_struc(n)%soil_tbl_name,      &
                                  NOAHMP36_struc(n)%gen_tbl_name,       &
                                  NOAHMP36_struc(n)%landuse_scheme_name,& 
                                  NOAHMP36_struc(n)%soil_scheme_name)
        ! SY: Begin for enabling OPTUE 
        do t = 1, LIS_rc%npatch(n, mtype)
            
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row

            ! SY: Begin lines following those in REDPRM

            ! SY: Begin SOIL PARAMETERS
            IF (NOAHMP36_struc(n)%noahmp36(t)%soiltype .gt. SLCATS) THEN
               write(LIS_logunit, *) 'SOILTYP must be less than SLCATS:'
               write(LIS_logunit, '("t = ", I6, "; SOILTYP = ", I6, ";    SLCATS = ", I6)') &
                         t, NOAHMP36_struc(n)%noahmp36(t)%soiltype, SLCATS
               write(LIS_logunit, *) 'NoahMP36_setup: Error: too many input soil types'
               write(LIS_logunit, *) 'program stopping ...'
               call LIS_endrun
            END IF
            NOAHMP36_struc(n)%noahmp36(t)%csoil = CSOIL_DATA
            NOAHMP36_struc(n)%noahmp36(t)%bexp = BB(NOAHMP36_struc(n)%noahmp36(t)%soiltype)
            NOAHMP36_struc(n)%noahmp36(t)%dksat = SATDK(NOAHMP36_struc(n)%noahmp36(t)%soiltype)
            NOAHMP36_struc(n)%noahmp36(t)%dwsat = SATDW(NOAHMP36_struc(n)%noahmp36(t)%soiltype)
            NOAHMP36_struc(n)%noahmp36(t)%psisat = SATPSI(NOAHMP36_struc(n)%noahmp36(t)%soiltype)
            NOAHMP36_struc(n)%noahmp36(t)%quartz = QTZ(NOAHMP36_struc(n)%noahmp36(t)%soiltype)
            NOAHMP36_struc(n)%noahmp36(t)%smcmax = MAXSMC(NOAHMP36_struc(n)%noahmp36(t)%soiltype)
            NOAHMP36_struc(n)%noahmp36(t)%smcref = REFSMC(NOAHMP36_struc(n)%noahmp36(t)%soiltype)
            NOAHMP36_struc(n)%noahmp36(t)%smcwlt = WLTSMC(NOAHMP36_struc(n)%noahmp36(t)%soiltype)
            ! SY: End SOIL PARAMETERS

            ! SY: Begin SOIL PARAMETER CONSTRAINT
            NOAHMP36_struc(n)%noahmp36(t)%smcdry = DRYSMC(NOAHMP36_struc(n)%noahmp36(t)%soiltype)
            ! SY: End SOIL PARAMETER CONSTRAINT

            ! SY: Begin UNIVERSAL PARAMETERS
            NOAHMP36_struc(n)%noahmp36(t)%czil = CZIL_DATA
            NOAHMP36_struc(n)%noahmp36(t)%frzk = FRZK_DATA
            NOAHMP36_struc(n)%noahmp36(t)%refdk = REFDK_DATA
            NOAHMP36_struc(n)%noahmp36(t)%refkdt = REFKDT_DATA
            NOAHMP36_struc(n)%noahmp36(t)%slope = SLOPE_DATA(NOAHMP36_struc(n)%noahmp36(t)%slopetype)
            ! SY: End UNIVERSAL PARAMETERS

            ! SY: Begin VEGETATION PARAMETERS
            IF (NOAHMP36_struc(n)%noahmp36(t)%vegetype .gt. LUCATS) THEN
               write(LIS_logunit, *) 'VEGTYP must be less than LUCATS:'
               write(LIS_logunit, '("t = ", I6, "; VEGTYP = ", I6, ";    LUCATS = ", I6)') &
                         t, NOAHMP36_struc(n)%noahmp36(t)%vegetype, LUCATS
               write(LIS_logunit, *) 'NoahMP36_setup: Error: too many input landuse types'
               write(LIS_logunit, *) 'program stopping ...'
               call LIS_endrun
            END IF
            NOAHMP36_struc(n)%noahmp36(t)%topt = TOPT_DATA
            NOAHMP36_struc(n)%noahmp36(t)%rgl = RGLTBL(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%rsmax = RSMAX_DATA
            NOAHMP36_struc(n)%noahmp36(t)%rsmin = RSTBL(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%hs = HSTBL(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%nroot = & 
                 NROTBL(NOAHMP36_struc(n)%noahmp36(t)%vegetype)

            IF (NROTBL(NOAHMP36_struc(n)%noahmp36(t)%vegetype) .gt. &
                NOAHMP36_struc(n)%nsoil) THEN
               WRITE (LIS_logunit,*) 'Warning: too many root layers'
               write (LIS_logunit,*) 'NROOT = ', NROTBL(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
               write (LIS_logunit,*) 'NSOIL = ', NOAHMP36_struc(n)%nsoil
               write(LIS_logunit, *) 'program stopping ...'
               call LIS_endrun
            END IF
            ! SY: End VEGETATION PARAMETERS

            ! SY: End lines following those in REDPRM

        enddo ! do t = 1, LIS_rc%npatch(n, mtype)
        ! SY: End for enabling OPTUE 

        call read_mp_veg_parameters(trim(NOAHMP36_struc(n)%noahmp_tbl_name), &
                                    trim(NOAHMP36_struc(n)%landuse_scheme_name))

        ! SY: Begin for enabling OPTUE 
        do t = 1, LIS_rc%npatch(n, mtype)
            
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row

            ! SY: Begin lines following those in read_mp_veg_parameters
            NOAHMP36_struc(n)%noahmp36(t)%CH2OP = CH2OP(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%DLEAF = DLEAF(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%Z0MVT = Z0MVT(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%HVT = HVT(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%HVB = HVB(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%RC = RC(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%RHOL1 = RHOL(NOAHMP36_struc(n)%noahmp36(t)%vegetype,1)
            NOAHMP36_struc(n)%noahmp36(t)%RHOL2 = RHOL(NOAHMP36_struc(n)%noahmp36(t)%vegetype,2)
            NOAHMP36_struc(n)%noahmp36(t)%RHOS1 = RHOS(NOAHMP36_struc(n)%noahmp36(t)%vegetype,1)
            NOAHMP36_struc(n)%noahmp36(t)%RHOS2 = RHOS(NOAHMP36_struc(n)%noahmp36(t)%vegetype,2)
            NOAHMP36_struc(n)%noahmp36(t)%TAUL1 = TAUL(NOAHMP36_struc(n)%noahmp36(t)%vegetype,1)
            NOAHMP36_struc(n)%noahmp36(t)%TAUL2 = TAUL(NOAHMP36_struc(n)%noahmp36(t)%vegetype,2)
            NOAHMP36_struc(n)%noahmp36(t)%TAUS1 = TAUS(NOAHMP36_struc(n)%noahmp36(t)%vegetype,1)
            NOAHMP36_struc(n)%noahmp36(t)%TAUS2 = TAUS(NOAHMP36_struc(n)%noahmp36(t)%vegetype,2)
            NOAHMP36_struc(n)%noahmp36(t)%XL = XL(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%CWPVT = CWPVT(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%C3PSN = C3PSN(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%KC25 = KC25(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%AKC = AKC(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%KO25 = KO25(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%AKO = AKO(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%AVCMX = AVCMX(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%AQE = AQE(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%LTOVRC = LTOVRC(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%DILEFC = DILEFC(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%DILEFW = DILEFW(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%RMF25 = RMF25(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%SLA = SLA(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%FRAGR = FRAGR(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%TMIN = TMIN(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%VCMX25 = VCMX25(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%TDLEF = TDLEF(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%BP = BP(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%MP = MP(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%QE25 = QE25(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%RMS25 = RMS25(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%RMR25 = RMR25(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%ARM = ARM(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%FOLNMX = FOLNMX(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%WDPOOL = WDPOOL(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%WRRAT = WRRAT(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            NOAHMP36_struc(n)%noahmp36(t)%MRP = MRP(NOAHMP36_struc(n)%noahmp36(t)%vegetype)
            ! SY: End lines following those in read_mp_veg_parameters

        enddo ! do t = 1, LIS_rc%npatch(n, mtype)
        ! SY: End for enabling OPTUE 
    enddo

end subroutine NoahMP36_setup

!BOP
!
! !ROUTINE: NOAHMP36_read_MULTILEVEL_param
!  \label{NOAHMP36_read_MULTILEVEL_param}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification for read_laiclimo
!  30 Oct  2013: Shugong Wang; Generalization for reading MULTILEVEL spatial parameter
!
! !INTERFACE:
subroutine NOAHMP36_read_MULTILEVEL_param(n, ncvar_name, level, placeholder)
! !USES:
    use netcdf
    use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_localPet,   &   
                            LIS_ews_halo_ind, LIS_ewe_halo_ind, &
                            LIS_nss_halo_ind, LIS_nse_halo_ind   
    use LIS_logMod,  only : LIS_logunit, LIS_verify, LIS_endrun
    use LIS_fileIOMod, only: LIS_read_param
    implicit none
! !ARGUMENTS: 
    integer, intent(in)          :: n
    integer, intent(in)          :: level
    character(len=*), intent(in) :: ncvar_name 
    real, intent(out)            :: placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n))
! !DESCRIPTION:
!  This subroutine reads MULTILEVEL parameters from the LIS
!  NetCDF parameter data file
!  
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[level]
!    level index (month, quarter, soil layer, snow layer) of the data to be read
!   \item[array]
!    array containing returned values
!   \end{description}
!
!EOP      

    integer       :: ios1
    integer       :: ios, nid, param_ID, nc_ID, nr_ID, dimids(3)
    integer       :: nc, nr, t, nlevel, k
    real, pointer :: level_data(:, :, :)
    logical       :: file_exists

    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then
        write(LIS_logunit, *) 'Reading '//trim(ncvar_name)//' map for level ', level

        ! open NetCDF parameter file
        ios = nf90_open(path=trim(LIS_rc%paramfile(n)), mode=NF90_NOWRITE, ncid=nid)
        call LIS_verify(ios, 'Error in nf90_open in NOAHMP36_read_MULTILEVEL_param')

        ! inquire the ID of east-west dimension
        ios = nf90_inq_dimid(nid, 'east_west', nc_ID)
        call LIS_verify(ios, 'Error in nf90_inq_dimid in NOAHMP36_read_MULTILEVEL_param')

        ! inquire the ID of north-south dimension
        ios = nf90_inq_dimid(nid, 'north_south', nr_ID)
        call LIS_verify(ios, 'Error in nf90_inq_dimid in NOAHMP36_read_MULTILEVEL_param')

        ! inquire the length of east-west dimension
        ios = nf90_inquire_dimension(nid, nc_ID, len=nc)
        call LIS_verify(ios, 'Error in nf90_inquire_dimension in NOAHMP36_read_MULTILEVEL_param')

        ! inquire the length of north-south dimension
        ios = nf90_inquire_dimension(nid, nr_ID, len=nr)
        call LIS_verify(ios, 'Error in nf90_inquire_dimension in NOAHMP36_read_MULTILEVEL_param')

        ! inquire the ID of parameter. 
        ios = nf90_inq_varid(nid, Trim(ncvar_name), param_ID)
        call LIS_verify(ios, trim(ncvar_name)//' field not found in the LIS param file')

        ! inquire the IDs of all dimensions. The third dimension is the level dimension
        ios = nf90_inquire_variable(nid, param_ID, dimids = dimids)
        call LIS_verify(ios, trim(ncvar_name)//' failed to inquire dimensions')

        ! inquire the length of the level dimension
        ios = nf90_inquire_dimension(nid, dimids(3), len=nlevel)
        call LIS_verify(ios, trim(ncvar_name)//' failed to inquire the length of the 3rd dimension')

        ! allocate memory
        allocate(level_data (LIS_rc%gnc(n), LIS_rc%gnr(n), nlevel))

        ! inquire the variable ID of parameter 
        ios = nf90_inq_varid(nid, trim(ncvar_name), param_ID)
        call LIS_verify(ios, trim(ncvar_name)//' field not found in the LIS param file')

        ! read parameter 
        ios = nf90_get_var(nid, param_ID, level_data)
        call LIS_verify(ios, 'Error in nf90_get_var in NOAHMP36_read_MULTILEVEL_param')

        ! close netcdf file 
        ios = nf90_close(nid)
        call LIS_verify(ios, 'Error in nf90_close in NOAHMP36_read_MULTILEVEL_param')

        ! grab parameter at specific level
        placeholder(:, :) = & 
             level_data(LIS_ews_halo_ind(n, LIS_localPet+1):LIS_ewe_halo_ind(n, LIS_localPet+1), &
                        LIS_nss_halo_ind(n, LIS_localPet+1):LIS_nse_halo_ind(n, LIS_localPet+1), level)

        ! free memory 
        deallocate(level_data)

    else
        write(LIS_logunit, *) 'MULTILEVEL parameter data file: ', LIS_rc%paramfile(n), ' does not exist'
        write(LIS_logunit, *) 'program stopping ...'
        call LIS_endrun
    endif
 end subroutine NOAHMP36_read_MULTILEVEL_param
                                          

