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
! !ROUTINE: RDHM356_setup
! \label{RDHM356_setup}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   3/4/14: Shugong Wang; initial implementation for LIS 7 and RDHM356
!
! !INTERFACE:
subroutine RDHM356_setup()
! !USES:
    use LIS_logMod,    only: LIS_logunit, LIS_verify, LIS_endrun
    use LIS_fileIOMod, only: LIS_read_param
    use LIS_coreMod,   only: LIS_rc, LIS_surface
    use RDHM356_lsmMod

!
! !DESCRIPTION:
!
!  This routine is the entry point to set up the parameters
!  required for RDHM356.  These include:  \newline
!    SoilAlb      - snow free ALBEDO (default value 0.15) [-] \newline
!    SnowAlb      - snow ALBEDO (default value 0.7) [-] \newline
!    SOILTYP      - Soil type [-] \newline
!    VEGETYP      - Vegetation type [-] \newline
!    UZTWM        - upper zone tension water maximum storage [mm] \newline
!    UZFWM        - upper zone free water maximum storage [mm] \newline
!    UZK          - upper zone free water latent depletion rate [day$^{-1}$] \newline
!    PCTIM        - impervious fraction of the watershad area [-] \newline
!    ADIMP        - additional impervious area [-] \newline
!    RIVA         - riparian vegetation area [-] \newline
!    ZPERC        - maximum percolation rate [-] \newline
!    REXP         - exponent of the percolation equation (percolation parameter) [-] \newline
!    LZTWM        - lower zone tension water maximum storage [mm] \newline
!    LZFSM        - lower zone supplemental free water (fast) maximum storage [mm] \newline
!    LZFPM        - lower zone primary free water (slow) maximum storage [mm] \newline
!    LZSK         - lower zone supplemental free water depletion rate [day$^{-1}$] \newline
!    LZPK         - lower zone primary free water depletion rate [day$^{-1}$] \newline
!    PFREE        - fraction percolation from upper to lower free water storage [day$^{-1}$] \newline
!    SIDE         - ratio of deep recharge to channel base flow [-] \newline
!    RSERV        - fraction of lower zone free water not transferable to tension water [-] \newline
!    EFC          - fraction of forest cover [-] \newline
!    TBOT         - bottom boundary soil temperature [C] \newline
!    RSMAX        - maximum residual porosity (usually = 0.58) [-] \newline
!    CKSL         - ratio of frozen to non-frozen surface (increase in frozen ground contact, usually = 8 s/m) [s/m] \newline
!    ZBOT         - lower boundary depth (negative value, usually = -2.5 m) [m] \newline
!    vegRCMIN     - minimal stomatal resistance table for SACHTET, 14 values [s/m] \newline
!    climRCMIN    - climate dependent miminal stomatal resistance for SACHTET, 14 values [s/m] \newline
!    RGL          - solar radiation threshold table for SACHTET, 14 values [W m-2] \newline
!    HS           - vapor pressure resistance factor table for SACHTET, 14 values [-] \newline
!    LAI          - leaf area index table for SACHTET, 14 values [-] \newline
!    D50          - the depth (cm) table at which 50% roots are allocated for SACHTET, 14 values [cm] \newline
!    CROOT        - root distribution parameter table for SACHTET, 14 values [-] \newline
!    Z0           - roughness coefficient of surface [m] \newline
!    CLAY         - clay content for SACHTET, 12 values [-] \newline
!    SAND         - sand content for sACHTET, 12 values [-] \newline
!    SATDK        - saturated hydraulic conductivityfor SACHTET, 12 values [m s-1] \newline
!    ALON         - logitude [-] \newline
!    ALAT         - latitude [-] \newline
!    SCF          - snow fall correction factor [-] \newline
!    MFMAX        - maximum melt factor [mm/(6hrC)] \newline
!    MFMIN        - minimum melt factor [mm/(6hrC)] \newline
!    NMF          - maximum negative melt factor [mm/(6hrC)] \newline
!    UADJ         - the average wind function during rain-on-snow periods [mm/mb] \newline
!    SI           - areal water-equivalent above which 100 percent areal snow cover [mm] \newline
!    MBASE        - base temperature for non-rain melt factor [C] \newline
!    PXTEMP       - temperature which spereates rain from snow [C] \newline
!    PLWHC        - maximum amount of liquid-water held against gravity drainage [-] \newline
!    TIPM         - antecedent snow temperature index parameter [-] \newline
!    GM           - daily ground melt [mm/day] \newline
!    ELEV         - elevation [m] \newline
!    LAEC         - snow-rain split temperature [C] \newline
! 
!  The routines invoked are:
!  \begin{description}
!  \item[LIS\_read\_param](\ref{LIS_read_param}) \newline
!    retrieves LIS parameter data from NetCDF file
!  \item[RDHM356\_read\_multilevel\_param](\ref{RDHM356_read_multilevel_param}) \newline
!    retrieves multilevel spatial parameter from NetCDF file
!  \end{description}
!EOP
    implicit none
    integer           :: mtype
    integer           :: t, k, n
    integer           :: col, row
    real, allocatable :: placeholder(:,:)
    
    mtype = LIS_rc%lsm_index
    
    do n=1, LIS_rc%nnest
        ! allocate memory for place holder
        allocate(placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n)))
        !------------------------------------!
        ! reading spatial spatial parameters !
        !------------------------------------!
        ! read: SoilAlb
        write(LIS_logunit,*) "RDHM356: reading parameter SOILALB from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_SoilAlb), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%soilalb = placeholder(col, row)
        enddo 

        ! read: SnowAlb
        write(LIS_logunit,*) "RDHM356: reading parameter SNOWALB from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_SnowAlb), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%snowalb = placeholder(col, row)
        enddo 
        
        ! SOILTYP takes value from the LIS built-in parameter soilt
        !TODO: convert soil texture into soil types according to scheme
        if(LIS_rc%usetexturemap(n) .ne. 'none') then
            write(LIS_logunit,*) "RDHM356: retrieve parameter SOILTYP from LIS"
            do t=1, LIS_rc%npatch(n, mtype)
                RDHM356_struc(n)%rdhm356(t)%SOILTYP= LIS_surface(n, mtype)%tile(t)%soilt
            enddo
        else 
            ! read: SOILTYP
            write(LIS_logunit,*) "RDHM356: reading parameter SOILTYP from ", trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_SOILTYP), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                RDHM356_struc(n)%rdhm356(t)%soiltyp = placeholder(col, row)
            enddo 
        endif
        
        ! VEGETYP takes value from the LIS built-in parameter vegt
        !TODO: convert vegetation data source into vegetation types
        if(LIS_rc%uselcmap(n) .ne. 'none') then
            write(LIS_logunit,*) "RDHM356: retrieve parameter VEGETYP from LIS"
            do t=1, LIS_rc%npatch(n, mtype)
                RDHM356_struc(n)%rdhm356(t)%VEGETYP= LIS_surface(n, mtype)%tile(t)%vegt
            enddo
        else 
            ! read: VEGETYP
            write(LIS_logunit,*) "RDHM356: reading parameter VEGETYP from ", trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_VEGETYP), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                RDHM356_struc(n)%rdhm356(t)%vegetyp = placeholder(col, row)
            enddo 
        endif

        ! read: UZTWM
        write(LIS_logunit,*) "RDHM356: reading parameter UZTWM from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_UZTWM), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%uztwm = placeholder(col, row)
        enddo 

        ! read: UZFWM
        write(LIS_logunit,*) "RDHM356: reading parameter UZFWM from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_UZFWM), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%uzfwm = placeholder(col, row)
        enddo 

        ! read: UZK
        write(LIS_logunit,*) "RDHM356: reading parameter UZK from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_UZK), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%uzk = placeholder(col, row)
        enddo 

        ! read: PCTIM
        write(LIS_logunit,*) "RDHM356: reading parameter PCTIM from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_PCTIM), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%pctim = placeholder(col, row)
        enddo 

        ! read: ADIMP
        write(LIS_logunit,*) "RDHM356: reading parameter ADIMP from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_ADIMP), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%adimp = placeholder(col, row)
        enddo 

        ! read: RIVA
        write(LIS_logunit,*) "RDHM356: reading parameter RIVA from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_RIVA), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%riva = placeholder(col, row)
        enddo 

        ! read: ZPERC
        write(LIS_logunit,*) "RDHM356: reading parameter ZPERC from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_ZPERC), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%zperc = placeholder(col, row)
        enddo 

        ! read: REXP
        write(LIS_logunit,*) "RDHM356: reading parameter REXP from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_REXP), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%rexp = placeholder(col, row)
        enddo 

        ! read: LZTWM
        write(LIS_logunit,*) "RDHM356: reading parameter LZTWM from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_LZTWM), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%lztwm = placeholder(col, row)
        enddo 

        ! read: LZFSM
        write(LIS_logunit,*) "RDHM356: reading parameter LZFSM from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_LZFSM), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%lzfsm = placeholder(col, row)
        enddo 

        ! read: LZFPM
        write(LIS_logunit,*) "RDHM356: reading parameter LZFPM from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_LZFPM), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%lzfpm = placeholder(col, row)
        enddo 

        ! read: LZSK
        write(LIS_logunit,*) "RDHM356: reading parameter LZSK from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_LZSK), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%lzsk = placeholder(col, row)
        enddo 

        ! read: LZPK
        write(LIS_logunit,*) "RDHM356: reading parameter LZPK from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_LZPK), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%lzpk = placeholder(col, row)
        enddo 

        ! read: PFREE
        write(LIS_logunit,*) "RDHM356: reading parameter PFREE from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_PFREE), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%pfree = placeholder(col, row)
        enddo 

        ! read: SIDE
        write(LIS_logunit,*) "RDHM356: reading parameter SIDE from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_SIDE), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%side = placeholder(col, row)
        enddo 

        ! read: RSERV
        write(LIS_logunit,*) "RDHM356: reading parameter RSERV from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_RSERV), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%rserv = placeholder(col, row)
        enddo 

        ! read: EFC
        write(LIS_logunit,*) "RDHM356: reading parameter EFC from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_EFC), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%efc = placeholder(col, row)
        enddo 

        ! read: TBOT
        write(LIS_logunit,*) "RDHM356: reading parameter TBOT from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_TBOT), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%tbot = placeholder(col, row)
        enddo 

        ! read: RSMAX
        write(LIS_logunit,*) "RDHM356: reading parameter RSMAX from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_RSMAX), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%rsmax = placeholder(col, row)
        enddo 

        ! read: CKSL
        write(LIS_logunit,*) "RDHM356: reading parameter CKSL from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_CKSL), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%cksl = placeholder(col, row)
        enddo 

        ! read: ZBOT
        write(LIS_logunit,*) "RDHM356: reading parameter ZBOT from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_ZBOT), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            ! zbot is a negtive number in RDHM 
            RDHM356_struc(n)%rdhm356(t)%zbot = -1.0*placeholder(col, row)
        enddo 

        ! read: vegRCMIN
        write(LIS_logunit,*) "RDHM356: reading parameter VEGRCMIN from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_vegRCMIN), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%vegrcmin = placeholder(col, row)
        enddo 

        ! read: climRCMIN
        write(LIS_logunit,*) "RDHM356: reading parameter CLIMRCMIN from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_climRCMIN), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%climrcmin = placeholder(col, row)
        enddo 

        ! read: RGL
        write(LIS_logunit,*) "RDHM356: reading parameter RGL from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_RGL), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%rgl = placeholder(col, row)
        enddo 

        ! read: HS
        write(LIS_logunit,*) "RDHM356: reading parameter HS from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_HS), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%hs = placeholder(col, row)
        enddo 

        ! read: LAI
        write(LIS_logunit,*) "RDHM356: reading parameter LAI from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_LAI), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%lai = placeholder(col, row)
        enddo 

        ! read: D50
        write(LIS_logunit,*) "RDHM356: reading parameter D50 from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_D50), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%d50 = placeholder(col, row)
        enddo 

        ! read: CROOT
        write(LIS_logunit,*) "RDHM356: reading parameter CROOT from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_CROOT), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%croot = placeholder(col, row)
        enddo 

        ! read: Z0
        write(LIS_logunit,*) "RDHM356: reading parameter Z0 from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_Z0), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%z0 = placeholder(col, row)
        enddo 

        ! read: CLAY
        write(LIS_logunit,*) "RDHM356: reading parameter CLAY from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_CLAY), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%clay = placeholder(col, row)
        enddo 

        ! read: SAND
        write(LIS_logunit,*) "RDHM356: reading parameter SAND from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_SAND), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%sand = placeholder(col, row)
        enddo 

        ! read: SATDK
        write(LIS_logunit,*) "RDHM356: reading parameter SATDK from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_SATDK), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%satdk = placeholder(col, row)*1E-6 
        enddo 

        ! read: ALON
        write(LIS_logunit,*) "RDHM356: reading parameter ALON from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_ALON), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%alon = placeholder(col, row)
        enddo 

        ! read: ALAT
        write(LIS_logunit,*) "RDHM356: reading parameter ALAT from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_ALAT), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%alat = placeholder(col, row)
        enddo 

        ! read: SCF
        write(LIS_logunit,*) "RDHM356: reading parameter SCF from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_SCF), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%scf = placeholder(col, row)
        enddo 

        ! read: MFMAX
        write(LIS_logunit,*) "RDHM356: reading parameter MFMAX from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_MFMAX), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%mfmax = placeholder(col, row)
        enddo 

        ! read: MFMIN
        write(LIS_logunit,*) "RDHM356: reading parameter MFMIN from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_MFMIN), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%mfmin = placeholder(col, row)
        enddo 

        ! read: NMF
        write(LIS_logunit,*) "RDHM356: reading parameter NMF from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_NMF), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%nmf = placeholder(col, row)
        enddo 

        ! read: UADJ
        write(LIS_logunit,*) "RDHM356: reading parameter UADJ from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_UADJ), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%uadj = placeholder(col, row)
        enddo 

        ! read: SI
        write(LIS_logunit,*) "RDHM356: reading parameter SI from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_SI), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%si = placeholder(col, row)
        enddo 

        ! read: MBASE
        write(LIS_logunit,*) "RDHM356: reading parameter MBASE from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_MBASE), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%mbase = placeholder(col, row)
        enddo 

        ! read: PXTEMP
        write(LIS_logunit,*) "RDHM356: reading parameter PXTEMP from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_PXTEMP), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%pxtemp = placeholder(col, row)
        enddo 

        ! read: PLWHC
        write(LIS_logunit,*) "RDHM356: reading parameter PLWHC from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_PLWHC), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%plwhc = placeholder(col, row)
        enddo 

        ! read: TIPM
        write(LIS_logunit,*) "RDHM356: reading parameter TIPM from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_TIPM), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%tipm = placeholder(col, row)
        enddo 

        ! read: GM
        write(LIS_logunit,*) "RDHM356: reading parameter GM from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_GM), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%gm = placeholder(col, row)
        enddo 

        ! ELEV takes value from the LIS built-in parameter elev
        if(LIS_rc%useelevationmap(n) .ne. 'none') then
            write(LIS_logunit,*) "RDHM356: retrieve parameter ELEV from LIS"
            do t=1, LIS_rc%npatch(n, mtype)
                RDHM356_struc(n)%rdhm356(t)%ELEV= LIS_surface(n, mtype)%tile(t)%elev
            enddo
        else 
            ! read: ELEV
            write(LIS_logunit,*) "RDHM356: reading parameter ELEV from ", trim(LIS_rc%paramfile(n))
            call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_ELEV), placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                RDHM356_struc(n)%rdhm356(t)%elev = placeholder(col, row)
            enddo 
        endif
        ! read: LAEC
        write(LIS_logunit,*) "RDHM356: reading parameter LAEC from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(RDHM356_struc(n)%LDT_ncvar_LAEC), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%laec = placeholder(col, row)
        enddo 

        !----------------------------------------------!
        ! multilevel reading spatial spatial parameters !
        !----------------------------------------------!
        ! read: PET_MON
        write(LIS_logunit,*) "RDHM356: reading parameter PET_MON from ", trim(LIS_rc%paramfile(n))
        do k = 1, 12
            call RDHM356_read_multilevel_param(n, RDHM356_struc(n)%LDT_ncvar_PET_MON, k, placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                RDHM356_struc(n)%rdhm356(t)%pet_mon(k) = placeholder(col, row)
            enddo 
        enddo 

        ! read: GRN_MON
        write(LIS_logunit,*) "RDHM356: reading parameter GRN_MON from ", trim(LIS_rc%paramfile(n))
        do k = 1, 12
            call RDHM356_read_multilevel_param(n, RDHM356_struc(n)%LDT_ncvar_GRN_MON, k, placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                RDHM356_struc(n)%rdhm356(t)%grn_mon(k) = placeholder(col, row)
            enddo 
        enddo 

        ! read: ADC
        write(LIS_logunit,*) "RDHM356: reading parameter ADC from ", trim(LIS_rc%paramfile(n))
        do k = 1, 11
            call RDHM356_read_multilevel_param(n, RDHM356_struc(n)%LDT_ncvar_ADC, k, placeholder)
            do t = 1, LIS_rc%npatch(n, mtype)
                col = LIS_surface(n, mtype)%tile(t)%col
                row = LIS_surface(n, mtype)%tile(t)%row
                RDHM356_struc(n)%rdhm356(t)%adc(k) = placeholder(col, row)
            enddo 
        enddo 

        deallocate(placeholder)
    enddo
end subroutine RDHM356_setup

!BOP
!
! !ROUTINE: RDHM356_read_multilevel_param
!  \label{RDHM356_read_multilevel_param}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification for read_laiclimo
!  30 Oct  2013: Shugong Wang; Generalization for reading multilevel spatial parameter
!
! !INTERFACE:
subroutine RDHM356_read_multilevel_param(n, ncvar_name, level, placeholder)
! !USES:
    use esmf
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
!  This subroutine reads multilevel parameters from the LIS
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
    integer       :: nc, nr, t, nlevel, k, rc 
    real, allocatable :: level_data(:, :, :)
    logical       :: file_exists
    
    type(ESMF_VM) :: vm
    ! initalize ESMF framework  
    call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, &
                         defaultlogfilename="TimeIntervalEx.Log", &
                         logkindflag=ESMF_LOGKIND_MULTI, rc=rc)

    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then
        write(LIS_logunit, *) 'Reading '//trim(ncvar_name)//' map for level ', level
    
        ios = nf90_open(path=trim(LIS_rc%paramfile(n)), mode=NF90_NOWRITE, ncid=nid)
        call LIS_verify(ios, 'Error in nf90_open in RDHM356_read_multilevel_param')

        ios = nf90_inq_dimid(nid, 'east_west', nc_ID)
        call LIS_verify(ios, 'Error in nf90_inq_dimid in RDHM356_read_multilevel_param')

        ios = nf90_inq_dimid(nid, 'north_south', nr_ID)
        call LIS_verify(ios, 'Error in nf90_inq_dimid in RDHM356_read_multilevel_param')

        ios = nf90_inquire_dimension(nid, nc_ID, len=nc)
        call LIS_verify(ios, 'Error in nf90_inquire_dimension in RDHM356_read_multilevel_param')

        ios = nf90_inquire_dimension(nid, nr_ID, len=nr)
        call LIS_verify(ios, 'Error in nf90_inquire_dimension in RDHM356_read_multilevel_param')

        ios = nf90_inq_varid(nid, trim(ncvar_name), param_ID)
        call LIS_verify(ios, trim(ncvar_name)//' field not found in the LIS param file')
        
        ios = nf90_inquire_variable(nid, param_ID, dimids=dimids)
        call LIS_verify(ios, trim(ncvar_name)//' failed to inquire dimensions')

        ios = nf90_inquire_dimension(nid, dimids(3), len=nlevel)
        call LIS_verify(ios, trim(ncvar_name)//' failed to inquire the length of the 3rd dimension')

        allocate(level_data (LIS_rc%gnc(n), LIS_rc%gnr(n), nlevel))

        ios = nf90_get_var(nid, param_ID, level_data)
        call LIS_verify(ios, 'Error in nf90_get_var in RDHM356_read_multilevel_param')

        ios = nf90_close(nid)
        call LIS_verify(ios, 'Error in nf90_close in RDHM356_read_multilevel_param')

        placeholder(:, :) = & 
             level_data(LIS_ews_halo_ind(n, LIS_localPet+1):LIS_ewe_halo_ind(n, LIS_localPet+1), &
                        LIS_nss_halo_ind(n, LIS_localPet+1):LIS_nse_halo_ind(n, LIS_localPet+1), level)

        deallocate(level_data)

    else
        write(LIS_logunit, *) 'multilevel parameter data file: ', LIS_rc%paramfile(n), ' does not exist'
        write(LIS_logunit, *) 'program stopping ...'
        call LIS_endrun
    endif
 end subroutine RDHM356_read_multilevel_param
                                          

