!------------------------------------------------------------------------------------------!
MODULE NOAHMP_VEG_PARAMETERS

    IMPLICIT NONE

    INTEGER, PARAMETER :: MAX_VEG_PARAMS = 33
    INTEGER, PARAMETER :: MVT   = 27
    INTEGER, PARAMETER :: MBAND = 2

    INTEGER, PRIVATE :: ISURBAN
    INTEGER :: ISWATER
    INTEGER :: ISBARREN
    INTEGER :: ISSNOW
    INTEGER :: EBLFOREST

    REAL :: CH2OP(MVT)       !maximum intercepted h2o per unit lai+sai (mm)
    REAL :: DLEAF(MVT)       !characteristic leaf dimension (m)
    REAL :: Z0MVT(MVT)       !momentum roughness length (m)
    REAL :: HVT(MVT)         !top of canopy (m)
    REAL :: HVB(MVT)         !bottom of canopy (m)
    REAL :: DEN(MVT)         !tree density (no. of trunks per m2)
    REAL :: RC(MVT)          !tree crown radius (m)
    REAL :: SAIM(MVT,12)     !monthly stem area index, one-sided
    REAL :: LAIM(MVT,12)     !monthly leaf area index, one-sided
    REAL :: SLA(MVT)         !single-side leaf area per Kg [m2/kg]
    REAL :: DILEFC(MVT)      !coeficient for leaf stress death [1/s]
    REAL :: DILEFW(MVT)      !coeficient for leaf stress death [1/s]
    REAL :: FRAGR(MVT)       !fraction of growth respiration  !original was 0.3 
    REAL :: LTOVRC(MVT)      !leaf turnover [1/s]

    REAL :: C3PSN(MVT)       !photosynthetic pathway: 0. = c4, 1. = c3
    REAL :: KC25(MVT)        !co2 michaelis-menten constant at 25c (pa)
    REAL :: AKC(MVT)         !q10 for kc25
    REAL :: KO25(MVT)        !o2 michaelis-menten constant at 25c (pa)
    REAL :: AKO(MVT)         !q10 for ko25
    REAL :: VCMX25(MVT)      !maximum rate of carboxylation at 25c (umol co2/m**2/s)
    REAL :: AVCMX(MVT)       !q10 for vcmx25
    REAL :: BP(MVT)          !minimum leaf conductance (umol/m**2/s)
    REAL :: MP(MVT)          !slope of conductance-to-photosynthesis relationship
    REAL :: QE25(MVT)        !quantum efficiency at 25c (umol co2 / umol photon)
    REAL :: AQE(MVT)         !q10 for qe25
    REAL :: RMF25(MVT)       !leaf maintenance respiration at 25c (umol co2/m**2/s)
    REAL :: RMS25(MVT)       !stem maintenance respiration at 25c (umol co2/kg bio/s)
    REAL :: RMR25(MVT)       !root maintenance respiration at 25c (umol co2/kg bio/s)
    REAL :: ARM(MVT)         !q10 for maintenance respiration
    REAL :: FOLNMX(MVT)      !foliage nitrogen concentration when f(n)=1 (%)
    REAL :: TMIN(MVT)        !minimum temperature for photosynthesis (k)

    REAL :: XL(MVT)          !leaf/stem orientation index
    REAL :: RHOL(MVT,MBAND)  !leaf reflectance: 1=vis, 2=nir
    REAL :: RHOS(MVT,MBAND)  !stem reflectance: 1=vis, 2=nir
    REAL :: TAUL(MVT,MBAND)  !leaf transmittance: 1=vis, 2=nir
    REAL :: TAUS(MVT,MBAND)  !stem transmittance: 1=vis, 2=nir

    REAL :: MRP(MVT)         !microbial respiration parameter (umol co2 /kg c/ s)
    REAL :: CWPVT(MVT)       !empirical canopy wind parameter

    REAL :: WRRAT(MVT)       !wood to non-wood ratio
    REAL :: WDPOOL(MVT)      !wood pool (switch 1 or 0) depending on woody or not [-]
    REAL :: TDLEF(MVT)       !characteristic T for leaf freezing [K]

    INTEGER :: IK,IM
    REAL :: TMP10(MVT*MBAND)
    REAL :: TMP11(MVT*MBAND)
    REAL :: TMP12(MVT*MBAND)
    REAL :: TMP13(MVT*MBAND)
    REAL :: TMP14(MVT*12)
    REAL :: TMP15(MVT*12)
    REAL :: TMP16(MVT*5)

    real slarea(MVT)
    real eps(MVT,5)

CONTAINS
  subroutine read_mp_veg_parameters(FILENAME_VEGTABLE,DATASET_IDENTIFIER)
    implicit none
    character(len=*), intent(in) :: FILENAME_VEGTABLE
    character(len=*), intent(in) :: DATASET_IDENTIFIER
    integer :: ierr

    ! Temporary arrays used in reshaping namelist arrays
    REAL :: TMP10(MVT*MBAND)
    REAL :: TMP11(MVT*MBAND)
    REAL :: TMP12(MVT*MBAND)
    REAL :: TMP13(MVT*MBAND)
    REAL :: TMP14(MVT*12)
    REAL :: TMP15(MVT*12)
    REAL :: TMP16(MVT*5)

    integer :: NVEG
    character(len=256) :: VEG_DATASET_DESCRIPTION

    NAMELIST / noah_mp_usgs_veg_categories / VEG_DATASET_DESCRIPTION, NVEG
    NAMELIST / noah_mp_usgs_parameters / ISURBAN, ISWATER, ISBARREN, ISSNOW, EBLFOREST, &
         CH2OP, DLEAF, Z0MVT, HVT, HVB, DEN, RC, RHOL,  RHOS, TAUL, TAUS, XL, CWPVT, C3PSN, KC25, AKC, KO25, AKO, AVCMX, AQE, &
         LTOVRC,  DILEFC,  DILEFW,  RMF25 ,  SLA   ,  FRAGR ,  TMIN  ,  VCMX25,  TDLEF ,  BP, MP, QE25, RMS25, RMR25, ARM, FOLNMX, WDPOOL, WRRAT, MRP,   &
         SAIM,  LAIM,  SLAREA, EPS

    NAMELIST / noah_mp_modis_veg_categories / VEG_DATASET_DESCRIPTION, NVEG
    NAMELIST / noah_mp_modis_parameters / ISURBAN, ISWATER, ISBARREN, ISSNOW, EBLFOREST, &
         CH2OP, DLEAF, Z0MVT, HVT, HVB, DEN, RC, RHOL,  RHOS, TAUL, TAUS, XL, CWPVT, C3PSN, KC25, AKC, KO25, AKO, AVCMX, AQE, &
         LTOVRC,  DILEFC,  DILEFW,  RMF25 ,  SLA   ,  FRAGR ,  TMIN  ,  VCMX25,  TDLEF ,  BP, MP, QE25, RMS25, RMR25, ARM, FOLNMX, WDPOOL, WRRAT, MRP,   &
         SAIM,  LAIM,  SLAREA, EPS

    ! MPC change: enable use of alternative veg tables
    !  - in this case, tables using attributes from other models used in the PLUMBER experiment

    NAMELIST / noah_mp_plumberCABLE_veg_categories / VEG_DATASET_DESCRIPTION, NVEG
    NAMELIST / noah_mp_plumberCABLE_parameters / ISURBAN, ISWATER, ISBARREN, ISSNOW, EBLFOREST, &
         CH2OP, DLEAF, Z0MVT, HVT, HVB, DEN, RC, RHOL,  RHOS, TAUL, TAUS, XL, CWPVT, C3PSN, KC25, AKC, KO25, AKO, AVCMX, AQE, &
         LTOVRC,  DILEFC,  DILEFW,  RMF25 ,  SLA   ,  FRAGR ,  TMIN  ,  VCMX25,  TDLEF ,  BP, MP, QE25, RMS25, RMR25, ARM, FOLNMX, WDPOOL, WRRAT, MRP,   &
         SAIM,  LAIM,  SLAREA, EPS

    NAMELIST / noah_mp_plumberCHTESSEL_veg_categories / VEG_DATASET_DESCRIPTION, NVEG
    NAMELIST / noah_mp_plumberCHTESSEL_parameters / ISURBAN, ISWATER, ISBARREN, ISSNOW, EBLFOREST, &
         CH2OP, DLEAF, Z0MVT, HVT, HVB, DEN, RC, RHOL,  RHOS, TAUL, TAUS, XL, CWPVT, C3PSN, KC25, AKC, KO25, AKO, AVCMX, AQE, &
         LTOVRC,  DILEFC,  DILEFW,  RMF25 ,  SLA   ,  FRAGR ,  TMIN  ,  VCMX25,  TDLEF ,  BP, MP, QE25, RMS25, RMR25, ARM, FOLNMX, WDPOOL, WRRAT, MRP,   &
         SAIM,  LAIM,  SLAREA, EPS

    NAMELIST / noah_mp_plumberSUMMA_veg_categories / VEG_DATASET_DESCRIPTION, NVEG
    NAMELIST / noah_mp_plumberSUMMA_parameters / ISURBAN, ISWATER, ISBARREN, ISSNOW, EBLFOREST, &
         CH2OP, DLEAF, Z0MVT, HVT, HVB, DEN, RC, RHOL,  RHOS, TAUL, TAUS, XL, CWPVT, C3PSN, KC25, AKC, KO25, AKO, AVCMX, AQE, &
         LTOVRC,  DILEFC,  DILEFW,  RMF25 ,  SLA   ,  FRAGR ,  TMIN  ,  VCMX25,  TDLEF ,  BP, MP, QE25, RMS25, RMR25, ARM, FOLNMX, WDPOOL, WRRAT, MRP,   &
         SAIM,  LAIM,  SLAREA, EPS

    ! Initialize our variables to bad values, so that if the namelist read fails, we come to a screeching halt as soon as we try to use anything.
    CH2OP  = -1.E36
    DLEAF  = -1.E36
    Z0MVT  = -1.E36
    HVT    = -1.E36
    HVB    = -1.E36
    DEN    = -1.E36
    RC     = -1.E36
    RHOL   = -1.E36
    RHOS   = -1.E36
    TAUL   = -1.E36
    TAUS   = -1.E36
    XL     = -1.E36
    CWPVT  = -1.E36
    C3PSN  = -1.E36
    KC25   = -1.E36
    AKC    = -1.E36
    KO25   = -1.E36
    AKO    = -1.E36
    AVCMX  = -1.E36
    AQE    = -1.E36
    LTOVRC = -1.E36
    DILEFC = -1.E36
    DILEFW = -1.E36
    RMF25  = -1.E36
    SLA    = -1.E36
    FRAGR  = -1.E36
    TMIN   = -1.E36
    VCMX25 = -1.E36
    TDLEF  = -1.E36
    BP     = -1.E36
    MP     = -1.E36
    QE25   = -1.E36
    RMS25  = -1.E36
    RMR25  = -1.E36
    ARM    = -1.E36
    FOLNMX = -1.E36
    WDPOOL = -1.E36
    WRRAT  = -1.E36
    MRP    = -1.E36
    SAIM   = -1.E36
    LAIM   = -1.E36
    SLAREA = -1.E36
    EPS    = -1.E36

    open(15, file=trim(FILENAME_VEGTABLE), status='old', form='formatted', action='read', iostat=ierr)
    if (ierr /= 0) then
       write(*,'("****** Error ******************************************************")')
       write(*,'("Cannot find file MPTABLE.TBL")')
       write(*,'("STOP")')
       write(*,'("*******************************************************************")')
       call wrf_error_fatal("STOP in Noah-MP read_mp_veg_parameters")
    endif

    if ( trim(DATASET_IDENTIFIER) == "USGS" ) then
       read(15,noah_mp_usgs_veg_categories)
       read(15,noah_mp_usgs_parameters)
    else if ( trim(DATASET_IDENTIFIER) == "MODIFIED_IGBP_MODIS_NOAH" ) then
       read(15,noah_mp_modis_veg_categories)
       read(15,noah_mp_modis_parameters)
    else if ( trim(DATASET_IDENTIFIER) == "plumberCABLE" ) then
       read(15,noah_mp_plumberCABLE_veg_categories)
       read(15,noah_mp_plumberCABLE_parameters)
    else if ( trim(DATASET_IDENTIFIER) == "plumberCHTESSEL" ) then
       read(15,noah_mp_plumberCHTESSEL_veg_categories)
       read(15,noah_mp_plumberCHTESSEL_parameters)
    else if ( trim(DATASET_IDENTIFIER) == "plumberSUMMA" ) then
       read(15,noah_mp_plumberSUMMA_veg_categories)
       read(15,noah_mp_plumberSUMMA_parameters)
    else
       write(*,'("Unrecognized DATASET_IDENTIFIER in subroutine READ_MP_VEG_PARAMETERS")')
       write(*,'("DATASET_IDENTIFIER = ''", A, "''")') trim(DATASET_IDENTIFIER)
       call wrf_error_fatal("STOP in Noah-MP read_mp_veg_parameters")
    endif
    close(15)

    ! Problem.  Namelist reading of 2-d arrays doesn't work well when the arrays are declared with larger dimension than the
    ! variables in the provided namelist.  So we need to reshape the 2-d arrays after we've read them.

    if ( MVT > NVEG ) then

       ! 
       ! Reshape the 2-d arrays:
       ! 

       TMP10 = reshape( RHOL, (/ MVT*size(RHOL,2) /))
       TMP11 = reshape( RHOS, (/ MVT*size(RHOS,2) /))
       TMP12 = reshape( TAUL, (/ MVT*size(TAUL,2) /))
       TMP13 = reshape( TAUS, (/ MVT*size(TAUS,2) /))
       TMP14 = reshape( SAIM, (/ MVT*size(SAIM,2) /))
       TMP15 = reshape( LAIM, (/ MVT*size(LAIM,2) /))
       TMP16 = reshape( EPS,  (/ MVT*size(EPS ,2) /))

       RHOL(1:NVEG,:) = reshape( TMP10, (/ NVEG, size(RHOL,2) /))
       RHOS(1:NVEG,:) = reshape( TMP11, (/ NVEG, size(RHOS,2) /))
       TAUL(1:NVEG,:) = reshape( TMP12, (/ NVEG, size(TAUL,2) /))
       TAUS(1:NVEG,:) = reshape( TMP13, (/ NVEG, size(TAUS,2) /))
       SAIM(1:NVEG,:) = reshape( TMP14, (/ NVEG, size(SAIM,2) /))
       LAIM(1:NVEG,:) = reshape( TMP15, (/ NVEG, size(LAIM,2) /))
       EPS(1:NVEG,:)  = reshape( TMP16, (/ NVEG, size(EPS,2)  /))

       RHOL(NVEG+1:MVT,:) = -1.E36
       RHOS(NVEG+1:MVT,:) = -1.E36
       TAUL(NVEG+1:MVT,:) = -1.E36
       TAUS(NVEG+1:MVT,:) = -1.E36
       SAIM(NVEG+1:MVT,:) = -1.E36
       LAIM(NVEG+1:MVT,:) = -1.E36
       EPS( NVEG+1:MVT,:) = -1.E36
    endif

  end subroutine read_mp_veg_parameters

END MODULE NOAHMP_VEG_PARAMETERS
! ==================================================================================================
