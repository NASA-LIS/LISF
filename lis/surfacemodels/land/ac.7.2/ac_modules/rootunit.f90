module ac_rootunit

use ac_global, only: undef_int, &
                     ActualRootingDepth, &
                     GetCompartment_Layer, &
                     GetCompartment_Theta, &
                     GetCompartment_Thickness, &
                     Getcrop_KsShapeFactorStomata, &
                     GetCrop_pdef, &
                     GetCrop_SmaxTop, &
                     GetCrop_SmaxBot, &
                     GetNrCompartments, &
                     GetSimulation_Germinate, &
                     GetSimulation_SCor, &
                     GetSimulParam_KsShapeFactorRoot, &
                     GetSoil_RootMax, &
                     GetSoilLayer_FC, &
                     GetSoilLayer_WP, &
                     GetSumWaBal_Tpot, &
                     GetSumWaBal_Tact, &
                     KsAny, &
                     SetSimulation_SCor, &
                     SetSoil_RootMax
use ac_kinds, only: int32, &
                    int16, &
                    int8, &
                    intEnum, &
                    sp
use ac_utils, only: roundc
implicit none


contains


real(sp) function AdjustedRootingDepth(&
        CCAct, CCpot, Tpot, Tact, StressLeaf, StressSenescence, DAP, L0,&
        LZmax, L1234, GDDL0, GDDLZmax, GDDL1234, SumGDDPrev, SumGDD, Zmin,&
        Zmax, Ziprev, ShapeFactor, TypeDays)
    real(sp), intent(in) :: CCAct
    real(sp), intent(in) :: CCpot
    real(sp), intent(in) :: Tpot
    real(sp), intent(in) :: Tact
    real(sp), intent(in) :: StressLeaf
    real(sp), intent(in) :: StressSenescence
    integer(int32), intent(in) :: DAP
    integer(int32), intent(in) :: L0
    integer(int32), intent(in) :: LZmax
    integer(int32), intent(in) :: L1234
    integer(int32), intent(in) :: GDDL0
    integer(int32), intent(in) :: GDDLZmax
    integer(int32), intent(in) :: GDDL1234
    real(sp), intent(in) :: SumGDDPrev
    real(sp), intent(in) :: SumGDD
    real(sp), intent(in) :: Zmin
    real(sp), intent(in) :: Zmax
    real(sp), intent(in) :: Ziprev
    integer(int8), intent(in) :: ShapeFactor
    integer(intEnum), intent(in) :: TypeDays

    real(sp) :: Zi, ZiUnlimM1, ZiUnlim, dZ, ZiTest, Zsoil, ThetaTreshold
    real(sp) :: TAWcompi, Wrel, pZexp, Zlimit, ZiMax, KsShapeFactorRoot
    integer(int32) :: compi, layer

    if (roundc(Ziprev, mold=1) == undef_int) then
        Zi = ActualRootingDepth(DAP, L0, LZmax, L1234, GDDL0, GDDLZmax,&
                                SumGDD, Zmin, Zmax, ShapeFactor, TypeDays)
    else
        ! 1. maximum rooting depth (ZiMax) that could have been reached at
        !    time t
        ! -- 1.1 Undo effect of restrictive soil layer(s)
        if (roundc(GetSoil_RootMax()*1000.0_sp, mold=1)&
                < roundc(Zmax*1000.0_sp, mold=1)) then
            Zlimit = GetSoil_RootMax()
            call SetSoil_RootMax(real(Zmax, kind=sp))
        else
            Zlimit = Zmax
        end if

        ! -- 1.2 Calculate ZiMax
        ZiMax = ActualRootingDepth(DAP, L0, LZmax, L1234, GDDL0, GDDLZmax,&
                                   SumGDD, Zmin, Zmax, ShapeFactor, TypeDays)
        ! -- 1.3 Restore effect of restrive soil layer(s)
        call SetSoil_RootMax(real(Zlimit, kind=sp))

        ! 2. increase (dZ) at time t
        ZiUnlimM1 = ActualRootingDepth(DAP-1, L0, LZmax, L1234, GDDL0, GDDLZmax,&
                                       SumGDDPrev, Zmin, Zmax, ShapeFactor, TypeDays)
        ZiUnlim = ActualRootingDepth(DAP, L0, LZmax, L1234, GDDL0, GDDLZmax,&
                                     SumGDD, Zmin, Zmax, ShapeFactor, TypeDays)
        dZ = ZiUnlim - ZiUnlimM1

        ! 3. corrections of dZ
        ! -- 3.1 correction for restrictive soil layer is already considered
        !    in ActualRootingDepth

        ! -- 3.2 correction for stomatal closure
        if ((Tpot > 0.0_sp) .and. (Tact < Tpot)&
                .and. (GetSimulParam_KsShapeFactorRoot() /= undef_int)) then
            if (GetSimulParam_KsShapeFactorRoot() >= 0) then
                dZ = dZ * (Tact/Tpot)   ! linear
            else
                KsShapeFactorRoot = real(GetSimulParam_KsShapeFactorRoot(),&
                                         kind=sp)
                dZ = dZ * (exp((Tact/Tpot)*KsShapeFactorRoot)-1.0_sp) &
                        / (exp(KsShapeFactorRoot)-1.0_sp) ! exponential
            end if
        end if

        ! -- 3.2 correction for dry soil at expansion front of actual root
        !        zone
        if (dZ > 0.001_sp) then
            ! soil water depletion threshold for root deepening
            pZexp = GetCrop_pdef() + (1-GetCrop_pdef())/2.0_sp
            ! restrictive soil layer is considered by ActualRootingDepth
            ZiTest = Ziprev + dZ
            compi = 0
            Zsoil = 0.0_sp
            do while ((Zsoil < ZiTest) .and. (compi < GetNrCompartments()))
                compi = compi + 1
                Zsoil = Zsoil + GetCompartment_Thickness(compi)
            end do
            layer = GetCompartment_Layer(compi)
            TAWcompi = GetSoilLayer_FC(layer)/100.0_sp &
                        - GetSoilLayer_WP(layer)/100.0_sp
            ThetaTreshold = GetSoilLayer_FC(layer)/100.0_sp - pZexp * TAWcompi
            if (GetCompartment_Theta(compi) < ThetaTreshold) then
                ! expansion is limited due to soil water content at
                ! expansion front
                if (GetCompartment_Theta(compi) &
                        <= GetSoilLayer_WP(layer)/100.0_sp) then
                    dZ = 0.0_sp
                else
                    Wrel = (GetSoilLayer_FC(layer)/100.0_sp - &
                                GetCompartment_theta(compi))/TAWcompi
                    dZ = dZ * KsAny(Wrel, pZexp, 1.0_sp,&
                                    GetCrop_KsShapeFactorStomata())
                end if
            end if
        end if

        ! -- 3.3 correction for early senescence
        if ((CCact <= epsilon(0.0_sp)) .and. (CCpot > 50.0_sp)) then
            dZ = 0.0_sp
        end if

        ! -- 3.4 correction for no germination
        if (.not. GetSimulation_Germinate()) then
            dZ = 0.0_sp
        end if

        ! 4. actual rooting depth (Zi)
        Zi = Ziprev + dZ

        ! 5. Correction for root density if root deepening is restricted
        !    (dry soil and/or restricitive layers)
        if (roundc(Zi*1000, mold=1) < roundc(ZiMax*1000, mold=1)) then
            ! Total extraction in restricted root zone (Zi) and max root
            ! zone (ZiMax) should be identical
            call SetSimulation_SCor(real((2*(ZiMax/Zi)&
                          *((GetCrop_SmaxTop()+GetCrop_SmaxBot())/2.0_sp)&
                          - GetCrop_SmaxTop())/GetCrop_SmaxBot(), kind=sp))
            ! consider part of the restricted deepening due to water stress
            ! (= less roots)
            if (GetSumWaBal_Tpot() > 0.0_sp) then
                call SetSimulation_SCor(real(GetSimulation_SCor()&
                      * (GetSumWaBal_Tact()/GetSumWaBal_Tpot()), kind=sp))
                if (GetSimulation_SCor() < 1.0_sp) then
                    call SetSimulation_SCor(1.0_sp)
                end if
            end if
        else
            call SetSimulation_SCor(1.0_sp)
        end if
    end if ! (roundc(Ziprev, mold=1) == undef_int)
    AdjustedRootingDepth = Zi
end function AdjustedRootingDepth

end module ac_rootunit
