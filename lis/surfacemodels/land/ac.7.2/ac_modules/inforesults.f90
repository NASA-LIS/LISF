module ac_inforesults

use ac_global, only:    typeObsSim_ObsSimCC, &
                        typeObsSim_ObsSimB, &
                        typeObsSim_ObsSimSWC, &
                        undef_int, &
                        GetPathNameSimul, &
                        GetSimulation_FromDayNr, &
                        typeProject_TypePRM, &
                        NameMonth
use ac_kinds, only: dp, &
                    int8, &
                    int16, &
                    int32, &
                    intEnum
use ac_utils, only: GetAquaCropDescriptionWithTimeStamp, &
                    roundc
use iso_fortran_env, only: iostat_end
implicit none


type rep_EventObsSim
    real(dp) :: Obsi
        !! Undocumented
    real(dp) :: StdObsi
        !! Undocumented
    real(dp) :: Simi
        !! Undocumented
    integer(int8) :: DDi
        !! Undocumented
    integer(int8) :: MMi
        !! Undocumented
    integer(int32) :: YYYYi
        !! Undocumented
end type rep_EventObsSim


contains


subroutine StatisticAnalysis(TypeObsSim, RangeObsMin, RangeObsMax, StrNr, &
                             Nobs, ObsAver, SimAver, PearsonCoeff, RMSE, &
                             NRMSE, NScoeff, IndexAg, ArrayObsSim)
    integer(intEnum), intent(in) :: TypeObsSim
    integer(int32), intent(in) :: RangeObsMin
    integer(int32), intent(in) :: RangeObsMax
    character(len=*), intent(in) :: StrNr
    integer(int32), intent(inout) :: Nobs
    real(dp), intent(inout) :: ObsAver
    real(dp), intent(inout) :: SimAver
    real(dp), intent(inout) :: PearsonCoeff
    real(dp), intent(inout) :: RMSE
    real(dp), intent(inout) :: NRMSE
    real(dp), intent(inout) :: NScoeff
    real(dp), intent(inout) :: IndexAg
    type(rep_EventObsSim), dimension(100), intent(inout) :: ArrayObsSim

    integer(int32) :: Nri
    real(dp) :: dDeNom, R2tel, SumSqrDobs, SumSqrDsim, SumSqrDiv, DeNom

    ! get data
    call GetObsSim(RangeObsMin, RangeObsMax, StrNr, Nobs, ArrayObsSim, &
                   ObsAver, SimAver)
    ! statistical evaluation
    if (Nobs > 1) then
        R2tel = 0._dp
        SumSqrDobs = 0._dp
        SumSqrDsim = 0._dp
        SumSqrDiv = 0._dp
        dDeNom = 0._dp
        do Nri = 1, Nobs
            R2tel = R2tel + (ArrayObsSim(Nri)%Obsi-ObsAver) &
                            * (ArrayObsSim(Nri)%Simi-SimAver)
            SumSqrDobs = SumSqrDobs + (ArrayObsSim(Nri)%Obsi-ObsAver)**2
            SumSqrDsim = SumSqrDsim + (ArrayObsSim(Nri)%Simi-SimAver)**2
            SumSqrDiv = SumSqrDiv &
                        + (ArrayObsSim(Nri)%Obsi-ArrayObsSim(Nri)%Simi)**2
            dDeNom = dDeNom + (abs(ArrayObsSim(Nri)%Simi-ObsAver) &
                               + abs(ArrayObsSim(Nri)%Obsi-ObsAver))**2
        end do
        ! R2
        DeNom = sqrt(SumSqrDobs*SumSqrDsim)
        if (DeNom > 0._dp) then
            PearsonCoeff = R2tel/DeNom
        else
            PearsonCoeff = real(undef_int, kind=dp)
        end if
        ! RMSE
        RMSE = sqrt(SumSqrDiv/NObs)
        ! NRMSE
        if (ObsAver > 0._dp) then
            NRMSE = 100._dp*(RMSE/ObsAver)
        else
            NRMSE = real(undef_int, kind=dp)
        end if
        ! Nash-Sutcliffe coefficient (EF)
        if (SumSqrDobs > 0._dp) then
            NScoeff = 1._dp - (SumSqrDiv/SumSqrDobs)
        else
            NScoeff = real(undef_int, kind=dp)
        end if
        ! Index of agreement (d)
        if (dDeNom > 0._dp) then
            IndexAg = 1._dp - (SumSqrDiv/dDeNom)
        else
            IndexAg = real(undef_int, kind=dp)
        end if
    else
        ObsAver = real(undef_int, kind=dp)
        SimAver = real(undef_int, kind=dp)
        PearsonCoeff = real(undef_int, kind=dp)
        RMSE = real(undef_int, kind=dp)
        NRMSE = real(undef_int, kind=dp)
        NScoeff = real(undef_int, kind=dp)
        IndexAg = real(undef_int, kind=dp)
    end if


    contains


    subroutine GetObsSim(RangeObsMin, RangeObsMax, StrNr, Nobs, &
                         ArrayObsSim, ObsAver, SimAver)
        integer(int32), intent(in) :: RangeObsMin
        integer(int32), intent(in) :: RangeObsMax
        character(len=*), intent(in) :: StrNr
        integer(int32), intent(inout) :: Nobs
        type(rep_EventObsSim), dimension(100), intent(inout) :: ArrayObsSim
        real(dp), intent(inout) :: ObsAver
        real(dp), intent(inout) :: SimAver

        character(len=:), allocatable :: OutputName
        character(len=1024) :: buffer
        integer :: fOut
        integer(int8) :: Dayi, Monthi
        integer(int32) :: SkipLines, NCobs, Yeari
        integer(int32) :: i, status
        real(dp) :: VarObsi, VarSimi, VarStdi
        real(dp), dimension(:), allocatable :: dummy_array

        Nobs = 0
        ObsAver = 0._dp
        SimAver = 0._dp

        ! open file
        OutputName = GetPathNameSimul() // 'EvalData' // trim(StrNr) // '.OUT'
        open(newunit=fOut, file=trim(OutputName), status='old', action='read')
        read(fOut, *) ! AquaCrop Version - Date and Time
        read(fOut, *) ! title
        read(fOut, *) !
        read(fOut, *) ! list of variables
        read(fOut, *) ! units

        ! find first day
        SkipLines = RangeObsMin - GetSimulation_FromDayNr()
        do i = 1, SkipLines
            read(fOut, *)
        end do

        ! get Sim and Obs in range
        select case (TypeObsSim)
        case(typeObsSim_ObsSimCC)
            NCobs = 2
        case (typeObsSim_ObsSimB)
            NCobs = 5
        case (typeObsSim_ObsSimSWC)
            NCobs = 8
        case default
            NCobs = 5
        end select

        allocate(dummy_array(NCobs))

        do i = RangeObsMin, RangeObsMax
            read(fOut, '(a)', iostat=status) buffer

            if (status /= iostat_end) then
                read(buffer, *) Dayi, Monthi, Yeari, dummy_array, &
                                VarSimi, VarObsi, VarStdi

                if ((roundc(VarObsi, mold=1) /= undef_int) &
                                .and. (Nobs < 100)) then
                    Nobs = Nobs + 1
                    ArrayObsSim(Nobs)%DDi = Dayi
                    ArrayObsSim(Nobs)%MMi = Monthi
                    ArrayObsSim(Nobs)%YYYYi = Yeari
                    ArrayObsSim(Nobs)%Simi = VarSimi
                    ArrayObsSim(Nobs)%Obsi = VarObsi
                    ArrayObsSim(Nobs)%StdObsi = VarStdi
                    SimAver = SimAver + VarSimi
                    ObsAver = ObsAver + VarObsi
                end if
            end if
        end do

        ! close file
        close(fOut)

        ! calculate averages
        if (Nobs > 0) then
            ObsAver = ObsAver/real(Nobs, kind=dp)
            SimAver = SimAver/real(Nobs, kind=dp)
        end if
    end subroutine GetObsSim
end subroutine StatisticAnalysis


subroutine WriteAssessmentSimulation(StrNr, totalnameEvalStat, &
                                     TheProjectType, RangeMin, RangeMax)
    character(len=*), intent(in) :: StrNr
    character(len=*), intent(in) :: totalnameEvalStat
    integer(intEnum), intent(in) :: TheProjectType
    integer(int32), intent(in) :: RangeMin
    integer(int32), intent(in) :: RangeMax

    integer :: fAssm
    integer(intEnum) :: TypeObsSim
    integer(int32) :: Nobs, Nri
    real(dp) :: ObsAver, SimAver, PearsonCoeff, RMSE, NRMSE, &
                NScoeff, IndexAg
    type(rep_EventObsSim), dimension(100) :: ArrayObsSim
    character(len=:), allocatable :: YearString

    ! 1. Open file for assessment
    open(newunit=fAssm, file=trim(totalnameEvalStat), status='replace', action='write')
    write(fAssm, '(a)') GetAquaCropDescriptionWithTimeStamp()
    write(fAssm, '(a)') 'Evaluation of simulation results - Statistics'
    if (TheProjectType == typeProject_TypePRM) then
        write(fAssm, '(2a)') '** Run number:', StrNr
    end if
    write(fAssm, '(a)')

    close(fAssm)

    ! 2. Run analysis

    ! 2.1 Canopy Cover
    TypeObsSim = typeObsSim_ObsSimCC
    call StatisticAnalysis(TypeObsSim, RangeMin, RangeMax, StrNr, Nobs, &
                           ObsAver, SimAver, PearsonCoeff, RMSE, NRMSE, &
                           NScoeff, IndexAg, ArrayObsSim)

    open(newunit=fAssm, file=trim(totalnameEvalStat), status='old', &
         position='append', action='write')

    write(fAssm, '(a)')
    write(fAssm, '(a)') &
    '  ASSESSMENT OF CANOPY COVER --------------------------------------'
    if (Nobs > 1) then
        write(fAssm, '(a)') &
        '              --------- Canopy Cover (%) ---------'
        write(fAssm, '(a)') &
        '    Nr        Observed    +/- St Dev     Simulated    Date'
        write(fAssm, '(a)') &
        '  ----------------------------------------------------------------'
        do Nri = 1, Nobs
            if (ArrayObsSim(Nri)%YYYYi <= 1901) then
                YearString = ''
            else
                YearString = '    '
                write(YearString, '(i4)') ArrayObsSim(Nri)%YYYYi
            end if
            write(fAssm, '(i6, 3f14.1, a, i2, 4a)') &
                            Nri, ArrayObsSim(Nri)%Obsi, &
                            ArrayObsSim(Nri)%StdObsi, &
                            ArrayObsSim(Nri)%Simi, '      ', &
                            ArrayObsSim(Nri)%DDi, ' ', &
                            trim(NameMonth(ArrayObsSim(Nri)%MMi)), &
                            ' ', trim(YearString)
        end do
        write(fAssm, '(a)')
        write(fAssm, '(a, i5)') &
        '  Valid observations/simulations sets (n) ....... : ', Nobs
        write(fAssm, '(a, f7.1, a)') &
        '  Average of observed Canopy Cover .............. : ', ObsAver, '   %'
        write(fAssm, '(a, f7.1, a)') &
        '  Average of simulated Canopy Cover ............. : ', SimAver, '   %'
        write(fAssm, '(a)')
        write(fAssm, '(a, f8.2)') &
        '  Pearson Correlation Coefficient (r) ........... : ', PearsonCoeff
        write(fAssm, '(a, f7.1, a)') &
        '  Root mean square error (RMSE) ................. : ', RMSE, '   % CC'
        write(fAssm, '(a, f7.1, a)') &
        '  Normalized root mean square error  CV(RMSE).... : ', NRMSE, '   %'
        write(fAssm, '(a, f8.2)') &
        '  Nash-Sutcliffe model efficiency coefficient (EF): ', NScoeff
        write(fAssm, '(a, f8.2, a)') &
        '  Willmotts index of agreement (d) .............. : ', IndexAg
    else
        write(fAssm, '(a)') &
        '  No statistic analysis (insufficient data)'
    end if
    write(fAssm, '(a)') &
    '  ----------------------------------------------------------------'

    ! 2.2 Biomass production
    TypeObsSim = typeObsSim_ObsSimB
    call StatisticAnalysis(TypeObsSim, RangeMin, RangeMax, StrNr, Nobs, &
                           ObsAver, SimAver, PearsonCoeff, RMSE, NRMSE, &
                           NScoeff, IndexAg, ArrayObsSim)
    write(fAssm, '(a)')
    write(fAssm, '(a)')
    write(fAssm, '(a)') &
    '  ASSESSMENT OF BIOMASS PRODUCTION --------------------------------'
    if (Nobs > 1) then
        write(fAssm, '(a)') &
        '              --------- Biomass (ton/ha) ---------'
        write(fAssm, '(a)') &
        '    Nr        Observed    +/- St Dev     Simulated    Date'
        write(fAssm, '(a)') &
        '  ----------------------------------------------------------------'
        do Nri = 1, Nobs
            if (ArrayObsSim(Nri)%YYYYi <= 1901) then
                YearString = ''
            else
                YearString = '    '
                write(YearString, '(i4)') ArrayObsSim(Nri)%YYYYi
            end if
            write(fAssm, '(i6, f16.3, 2f14.3, a, i2, 4a)') &
                            Nri, ArrayObsSim(Nri)%Obsi, &
                            ArrayObsSim(Nri)%StdObsi, &
                            ArrayObsSim(Nri)%Simi, '      ', &
                            ArrayObsSim(Nri)%DDi, ' ', &
                            trim(NameMonth(ArrayObsSim(Nri)%MMi)), &
                            ' ', trim(YearString)
        end do
        write(fAssm, '(a)')
        write(fAssm, '(a, i5)') &
        '  Valid observations/simulations sets (n) ....... : ', Nobs
        write(fAssm, '(a, f9.3, a)') &
        '  Average of observed Biomass production ........ : ', ObsAver, '   ton/ha'
        write(fAssm, '(a, f9.3, a)') &
        '  Average of simulated Biomass production ....... : ', SimAver, '   ton/ha'
        write(fAssm, '(a)')
        write(fAssm, '(a, f8.2)') &
        '  Pearson Correlation Coefficient (r) ........... : ', PearsonCoeff
        write(fAssm, '(a, f9.3, a)') &
        '  Root mean square error (RMSE) ................. : ', RMSE, '   ton/ha'
        write(fAssm, '(a, f7.1, a)') &
        '  Normalized root mean square error  CV(RMSE).... : ', NRMSE, '   %'
        write(fAssm, '(a, f8.2)') &
        '  Nash-Sutcliffe model efficiency coefficient (EF): ', NScoeff
        write(fAssm, '(a, f8.2)') &
        '  Willmotts index of agreement (d) .............. : ', IndexAg
    else
        write(fAssm, '(a)') '  No statistic analysis (insufficient data)'
    end if
    write(fAssm, '(a)') &
    '  ----------------------------------------------------------------'

    ! 2.3 Soil Water Content
    TypeObsSim = typeObsSim_ObsSimSWC
    call StatisticAnalysis(TypeObsSim, RangeMin, RangeMax, StrNr, Nobs, &
                           ObsAver, SimAver, PearsonCoeff, RMSE, NRMSE, &
                           NScoeff, IndexAg, ArrayObsSim)
    write(fAssm, '(a)')
    write(fAssm, '(a)')
    write(fAssm, '(a)') &
    '  ASSESSMENT OF SOIL WATER CONTENT --------------------------------'
    if (Nobs > 1) then
        write(fAssm, '(a)') &
        '              ------ Soil water content (mm) -----'
        write(fAssm, '(a)') &
        '    Nr        Observed    +/- St Dev     Simulated    Date'
        write(fAssm, '(a)') &
        '  ----------------------------------------------------------------'
        do Nri = 1, Nobs
            if (ArrayObsSim(Nri)%YYYYi <= 1901) then
                YearString = ''
            else
                YearString = '    '
                write(YearString, '(i4)') ArrayObsSim(Nri)%YYYYi
            end if
            write(fAssm, '(i6, 3f14.1, a, i2, 4a)') &
                            Nri, ArrayObsSim(Nri)%Obsi, &
                            ArrayObsSim(Nri)%StdObsi, &
                            ArrayObsSim(Nri)%Simi, '      ', &
                            ArrayObsSim(Nri)%DDi, ' ', &
                            trim(NameMonth(ArrayObsSim(Nri)%MMi)), &
                            ' ', trim(YearString)
        end do
        write(fAssm, '(a)')
        write(fAssm, '(a, i5)') &
        '  Valid observations/simulations sets (n) ....... : ', Nobs
        write(fAssm, '(a, f7.1, a)') &
        '  Average of observed Soil water content ........ : ', ObsAver, '   mm'
        write(fAssm, '(a, f7.1, a)') &
        '  Average of simulated Soil water content ....... : ', SimAver, '   mm'
        write(fAssm, '(a)')
        write(fAssm, '(a, f8.2)') &
        '  Pearson Correlation Coefficient (r) ........... : ', PearsonCoeff
        write(fAssm, '(a, f7.1, a)') &
        '  Root mean square error (RMSE) ................. : ', RMSE, '   mm'
        write(fAssm, '(a, f7.1, a)') &
        '  Normalized root mean square error  CV(RMSE).... : ', NRMSE, '   %'
        write(fAssm, '(a, f8.2)') &
        '  Nash-Sutcliffe model efficiency coefficient (EF): ', NScoeff
        write(fAssm, '(a, f8.2, a)') &
        '  Willmotts index of agreement (d) .............. : ', IndexAg
    else
        write(fAssm, '(a)') '  No statistic analysis (insufficient data)'
    end if
    write(fAssm, '(a)') &
    '  ----------------------------------------------------------------'
    write(fAssm, '(a)')

    ! 3. Close file for assessment
    close(fAssm)
end subroutine WriteAssessmentSimulation

end module ac_inforesults
