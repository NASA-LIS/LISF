module ac_climprocessing

use ac_global, only:    DaysInMonth, &
                        DetermineDate, &
                        DetermineDayNr, &
                        GetEToFilefull, &
                        GetEToRecord_FromD, &
                        GetEToRecord_FromM, &
                        GetEToRecord_FromY, &
                        GetEToRecord_NrObs, &
                        GetEToRecord_ToD, &
                        GetEToRecord_ToM, &
                        GetEToRecord_ToY, &
                        GetRainFilefull, &
                        GetRainRecord_FromD, &
                        GetRainRecord_FromM, &
                        GetRainRecord_FromY, &
                        GetRainRecord_NrObs, &
                        GetRainRecord_ToM, &
                        GetRainRecord_ToY, &
                        LeapYear, &
                        rep_DayEventDbl
use ac_kinds, only: dp, &
                    int8, &
                    int32, &
                    sp
implicit none


contains


subroutine AdjustDecadeMONTHandYEAR(DecFile, Mfile, Yfile)
    integer(int32), intent(inout) :: DecFile
    integer(int32), intent(inout) :: Mfile
    integer(int32), intent(inout) :: Yfile

    DecFile = 1
    Mfile = Mfile + 1
    if (Mfile > 12) then
        Mfile = 1
        YFile = Yfile + 1
    end if
end subroutine AdjustDecadeMONTHandYEAR


subroutine AdjustMONTHandYEAR(Mfile, Yfile)
    integer(int32), intent(inout) :: Mfile
    integer(int32), intent(inout) :: Yfile

    Mfile = Mfile - 12
    YFile = Yfile + 1
end subroutine AdjustMONTHandYEAR


subroutine GetParameters(C1, C2, C3, UL, LL, Mid)
    real(dp), intent(in) :: C1
    real(dp), intent(in) :: C2
    real(dp), intent(in) :: C3
    real(dp), intent(inout) :: UL
    real(dp), intent(inout) :: LL
    real(dp), intent(inout) :: Mid

    UL = (C1+C2)/2.0_dp
    LL = (C2+C3)/2.0_dp
    Mid = 2.0_dp*C2 - (UL+LL)/2.0_dp
    ! --previous decade-->/UL/....... Mid ......../LL/<--next decade--
end subroutine GetParameters


subroutine GetInterpolationParameters(C1, C2, C3, aOver3, bOver2, c)
    real(dp), intent(in) :: C1
    real(dp), intent(in) :: C2
    real(dp), intent(in) :: C3
    real(dp), intent(inout) :: aOver3
    real(dp), intent(inout) :: bOver2
    real(dp), intent(inout) :: c

    ! n1=n2=n3=30 --> better parabola
    aOver3 = (C1-2*C2+C3)/(6*30*30*30)
    bOver2 = (-6*C1+9*C2-3*C3)/(6*30*30)
    c = (11*C1-7*C2+2*C3)/(6*30)

end subroutine GetInterpolationParameters


subroutine GetMonthlyEToDataSet(DayNri, EToDataSet)
    integer(int32), intent(in) :: DayNri
    type(rep_DayEventDbl), dimension(31), intent(inout) :: EToDataSet

    integer(int32) :: Dayi, Monthi, Yeari, DayN
    integer(int32) :: DNR
    integer(int32) :: X1, X2, X3, t1, t2
    real(dp) :: C1, C2, C3
    real(dp) :: aOver3, bOver2, c

    ! GetMonthlyEToDataSet
    call DetermineDate(DayNri, Dayi, Monthi, Yeari)
    call GetSetofThreeMonths(Monthi, Yeari, C1, C2, C3, X1, X2, X3, t1)

    Dayi = 1
    call DetermineDayNr(Dayi, Monthi, Yeari, DNR)
    DayN = DaysInMonth(Monthi)
    if ((Monthi == 2) .and. LeapYear(Yeari)) then
        DayN = DayN + 1
    end if
    call GetInterpolationParameters(C1, C2, C3, aOver3, bOver2, c)
    do Dayi = 1, DayN
        t2 = t1 + 1
        EToDataSet(Dayi)%DayNr = DNR+Dayi-1
        EToDataSet(Dayi)%Param = aOver3*(t2*t2*t2-t1*t1*t1) &
                                 + bOver2*(t2*t2-t1*t1) &
                                 + c*(t2-t1)
        if (EToDataSet(Dayi)%Param < 0) then
            EToDataSet(Dayi)%Param = 0
        end if
        t1 = t2
    end do
    do Dayi = (DayN+1), 31
        EToDataSet(Dayi)%DayNr = DNR+DayN-1
        EToDataSet(Dayi)%Param = 0._dp
    end do


    contains


    subroutine GetSetofThreeMonths(Monthi, Yeari, C1, C2, C3, X1, X2, X3, t1)
        integer(int32), intent(in) :: Monthi
        integer(int32), intent(in) :: Yeari
        real(dp), intent(inout) :: C1
        real(dp), intent(inout) :: C2
        real(dp), intent(inout) :: C3
        integer(int32), intent(inout) :: X1
        integer(int32), intent(inout) :: X2
        integer(int32), intent(inout) :: X3
        integer(int32), intent(inout) :: t1

        integer :: fETo
        integer(int32) :: Mfile, Yfile, Nri, Obsi
        integer(int32), parameter :: ni=30
        logical :: OK3

        ! 1. Prepare record
        open(newunit=fETo, file=trim(GetEToFilefull()), status='old', &
                                                        action='read')
        read(fETo, *) ! description
        read(fETo, *) ! time step
        read(fETo, *) ! day
        read(fETo, *) ! month
        read(fETo, *) ! year
        read(fETo, *)
        read(fETo, *)
        read(fETo, *)

        Mfile = GetEToRecord_FromM()
        if (GetEToRecord_FromY() == 1901) then
            Yfile = Yeari
        else
            Yfile = GetEToRecord_FromY()
        end if
        OK3 = .false.

        ! 2. IF 3 or less records
        if (GetEToRecord_NrObs() <= 3) then
            read(fETo, *) C1
            C1 = C1 * ni
            X1 = ni
            select case (GetEToRecord_NrObs())
                case(1)
                    t1 = X1
                    X2 = X1 + ni
                    C2 = C1
                    X3 = X2 + ni
                    C3 = C1
                case(2)
                    t1 = X1
                    Mfile = Mfile + 1
                    if (Mfile > 12) then
                        call AdjustMONTHandYEAR(Mfile, Yfile)
                    end if
                    read(fETo, *) C3
                    C3 = C3 * ni
                    if (Monthi == Mfile) then
                        C2 = C3
                        X2 = X1 + ni
                        X3 = X2 + ni
                    else
                        C2 = C1
                        X2 = X1 + ni
                        X3 = X2 + ni
                    end if
                case(3)
                    if (Monthi == Mfile) then
                        t1 = 0
                    end if
                    Mfile = Mfile + 1
                    if (Mfile > 12) then
                        call AdjustMONTHandYEAR(Mfile, Yfile)
                    end if
                    read(fETo, *) C2
                    C2 = C2 * ni
                    X2 = X1 + ni
                    if (Monthi == Mfile) then
                        t1 = X1
                    end if
                    Mfile = Mfile + 1
                    if (Mfile > 12) then
                        call AdjustMONTHandYEAR(Mfile, Yfile)
                    end if
                    read(fETo, *) C3
                    C3 = C3 * ni
                    X3 = X2 + ni
                    if (Monthi == Mfile) then
                        t1 = X2
                    end if
            end select
            OK3 = .true.
        end if

        ! 3. If first observation
        if ((.not. OK3) .and. ((Monthi == Mfile) .and. (Yeari == Yfile))) then
            t1 = 0
            read(fETo, *) C1
            C1 = C1 * ni
            X1 = ni
            Mfile = Mfile + 1
            if (Mfile > 12) then
                call AdjustMONTHandYEAR(Mfile, Yfile)
            end if
            read(fETo, *) C2
            C2 = C2 * ni
            X2 = X1 + ni
            Mfile = Mfile + 1
            if (Mfile > 12) then
                call AdjustMONTHandYEAR(Mfile, Yfile)
            end if
            read(fETo, *) C3
            C3 = C3 * ni
            X3 = X2 + ni
            OK3 = .true.
        end if

        ! 4. If last observation
        if ((.not. OK3) .and. (Monthi == GetEToRecord_ToM())) then
            if ((GetEToRecord_FromY() == 1901) &
                        .or. (Yeari == GetEToRecord_ToY())) then
                do Nri = 1, (GetEToRecord_NrObs()-3)
                    read(fETo, *)
                    Mfile = Mfile + 1
                    if (Mfile > 12) then
                        call AdjustMONTHandYEAR(Mfile, Yfile)
                    end if
                end do
                read(fETo, *) C1
                C1 = C1 * ni
                X1 = ni
                Mfile = Mfile + 1
                if (Mfile > 12) then
                    call AdjustMONTHandYEAR(Mfile, Yfile)
                end if
                read(fETo, *) C2
                C2 = C2 * ni
                X2 = X1 + ni
                t1 = X2
                Mfile = Mfile + 1
                if (Mfile > 12) then
                    call AdjustMONTHandYEAR(Mfile, Yfile)
                end if
                read(fETo, *) C3
                C3 = C3 * ni
                X3 = X2 + ni
                OK3 = .true.
            end if
        end if

        ! 5. IF not previous cases
        if (.not. OK3) then
            Obsi = 1
            loop: do
                if ((Monthi == Mfile) .and. (Yeari == Yfile)) then
                    OK3 = .true.
                else
                    Mfile = Mfile + 1
                    if (Mfile > 12) then
                        call AdjustMONTHandYEAR(Mfile, Yfile)
                    end if
                    Obsi = Obsi + 1
                end if
                if (OK3) exit loop
            end do loop
            Mfile = GetEToRecord_FromM()
            do Nri = 1, (Obsi-2)
                read(fETo, *)
                Mfile = Mfile + 1
                if (Mfile > 12) then
                    call AdjustMONTHandYEAR(Mfile, Yfile)
                end if
            end do
            read(fETo, *) C1
            C1 = C1 * ni
            X1 = ni
            t1 = X1
            Mfile = Mfile + 1
            if (Mfile > 12) then
                call AdjustMONTHandYEAR(Mfile, Yfile)
            end if
            read(fETo, *) C2
            C2 = C2 * ni
            X2 = X1 + ni
            Mfile = Mfile + 1
            if (Mfile > 12) then
                call AdjustMONTHandYEAR(Mfile, Yfile)
            end if
            read(fETo, *) C3
            C3 = C3 * ni
            X3 = X2 + ni
        end if

        close(fETo)
    end subroutine GetSetofThreeMonths
end subroutine GetMonthlyEToDataSet


subroutine GetDecadeEToDataSet(DayNri, EToDataSet)
    integer(int32), intent(in) :: DayNri
    type(rep_DayEventDbl), dimension(31), intent(inout) :: EToDataSet

    integer(int32) :: Nri, ni, Dayi, Deci, Monthi, Yeari, DayN
    integer(int32) :: DNR
    real(dp) :: C1, C2, C3
    real(dp) :: Ul, LL, Mid

    call DetermineDate(DayNri, Dayi, Monthi, Yeari)
    if (Dayi > 20) then
        Deci = 3
        Dayi = 21
        DayN = DaysInMonth(Monthi)
        if ((Monthi == 2) .and. LeapYear(Yeari)) then
            DayN = DayN + 1
        end if
        ni = DayN - Dayi + 1
    elseif (Dayi > 10) then
        Deci = 2
        Dayi = 11
        DayN = 20
        ni = 10
    else
        Deci = 1
        Dayi = 1
        DayN = 10
        ni = 10
    end if
    call GetSetofThree(DayN, Deci, Monthi, Yeari, C1, C2, C3)
    call DetermineDayNr(Dayi, Monthi, Yeari, DNR)
    if (abs(C2) < epsilon(0._dp)) then
        do Nri = 1, ni
            EToDataSet(Nri)%DayNr = DNR+Nri-1
            EToDataSet(Nri)%Param = 0._dp
        end do
    else
        call GetParameters(C1, C2, C3, UL, LL, Mid)
        do Nri = 1, ni
            EToDataSet(Nri)%DayNr = DNR+Nri-1
            if (Nri <= (ni/2._dp+0.01)) then
                EToDataSet(Nri)%Param = (2._dp*UL + (Mid-UL)*(2._dp*Nri-1._dp) &
                                        / (ni/2._dp))/2._dp
            else
                if (((ni == 11) .or. (ni == 9)) .and. (Nri < (ni+1.01)/2)) then
                    EToDataSet(Nri)%Param = Mid
                else
                    EToDataSet(Nri)%Param = (2._dp*Mid &
                                             + (LL-Mid) &
                                             * (2._dp*Nri &
                                                -(ni+1._dp))/(ni/2._dp))/2._dp
                end if
            end if
            if (EToDataSet(Nri)%Param < 0._dp) then
                EToDataSet(Nri)%Param = 0._dp
            end if
        end do
    end if

    do Nri = (ni+1), 31
        EToDataSet(Nri)%DayNr = DNR+ni-1
        EToDataSet(Nri)%Param = 0._dp
    end do


    contains


    subroutine GetSetofThree(DayN, Deci, Monthi, Yeari, C1, C2, C3)
        integer(int32), intent(in) :: DayN
        integer(int32), intent(in) :: Deci
        integer(int32), intent(in) :: Monthi
        integer(int32), intent(in) :: Yeari
        real(dp), intent(inout) :: C1
        real(dp), intent(inout) :: C2
        real(dp), intent(inout) :: C3

        integer :: fETo
        integer(int32) :: DecFile, Mfile, Yfile, Nri, Obsi
        logical :: OK3


        !! 1 = previous decade, 2 = Actual decade, 3 = Next decade;
        open(newunit=fETo, file=trim(GetEToFileFull()), status='old', &
                                                        action='read')
        read(fETo, *) ! description
        read(fETo, *) ! time step
        read(fETo, *) ! day
        read(fETo, *) ! month
        read(fETo, *) ! year
        read(fETo, *)
        read(fETo, *)
        read(fETo, *)

        if (GetEToRecord_FromD() > 20) then
            DecFile = 3
        elseif (GetEToRecord_FromD() > 10) then
            DecFile = 2
        else
            DecFile = 1
        end if
        Mfile = GetEToRecord_FromM()
        if (GetEToRecord_FromY() == 1901) then
            Yfile = Yeari
        else
            Yfile = GetEToRecord_FromY()
        end if
        OK3 = .false.

        if (GetEToRecord_NrObs() <= 2) then
            read(fETo, *) C1
            select case (GetEToRecord_NrObs())
                case(1)
                    C2 = C1
                    C3 = C1
                case(2)
                    DecFile = DecFile + 1
                    if (DecFile > 3) then
                        call AdjustDecadeMONTHandYEAR(DecFile, Mfile, Yfile)
                    end if
                    read(fETo, *) C3
                    if (Deci == DecFile) then
                        C2 = C3
                        C3 = C2+(C2-C1)/4._dp
                    else
                        C2 = C1
                        C1 = C2 + (C2-C3)/4._dp
                    end if
            end select
            OK3 = .true.
        end if

        if ((.not. OK3) .and. ((Deci == DecFile) &
                        .and. (Monthi == Mfile) &
                        .and. (Yeari == Yfile))) then
            read(fETo, *) C1
            C2 = C1
            read(fETo, *) C3
            C1 = C2 + (C2-C3)/4._dp
            OK3 = .true.
        end if

        if ((.not. OK3) .and. ((DayN == GetEToRecord_ToD()) &
                        .and. (Monthi == GetEToRecord_ToM()))) then
            if ((GetEToRecord_FromY() == 1901) &
                        .or. (Yeari == GetEToRecord_ToY())) then
                do Nri = 1, (GetEToRecord_NrObs()-2)
                    read(fETo, *)
                end do
                read(fETo, *) C1
                read(fETo, *) C2
                C3 = C2+(C2-C1)/4._dp
                OK3 = .true.
            end if
        end if

        if (.not. OK3) then
            Obsi = 1
            loop: do
                if ((Deci == DecFile) .and. (Monthi == Mfile) &
                                      .and. (Yeari == Yfile)) then
                    OK3 = .true.
                else
                    DecFile = DecFile + 1
                    if (DecFile > 3) then
                        call AdjustDecadeMONTHandYEAR(DecFile, Mfile, Yfile)
                    end if
                    Obsi = Obsi + 1
                end if
                if (OK3) exit loop
            end do loop
            if (GetEToRecord_FromD() > 20) then
                DecFile = 3
            elseif (GetEToRecord_FromD() > 10) then
                DecFile = 2
            else
                DecFile = 1
            end if
            do Nri = 1, (Obsi-2)
                read(fETo, *)
            end do
            read(fETo, *) C1
            read(fETo, *) C2
            read(fETo, *) C3
        end if
        close(fETo)
    end subroutine GetSetofThree
end subroutine GetDecadeEToDataSet


subroutine GetDecadeRainDataSet(DayNri, RainDataSet)
    integer(int32), intent(in) :: DayNri
    type(rep_DayEventDbl), dimension(31), intent(inout) :: RainDataSet

    integer(int32) :: Nri, Day1, Deci, Monthi, Yeari, ni, DecFile, Mfile, Yfile
    integer(int32) :: DNR
    integer :: fRain
    logical :: OKRain
    real(dp) :: C

    call DetermineDate(DayNri, Day1, Monthi, Yeari)

    ! 0. Set Monthly Parameters

    ! 1. Which decade ?
    if (Day1 > 20) then
        Deci = 3
        Day1 = 21
        ni = DaysInMonth(Monthi) - Day1 + 1
        if ((Monthi == 2) .and. LeapYear(Yeari)) then
            ni = ni + 1
        end if
    elseif (Day1 > 10) then
        Deci = 2
        Day1 = 11
        ni = 10
    else
        Deci = 1
        Day1 = 1
        ni = 10
    end if

    ! 2. Load datafile
    open(newunit=fRain, file=trim(GetRainfilefull()), status='old', &
                                                      action='read')
    read(fRain, *) ! description
    read(fRain, *) ! time step
    read(fRain, *) ! day
    read(fRain, *) ! month
    read(fRain, *) ! year
    read(fRain, *)
    read(fRain, *)
    read(fRain, *)
    if (GetRainRecord_FromD() > 20) then
        DecFile = 3
    elseif (GetRainRecord_FromD() > 10) then
        DecFile = 2
    else
        DecFile = 1
    end if
    Mfile = GetRainRecord_FromM()
    if (GetRainRecord_FromY() == 1901) then
        Yfile = Yeari
    else
        Yfile = GetRainRecord_FromY()
    end if

    ! 3. Find decade
    OKRain = .false.
    C = 999._dp
    loop: do
        if ((Deci == DecFile) .and. (Monthi == Mfile) &
                              .and. (Yeari == Yfile)) then
            read(fRain, *) C
            OKRain = .true.
        else
            read(fRain, *)
            DecFile = DecFile + 1
            if (DecFile > 3) then
                call AdjustDecadeMONTHandYEAR(DecFile, Mfile, Yfile)
            end if
        end if
        if (OKRain) exit loop
    end do loop
    close(fRain)

    ! 4. Process data
    call DetermineDayNr(Day1, Monthi, Yeari, DNR)
    do Nri = 1, ni
        RainDataSet(Nri)%DayNr = DNR+Nri-1
        RainDataSet(Nri)%Param = C/ni
    end do
    do Nri = (ni+1), 31
        RainDataSet(Nri)%DayNr = DNR+ni-1
        RainDataSet(Nri)%Param = 0._dp
    end do
end subroutine GetDecadeRainDataSet


subroutine GetMonthlyRainDataSet(DayNri, RainDataSet)
    integer(int32), intent(in) :: DayNri
    type(rep_DayEventDbl), dimension(31), intent(inout) :: RainDataSet

    integer(int32) :: Dayi, DayN, Monthi, Yeari
    real(dp) :: C1, C2, C3, RainDec1, RainDec2, RainDec3
    integer(int32) :: DNR

    call DetermineDate(DayNri, Dayi, Monthi, Yeari)

    ! Set Monthly Parameters

    call GetSetofThreeMonths(Monthi, Yeari, C1, C2, C3)

    Dayi = 1
    call DetermineDayNr(Dayi, Monthi, Yeari, DNR)
    DayN = DaysInMonth(Monthi)
    if ((Monthi == 2) .and. LeapYear(Yeari)) then
        DayN = DayN + 1
    end if
    if (C2 > epsilon(0._dp)) then
        RainDec1 = (5._dp*C1 + 26._dp*C2 - 4._dp*C3)/(27._dp*3._dp) ! mm/dec
        RainDec2 = (-C1 + 29._dp*C2 - C3)/(27._dp*3._dp)
        RainDec3 = (-4._dp*C1 + 26._dp*C2 + 5._dp*C3)/(27._dp*3._dp)
        do Dayi = 1, 10
            RainDataSet(Dayi)%DayNr = DNR+Dayi-1
            RainDataSet(Dayi)%Param = RainDec1/10._dp
            if (RainDataSet(Dayi)%Param < epsilon(0._dp)) then
                RainDataSet(Dayi)%Param = 0._dp
            end if
        end do
        do Dayi = 11, 20
            RainDataSet(Dayi)%DayNr = DNR+Dayi-1
            RainDataSet(Dayi)%Param = RainDec2/10._dp
            if (RainDataSet(Dayi)%Param < epsilon(0._dp)) then
                RainDataSet(Dayi)%Param = 0._dp
            end if
        end do
        do Dayi = 21, DayN
            RainDataSet(Dayi)%DayNr = DNR+Dayi-1
            RainDataSet(Dayi)%Param = RainDec3/(DayN-21._dp+1._dp)
            if (RainDataSet(Dayi)%Param < epsilon(0._dp)) then
                RainDataSet(Dayi)%Param = 0._dp
            end if
        end do
    else
        do Dayi = 1, DayN
        RainDataSet(Dayi)%DayNr = DNR+Dayi-1
        RainDataSet(Dayi)%Param = 0._dp
        end do
    end if

    do Dayi = (DayN+1), 31
        RainDataSet(Dayi)%DayNr = DNR+DayN-1
        RainDataSet(Dayi)%Param = 0._dp
    end do


    contains


    subroutine GetSetofThreeMonths(Monthi, Yeari, C1, C2, C3)
        integer(int32), intent(in) :: Monthi
        integer(int32), intent(in) :: Yeari
        real(dp), intent(inout) :: C1
        real(dp), intent(inout) :: C2
        real(dp), intent(inout) :: C3

        integer :: fRain
        integer(int32) :: Mfile, Yfile, Nri, Obsi
        logical :: OK3

        ! 1. Prepare record
        open(newunit=fRain, file=trim(GetRainFilefull()), status='old', &
                                                          action='read')
        read(fRain, *) ! description
        read(fRain, *) ! time step
        read(fRain, *) ! day
        read(fRain, *) ! month
        read(fRain, *) ! year
        read(fRain, *)
        read(fRain, *)
        read(fRain, *)
        Mfile = GetRainRecord_FromM()
        if (GetRainRecord_FromY() == 1901) then
            Yfile = Yeari
        else
            Yfile = GetRainRecord_FromY()
        end if
        OK3 = .false.

        ! 2. IF 2 or less records
        if (GetRainRecord_NrObs() <= 2) then
            read(fRain, *) C1
            select case (GetRainRecord_NrObs())
                case(1)
                    C2 = C1
                    C3 = C1
                case(2)
                    Mfile = Mfile + 1
                    if (Mfile > 12) then
                        call AdjustMONTHandYEAR(Mfile, Yfile)
                    end if
                    read(fRain, *) C3
                    if (Monthi == Mfile) then
                        C2 = C3
                    else
                        C2 = C1
                    end if
            end select
            OK3 = .true.
        end if

        ! 3. If first observation
        if ((.not. OK3) .and. ((Monthi == Mfile) &
                        .and. (Yeari == Yfile))) then
            read(fRain, *) C1
            C2 = C1
            read(fRain, *) C3
            OK3 = .true.
        end if

        ! 4. If last observation
        if ((.not. OK3) .and. (Monthi == GetRainRecord_ToM())) then
            if ((GetRainRecord_FromY() == 1901) .or. (Yeari == GetRainRecord_ToY())) then
                do Nri = 1, (GetRainRecord_NrObs()-2)
                    read(fRain, *)
                end do
                read(fRain, *) C1
                read(fRain, *) C2
                C3 = C2
                OK3 = .true.
            end if
        end if

        ! 5. IF not previous cases
        if (.not. OK3) then
            Obsi = 1
            loop: do
                if ((Monthi == Mfile) .and. (Yeari == Yfile)) then
                    OK3 = .true.
                else
                    Mfile = Mfile + 1
                    if (Mfile > 12) then
                        call AdjustMONTHandYEAR(Mfile, Yfile)
                    end if
                    Obsi = Obsi + 1
                end if
                if (OK3) exit loop
            end do loop
            Mfile = GetRainRecord_FromM()
            do Nri = 1, (Obsi-2)
                read(fRain, *)
                Mfile = Mfile + 1
                if (Mfile > 12) then
                    call AdjustMONTHandYEAR(Mfile, Yfile)
                end if
            end do
            read(fRain, *) C1
            read(fRain, *) C2
            read(fRain, *) C3
        end if
        close(fRain)
    end subroutine GetSetofThreeMonths
end subroutine GetMonthlyRainDataSet

end module ac_climprocessing
