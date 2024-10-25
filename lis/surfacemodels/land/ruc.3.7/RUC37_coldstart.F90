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
! !ROUTINE: RUC37_coldstart
! \label{RUC37_coldstart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   1/15/15: Shugong Wang; initial implementation for LIS 7 and RUC37
!
! !INTERFACE:
subroutine RUC37_coldstart(mtype)
! !USES:
    use LIS_coreMod, only: LIS_rc
    use LIS_logMod, only: LIS_logunit
    use LIS_timeMgrMod, only: LIS_date2time
    use RUC37_lsmMod
   
!
! !DESCRIPTION:
!
!  This routine initializes the RUC37 state variables with some
!  predefined values constantly for the entire domain. 
!
!EOP
 
    implicit none

    integer :: mtype
    integer :: t, l, n, i
    integer :: c, r
    do n=1, LIS_rc%nnest
        if (trim(LIS_rc%startcode) .eq. "coldstart") then
            write(LIS_logunit,*) "MSG: RUC37_coldstart -- cold-starting RUC37"
            do t=1, LIS_rc%npatch(n,mtype)
                RUC37_struc(n)%ruc37(t)%emiss = RUC37_struc(n)%init_emiss
                RUC37_struc(n)%ruc37(t)%ch = RUC37_struc(n)%init_ch
                RUC37_struc(n)%ruc37(t)%cm = RUC37_struc(n)%init_cm
                RUC37_struc(n)%ruc37(t)%sneqv = RUC37_struc(n)%init_sneqv
                RUC37_struc(n)%ruc37(t)%snowh = RUC37_struc(n)%init_snowh
                RUC37_struc(n)%ruc37(t)%snowc = RUC37_struc(n)%init_snowc
                RUC37_struc(n)%ruc37(t)%canwat = RUC37_struc(n)%init_canwat
                RUC37_struc(n)%ruc37(t)%alb = RUC37_struc(n)%init_alb
                do l=1, RUC37_struc(n)%nsoil
                    RUC37_struc(n)%ruc37(t)%smc(l) = RUC37_struc(n)%init_smc(l)
                enddo
                do l=1, RUC37_struc(n)%nsoil
                    RUC37_struc(n)%ruc37(t)%sho(l) = RUC37_struc(n)%init_sho(l)
                enddo
                do l=1, RUC37_struc(n)%nsoil
                    RUC37_struc(n)%ruc37(t)%stc(l) = RUC37_struc(n)%init_stc(l)
                enddo
                do l=1, RUC37_struc(n)%nsoil
                    RUC37_struc(n)%ruc37(t)%smfr(l) = RUC37_struc(n)%init_smfr(l)
                enddo
                do l=1, RUC37_struc(n)%nsoil
                    RUC37_struc(n)%ruc37(t)%keepfr(l) = RUC37_struc(n)%init_keepfr(l)
                enddo

                RUC37_struc(n)%ruc37(t)%tskin = RUC37_struc(n)%init_tskin
                RUC37_struc(n)%ruc37(t)%qvg = RUC37_struc(n)%init_qvg
                RUC37_struc(n)%ruc37(t)%qsfc = RUC37_struc(n)%init_qsfc
                RUC37_struc(n)%ruc37(t)%qcg = RUC37_struc(n)%init_qcg
                RUC37_struc(n)%ruc37(t)%qsg = RUC37_struc(n)%init_qsg
                RUC37_struc(n)%ruc37(t)%snt75cm = RUC37_struc(n)%init_snt75cm
                RUC37_struc(n)%ruc37(t)%tsnav = RUC37_struc(n)%init_tsnav
                RUC37_struc(n)%ruc37(t)%soilm = RUC37_struc(n)%init_soilm
                RUC37_struc(n)%ruc37(t)%smroot = RUC37_struc(n)%init_smroot
            enddo
        endif
    
        LIS_rc%yr = LIS_rc%syr
        LIS_rc%mo = LIS_rc%smo
        LIS_rc%da = LIS_rc%sda
        LIS_rc%hr = LIS_rc%shr
        LIS_rc%mn = LIS_rc%smn
        LIS_rc%ss = LIS_rc%sss
        
        call LIS_date2time(LIS_rc%time, LIS_rc%doy, LIS_rc%gmt, LIS_rc%yr,      &
                           LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn, LIS_rc%ss)
        write(LIS_logunit,*) "MSG: RUC37_coldstart -- ",     &
                             "Using the specified start time ", LIS_rc%time
    enddo
end subroutine RUC37_coldstart
