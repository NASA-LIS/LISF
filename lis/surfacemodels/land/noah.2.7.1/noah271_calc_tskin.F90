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
! !ROUTINE: noah271_calc_tskin
! \label{noah271_calc_tskin}
! 
! !REVISION HISTORY: 
!  10 Apr 2008: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine noah271_calc_tskin (n, stc,t1)
! !USES: 
  use LIS_coreMod,      only : LIS_rc
  use LIS_constantsMod, only : LIS_CONST_TKFRZ
  use LIS_vegDataMod,     only : LIS_gfrac
  use noah271_lsmMod,      only : noah271_struc
  use module_sf_noah271lsm, only : csnow, tdfcnd, penman

! !ARGUMENTS: 
  implicit none

  real, parameter         :: sbeta = -2.0
  real, parameter         :: sigma = 5.67E-8
  real, parameter         :: cpice = 2.106E+3
  real, parameter         :: cph2o = 4.218E+3
  real, parameter         :: LSUBC = 2.501000E+6
  real, parameter         :: LSUBS = 2.83E+6
  integer, intent(in)     :: n
  real, intent(in)        :: stc(LIS_rc%npatch(n,LIS_rc%lsm_index))
!  real, intent(out)       :: t1(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                    :: t1(LIS_rc%npatch(n,LIS_rc%lsm_index))
!
! !DESCRIPTION: 
!   Diagnostic routine to update noah variables based on an updated 
!   T1 value
! 
!EOP
  integer                 :: t, kz
  integer                 :: nsoil
  real                    :: zsoil(4)
  real                    :: sfctmp
  real                    :: fdown
  real                    :: t24
  real                    :: df1
  real                    :: yynum
  real                    :: esat
  real                    :: q2sat
  real                    :: q2
  real                    :: th2
  real                    :: t2v
  real                    :: ssoil
  real                    :: etp
  logical*1               :: snowng
  logical*1               :: frzgra
  real                    :: ffrozp
  real                    :: flx1
  real                    :: flx2
  real                    :: dqsdt2
  real                    :: dqsdt
  real                    :: rr
  real                    :: yy
  real                    :: zz1
  real                    :: beta
  real                    :: rch
  real                    :: e
  real                    :: epsca
  real                    :: t12a
  real                    :: etanrg
  real                    :: df1a
  real                    :: df1h
  real                    :: t12
  real                    :: t12b
  real                    :: denom
  real                    :: frcsoi
  real                    :: frcsno
  real                    :: dtot
  real                    :: dsoil
  real                    :: sncond
  real                    :: sndens
!  real                    :: csnow

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
!     if(t.eq.641) then 
!        print*, 'bef upd ',t, stc(t), t1(t)
!     endif
     nsoil = noah271_struc(n)%nslay
     sfctmp = noah271_struc(n)%noah(t)%tair

     snowng = .false.
     ffrozp = 0.0
     if(noah271_struc(n)%noah(t)%rainf.gt.0.0) then 
        if(sfctmp.le.LIS_CONST_TKFRZ) then 
           ffrozp = 1.0
        endif
     endif
     
     if(noah271_struc(n)%noah(t)%rainf.gt.1e-20) then 
        if(ffrozp.gt.0.5) then 
           snowng = .true. 
        else
           if(t1(t).le.LIS_CONST_TKFRZ) frzgra = .true.
        endif
     endif
     if(noah271_struc(n)%noah(t)%sneqv.eq.0.0) then 
        sncond = 1.0       
     else
        sndens = noah271_struc(n)%noah(t)%sneqv/noah271_struc(n)%noah(t)%snowh
!        sncond = csnow(sndens)
        call csnow(sndens, sncond)
     endif

     if(noah271_struc(n)%noah(t)%vegt.eq.LIS_rc%waterclass) then 
        do kz = 1,nsoil
           zsoil(kz) = -3.0*float(kz)/float(nsoil)
        enddo
     else
        zsoil(1) = -noah271_struc(n)%lyrthk(1)
        do kz = 2, nsoil
           zsoil(kz) = -noah271_struc(n)%lyrthk(kz)+zsoil(kz-1)
        enddo
     endif

!nopac
     fdown = noah271_struc(n)%noah(t)%swdown * &
          (1- noah271_struc(n)%noah(t)%albedo) + &
          noah271_struc(n)%noah(t)%lwdown
     t24    = sfctmp * sfctmp * sfctmp* sfctmp

     call tdfcnd(df1,noah271_struc(n)%noah(t)%smc(1),&
          noah271_struc(n)%noah(t)%quartz,&
          noah271_struc(n)%noah(t)%smcmax,&
          noah271_struc(n)%noah(t)%sh2o(1))
     df1 = df1* exp(sbeta*LIS_gfrac(n)%greenness(t))
     dqsdt2 = dqsdt(sfctmp, noah271_struc(n)%noah(t)%psurf)
     q2 = noah271_struc(n)%noah(t)%qair
     esat = e(sfctmp)
     q2sat = 0.622 * esat / (noah271_struc(n)%noah(t)%psurf&
          - (1-0.622)*esat)
     if(q2.ge.q2sat) q2 = q2sat * 0.99

     th2 = sfctmp + 0.0098*noah271_struc(n)%noah(t)%z
     t2v = sfctmp * (1.0 + 0.61 * Q2)     
     yynum = fdown - sigma * t24

!snow-based update. 
     dsoil = -(0.5*zsoil(1))
     if(noah271_struc(n)%noah(t)%sneqv.eq.0) then 
        ssoil = df1*(t1(t)-stc(t))/dsoil
     else
        dtot = noah271_struc(n)%noah(t)%snowh + dsoil
        frcsno = noah271_struc(n)%noah(t)%snowh/dtot
        frcsoi = dsoil/dtot
        df1h = (sncond*df1)/(frcsoi*sncond+frcsno*df1)
        df1a = (frcsno*sncond + frcsoi *df1)
        df1 = df1a * noah271_struc(n)%noah(t)%sca + &
             df1*(1.0 - noah271_struc(n)%noah(t)%sca)
        ssoil = df1*(t1(t)-stc(t))/dtot
     endif

     call penman(sfctmp, noah271_struc(n)%noah(t)%psurf, &
          noah271_struc(n)%noah(t)%ch, t2v, th2, &
          noah271_struc(n)%noah(t)%rainf, fdown, t24, ssoil, &
          q2, q2sat, etp, rch, epsca, rr, snowng, frzgra, dqsdt2, flx2,&
          noah271_struc(n)%noah(t)%emiss)
! ----------------------------------------------------------------------
! BASED ON ETP AND E VALUES, DETERMINE BETA
! ----------------------------------------------------------------------
     if(noah271_struc(n)%noah(t)%sneqv.eq.0) then 
        beta = noah271_struc(n)%noah(t)%beta
        yy = sfctmp + (yynum/rch + th2 - sfctmp - beta*epsca)/ rr
        zz1 = df1 /(-0.5 * zsoil(1)*rch*rr)+ 1.0
        t1(t) = (yy + (zz1 - 1.0) * stc(t)) /zz1
     else
!computation of etanrg?
!        etp1 = etp*0.001
!        if(etp.lt.0.0) then 
!           dew = -etp1
!           esnow2 = etp1*LIS_rc%ts
!           etanrg = etp* ((1-noah271_struc(n)%noah(t)%sca)*lsubc + & 
!                noah271_struc(n)%noah(t)%sca*lsubs)
!        else
!           print*, 'stopping in noah271_calc_tskin'
!           print*, t, noah271_struc(n)%noah(t)%sca
!           stop
!        endif

        etanrg = noah271_struc(n)%noah(t)%etanrg
        flx1 = 0.0
        if(snowng) then 
           flx1 = cpice * noah271_struc(n)%noah(t)%rainf * (t1(t)-sfctmp)
        else
           if(noah271_struc(n)%noah(t)%rainf.gt.0.0) &
                flx1 = cph2o* noah271_struc(n)%noah(t)%rainf * (t1(t)-sfctmp)
        endif
        dtot = noah271_struc(n)%noah(t)%snowh + dsoil
        denom = 1.0 + df1 /(dtot * rr * rch)
        t12a = ((fdown - flx1 -flx2 - & 
             (0.95 * noah271_struc(n)%noah(t)%sca + &
             (1.0 - noah271_struc(n)%noah(t)%sca))*sigma * t24)/rch &
             + th2 - sfctmp - etanrg /rch)/rr
        t12b = df1*stc(t)/(dtot * rr * rch)
        t12 = (sfctmp + t12a + t12b) /denom
        if(t12.le.LIS_CONST_TKFRZ) then 
           t1(t) = t12
        else
           t1(t) = LIS_CONST_TKFRZ*noah271_struc(n)%noah(t)%sca**2.0 + &
                t12 * (1.0 - noah271_struc(n)%noah(t)%sca*2.0)
        endif
     endif
!     if(t.eq.641) then 
!        print*, 'aft upd ',t, stc(t), t1(t)
!     endif
  enddo

end subroutine noah271_calc_tskin

