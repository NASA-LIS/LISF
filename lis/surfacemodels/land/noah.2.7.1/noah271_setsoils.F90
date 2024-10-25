!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!#include "LIS_misc.h"
!BOP
!
! !ROUTINE: noah271_setsoils
!  \label{noah271_setsoils}
!
! !REVISION HISTORY:
!  28 Apr 2002: Kristi Arsenault;  Added Noah2.5 LSM, Initial Code
!  13 Oct 2003: Sujay Kumar; Domain independent modifications
!  15 Sep 2005: Matthew Garcia; update for ARMS project (Windows capability)
!  04 Oct 2005: Matthew Garcia; additions for Noah Unified
!   9 Jul 2008: Chris Franks/AFWA; Save wilting point for calculation of relative soil moisture
!  27 Oct 2010: David Mocko, changes for Noah2.7.1 in LIS6.1
!
! !INTERFACE:
subroutine noah271_setsoils(mtype)
! !USES:
  use LIS_coreMod
  use LIS_soilsMod,  only : LIS_soils
  use LIS_logMod,    only : LIS_logunit, LIS_endrun
  use noah271_lsmMod     

!
! !DESCRIPTION:
!  This subroutine sets the soil parameters in Noah2.7.1.  Noah2.7.1 uses a
!  lookup table based on the input soil texture to derive these parameters. 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[LIS\_mapSoilType](\ref{LIS_mapSoilType}) \newline
!     method to derive the soil texture type from the sand, silt, and 
!     clay fractions.
!   \end{description}
!EOP      
  implicit none
  integer :: mtype
  integer :: jj,i,n
  real, allocatable :: basicset(:,:)
!  real :: maxdt,mindz,tiledt,psis,ks,b
!  integer :: newdt,idiv

  real  :: sand, silt, clay
  real  :: sandf, siltf, clayf
  real  :: poros, bexp, psisat1, dksat1
  real  :: soildepths(4)
  real  :: czil
  real, allocatable :: temp(:,:)
  integer       :: c,r

  do n=1,LIS_rc%nnest

     allocate(basicset(noah271_struc(n)%nstxts,noah271_struc(n)%nsoilp))
     write(unit=LIS_logunit,fmt=*)                                     &
                           'MSG: noah271_setsoils -- reading soil files'
     
     if(trim(LIS_rc%usetexturemap(n)).eq."none") then
        do i=1,LIS_rc%npatch(n,mtype)
           call LIS_mapSoilType(noah271_struc(n)%soilscheme,&
                LIS_surface(n,mtype)%tile(i)%sand, &
                LIS_surface(n,mtype)%tile(i)%clay, &
                LIS_surface(n,mtype)%tile(i)%silt,&
                LIS_surface(n,mtype)%tile(i)%soilt)
        enddo
     endif

     do i=1,LIS_rc%npatch(n,mtype)
        noah271_struc(n)%noah(i)%soiltype = &
             LIS_surface(n,mtype)%tile(i)%soilt
        if ( noah271_struc(n)%noah(i)%vegt .eq. LIS_rc%waterclass ) then
           if ( noah271_struc(n)%soilscheme == 2 ) then ! statsgo
              noah271_struc(n)%noah(i)%soiltype = 14 ! water
             LIS_surface(n,mtype)%tile(i)%soilt = 14
           else
              noah271_struc(n)%noah(i)%soiltype = 1 ! sand
              LIS_surface(n,mtype)%tile(i)%soilt = 1
           endif
        else
           if ( noah271_struc(n)%noah(i)%soiltype .eq. -9999 .or. &
                noah271_struc(n)%noah(i)%soiltype .eq. 0 .or. &
                noah271_struc(n)%noah(i)%soiltype .eq. 14 ) then
              if ( noah271_struc(n)%soilscheme == 2 ) then ! statsgo
                 noah271_struc(n)%noah(i)%soiltype = 1 ! sand
                 LIS_surface(n,mtype)%tile(i)%soilt = 1
              else ! zobler
                 noah271_struc(n)%noah(i)%soiltype = 1 ! loamy sand
                 LIS_surface(n,mtype)%tile(i)%soilt = 1

              endif
           endif
        endif
     enddo
!Hardcoded for STATSGO classes only for soil types found to be water
! at a land point

!-----------------------------------------------------------------------
! Read in the Noah2.7.1 Soil Parameter File
!-----------------------------------------------------------------------
     write(LIS_logunit,*) 'Reading Noah2.7.1 soil parameter file: ',   &
                                            trim(noah271_struc(n)%sfile)
     open(unit=18,file=noah271_struc(n)%sfile,status='old',            &
                                                    access='sequential')
     read(18,*)
! Read soil depths from parameter file but do nothing with them.
! Soil layers/depths are now set in lis.config file instead.  - D. Mocko
     read(18,*) soildepths
     do i=1,noah271_struc(n)%nsoilp
        read(18,*)
        read(18,*)(basicset(jj,i),jj=1,noah271_struc(n)%nstxts)
     enddo
     !read the czil value. 
     read(18,*) 
     read(18,*) czil
     close(18)
     
!-----------------------------------------------------------------------
! Assign SOIL Parameters to each tile based on the
! type of Zobler soil class present in that tile.
!-----------------------------------------------------------------------
! Set tile soil porosity
     if (trim(LIS_rc%useporositymap(n)).eq."none") then ! default, from look-up table
        do i = 1,LIS_rc%npatch(n,mtype)
           noah271_struc(n)%noah(i)%smcmax = basicset(noah271_struc(n)%noah(i)%soiltype,1)
           noah271_struc(n)%noah(i)%smcwlt = basicset(noah271_struc(n)%noah(i)%soiltype,7)
        end do
     else
        do i = 1,LIS_rc%npatch(n,mtype)
           noah271_struc(n)%noah(i)%smcmax = &
                LIS_soils(n)%porosity(LIS_surface(n,mtype)%tile(i)%col,&
                LIS_surface(n,mtype)%tile(i)%row,1)
        end do
     end if
     
! Set tile soil Psi-sat
     if (trim(LIS_rc%usepsisatmap(n)).eq."none") then ! default, from look-up table
        do i = 1,LIS_rc%npatch(n,mtype)
           noah271_struc(n)%noah(i)%psisat = &
                basicset(noah271_struc(n)%noah(i)%soiltype,2)
        end do
     else
        do i = 1,LIS_rc%npatch(n,mtype)
           noah271_struc(n)%noah(i)%psisat = &
                LIS_soils(n)%psisat(LIS_surface(n,mtype)%tile(i)%col,&
                LIS_surface(n,mtype)%tile(i)%row)
        end do
     end if

! Set tile soil K-sat
     if (trim(LIS_rc%useksatmap(n)).eq."none") then ! default, from look-up table
        do i = 1,LIS_rc%npatch(n,mtype)
           noah271_struc(n)%noah(i)%dksat = &
                basicset(noah271_struc(n)%noah(i)%soiltype,3)
        end do
     else
        do i = 1,LIS_rc%npatch(n,mtype)
           noah271_struc(n)%noah(i)%dksat = &
                LIS_soils(n)%ksat(LIS_surface(n,mtype)%tile(i)%col,&
                LIS_surface(n,mtype)%tile(i)%row)
        end do
     end if
     
! Set tile soil b-parameter
     if (trim(LIS_rc%usebexpmap(n)).eq."none") then ! default, from look-up table
        do i = 1,LIS_rc%npatch(n,mtype)
           noah271_struc(n)%noah(i)%bexp = &
                basicset(noah271_struc(n)%noah(i)%soiltype,4)
        end do
     else
        do i = 1,LIS_rc%npatch(n,mtype)
           noah271_struc(n)%noah(i)%bexp = &
                LIS_soils(n)%bexp(LIS_surface(n,mtype)%tile(i)%col,&
                LIS_surface(n,mtype)%tile(i)%row)
        end do
     end if
     
! Set tile soil Qz content
     if (trim(LIS_rc%usequartzmap(n)).eq."none") then ! default
        do i = 1,LIS_rc%npatch(n,mtype)
           noah271_struc(n)%noah(i)%quartz = &
                basicset(noah271_struc(n)%noah(i)%soiltype,5)
        end do
     else 
        do i = 1,LIS_rc%npatch(n,mtype)
           noah271_struc(n)%noah(i)%quartz = &
                LIS_soils(n)%quartz(LIS_surface(n,mtype)%tile(i)%col,&
                LIS_surface(n,mtype)%tile(i)%row)
        end do
     end if
#if 0 
!temporary ... 
     allocate(temp(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     temp = -9999.0
     do i=1,LIS_rc%npatch(n,mtype)
        c = LIS_surface(n,mtype)%tile(i)%col
        r = LIS_surface(n,mtype)%tile(i)%row
        temp(c,r) = noah271_struc(n)%noah(i)%smcmax
     enddo
     open(100,file='smcmax_def.1gd4r',form='unformatted',              &
          access='direct',recl=LIS_rc%lnc(n)*LIS_rc%lnr(n)*4)
     write(100,rec=1) temp
     close(100)


     temp = -9999.0
     do i=1,LIS_rc%npatch(n,mtype)
        c = LIS_surface(n,mtype)%tile(i)%col
        r = LIS_surface(n,mtype)%tile(i)%row
        temp(c,r) = noah271_struc(n)%noah(i)%psisat
     enddo
     open(100,file='psisat_def.1gd4r',form='unformatted',              &
          access='direct',recl=LIS_rc%lnc(n)*LIS_rc%lnr(n)*4)
     write(100,rec=1) temp
     close(100)

     temp = -9999.0
     do i=1,LIS_rc%npatch(n,mtype)
        c = LIS_surface(n,mtype)%tile(i)%col
        r = LIS_surface(n,mtype)%tile(i)%row
        temp(c,r) = noah271_struc(n)%noah(i)%dksat
     enddo
     open(100,file='dksat_def.1gd4r',form='unformatted',               &
          access='direct',recl=LIS_rc%lnc(n)*LIS_rc%lnr(n)*4)
     write(100,rec=1) temp
     close(100)

     temp = -9999.0
     do i=1,LIS_rc%npatch(n,mtype)
        c = LIS_surface(n,mtype)%tile(i)%col
        r = LIS_surface(n,mtype)%tile(i)%row
        temp(c,r) = noah271_struc(n)%noah(i)%dwsat
     enddo
     open(100,file='dwsat_def.1gd4r',form='unformatted',               &
          access='direct',recl=LIS_rc%lnc(n)*LIS_rc%lnr(n)*4)
     write(100,rec=1) temp
     close(100)
     
     temp = -9999.0
     do i=1,LIS_rc%npatch(n,mtype)
        c = LIS_surface(n,mtype)%tile(i)%col
        r = LIS_surface(n,mtype)%tile(i)%row
        temp(c,r) = LIS_soils(n)%texture(c,r)
     enddo
     open(100,file='texture_def.1gd4r',form='unformatted',             &
          access='direct',recl=LIS_rc%lnc(n)*LIS_rc%lnr(n)*4)
     write(100,rec=1) temp
     close(100)
     stop
#endif
     
! standard values: 
     do i = 1,LIS_rc%npatch(n,mtype)
        noah271_struc(n)%noah(i)%sbeta = -2.0
        noah271_struc(n)%noah(i)%czil = czil
        noah271_struc(n)%noah(i)%fxexp = 2.0
        noah271_struc(n)%noah(i)%cfactr = 0.5
        noah271_struc(n)%noah(i)%cmcmax = 0.5E-3
        noah271_struc(n)%noah(i)%rsmax = 5000.0
        noah271_struc(n)%noah(i)%topt = 298.0
        noah271_struc(n)%noah(i)%refdk = 2.0E-6
        noah271_struc(n)%noah(i)%refkdt = 3.0
        noah271_struc(n)%noah(i)%csoil = 2.00E6
        noah271_struc(n)%noah(i)%frzk = 0.15
        noah271_struc(n)%noah(i)%emiss = 1.0
     end do

! overwrite default values with PTF 
     if(noah271_struc(n)%useptf.eq.1) then !use cosby et al.(1984) PTFs
        do i = 1,LIS_rc%npatch(n,mtype)
           
           noah271_struc(n)%noah(i)%smcmax = 0.489 - &
                0.00126*LIS_surface(n,mtype)%tile(i)%sand*100.0

           noah271_struc(n)%noah(i)%psisat = (10.0*(10.0**(1.88-0.0131*&
                LIS_surface(n,mtype)%tile(i)%sand*100.0)))/1000.0

           noah271_struc(n)%noah(i)%dksat =  (0.0070556*(10.0**(-0.884+0.0153*&
                LIS_surface(n,mtype)%tile(i)%sand)))/1000.0

           noah271_struc(n)%noah(i)%bexp = 2.91 + (0.159*&
                LIS_surface(n,mtype)%tile(i)%clay*100.0)

        end do
     elseif(noah271_struc(n)%useptf.eq.2) then !use Rawls & Brakensiek (1985)
        do i=1, LIS_rc%npatch(n,mtype)
           sandf = LIS_surface(n,mtype)%tile(i)%sand
           clayf = LIS_surface(n,mtype)%tile(i)%clay
           siltf = LIS_surface(n,mtype)%tile(i)%silt

           sand = sandf*100.0
           clay = clayf*100.0
           silt = siltf*100.0
           poros = noah271_struc(n)%noah(i)%smcmax

           dksat1 =  &
                exp((19.52348 * poros) - 8.96847 - (0.028212 * clay)     &
                + (0.00018107 * sand**2.0) - (0.0094125 * clay**2.0)     &
                - (8.395215 * poros**2.0) + (0.077718 * sand * poros)    &
                - (0.00298 * sand**2.0 * poros**2.0) - (0.019492         &
                * clay**2.0 * poros**2.0) + (0.0000173 * sand**2.0       &
                * clay) + (0.02733 * clay**2.0 * poros) + (0.001434      &
                * sand**2.0 * poros) - (0.0000035 * clay**2.0 * sand))

           noah271_struc(n)%noah(i)%dksat = dksat1*10 /3600.0 /1000.0
           
           poros = noah271_struc(n)%noah(i)%smcmax *0.9
           
           psisat1 = exp(6.5309 - (7.32561 * poros) + (0.001583 * clay**2.0) &
                + (3.809479 * poros**2.0) + (0.000344 * sand * clay)   &
                - (0.049837 * sand * poros) + (0.001608 * sand**2.0    &
                * poros**2.0) + (0.001602 * clay**2.0 * poros**2.0)    &
                - (0.0000136 * sand**2.0 * clay) - (0.0003479          &
                * clay**2.0 * poros) - (0.000799 * sand**2.0 * poros))

           noah271_struc(n)%noah(i)%psisat = psisat1/100.0

           bexp = exp(-0.7842831 + (0.0177544 * sand) - (1.062498 * poros)  &
              - (0.00005304 * sand**2.0) - (0.00273493 * clay**2.0)    &
              + (1.11134946 * poros**2.0) - (0.03088295 * sand         &
              * poros) + (0.00026587 * sand**2.0 * poros**2.0)         &
              - (0.00610522 * clay**2.0 * poros**2.0) - (0.00000235    &
              * sand**2.0 * clay) + (0.00798746 * clay**2.0 * poros)   &
              - (0.00674491 * poros**2.0 * clay))

           noah271_struc(n)%noah(i)%bexp = 1.0 / BEXP
        enddo
     endif
     
! site5-optimized values - GA + likelihood function
!     do i=1, LIS_rc%npatch(n,mtype)
!        noah271_struc(n)%noah(i)%smcmax = 0.322
!        noah271_struc(n)%noah(i)%psisat = 0.1436
!        noah271_struc(n)%noah(i)%dksat = 4.833e-06
!        noah271_struc(n)%noah(i)%bexp = 4.494
!     enddo

! site1-cosby optimized values GA
!     do i=1, LIS_rc%npatch(n,mtype)
!        noah271_struc(n)%noah(i)%smcmax = 0.187674
!        noah271_struc(n)%noah(i)%psisat = 0.0101895
!        noah271_struc(n)%noah(i)%dksat = 5.00595e-07
!        noah271_struc(n)%noah(i)%bexp = 7.44191
!        noah271_struc(n)%noah(i)%quartz = 0.896606
! site 5 cosby optimized values
!        noah271_struc(n)%noah(i)%smcmax = 0.231382
!        noah271_struc(n)%noah(i)%psisat = 0.0105896
!        noah271_struc(n)%noah(i)%dksat = 5.03571e-07
!        noah271_struc(n)%noah(i)%bexp = 3.385806
!        noah271_struc(n)%noah(i)%quartz = 0.120313
!     enddo

! site1-cosby optimized values LM
!     do i=1, LIS_rc%npatch(n,mtype)
!        noah271_struc(n)%noah(i)%smcmax = 0.197421
!        noah271_struc(n)%noah(i)%psisat = 0.011878
!        noah271_struc(n)%noah(i)%dksat = 5.13e-07
!        noah271_struc(n)%noah(i)%bexp = 3.138477
!        noah271_struc(n)%noah(i)%quartz = 0.899985
! site 5 cosby optimized values
!        noah271_struc(n)%noah(i)%smcmax = 0.234308
!        noah271_struc(n)%noah(i)%psisat = 0.0903
!        noah271_struc(n)%noah(i)%dksat = 5.12e-07
!        noah271_struc(n)%noah(i)%bexp = 7.963148
!        noah271_struc(n)%noah(i)%quartz = 0.102386
!     enddo
! site1-cosby optimized values SCE
!     do i=1, LIS_rc%npatch(n,mtype)
!        noah271_struc(n)%noah(i)%smcmax = 0.195243
!        noah271_struc(n)%noah(i)%psisat = 1.16E-2
!        noah271_struc(n)%noah(i)%dksat = 5.08e-07
!        noah271_struc(n)%noah(i)%bexp = 4.52
!        noah271_struc(n)%noah(i)%quartz = 0.554
! site 5 cosby optimized values
!        noah271_struc(n)%noah(i)%smcmax = 0.242
!        noah271_struc(n)%noah(i)%psisat = 5.62E-2
!        noah271_struc(n)%noah(i)%dksat = 7.09E-7
!        noah271_struc(n)%noah(i)%bexp = 5.146213
!        noah271_struc(n)%noah(i)%quartz = 0.552939
!     enddo


!  call soils_finalize()

!-----------------------------------------------------------------------
! Check numerical stability at saturation on basis of timestep and 
! soil parameters.  From formula derived by Yihua Wu, NASA/GSFC HSB
!-----------------------------------------------------------------------
!This routine as of now is "unsafe" in a parallel mode since each processor
!can change the timestep. There needs to be some gather/scatter...
#if 0 
     maxdt = LIS_rc%ts
     mindz = noahattrib%lyrthk(1) 
     do i = 2,noah271_struc(n)%nslay ! soil layer loop
        mindz = min(mindz,noahattrib%lyrthk(i))
     end do ! i
     do i = 1,LIS_rc%npatch(n,mtype)          ! tile loop
! NOTE: depends on *positive* Psi-sat values, as in Noah2.7.1 LSM
        psis = noah271_struc(n)%noah(i)%psisat
        ks = noah271_struc(n)%noah(i)%dksat
        b = noah271_struc(n)%noah(i)%bexp
        tiledt = 2 * mindz**2 / (b * psis * ks)
        maxdt = min(maxdt,tiledt) 
     end do ! i
     if (LIS_rc%ts.gt.maxdt) then
        write(LIS_logunit,*) 'WRN: noah271_setsoils -- timestep may',  &
                             'be too large for numerical stability'
        write(LIS_logunit,*) 'WRN: noah271_setsoils -- soil parameters',&
                             'suggest dt .le.',maxdt
        newdt = 3600
        idiv = 1
        do while (newdt.gt.maxdt)
           if (mod(3600,idiv).eq.0) then
              newdt = 3600 / idiv
           endif
           idiv = idiv + 1
           if (idiv.gt.60) exit
        enddo
        if (newdt.gt.maxdt) then
           write(LIS_logunit,*) 'MSG: noah271_setsoils -- new timestep',&
                                'is less than 1 minute -- stopping'
           call LIS_endrun
        else
           write(LIS_logunit,*) 'MSG: noah271_setsoils -- setting',    &
                                'timestep to ',newdt,' sec'
           LIS_rc%ts = newdt
        end if
     else
        write(LIS_logunit,*) 'MSG: noah271_setsoils -- timestep OK',   &
                      'for soil hydraulic stability (',LIS_rc%ts,' sec)'
     end if
#endif
     deallocate(basicset)
  enddo
  return

end subroutine noah271_setsoils
