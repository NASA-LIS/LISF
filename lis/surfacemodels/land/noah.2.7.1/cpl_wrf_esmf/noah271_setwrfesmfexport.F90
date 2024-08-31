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
! !ROUTINE: noah271_setwrfesmfexport.F90
!
! !DESCRIPTION:  
!  Defines the export states from Noah to WRF in coupled mode
!
! !REVISION HISTORY:
! 02 Dec 2003; Sujay Kumar, Initial Version
! 17 Nov 2008; Sujay Kumar, Modified for the ESMF coupled version
! 
! !INTERFACE:
subroutine noah271_setwrfesmfexport(n, LISWRF_ExpState)
! !USES:
  use ESMF
  use LIS_coreMod,    only : LIS_rc, LIS_domain
  use LIS_historyMod, only : LIS_tile2grid
  use LIS_logMod,     only : LIS_verify
!  use module_lis_component, only : LISWRF_ExpState
  use noah271_lsmMod
  
  implicit none
  integer, intent(in) :: n 
  type(ESMF_State)    :: LISWRF_ExpState
!EOP

  integer          :: c,r, t,gid
  type(ESMF_Field) :: qleField, qhField, qgField, avgsurftField, albedoField, &
       emissField, z0Field, xiceField, qsfcField, cqs2Field, chs2Field, &
       t2Field, th2Field,q2Field, psfcField,qfxField,chField
  real, allocatable    :: qle(:,:), qh(:,:), qg(:,:), avgsurft(:,:), albedo(:,:),&
       emiss(:,:), z0(:,:), xice(:,:), qsfc(:,:), cqs2(:,:), chs2(:,:), &
       t2(:,:),q2(:,:),th2(:,:),psfc(:,:),qfx(:,:),ch(:,:)
  integer          :: rc
  
  integer :: nobjs
  
  call ESMF_StateGet(LISWRF_expState, 'Latent Heat Flux', qlefield, rc=rc)
  call LIS_verify(rc,'StateGet:Qle in noah271_setwrfesmfexport failed') 
  call ESMF_FieldGet(qlefield,localDE=0,farrayPtr=qle,rc=rc)
  call LIS_verify(rc,'FieldGet:qle in noah271_setwrfesmfexport failed')


  call ESMF_StateGet(LISWRF_expState, 'Sensible Heat Flux', qhfield, rc=rc)
  call LIS_verify(rc,'StateGet:Qh in noah271_setwrfesmfexport failed') 
  call ESMF_FieldGet(qhfield,localDE=0,farrayPtr=qh,rc=rc)
  call LIS_verify(rc,'FieldGet:qh in noah271_setwrfesmfexport failed')


  call ESMF_StateGet(LISWRF_expState, 'Ground Heat Flux', qgfield, rc=rc)
  call LIS_verify(rc,'StateGet:Qg in noah271_setwrfesmfexport failed') 
  call ESMF_FieldGet(qgfield,localDE=0,farrayPtr=qg,rc=rc)
  call LIS_verify(rc,'FieldGet:qg in noah271_setwrfesmfexport failed')


  call ESMF_StateGet(LISWRF_expState, 'Surface Skin Temperature', avgsurftfield, rc=rc)
  call LIS_verify(rc,'StateGet:Avgsurft in noah271_setwrfesmfexport failed') 
  call ESMF_FieldGet(avgsurftfield,localDE=0,farrayPtr=avgsurft,rc=rc)
  call LIS_verify(rc,'FieldGet:avgsurft in noah271_setwrfesmfexport failed')


  call ESMF_StateGet(LISWRF_expState, 'Surface Albedo', albedofield, rc=rc)
  call LIS_verify(rc,'StateGet:Albedo in noah271_setwrfesmfexport failed') 
  call ESMF_FieldGet(albedofield,localDE=0,farrayPtr=albedo,rc=rc)
  call LIS_verify(rc,'FieldGet:albedo in noah271_setwrfesmfexport failed')


  call ESMF_StateGet(LISWRF_expState, 'Surface Emissivity', emissfield, rc=rc)
  call LIS_verify(rc,'StateGet:Emiss in noah271_setwrfesmfexport failed') 
  call ESMF_FieldGet(emissfield,localDE=0,farrayPtr=emiss,rc=rc)
  call LIS_verify(rc,'FieldGet:emiss in noah271_setwrfesmfexport failed')


  call ESMF_StateGet(LISWRF_expState, 'Surface Roughness', z0field, rc=rc)
  call LIS_verify(rc,'StateGet:Z0 in noah271_setwrfesmfexport failed') 
  call ESMF_FieldGet(z0field,localDE=0,farrayPtr=z0,rc=rc)
  call LIS_verify(rc,'FieldGet:z0 in noah271_setwrfesmfexport failed')

  call ESMF_StateGet(LISWRF_expState, 'Sea Ice Mask', xicefield, rc=rc)
  call LIS_verify(rc,'StateGet:XICE in noah271_setwrfesmfexport failed') 
  call ESMF_FieldGet(xicefield,localDE=0,farrayPtr=xice,rc=rc)
  call LIS_verify(rc,'FieldGet:xice in noah271_setwrfesmfexport failed')

  call ESMF_StateGet(LISWRF_expState, 'Surface Specific Humidity', qsfcfield, rc=rc)
  call LIS_verify(rc,'StateGet:QSFC in noah271_setwrfesmfexport failed') 
  call ESMF_FieldGet(qsfcfield,localDE=0,farrayPtr=qsfc,rc=rc)
  call LIS_verify(rc,'FieldGet:qsfc in noah271_setwrfesmfexport failed')

  call ESMF_StateGet(LISWRF_expState, '2m Surface Exchange Coefficient for Heat', &
       chs2field, rc=rc)
  call LIS_verify(rc,'StateGet:CHS2 in noah271_setwrfesmfexport failed') 
  call ESMF_FieldGet(chs2field,localDE=0,farrayPtr=chs2,rc=rc)
  call LIS_verify(rc,'FieldGet:chs2 in noah271_setwrfesmfexport failed')

  call ESMF_StateGet(LISWRF_expState, '2m Surface Exchange Coefficient for Moisture', &
       cqs2field, rc=rc)
  call LIS_verify(rc,'StateGet:CQS2 in noah271_setwrfesmfexport failed') 
  call ESMF_FieldGet(cqs2field,localDE=0,farrayPtr=cqs2,rc=rc)
  call LIS_verify(rc,'FieldGet:cqs2 in noah271_setwrfesmfexport failed')


  call ESMF_StateGet(LISWRF_expState, '2m Air Temperature', &
       t2field, rc=rc)
  call LIS_verify(rc,'StateGet:T2 in noah271_setwrfesmfexport failed') 
  call ESMF_FieldGet(t2field,localDE=0,farrayPtr=t2,rc=rc)
  call LIS_verify(rc,'FieldGet:t2 in noah271_setwrfesmfexport failed')

  call ESMF_StateGet(LISWRF_expState, '2m Specific Humidity', &
       q2field, rc=rc)
  call LIS_verify(rc,'StateGet:Q2 in noah271_setwrfesmfexport failed') 
  call ESMF_FieldGet(q2field,localDE=0,farrayPtr=q2,rc=rc)
  call LIS_verify(rc,'FieldGet:q2 in noah271_setwrfesmfexport failed')

  call ESMF_StateGet(LISWRF_expState, '2m Potential Temperature', &
       th2field, rc=rc)
  call LIS_verify(rc,'StateGet:TH2 in noah271_setwrfesmfexport failed') 
  call ESMF_FieldGet(th2field,localDE=0,farrayPtr=th2,rc=rc)
  call LIS_verify(rc,'FieldGet:th2 in noah271_setwrfesmfexport failed')

  call ESMF_StateGet(LISWRF_expState, 'Surface Pressure', &
       psfcfield, rc=rc)
  call LIS_verify(rc,'StateGet:PSFC in noah271_setwrfesmfexport failed') 
  call ESMF_FieldGet(psfcfield,localDE=0,farrayPtr=psfc,rc=rc)
  call LIS_verify(rc,'FieldGet:psfc in noah271_setwrfesmfexport failed')

  call ESMF_StateGet(LISWRF_expState, 'Upward Moisture Flux at the Surface', &
       qfxfield, rc=rc)
  call LIS_verify(rc,'StateGet:QFX in noah271_setwrfesmfexport failed') 
  call ESMF_FieldGet(qfxfield,localDE=0,farrayPtr=qfx,rc=rc)
  call LIS_verify(rc,'FieldGet:qfx in noah271_setwrfesmfexport failed')

  call ESMF_StateGet(LISWRF_expState, 'Surface Exchange Coefficient for Heat', &
       chfield, rc=rc)
  call LIS_verify(rc,'StateGet:CH in noah271_setwrfesmfexport failed') 
  call ESMF_FieldGet(chfield,localDE=0,farrayPtr=ch,rc=rc)
  call LIS_verify(rc,'FieldGet:ch in noah271_setwrfesmfexport failed')

  call LIS_tile2grid(n,qle,noah271_struc(n)%noah%qle)
  call LIS_tile2grid(n,qh,noah271_struc(n)%noah%qh)
  call LIS_tile2grid(n,qg,noah271_struc(n)%noah%qg)
  call LIS_tile2grid(n,avgsurft,noah271_struc(n)%noah%t1)
  call LIS_tile2grid(n,albedo,noah271_struc(n)%noah%albedo)
  call LIS_tile2grid(n,emiss,noah271_struc(n)%noah%emiss)
  call LIS_tile2grid(n,xice,noah271_struc(n)%noah%xice)
  call LIS_tile2grid(n,z0,noah271_struc(n)%noah%z0)
  call LIS_tile2grid(n,qsfc,noah271_struc(n)%noah%qsfc)
  call LIS_tile2grid(n,chs2,noah271_struc(n)%noah%chs2)
  call LIS_tile2grid(n,cqs2,noah271_struc(n)%noah%cqs2)
  call LIS_tile2grid(n,t2,noah271_struc(n)%noah%t2)
  call LIS_tile2grid(n,th2,noah271_struc(n)%noah%th2)
  call LIS_tile2grid(n,q2,noah271_struc(n)%noah%q2)
  call LIS_tile2grid(n,psfc,noah271_struc(n)%noah%psurf)
  call LIS_tile2grid(n,qfx,noah271_struc(n)%noah%eta_kinematic)
  call LIS_tile2grid(n,ch,noah271_struc(n)%noah%ch)

!  do r=1,LIS_rc%lnr(n)
!     do c=1,LIS_rc%lnc(n)
!        gid = LIS_domain(n)%gindex(c,r)
!        print*, 'NOAH ',z0(c,r)
!        print*, 'E ',qle(c,r), qh(c,r), qg(c,r),avgsurft(c,r), emiss(c,r), albedo(c,r), xice(c,r), qsfc(c,r), chs2(c,r), cqs2(c,r), t2(c,r), th2(c,r), q2(c,r)
!     enddo
!  enddo

#if 0 
! trying this for now: (assumes that ngrid is same as ntiles)
  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        gid = LIS_domain(n)%gindex(c,r)
        qle(c,r) = noah271_struc(n)%noah(gid)%qle
        qh(c,r) = noah271_struc(n)%noah(gid)%qh        
        qg(c,r) = noah271_struc(n)%noah(gid)%qg
        avgsurft(c,r) = noah271_struc(n)%noah(gid)%t1
        albedo(c,r) = noah271_struc(n)%noah(gid)%albedo
        emiss(c,r) = noah271_struc(n)%noah(gid)%emiss
        xice(c,r) = noah271_struc(n)%noah(gid)%xice
        t2(c,r) = noah271_struc(n)%noah(gid)%t2
        th2(c,r) = noah271_struc(n)%noah(gid)%th2
        q2(c,r) = noah271_struc(n)%noah(gid)%q2
!        print*, 'qg: ', c,r,qg(c,r)
     enddo
  enddo  
#endif
  print*, 'EXITING noah271_setwrfesmfexport '
end subroutine noah271_setwrfesmfexport
 









