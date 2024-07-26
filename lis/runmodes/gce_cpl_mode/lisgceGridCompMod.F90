!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module lisgceGridCompMod
  
  use ESMF
  use lisgceimport_module
  use lisgceexport_module

  implicit none

  PRIVATE
  
  public :: lisgce_alloc_states
  public :: lisgce_import, lisgce_export

  type(lisgceimport) :: lisgce_import
  type(lisgceexport) :: lisgce_export
  contains

    subroutine lisgce_alloc_states

      use LIS_coreMod, only : LIS_rc

      implicit none
      integer :: n
      
      n = 1
         ! Import states      
      allocate(lisgce_import%data_tair(LIS_rc%lnc(n),LIS_rc%lnr(n)))
      allocate(lisgce_import%data_qair(LIS_rc%lnc(n),LIS_rc%lnr(n)))
      allocate(lisgce_import%data_swd(LIS_rc%lnc(n),LIS_rc%lnr(n)))
      allocate(lisgce_import%data_lwd(LIS_rc%lnc(n),LIS_rc%lnr(n)))
      allocate(lisgce_import%data_psurf(LIS_rc%lnc(n),LIS_rc%lnr(n)))
      allocate(lisgce_import%data_prcp(LIS_rc%lnc(n),LIS_rc%lnr(n)))
      allocate(lisgce_import%data_u(LIS_rc%lnc(n),LIS_rc%lnr(n)))
      allocate(lisgce_import%data_v(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         !Export states
      allocate(lisgce_export%swnet(LIS_rc%lnc(n),LIS_rc%lnr(n)))
      allocate(lisgce_export%lwnet(LIS_rc%lnc(n),LIS_rc%lnr(n)))
      allocate(lisgce_export%qle(LIS_rc%lnc(n),LIS_rc%lnr(n)))
      allocate(lisgce_export%qh(LIS_rc%lnc(n),LIS_rc%lnr(n)))
      allocate(lisgce_export%qg(LIS_rc%lnc(n),LIS_rc%lnr(n)))
      allocate(lisgce_export%tau(LIS_rc%lnc(n),LIS_rc%lnr(n)))
      allocate(lisgce_export%tauu(LIS_rc%lnc(n),LIS_rc%lnr(n)))
      allocate(lisgce_export%tauv(LIS_rc%lnc(n),LIS_rc%lnr(n)))
      allocate(lisgce_export%albedo(LIS_rc%lnc(n),LIS_rc%lnr(n)))
      allocate(lisgce_export%avgsurft(LIS_rc%lnc(n),LIS_rc%lnr(n)))

      
#if 0 
  n = 1
  call ESMF_ArraySpecSet(arrayspec,rank=1,typekind=ESMF_TYPEKIND_R4)
  field_tair = ESMF_FieldCreate(grid = lisGrid(n),arrayspec=arrayspec,&
       name="Tair",rc=status)
  call verify(status)

  field_qair = ESMF_FieldCreate(grid = lisGrid(n),arrayspec=arrayspec,&
       name="Qair",rc=status)
  call verify(status)

  field_swdown = ESMF_FieldCreate(grid = lisGrid(n),arrayspec=arrayspec,&
       name="SWdown",rc=status)
  call verify(status)

  field_lwdown = ESMF_FieldCreate(grid = lisGrid(n),arrayspec=arrayspec,&
       name="LWdown",rc=status)
  call verify(status)

  field_uwind = ESMF_FieldCreate(grid = lisGrid(n),arrayspec=arrayspec,&
       name="Uwind",rc=status)
  call verify(status)

  field_vwind = ESMF_FieldCreate(grid = lisGrid(n),arrayspec=arrayspec,&
       name="Vwind",rc=status)
  call verify(status)

  field_psurf = ESMF_FieldCreate(grid = lisGrid(n),arrayspec=arrayspec,&
       name="Psurf",rc=status)
  call verify(status)

  field_rainf = ESMF_FieldCreate(grid = lisGrid(n),arrayspec=arrayspec,&
       name="Rainf",rc=status)
  call verify(status)

  call ESMF_StateAddField(LISimp, field_tair,rc)
  call ESMF_StateAddField(LISimp, field_qair,rc)
  call ESMF_StateAddField(LISimp, field_swdown,rc)
  call ESMF_StateAddField(LISimp, field_lwdown,rc)
  call ESMF_StateAddField(LISimp, field_uwind,rc)
  call ESMF_StateAddField(LISimp, field_vwind,rc)
  call ESMF_StateAddField(LISimp, field_psurf,rc)
  call ESMF_StateAddField(LISimp, field_rainf,rc)
#endif

    end subroutine lisgce_alloc_states
  end module lisgceGridCompMod
