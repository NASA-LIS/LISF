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
! !MODULE: pptEnsFcst_VariablesMod
! \label{pptEnsFcst_VariablesMod}
! 
! !DESCRIPTION: 
!  This module contains routines that initialize and read in forcing
!   variable fields and arrays for use in LIS-7 model run.
!   
! !REVISION HISTORY: 
!  27Sep2016 -- KR Arsenault;  Initial Specification
! 
! !INTERFACE:
module pptEnsFcst_VariablesMod
!
! !USES:
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LIS_FORC_AttributesMod

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: pptEnsFcst_Variables_init   ! initialize forcing variables
  public :: pptEnsFcst_Variables_read   ! Read NetCDF forcing variables
!  public :: pptEnsFcst_Variables_reset  ! deallocate/reinit any variable arrays/inputs
!  public :: pptEnsFcst_Variables_final  ! deallocate any variable arrays/inputs
!EOP

!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: ensfcstppt_struc
  public :: forcopts
!EOP

  type forcvars_type_dec

    real, allocatable :: airtmp(:,:)!,:)
    real, allocatable :: spechum(:,:)!,:)
    real, allocatable :: psurf(:,:)!,:)
    real, allocatable :: swdown(:,:)!,:)
    real, allocatable :: lwdown(:,:)!,:)
    real, allocatable :: uwind(:,:)!,:)
    real, allocatable :: vwind(:,:)!,:)
    real, allocatable :: rainf(:,:)!,:)
    real, allocatable :: cpcp(:,:)!,:)
    real, allocatable :: snowf(:,:)!,:)

  end type forcvars_type_dec
  type(forcvars_type_dec), allocatable :: ensfcstppt_struc(:)

  type forcopts_type_dec

  ! Forcing field Present Flag:
    logical       :: read_airtmp
    logical       :: read_spechum
    logical       :: read_psurf
    logical       :: read_swdown
    logical       :: read_lwdown
    logical       :: read_uwind
    logical       :: read_vwind
    logical       :: read_rainf
    logical       :: read_cpcp
    logical       :: read_snowf

  ! Forcing field Present Flag:
    integer       :: index_airtmp
    integer       :: index_spechum
    integer       :: index_psurf
    integer       :: index_swdown
    integer       :: index_lwdown
    integer       :: index_uwind
    integer       :: index_vwind
    integer       :: index_rainf
    integer       :: index_cpcp
    integer       :: index_snowf

  end type forcopts_type_dec
  type(forcopts_type_dec) :: forcopts

contains

!BOP
! !ROUTINE: pptEnsFcst_Variables_init 
!  \label{pptEnsFcst_Variables_init}
! 
! !INTERFACE: 
!
 subroutine pptEnsFcst_Variables_init( findex, ftn, &
                  inc, inr, tindex, num_vars )
!
! !USES: 
   use LIS_coreMod, only : LIS_rc, LIS_domain
   use LIS_logMod,  only : LIS_logunit, LIS_endrun
   use LIS_FORC_AttributesMod

   implicit none
!
! !ARGUMENTS: 
   integer, intent(in) :: findex       ! Forcing index
   integer, intent(in) :: ftn          ! Forcing file unit (netcdf-format)
   integer, intent(in) :: inc, inr     ! Input forcing cols, rows
   integer, intent(in) :: tindex       ! Index number of daily time point
   integer, intent(out):: num_vars     ! Number of forcing variables counted
!
! !DESCRIPTION: 
!  This routine allocates and initializes the forcing
!   variables needed for driving a LSM.
!EOP  
   integer        :: n, i
   integer        :: ios
   character(40)  :: varName
! _______________________________________________

   write(LIS_logunit,*) "[INFO] -- Allocating and Initializing Forcing Variables -- "
   allocate( ensfcstppt_struc(LIS_rc%nnest) )

!- Allocate and initialize all read-in forcing variable fields:

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
   num_vars = 0
   do i = 1, 50   ! Loop over potential variables
    ! Check variable id/Name:
      ios = nf90_inquire_variable(ftn, &
                 varid = i, &
                 name = varName )
      if( ios == 0 ) then
     !- Estimate number of forcing variables being readin from netcdf file:
        if( varName == "lat" .or. &
            varName == "lon" .or. &
            varName == "time" ) then
           cycle
        endif

      ! AIRTEMP:
        if( LIS_FORC_Tair%selectOpt == 1 ) then  ! Air temp. field selected
          if( varName == "T2M" ) then  ! Air temp. field present 
            write(LIS_logunit,*) &
             "[INFO] PPTEnsFcst variable being read in: ",trim(varName)
            num_vars = num_vars + 1
            do n=1, LIS_rc%nnest
               allocate( ensfcstppt_struc(n)%airtmp(&
                         LIS_rc%lnc(n),LIS_rc%lnr(n) ) )
               ensfcstppt_struc(n)%airtmp = LIS_rc%udef
            end do
            forcopts%read_airtmp = .true.
            forcopts%index_airtmp = num_vars ! EMK
          endif
        endif

      ! SPECHUM (Q2):
        if( LIS_FORC_Qair%selectOpt == 1 ) then  ! Spec. humidity field selected
          if( varName == "QV2M" ) then   ! field present 
            write(LIS_logunit,*) &
             "[INFO] PPTEnsFcst variable being read in: ",trim(varName)
            num_vars = num_vars + 1
            do n=1, LIS_rc%nnest
               allocate( ensfcstppt_struc(n)%spechum(&
                         LIS_rc%lnc(n),LIS_rc%lnr(n) ) )
               ensfcstppt_struc(n)%spechum = LIS_rc%udef
            end do
            forcopts%read_spechum = .true.
            forcopts%index_spechum = num_vars ! EMK
          endif
        endif

      ! Downward shortwave radiation (SWGDN):
        if( LIS_FORC_SWdown%selectOpt == 1 ) then   ! field selected
          if( varName == "SWGDN" ) then   ! field present 
            write(LIS_logunit,*) &
             "[INFO] PPTEnsFcst variable being read in: ",trim(varName)
            num_vars = num_vars + 1
            do n=1, LIS_rc%nnest
               allocate( ensfcstppt_struc(n)%swdown(&
                         LIS_rc%lnc(n),LIS_rc%lnr(n) ) )
               ensfcstppt_struc(n)%swdown = LIS_rc%udef
            end do
            forcopts%read_swdown = .true.
            forcopts%index_swdown = num_vars ! EMK
          endif
        endif

      ! Downward longwave radiation (LWGAB):
        if( LIS_FORC_LWdown%selectOpt == 1 ) then   ! field selected
          if( varName == "LWGAB" ) then   ! field present 
            write(LIS_logunit,*) &
             "[INFO] PPTEnsFcst variable being read in: ",trim(varName)
            num_vars = num_vars + 1
            do n=1, LIS_rc%nnest
               allocate( ensfcstppt_struc(n)%lwdown(&
                         LIS_rc%lnc(n),LIS_rc%lnr(n) ) )
               ensfcstppt_struc(n)%lwdown = LIS_rc%udef
            end do
            forcopts%read_lwdown = .true.
            forcopts%index_lwdown = num_vars ! EMK
          endif
        endif

      ! Wind E-dir (U10M):
        if( LIS_FORC_Wind_E%selectOpt == 1 ) then   ! field selected
          if( varName == "U10M" ) then   ! field present 
            write(LIS_logunit,*) &
             "[INFO] PPTEnsFcst variable being read in: ",trim(varName)
            num_vars = num_vars + 1
            do n=1, LIS_rc%nnest
               allocate( ensfcstppt_struc(n)%uwind(&
                         LIS_rc%lnc(n),LIS_rc%lnr(n) ) )
               ensfcstppt_struc(n)%uwind = LIS_rc%udef
            end do
            forcopts%read_uwind = .true.
            forcopts%index_uwind = num_vars ! EMK
          endif
        endif

      ! Wind - North direction (V10M):
        if( LIS_FORC_Wind_N%selectOpt == 1 ) then   ! field selected
          if( varName == "V10M" ) then   ! field present 
            write(LIS_logunit,*) &
             "[INFO] PPTEnsFcst variable being read in: ",trim(varName)
            num_vars = num_vars + 1
            do n=1, LIS_rc%nnest
               allocate( ensfcstppt_struc(n)%vwind(&
                         LIS_rc%lnc(n),LIS_rc%lnr(n) ) )
               ensfcstppt_struc(n)%vwind = LIS_rc%udef
            end do
            forcopts%read_vwind = .true.
            forcopts%index_vwind = num_vars ! EMK
          endif
        endif

      ! Surface Pressure (PS):
        if( LIS_FORC_Psurf%selectOpt == 1 ) then   ! field selected
          if( varName == "PS" ) then   ! field present 
            write(LIS_logunit,*) &
             "[INFO] PPTEnsFcst variable being read in: ",trim(varName)
            num_vars = num_vars + 1
            do n=1, LIS_rc%nnest
               allocate( ensfcstppt_struc(n)%psurf(&
                         LIS_rc%lnc(n),LIS_rc%lnr(n) ) )
               ensfcstppt_struc(n)%psurf = LIS_rc%udef
            end do
            forcopts%read_psurf = .true.
            forcopts%index_psurf = num_vars ! EMK
          endif
        endif

      ! RAINFALL:
        if( LIS_FORC_Rainf%selectOpt == 1 ) then  ! Rainfall field selected
          if( varName == "PRECTOT" ) then  ! Rainfall field present 
            write(LIS_logunit,*) &
             "[INFO] PPTEnsFcst variable being read in: ",trim(varName)
            num_vars = num_vars + 1
            do n=1, LIS_rc%nnest
               allocate( ensfcstppt_struc(n)%rainf(&
                         LIS_rc%lnc(n),LIS_rc%lnr(n) ) )
               ensfcstppt_struc(n)%rainf = LIS_rc%udef
            end do
            forcopts%read_rainf = .true.
            forcopts%index_rainf = num_vars ! EMK
          endif
        endif

      ! Convective Rainfall (PRECCON):
        if( LIS_FORC_CRainf%selectOpt == 1 ) then   ! field selected
          if( varName == "PRECCON" ) then   ! field present 
            write(LIS_logunit,*) &
             "[INFO] PPTEnsFcst variable being read in: ",trim(varName)
            num_vars = num_vars + 1
            do n=1, LIS_rc%nnest
               allocate( ensfcstppt_struc(n)%cpcp(&
                         LIS_rc%lnc(n),LIS_rc%lnr(n) ) )
               ensfcstppt_struc(n)%cpcp = LIS_rc%udef
            end do
            forcopts%read_cpcp = .true.
            forcopts%index_cpcp = num_vars ! EMK
          endif
        endif

      endif
   end do
#endif

 end subroutine pptEnsFcst_Variables_init

!BOP
! !ROUTINE: pptEnsFcst_Variables_read
!  \label{pptEnsFcst_Variables_read}
! 
! !INTERFACE: 
!
 subroutine pptEnsFcst_Variables_read(findex, filename, inc, inr, tindex )
!                       start_inc, start_inr, tindex)
!
! !USES: 
   use LIS_coreMod, only : LIS_rc, LIS_domain
   use LIS_logMod,  only : LIS_logunit, LIS_verify, LIS_endrun
   use pptEnsFcst_SpatialInterpMod

   implicit none
!
! !ARGUMENTS: 
   integer, intent(in) :: findex          ! Forcing index
   character(len=*), intent(in) :: filename ! Forcing filename path
   integer, intent(in) :: inc, inr        ! Input forcing cols, rows
!   integer, intent(in) :: start_inc, start_inr  ! Initial col / row points of subsetted domain 
   integer, intent(in) :: tindex          ! Index of daily time pt
!
! !DESCRIPTION: 
!  This routine reads in the forcing
!   variables needed for driving a LSM.
!EOP  
!
   integer    :: n       ! Nest index
   integer    :: ios     ! Input/output status
   integer    :: nid     ! Netcdf file unit ID 
   integer    :: varid   ! Netcdf file id
   integer    :: c,r
   real       :: input_var(inc,inr)
! _______________________________________________

  write(LIS_logunit,*)"[INFO] Reading, reprojecting ... Ensemble forecast fields"

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

!- Open and check PPTEnsFcst netcdf file:
    ios = nf90_open(path=filename,&
          mode=NF90_NOWRITE,ncId=nid)
    call LIS_verify(ios,'Error in nf90_open in pptEnsFcst_Variables_read')

!- Read in PPTEnsFcst netcdf file:

  ! Read-in and spatially reproject Air temperature:
    if( forcopts%read_airtmp ) then
      ios = nf90_inq_varid( nid, "T2M", varid )
      call LIS_verify(ios,'Error in nf90_inq_varid(T2M) in pptEnsFcst_Variables_read')

      input_var = LIS_rc%udef
      ios = nf90_get_var( nid, varid, input_var, &
! future subset option:
!              start=(/start_inc, start_inr, tindex/),&
              start=(/1, 1, tindex/),&
              count=(/inc, inr, 1/) )
      call LIS_verify(ios,'Error in nf90_get_var(T2M) in pptEnsFcst_Variables_read')
      do n = 1, LIS_rc%nnest
         call pptEnsFcst_interp_data( n, findex, inc, inr, tindex,&
                 input_var, .false., ensfcstppt_struc(n)%airtmp )
      end do
    endif

  ! Read-in and spatially reproject Spec Humidity:
    if( forcopts%read_spechum ) then
      ios = nf90_inq_varid( nid, "QV2M", varid )
      call LIS_verify(ios,'Error in nf90_inq_varid(QV2M) in pptEnsFcst_Variables_read')

      input_var = LIS_rc%udef
      ios = nf90_get_var( nid, varid, input_var, & 
!              start=(/start_inc, start_inr, tindex/),&
              start=(/1, 1, tindex/),&
              count=(/inc, inr, 1/) )
      call LIS_verify(ios,'Error in nf90_get_var(QV2M) in pptEnsFcst_Variables_read')
      do n = 1, LIS_rc%nnest
         call pptEnsFcst_interp_data( n, findex, inc, inr, tindex,&
                 input_var, .false., ensfcstppt_struc(n)%spechum )
      end do
    endif

  ! Read-in and spatially reproject SW Down:
    if( forcopts%read_swdown ) then
      ios = nf90_inq_varid( nid, "SWGDN", varid )
      call LIS_verify(ios,'Error in nf90_inq_varid(SWGDN) in pptEnsFcst_Variables_read')

      input_var = LIS_rc%udef
      ios = nf90_get_var( nid, varid, input_var, & 
!              start=(/start_inc, start_inr, tindex/),&
              start=(/1, 1, tindex/),&
              count=(/inc, inr, 1/) )
      call LIS_verify(ios,'Error in nf90_get_var(SWGDN) in pptEnsFcst_Variables_read')
      do n = 1, LIS_rc%nnest
         call pptEnsFcst_interp_data( n, findex, inc, inr, tindex,&
                 input_var, .false., ensfcstppt_struc(n)%swdown )
      end do
    endif

  ! Read-in and spatially reproject LWdown:
    if( forcopts%read_lwdown ) then
      ios = nf90_inq_varid( nid, "LWGAB", varid )
      call LIS_verify(ios,'Error in nf90_inq_varid(LWGAB) in pptEnsFcst_Variables_read')

      input_var = LIS_rc%udef
      ios = nf90_get_var( nid, varid, input_var, & 
!              start=(/start_inc, start_inr, tindex/),&
              start=(/1, 1, tindex/),&
              count=(/inc, inr, 1/) )
      call LIS_verify(ios,'Error in nf90_get_var(LWGAB) in pptEnsFcst_Variables_read')
      do n = 1, LIS_rc%nnest
         call pptEnsFcst_interp_data( n, findex, inc, inr, tindex,&
                 input_var, .false., ensfcstppt_struc(n)%lwdown )
      end do
    endif

  ! Read-in and spatially reproject U-wind:
    if( forcopts%read_uwind ) then
      ios = nf90_inq_varid( nid, "U10M", varid )
      call LIS_verify(ios,'Error in nf90_inq_varid(U10M) in pptEnsFcst_Variables_read')

      input_var = LIS_rc%udef
      ios = nf90_get_var( nid, varid, input_var, & 
!              start=(/start_inc, start_inr, tindex/),&
              start=(/1, 1, tindex/),&
              count=(/inc, inr, 1/) )
      call LIS_verify(ios,'Error in nf90_get_var(U10M) in pptEnsFcst_Variables_read')
      do n = 1, LIS_rc%nnest
         call pptEnsFcst_interp_data( n, findex, inc, inr, tindex,&
                 input_var, .false., ensfcstppt_struc(n)%uwind )
      end do
    endif

  ! Read-in and spatially reproject V-wind:
    if( forcopts%read_vwind ) then
      ios = nf90_inq_varid( nid, "V10M", varid )
      call LIS_verify(ios,'Error in nf90_inq_varid(V10M) in pptEnsFcst_Variables_read')

      input_var = LIS_rc%udef
      ios = nf90_get_var( nid, varid, input_var, & 
!              start=(/start_inc, start_inr, tindex/),&
              start=(/1, 1, tindex/),&
              count=(/inc, inr, 1/) )
      call LIS_verify(ios,'Error in nf90_get_var(V10M) in pptEnsFcst_Variables_read')
      do n = 1, LIS_rc%nnest
         call pptEnsFcst_interp_data( n, findex, inc, inr, tindex,&
                 input_var, .false., ensfcstppt_struc(n)%vwind )
      end do
    endif

  ! Read-in and spatially reproject PS:
    if( forcopts%read_psurf ) then
      ios = nf90_inq_varid( nid, "PS", varid )
      call LIS_verify(ios,'Error in nf90_inq_varid(PS) in pptEnsFcst_Variables_read')

      input_var = LIS_rc%udef
      ios = nf90_get_var( nid, varid, input_var, & 
!              start=(/start_inc, start_inr, tindex/),&
              start=(/1, 1, tindex/),&
              count=(/inc, inr, 1/) )
      call LIS_verify(ios,'Error in nf90_get_var(PS) in pptEnsFcst_Variables_read')
      do n = 1, LIS_rc%nnest
         call pptEnsFcst_interp_data( n, findex, inc, inr, tindex,&
                 input_var, .false., ensfcstppt_struc(n)%psurf )
      end do
    endif

  ! Read-in and spatially reproject Rainfall:
    if( forcopts%read_rainf ) then
      ios = nf90_inq_varid( nid, "PRECTOT", varid )
      call LIS_verify(ios,'Error in nf90_inq_varid(PRECTOT) in pptEnsFcst_Variables_read')

      input_var = LIS_rc%udef
      ios = nf90_get_var( nid, varid, input_var, & 
!              start=(/start_inc, start_inr, tindex/),&
              start=(/1, 1, tindex/),&
              count=(/inc, inr, 1/) )
      call LIS_verify(ios,'Error in nf90_get_var(PRECTOT) in pptEnsFcst_Variables_read')
      do n = 1, LIS_rc%nnest
         call pptEnsFcst_interp_data( n, findex, inc, inr, tindex,&
                 input_var, .true., ensfcstppt_struc(n)%rainf )
      end do   
    endif

  ! Read-in and spatially reproject Conv. Rainfall:
    if( forcopts%read_cpcp ) then
      ios = nf90_inq_varid( nid, "PRECCON", varid )
      call LIS_verify(ios,'Error in nf90_inq_varid(PRECCON) in pptEnsFcst_Variables_read')

      input_var = LIS_rc%udef
      ios = nf90_get_var( nid, varid, input_var, & 
!              start=(/start_inc, start_inr, tindex/),&
              start=(/1, 1, tindex/),&
              count=(/inc, inr, 1/) )
      call LIS_verify(ios,'Error in nf90_get_var(PRECCON) in pptEnsFcst_Variables_read')
      do n = 1, LIS_rc%nnest
         call pptEnsFcst_interp_data( n, findex, inc, inr, tindex,&
                 input_var, .true., ensfcstppt_struc(n)%cpcp )
      end do
    endif

!- Close netCDF file.
   ios=nf90_close(nid)
   call LIS_verify(ios,'Error in nf90_close in pptEnsFcst_Variables_read')
!-
#endif

 end subroutine pptEnsFcst_Variables_read

!
!  subroutine pptEnsFcst_Variables_reset( n, findex )
!
!  end subroutine pptEnsFcst_Variables_reset

!
!  subroutine pptEnsFcst_Variables_final( n, findex )
!
!  end subroutine pptEnsFcst_Variables_final

end module pptEnsFcst_VariablesMod
