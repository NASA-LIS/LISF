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
! !MODULE: climatology_VariablesMod
! \label{climatology_VariablesMod}
! 
! !DESCRIPTION: 
!  This module contains routines that initialize and read in forcing
!   variable fields and arrays for use in LIS-7 model run.
!   
! !REVISION HISTORY: 
!  27Sep2016 -- KR Arsenault;  Initial Specification
! 
! !INTERFACE:
module climatology_VariablesMod
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
  public :: climatology_Variables_init   ! initialize forcing variables
  public :: climatology_Variables_read   ! Read NetCDF forcing variables
!  public :: climatology_Variables_reset  ! deallocate/reinit any variable arrays/inputs
!  public :: climatology_Variables_final  ! deallocate any variable arrays/inputs
!EOP

!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: climvars_struc
  public :: forcopts
!EOP

  type forcvars_type_dec

    real, allocatable :: airtmp(:,:,:)
    real, allocatable :: spechum(:,:,:)
    real, allocatable :: psurf(:,:,:)
    real, allocatable :: swdown(:,:,:)
    real, allocatable :: lwdown(:,:,:)
    real, allocatable :: uwind(:,:,:)
    real, allocatable :: vwind(:,:,:)
    real, allocatable :: rainf(:,:,:)
    real, allocatable :: cpcp(:,:,:)
    real, allocatable :: snowf(:,:,:)

  end type forcvars_type_dec
  type(forcvars_type_dec), allocatable :: climvars_struc(:)

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
! !ROUTINE: climatology_Variables_init 
!  \label{climatology_Variables_init}
! 
! !INTERFACE: 
!
 subroutine climatology_Variables_init( findex, ftn, inc, inr, &
                        ntimes, num_vars )
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
   integer, intent(in) :: ntimes       ! Input number of daily time points
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
   allocate( climvars_struc(LIS_rc%nnest) )

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
          if( varName == "Tair" ) then  ! Air temp. field present 
            write(LIS_logunit,*) &
             "[INFO] LDT-climatology variable being read in: ",trim(varName)
            num_vars = num_vars + 1
            do n=1, LIS_rc%nnest
               allocate( climvars_struc(n)%airtmp(&
                         LIS_rc%lnc(n),LIS_rc%lnr(n),ntimes) )
               climvars_struc(n)%airtmp = LIS_rc%udef
            end do
            forcopts%read_airtmp = .true.
          endif
        endif

      ! SPECHUM (Q2):
        if( LIS_FORC_Qair%selectOpt == 1 ) then  ! Spec. humidity field selected
          if( varName == "Qair" ) then   ! field present 
            write(LIS_logunit,*) &
             "[INFO] LDT-climatology variable being read in: ",trim(varName)
            num_vars = num_vars + 1
            do n=1, LIS_rc%nnest
               allocate( climvars_struc(n)%spechum(&
                         LIS_rc%lnc(n),LIS_rc%lnr(n),ntimes) )
               climvars_struc(n)%spechum = LIS_rc%udef
            end do
            forcopts%read_spechum = .true.
          endif
        endif

      ! Downward shortwave radiation (SWdown):
        if( LIS_FORC_SWdown%selectOpt == 1 ) then   ! field selected
          if( varName == "SWdown" ) then   ! field present 
            write(LIS_logunit,*) &
             "[INFO] LDT-climatology variable being read in: ",trim(varName)
            num_vars = num_vars + 1
            do n=1, LIS_rc%nnest
               allocate( climvars_struc(n)%swdown(&
                         LIS_rc%lnc(n),LIS_rc%lnr(n),ntimes) )
               climvars_struc(n)%swdown = LIS_rc%udef
            end do
            forcopts%read_swdown = .true.
          endif
        endif

      ! Downward longwave radiation (LWdown):
        if( LIS_FORC_LWdown%selectOpt == 1 ) then   ! field selected
          if( varName == "LWdown" ) then   ! field present 
            write(LIS_logunit,*) &
             "[INFO] LDT-climatology variable being read in: ",trim(varName)
            num_vars = num_vars + 1
            do n=1, LIS_rc%nnest
               allocate( climvars_struc(n)%lwdown(&
                         LIS_rc%lnc(n),LIS_rc%lnr(n),ntimes) )
               climvars_struc(n)%lwdown = LIS_rc%udef
            end do
            forcopts%read_lwdown = .true.
          endif
        endif

      ! Wind E-dir (Uwind):
        if( LIS_FORC_Wind_E%selectOpt == 1 ) then   ! field selected
          if( varName == "Uwind" ) then   ! field present 
            write(LIS_logunit,*) &
             "[INFO] LDT-climatology variable being read in: ",trim(varName)
            num_vars = num_vars + 1
            do n=1, LIS_rc%nnest
               allocate( climvars_struc(n)%uwind(&
                         LIS_rc%lnc(n),LIS_rc%lnr(n),ntimes) )
               climvars_struc(n)%uwind = LIS_rc%udef
            end do
            forcopts%read_uwind = .true.
          endif
        endif

      ! Wind - North direction (Vwind):
        if( LIS_FORC_Wind_N%selectOpt == 1 ) then   ! field selected
          if( varName == "Vwind" ) then   ! field present 
            write(LIS_logunit,*) &
             "[INFO] LDT-climatology variable being read in: ",trim(varName)
            num_vars = num_vars + 1
            do n=1, LIS_rc%nnest
               allocate( climvars_struc(n)%vwind(&
                         LIS_rc%lnc(n),LIS_rc%lnr(n),ntimes) )
               climvars_struc(n)%vwind = LIS_rc%udef
            end do
            forcopts%read_vwind = .true.
          endif
        endif

      ! Surface Pressure (Psurf):
        if( LIS_FORC_Psurf%selectOpt == 1 ) then   ! field selected
          if( varName == "Psurf" ) then   ! field present 
            write(LIS_logunit,*) &
             "[INFO] LDT-climatology variable being read in: ",trim(varName)
            num_vars = num_vars + 1
            do n=1, LIS_rc%nnest
               allocate( climvars_struc(n)%psurf(&
                         LIS_rc%lnc(n),LIS_rc%lnr(n),ntimes) )
               climvars_struc(n)%psurf = LIS_rc%udef
            end do
            forcopts%read_psurf = .true.
          endif
        endif

      ! RAINFALL:
        if( LIS_FORC_Rainf%selectOpt == 1 ) then  ! Rainfall field selected
          if( varName == "Rainf" ) then  ! Rainfall field present 
            write(LIS_logunit,*) &
             "[INFO] LDT-climatology variable being read in: ",trim(varName)
            num_vars = num_vars + 1
            do n=1, LIS_rc%nnest
               allocate( climvars_struc(n)%rainf(&
                         LIS_rc%lnc(n),LIS_rc%lnr(n),ntimes) )
               climvars_struc(n)%rainf = LIS_rc%udef
            end do
            forcopts%read_rainf = .true.
          endif
        endif

      ! Convective Rainfall (CRainf):
        if( LIS_FORC_CRainf%selectOpt == 1 ) then   ! field selected
          if( varName == "CRainf" ) then   ! field present 
            write(LIS_logunit,*) &
             "[INFO] LDT-climatology variable being read in: ",trim(varName)
            num_vars = num_vars + 1
            do n=1, LIS_rc%nnest
               allocate( climvars_struc(n)%cpcp(&
                         LIS_rc%lnc(n),LIS_rc%lnr(n),ntimes) )
               climvars_struc(n)%cpcp = LIS_rc%udef
            end do
            forcopts%read_cpcp = .true.
          endif
        endif

      endif
   end do
#endif

 end subroutine climatology_Variables_init

!BOP
! !ROUTINE: climatology_Variables_read
!  \label{climatology_Variables_read}
! 
! !INTERFACE: 
!
 subroutine climatology_Variables_read(findex, filename, inc, inr, ntimes)
!
! !USES: 
   use LIS_coreMod, only : LIS_rc, LIS_domain
   use LIS_logMod,  only : LIS_logunit, LIS_verify, LIS_endrun
   use climatology_SpatialInterpMod

   implicit none
!
! !ARGUMENTS: 
   integer, intent(in) :: findex          ! Forcing index
   character(len=*), intent(in) :: filename ! Forcing filename path
   integer, intent(in) :: inc, inr        ! Input forcing cols, rows
   integer, intent(in) :: ntimes          ! Input number of daily time pts 
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
   real       :: input_var(inc,inr,ntimes)
! _______________________________________________

  write(LIS_logunit,*)"[INFO] Reading and reprojecting Climatology forcing fields "

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

!- Open and check LDT-climatology netcdf file:
    ios = nf90_open(path=filename,&
          mode=NF90_NOWRITE,ncId=nid)
    call LIS_verify(ios,'Error in nf90_open in climatology_Variables_read')

!- Read in LDT-climatology netcdf file:

  ! Read-in and spatially reproject Air temperature:
    if( forcopts%read_airtmp ) then
      ios = nf90_inq_varid( nid, "Tair", varid )
      call LIS_verify(ios,'Error in nf90_inq_varid(Tair) in climatology_Variables_read')
      forcopts%index_airtmp = varid-2    ! Subtract off lat,lon id nums

      input_var = LIS_rc%udef
      ios = nf90_get_var( nid, varid, input_var )
      call LIS_verify(ios,'Error in nf90_get_var(Tair) in climatology_Variables_read')
      do n = 1, LIS_rc%nnest
         call climatology_interp_data( n, findex, inc, inr, ntimes,&
                 input_var, .false., climvars_struc(n)%airtmp )
      end do
    endif

  ! Read-in and spatially reproject Spec Humidity:
    if( forcopts%read_spechum ) then
      ios = nf90_inq_varid( nid, "Qair", varid )
      call LIS_verify(ios,'Error in nf90_inq_varid(Qair) in climatology_Variables_read')
      forcopts%index_spechum = varid-2    ! Subtract off lat,lon id nums

      input_var = LIS_rc%udef
      ios = nf90_get_var( nid, varid, input_var )
      call LIS_verify(ios,'Error in nf90_get_var(Qair) in climatology_Variables_read')
      do n = 1, LIS_rc%nnest
         call climatology_interp_data( n, findex, inc, inr, ntimes,&
                 input_var, .false., climvars_struc(n)%spechum )
      end do
    endif

  ! Read-in and spatially reproject SW Down:
    if( forcopts%read_swdown ) then
      ios = nf90_inq_varid( nid, "SWdown", varid )
      call LIS_verify(ios,'Error in nf90_inq_varid(SWdown) in climatology_Variables_read')
      forcopts%index_swdown = varid-2    ! Subtract off lat,lon id nums

      input_var = LIS_rc%udef
      ios = nf90_get_var( nid, varid, input_var )
      call LIS_verify(ios,'Error in nf90_get_var(SWdown) in climatology_Variables_read')
      do n = 1, LIS_rc%nnest
         call climatology_interp_data( n, findex, inc, inr, ntimes,&
                 input_var, .false., climvars_struc(n)%swdown )
      end do
    endif

  ! Read-in and spatially reproject LWdown:
    if( forcopts%read_lwdown ) then
      ios = nf90_inq_varid( nid, "LWdown", varid )
      call LIS_verify(ios,'Error in nf90_inq_varid(LWdown) in climatology_Variables_read')
      forcopts%index_lwdown = varid-2    ! Subtract off lat,lon id nums

      input_var = LIS_rc%udef
      ios = nf90_get_var( nid, varid, input_var )
      call LIS_verify(ios,'Error in nf90_get_var(LWdown) in climatology_Variables_read')
      do n = 1, LIS_rc%nnest
         call climatology_interp_data( n, findex, inc, inr, ntimes,&
                 input_var, .false., climvars_struc(n)%lwdown )
      end do
    endif

  ! Read-in and spatially reproject U-wind:
    if( forcopts%read_uwind ) then
      ios = nf90_inq_varid( nid, "Uwind", varid )
      call LIS_verify(ios,'Error in nf90_inq_varid(Uwind) in climatology_Variables_read')
      forcopts%index_uwind = varid-2    ! Subtract off lat,lon id nums

      input_var = LIS_rc%udef
      ios = nf90_get_var( nid, varid, input_var )
      call LIS_verify(ios,'Error in nf90_get_var(Uwind) in climatology_Variables_read')
      do n = 1, LIS_rc%nnest
         call climatology_interp_data( n, findex, inc, inr, ntimes,&
                 input_var, .false., climvars_struc(n)%uwind )
      end do
    endif

  ! Read-in and spatially reproject V-wind:
    if( forcopts%read_vwind ) then
      ios = nf90_inq_varid( nid, "Vwind", varid )
      call LIS_verify(ios,'Error in nf90_inq_varid(Vwind) in climatology_Variables_read')
      forcopts%index_vwind = varid-2    ! Subtract off lat,lon id nums

      input_var = LIS_rc%udef
      ios = nf90_get_var( nid, varid, input_var )
      call LIS_verify(ios,'Error in nf90_get_var(Vwind) in climatology_Variables_read')
      do n = 1, LIS_rc%nnest
         call climatology_interp_data( n, findex, inc, inr, ntimes,&
                 input_var, .false., climvars_struc(n)%vwind )
      end do
    endif

  ! Read-in and spatially reproject Psurf:
    if( forcopts%read_psurf ) then
      ios = nf90_inq_varid( nid, "Psurf", varid )
      call LIS_verify(ios,'Error in nf90_inq_varid(Psurf) in climatology_Variables_read')
      forcopts%index_psurf = varid-2    ! Subtract off lat,lon id nums

      input_var = LIS_rc%udef
      ios = nf90_get_var( nid, varid, input_var )
      call LIS_verify(ios,'Error in nf90_get_var(Psurf) in climatology_Variables_read')
      do n = 1, LIS_rc%nnest
         call climatology_interp_data( n, findex, inc, inr, ntimes,&
                 input_var, .false., climvars_struc(n)%psurf )
      end do
    endif

  ! Read-in and spatially reproject Rainfall:
    if( forcopts%read_rainf ) then
      ios = nf90_inq_varid( nid, "Rainf", varid )
      call LIS_verify(ios,'Error in nf90_inq_varid(Rainf) in climatology_Variables_read')
      forcopts%index_rainf = varid-2    ! Subtract off lat,lon id nums

      input_var = LIS_rc%udef
      ios = nf90_get_var( nid, varid, input_var )
      call LIS_verify(ios,'Error in nf90_get_var(Rainf) in climatology_Variables_read')
      do n = 1, LIS_rc%nnest
         call climatology_interp_data( n, findex, inc, inr, ntimes,&
                 input_var, .true., climvars_struc(n)%rainf )
      end do   
    endif

  ! Read-in and spatially reproject Conv. Rainfall:
    if( forcopts%read_cpcp ) then
      ios = nf90_inq_varid( nid, "CRainf", varid )
      call LIS_verify(ios,'Error in nf90_inq_varid(CRainf) in climatology_Variables_read')
      forcopts%index_cpcp = varid-2    ! Subtract off lat,lon id nums

      input_var = LIS_rc%udef
      ios = nf90_get_var( nid, varid, input_var )
      call LIS_verify(ios,'Error in nf90_get_var(CRainf) in climatology_Variables_read')
      do n = 1, LIS_rc%nnest
         call climatology_interp_data( n, findex, inc, inr, ntimes,&
                 input_var, .true., climvars_struc(n)%cpcp )
      end do
    endif

!- Close netCDF file.
   ios=nf90_close(nid)
   call LIS_verify(ios,'Error in nf90_close in climatology_Variables_read')
!-
#endif

!   write(LIS_logunit,*) " Done reading and spatially reprojecting forcing fields - "

 end subroutine climatology_Variables_read

!
!  subroutine climatology_Variables_reset( n, findex )
!
!  end subroutine climatology_Variables_reset

!
!  subroutine climatology_Variables_final( n, findex )
!
!  end subroutine climatology_Variables_final

end module climatology_VariablesMod
