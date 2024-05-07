!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module narr_forcingMod
!BOP
! !MODULE: narr_forcingMod
! 
! !DESCRIPTION: 
!  
! !USES: 
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_NARR      !defines the native resolution of 
                                  !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: narr_struc
!EOP

  type, public :: narr_type_dec
     
     real             :: ts
     real*8           :: narrtime1,narrtime2
     character(len=LDT_CONST_PATH_LEN)    :: narrdir
     integer          :: nc,nr
     integer          :: mi
     integer          :: nlevels
     integer, allocatable :: n111(:)
     integer, allocatable :: n121(:)
     integer, allocatable :: n211(:)
     integer, allocatable :: n221(:)
     real,    allocatable :: w111(:)
     real,    allocatable :: w121(:)
     real,    allocatable :: w211(:)
     real,    allocatable :: w221(:)

     integer          :: findtime1, findtime2
  end type narr_type_dec
  
  type(narr_type_dec), allocatable :: narr_struc(:)
contains
  
!BOP
!
! !ROUTINE: init_narr
!  \label{init_narr}
!
! !REVISION HISTORY: 
! 30 APR 2009: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine init_narr(findex)
! !USES: 
    use LDT_coreMod,    only : LDT_rc
    use LDT_timeMgrMod, only : LDT_date2time, LDT_update_timestep
    use LDT_logMod,     only : LDT_logunit, LDT_endrun

    implicit none
    
    integer,  intent(in) :: findex
! 
! !DESCRIPTION: 
!
!EOP
    integer   :: n
    real      :: gridDesci(20)

    allocate(narr_struc(LDT_rc%nnest))

   write(LDT_logunit,fmt=*)"MSG: Initializing NARR forcing grid ... "

!    call readcrd_narr()

! Remove this warning once the met_nf variable is set.     
!    LDT_rc%met_nf(findex) = ??
    print*, 'set the number of forcing variables in init_narr'
    print*, 'stopping..'
    stop

    LDT_rc%met_ts(findex) = 3*3600
    LDT_rc%met_zterp(findex) = .false.

    narr_struc%nc = 349
    narr_struc%nr = 277
    LDT_rc%met_nc(findex) = narr_struc(1)%nc
    LDT_rc%met_nr(findex) = narr_struc(1)%nr

 !- NARR Grid description:
    LDT_rc%met_proj(findex)        = "lambert"
    LDT_rc%met_gridDesc(findex,1)  = 3
    LDT_rc%met_gridDesc(findex,2)  = narr_struc(1)%nc     ! 349
    LDT_rc%met_gridDesc(findex,3)  = narr_struc(1)%nr     ! 277
    LDT_rc%met_gridDesc(findex,4)  = 1.0      ! latitude of origin
    LDT_rc%met_gridDesc(findex,5)  = -145.5   ! longitude of origin
    LDT_rc%met_gridDesc(findex,6)  = 128
    LDT_rc%met_gridDesc(findex,7)  = 50.0     ! truelat1
    LDT_rc%met_gridDesc(findex,8)  = 32.463   ! grid spacing in meters
    LDT_rc%met_gridDesc(findex,9)  = 0.469
    LDT_rc%met_gridDesc(findex,10) = 50.0     ! truelat2
    LDT_rc%met_gridDesc(findex,11) = -107.00  ! standard longitude
    LDT_rc%met_gridDesc(findex,20) = 0.0       

    gridDesci(:) = LDT_rc%met_gridDesc(findex,:)

    narr_struc%mi = narr_struc%nc*narr_struc%nr

 !- If only processing parameters, then return to main routine calls ...
    if( LDT_rc%runmode == "LSM parameter processing" ) return


#if 0
    do n=1,LDT_rc%nnest
       narr_struc(n)%ts = 3*3600 
       call LDT_update_timestep(LDT_rc, n, narr_struc(n)%ts)

       if(trim(LDT_rc%met_gridtransform(findex)).eq."bilinear") then 

          allocate(narr_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(narr_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(narr_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(narr_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          allocate(narr_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(narr_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(narr_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(narr_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          
          call bilinear_interp_input(n, LDT_rc%met_gridDesc(findex,:),&
               narr_struc(n)%n111,&
               narr_struc(n)%n121,narr_struc(n)%n211,&
               narr_struc(n)%n221,narr_struc(n)%w111,&
               narr_struc(n)%w121,narr_struc(n)%w211,&
               narr_struc(n)%w221)

       endif
    enddo
#endif

  end subroutine init_narr

end module narr_forcingMod

