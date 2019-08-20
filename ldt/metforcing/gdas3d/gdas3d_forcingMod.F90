!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module gdas3d_forcingMod
!BOP
! !MODULE: gdas3d_forcingMod
! 
! !DESCRIPTION: 
!  
!
! !USES: 
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_GDAS3D      !defines the native resolution of 
                             !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: gdas3d_struc
!EOP

  type, public :: gdas3d_type_dec
     real             :: ts
     real*8           :: gdastime1,gdastime2
     character*100    :: gdasdir
     integer          :: nc,nr
     integer          :: nlayer
     integer, allocatable :: n111(:)
     integer, allocatable :: n121(:)
     integer, allocatable :: n211(:)
     integer, allocatable :: n221(:)
     real,    allocatable :: w111(:)
     real,    allocatable :: w121(:)
     real,    allocatable :: w211(:)
     real,    allocatable :: w221(:)
     integer          :: findtime1, findtime2
  end type gdas3d_type_dec
  
  type(gdas3d_type_dec), allocatable :: gdas3d_struc(:)
contains
  
!BOP
!
! !ROUTINE: init_GDAS3d
!  \label{init_GDAS3d}
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
  subroutine init_GDAS3d(findex)

! !USES: 
   use LDT_coreMod,    only : LDT_rc
   use LDT_timeMgrMod, only : LDT_date2time, LDT_update_timestep
   use LDT_logMod,     only : LDT_logunit,LDT_endrun

    implicit none
    
    integer,  intent(in) :: findex
! 
! !DESCRIPTION: 
!
!EOP
    integer  :: n

    allocate(gdas3d_struc(LDT_rc%nnest))

    write(LDT_logunit,*) "Reading the GDAS-3D forcing data"

!    call readcrd_gdas3d()

    gdas3d_struc%nc = 768
    gdas3d_struc%nr = 384
    LDT_rc%met_nc(findex) = gdas3d_struc(1)%nc
    LDT_rc%met_nr(findex) = gdas3d_struc(1)%nr

 !- GDAS 3D Grid description:
    LDT_rc%met_proj(findex)        = "gaussian"
    LDT_rc%met_gridDesc(findex,1)  = 4
    LDT_rc%met_gridDesc(findex,2)  = LDT_rc%met_nc(findex)
    LDT_rc%met_gridDesc(findex,3)  = LDT_rc%met_nr(findex)
    LDT_rc%met_gridDesc(findex,4)  = 89.642
    LDT_rc%met_gridDesc(findex,5)  = 0
    LDT_rc%met_gridDesc(findex,6)  = 128
    LDT_rc%met_gridDesc(findex,7)  = 89.642
    LDT_rc%met_gridDesc(findex,8)  = -0.469
    LDT_rc%met_gridDesc(findex,9)  = 0.469
    LDT_rc%met_gridDesc(findex,10) = 192
    LDT_rc%met_gridDesc(findex,20) = 0.0

!    LDT_rc%met_nf(findex) = ??
!remove this warning once the met_nf variable is set.     
    print*, 'set the number of forcing variables in init_gdas3d'
    print*, 'stopping..'
    stop

    LDT_rc%met_ts(findex) = 6*3600

    do n=1, LDT_rc%nnest
       gdas3d_struc(n)%ts = 6*3600 
       call LDT_update_timestep(LDT_rc, n, gdas3d_struc(n)%ts)

       gdas3d_struc(n)%findtime1 = 0 
       gdas3d_struc(n)%findtime2 = 0 

       if(trim(LDT_rc%met_gridtransform(findex)).eq."bilinear") then 

          allocate(gdas3d_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gdas3d_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gdas3d_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gdas3d_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          allocate(gdas3d_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gdas3d_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gdas3d_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(gdas3d_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          
          call bilinear_interp_input(n, LDT_rc%met_gridDesc(findex,:), &
               gdas3d_struc(n)%n111,&
               gdas3d_struc(n)%n121,gdas3d_struc(n)%n211,&
               gdas3d_struc(n)%n221,gdas3d_struc(n)%w111,&
               gdas3d_struc(n)%w121,gdas3d_struc(n)%w211,&
               gdas3d_struc(n)%w221)

       endif
    enddo

  end subroutine init_GDAS3d
end module gdas3d_forcingMod

