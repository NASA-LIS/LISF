!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! !MODULE: HYMAP3_managMod
!
! !DESCRIPTION
!
! !REVISION HISTORY:
! 1 May 2021: Augusto Getirana, Initial implementation
!17 Apr 2023: Augusto Getirana, Add option for prescribed monthly climatology
!
module HYMAP3_managMod
  !module containing routines for water management rules

  implicit none

  private

  public :: HYMAP3_manag_rules
  public :: HYMAP3_get_management_rules

contains
  ! ================================================
  subroutine HYMAP3_manag_rules(time, managtype, nmanagcoef, managqmax, &
       managact, managcoef, rivout, rivinf, rivinf1)

    use LIS_logMod

    implicit none

    real*8,  intent(in)    :: time         !current time [year]
    integer, intent(in)    :: managtype,nmanagcoef
    integer, intent(inout) :: managact
    real,    intent(in)    :: managcoef(nmanagcoef)
    real,    intent(inout) :: rivout
    real,    intent(inout) :: rivinf,rivinf1,managqmax

    real                   :: zdistmp
    integer                :: im

    if(managact==0.and.rivout>=managcoef(10))managact=1
    if(managact==1.and.rivout<managcoef(11))then
       managqmax=0
       managact=0
    endif
    if((managtype==1.or.managtype==2).and.managact==0)return

    !two polynomial equations (ax2 + bx + c) used as a function of
    !pre-defined min-max discharge
    zdistmp=0.
    if(managtype==1)then
       ! single linear equation (ax + b)
       if(rivout>=managcoef(7).and.rivout<managcoef(8))then
          zdistmp=managcoef(1)*(rivout**2)+managcoef(2)*rivout+managcoef(3)
       elseif(rivout>=managcoef(8).and.rivout<=managcoef(9))then
          zdistmp=managcoef(4)*(rivout**2)+managcoef(5)*rivout+managcoef(6)
       endif
       zdistmp=max(0.,min(zdistmp,rivinf))
       rivinf=max(0.,rivinf-zdistmp)
       rivinf1=rivinf1+zdistmp

       ! two polynomial equations (ax2 + bx + c) used as a function of
       ! simulated peak discharge
    elseif(managtype==2)then
       if(rivout>=managqmax)then
          zdistmp=min(max(managcoef(1)*(rivout**2)+ &
               managcoef(2)*rivout+managcoef(3),managcoef(7)),managcoef(8))
       else
          zdistmp=min(max(managcoef(4)*(rivout**2)+ &
               managcoef(5)*rivout+managcoef(6),managcoef(7)),managcoef(8))
       endif
       managqmax=rivout
       zdistmp=max(0.,min(zdistmp,rivinf))
       rivinf=max(0.,rivinf-zdistmp)
       rivinf1=rivinf1+zdistmp

    !prescribed monthly climatology of water use
    !here, managcoef corresponds to monthly water use amounts in m3/s
    !this option considers that the water use doesn't go to another
    !grid point
    elseif(managtype==3)then
       im=mod(int(time)/100,100)
       zdistmp=min(min(managcoef(im),rivinf),rivout)
       rivout=max(0.,rivout-zdistmp)
       rivinf=max(0.,rivinf-zdistmp)
    else
       write(LIS_logunit,*) &
            '[ERR] water management type not defined ', managtype
       call LIS_endrun()
    endif

  end subroutine HYMAP3_manag_rules
  !=============================================
  subroutine HYMAP3_get_management_rules(managheader, nx, ny, sindex, &
       nmanag, &
       nmanagcoef, managloc, managqmax, managact, managtype, managcoef)

    use LIS_logMod

    implicit none

    character(*), intent(in)  :: managheader
    integer,      intent(in)  :: nx,ny
    integer,      intent(in)  :: sindex(nx,ny)
    integer,      intent(in)  :: nmanag,nmanagcoef
    integer,      intent(out) :: managloc(nmanag,2)
    integer,      intent(out) :: managact(nmanag)
    integer,      intent(out) :: managtype(nmanag)
    real,         intent(out) :: managcoef(nmanag,nmanagcoef), &
         managqmax(nmanag)

    integer                   :: ii,ix,iy,jx,jy
    logical                   :: file_exists
    integer :: ftn

    inquire(file=managheader,exist=file_exists)
    if(file_exists) then
       ftn = LIS_getNextUnitNumber()
       write(LIS_logunit,*)&
            '[INFO] get management rules '//trim(managheader)
       open(ftn,file=trim(managheader), status='old')
       read(ftn,*)ii
       do ii=1,nmanag
          read(ftn,*,end=10)ix,iy,jx,jy,managtype(ii),managcoef(ii,:)
          managloc(ii,1)=sindex(ix,iy)
          managloc(ii,2)=sindex(jx,jy)
          managact(ii)=0
          managqmax(ii)=0.
       enddo
       close(ftn)
       call LIS_releaseUnitNumber(ftn)
    else
       write(LIS_logunit,*) &
            '[ERR] HYMAP3_routingMod: Failed to open file ' &
            //trim(managheader)
       call LIS_endrun()
    endif
    return
10  continue
    write(LIS_logunit,*) '[ERR] failed in getting management rules'
    write(LIS_logunit,*) '[ERR] check header file '//trim(managheader)
    call LIS_endrun()

  end subroutine HYMAP3_get_management_rules
  !=============================================
end module HYMAP3_managMod
