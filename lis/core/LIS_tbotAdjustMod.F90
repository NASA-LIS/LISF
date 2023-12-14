!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module LIS_tbotAdjustMod

!BOP
!
! !MODULE: LIS_tbotAdjustMod
!
! !DESCRIPTION:
!  The code in this file adjusts the deep soil temperature based on options
!  included in WRF.
!
! \subsubsection{Overview}
!  Options exist to adjust the initial deep soil temperature to account
!  for terrain influences, and to update the deep soil temperature as a
!  function of lagged skin temperature. Also includes routines to write and
!  read time lag information to and from a restart file.
!  \begin{description}
!   \item[LIS\_tmnUpdateTileDec] 
!    User defined type containing time lag information used for
!    dynamically adjusting deep soil temperature as function of skin 
!    temperature for a single tile.
!   \item[LIS\_tmnUpdateDec]
!    User defined type containing arrays of type LIS\_tmnUpdateTileDec,
!    storing data for all tiles.
!   \item[LIS\_tmnUpdate]
!    Array of type LIS\_tmnUpdate, storing data for all domains.
!  \end{description}
!
! !REVISION HISTORY:
! 3 June 2013: Eric Kemp; Initial implementation
!
!EOP

   implicit none

   private

   type LIS_tmnUpdateTileDec
      real,allocatable :: tlag(:) ! Time series of daily average skin temperature
      real         :: tdly    ! Accumulated daily mean skin temperature of
                              ! the current day
      real         :: tyr     ! Annual mean skin temperature of previous year 
                              ! (K)
      real         :: tyra    ! Accumulated skin temperature in the current
                              ! year (K)
      integer :: nyear        ! Accumulated time (days) during current UTC year
      real :: nday            ! Accumulated time (sec) during current UTC day
   end type LIS_tmnUpdateTileDec

   type LIS_tmnUpdateDec
      type(LIS_tmnUpdateTileDec),allocatable :: tile(:)
   end type LIS_tmnUpdateDec

  type(LIS_tmnUpdateDec), allocatable :: LIS_tmnUpdate(:)

!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------

   public :: LIS_createTmnUpdate
   public :: LIS_initTmnUpdateTile
   public :: LIS_tbotTimeUtil
   public :: LIS_updateTbot
   !public :: LIS_tbotTerrainAdjustment
   public :: LIS_writeTmnUpdateRestart
   public :: LIS_readTmnUpdateRestart
   
contains

!BOP
! !ROUTINE: LIS_createTmnUpdate
! \label{LIS_createTmnUpdate}
!
! !INTERFACE:
   subroutine LIS_createTmnUpdate()
! !USES
      use LIS_coreMod,   only : LIS_rc
      use LIS_logMod,    only : LIS_logunit, LIS_endrun
! !DESCRIPTION:
!  This routine allocates the LIS\_tmnUpdate array and component data
!  structures. Initial values of all member variables is zero.
!
!EOP
      ! Local variables
      integer :: n,i
      integer :: lagday

      ! Bail if option is turned off.
      if (LIS_rc%tbot_update_lag .eq. 0) return

      ! Sanity check value of lagdays. Must be positive
      lagday = LIS_rc%tbot_lagday
      if (lagday < 1) then
         write(LIS_logunit,*) '[ERR] tbot_lagday is non-positive!'
         call LIS_endrun()         
      end if

      ! Allocate array of data structures, 1 per domain
      n = LIS_rc%nnest
      allocate(LIS_tmnUpdate(n))

      ! Allocate the tiles per domain, allocate daily mean skin temperature 
      ! time series array, and set values to zero. Later more appropriate
      ! values will be inserted on a tile by tile basis.
      do n = 1, LIS_rc%nnest
         allocate(LIS_tmnUpdate(n)%tile(LIS_rc%npatch(n,LIS_rc%lsm_index)))
         do i = 1, LIS_rc%npatch(n,LIS_rc%lsm_index)
            allocate(LIS_tmnUpdate(n)%tile(i)%tlag(lagday))
            LIS_tmnUpdate(n)%tile(i)%tlag = 0 ! Later use input deep soil temp
            LIS_tmnUpdate(n)%tile(i)%tdly = 0
            LIS_tmnUpdate(n)%tile(i)%tyr = 0 ! Later use input deep soil temp
            LIS_tmnUpdate(n)%tile(i)%tyra = 0         
            LIS_tmnUpdate(n)%tile(i)%nday = 0
            LIS_tmnUpdate(n)%tile(i)%nyear = 0
         end do
      end do

   end subroutine LIS_createTmnUpdate

   !---------------------------------------------------------------------------

!BOP
! !ROUTINE: LIS_initTmnUpdateTile
! \label{LIS_initTmnUpdateTile}
!
! !INTERFACE:
   subroutine LIS_initTmnUpdateTile(n,t,tmn)
! !USES:
      use LIS_coreMod,   only : LIS_rc
! !ARGUMENTS:
      integer,intent(in) :: n ! Domain
      integer,intent(in) :: t ! Tile
      real,intent(in) :: tmn  ! Original deep soil temperature
! !DESCRIPTION:
!
! Sets annual mean skin temperature for previous year and time series of
! daily mean skin temperatures to climatological values on tile by tile
! basis.
!EOP

      ! Local variables
      integer :: i

      ! Bail if option is turned off.
      if (LIS_rc%tbot_update_lag .eq. 0) return

      ! For cold start, set previous year annual mean skin temperature and 
      ! time series of daily mean skin temperatures to input climatological
      ! values.
      LIS_tmnUpdate(n)%tile(t)%tyr = tmn 
      do i = 1,LIS_rc%tbot_lagday
         LIS_tmnUpdate(n)%tile(t)%tlag(i) = tmn
      end do

   end subroutine LIS_initTmnUpdateTile

!------------------------------------------------------------------------------

   subroutine LIS_tbotTimeUtil(julianDay,yr)

      ! Uses
      use ESMF
      use LIS_timeMgrMod, only : LIS_clock

      ! Arguments
      real,intent(inout) :: julianDay
      integer,intent(inout) :: yr

      ! Local variables
      integer                 :: mo, da, hr, mn, ss, Sn, Sd
      type(ESMF_Time)         :: currTime
      integer                 :: status
      real(ESMF_KIND_R8)      :: rsec
      real(ESMF_KIND_R8)      :: dayOfYear_r8

      ! Get time information. Partially based on WRF since ESMF 3.1.0rp2
      ! doesn't return fractional day of the year.                            
      call ESMF_ClockGet(LIS_clock,currTime=currTime,rc=status)
      call ESMF_TimeGet(currTime,yy=yr,s=ss,sN=Sn,sD=Sd,rc=status)
      ! 64-bit IEEE 754 has 52-bit mantisssa -- only need 25 bits to hold
      ! number of seconds in a year...
      rsec = real(ss, ESMF_KIND_R8)
      if (Sd .ne. 0) then
         rsec = rsec + ( real(Sn, ESMF_KIND_R8) / real(Sd, ESMF_KIND_R8) )
      end if
      dayOfYear_r8 = rsec / real(86400_ESMF_KIND_I8, ESMF_KIND_R8)
      julianDay = real(dayOfYear_r8)

   end subroutine LIS_tbotTimeUtil

!------------------------------------------------------------------------------

!BOP
! !ROUTINE: LIS_updateTbot
! \label{LIS_updateTbot}
!
! !INTERFACE:
   subroutine LIS_updateTbot(n,i,julian_in,yr,dt,tsk,tmn)
! !USES:
      use LIS_coreMod,   only : LIS_rc
      use LIS_logMod,    only : LIS_logunit
! !ARGUMENTS:
      integer,intent(in) :: n      ! Domain
      integer,intent(in) :: i      ! Tile in domain
      real,intent(in) :: julian_in ! Day of year
      integer,intent(in) :: yr     ! Current year
      real,intent(in) :: dt        ! Time step
      real,intent(in) :: tsk       ! Current skin temperature
      real,intent(inout) :: tmn    ! Modified deep soil temperature
! !DESCRIPTION:
! 
! Dynamically adjusts deep soil temperature as weighted average of previous
! year's annual mean skin temperature and mean of time series of recent
! daily mean skin temperatures. The length of the time series is set by
! lagday. Based on WRF.
!EOP
      ! Local variables
      real :: yrday
      real :: deltat
      real :: tprior
      real :: julian
      integer :: lagday
      integer :: m

      real, parameter :: INVDAY = 1./86400.
      real, parameter :: TCONST = 0.6

      ! Bail if option is turned off.
      if (LIS_rc%tbot_update_lag .eq. 0) return

      lagday = LIS_rc%tbot_lagday

      ! Set days in year. This is the code used in WRF for leap years
      yrday = 365.
      if (mod(yr,4).eq.0) yrday=366.

      ! Integrate skin temperature and time in current day.
      LIS_tmnUpdate(n)%tile(i)%tdly = LIS_tmnUpdate(n)%tile(i)%tdly + tsk*dt
      LIS_tmnUpdate(n)%tile(i)%nday = LIS_tmnUpdate(n)%tile(i)%nday + 1.*dt

      ! Update deep soil temperature
      ! If it is the end of the day, update variables
      deltat=(julian_in-nint(julian_in))*86400.
      if (abs(deltat) .le. dt*0.5) then ! Start of new day
         julian=(julian_in-1.)+(dt*INVDAY)
         tprior = 0.0
         do m = 1,lagday
            tprior = tprior + LIS_tmnUpdate(n)%tile(i)%tlag(m)
         end do
         tprior=tprior/float(lagday)
         tmn = TCONST*LIS_tmnUpdate(n)%tile(i)%tyr + &
               (1.-TCONST)*tprior
         if (tmn < 250.) then
            write(LIS_logunit,*) '[ERR] Cold tmn for domain/tile: ',n,i
            write(LIS_logunit,*) '[ERR] tmn = ',tmn
            write(LIS_logunit,*) '[ERR] tyr = ',LIS_tmnUpdate(n)%tile(i)%tyr
            write(LIS_logunit,*) '[ERR] tprior = ',tprior
         end if

         do m=1,lagday-1
            LIS_tmnUpdate(n)%tile(i)%tlag(m) = &
                 LIS_tmnUpdate(n)%tile(i)%tlag(m+1)            
         end do
         LIS_tmnUpdate(n)%tile(i)%tlag(lagday) = &
              LIS_tmnUpdate(n)%tile(i)%tdly / LIS_tmnUpdate(n)%tile(i)%nday
         LIS_tmnUpdate(n)%tile(i)%tdly = 0.0
         LIS_tmnUpdate(n)%tile(i)%nday = 0.

         ! Update tyr if it is the end of the year
         if ((yrday-julian).le.1.) then
            LIS_tmnUpdate(n)%tile(i)%tyr = &
                 LIS_tmnUpdate(n)%tile(i)%tyra / LIS_tmnUpdate(n)%tile(i)%nyear
            LIS_tmnUpdate(n)%tile(i)%tyra = 0.0
            LIS_tmnUpdate(n)%tile(i)%nyear = 0
         else
            LIS_tmnUpdate(n)%tile(i)%tyra = &
                 LIS_tmnUpdate(n)%tile(i)%tyra + &
                 LIS_tmnUpdate(n)%tile(i)%tlag(lagday)
            LIS_tmnUpdate(n)%tile(i)%nyear = &
                 LIS_tmnUpdate(n)%tile(i)%nyear + 1
         end if
      end if

      return
   end subroutine LIS_updateTbot

!------------------------------------------------------------------------------

!BOP
! !ROUTINE: LIS_writeTmnUpdateRestart
! \label{LIS_writeTmnUpdateRestart}
!
! !INTERFACE:
   subroutine LIS_writeTmnUpdateRestart(n,ftn,dimID,wformat)
! !USES:
      use LIS_coreMod,    only : LIS_rc, LIS_masterproc
      use LIS_historyMod, only : LIS_writevar_restart, LIS_writeHeader_restart
      use LIS_logMod,     only : LIS_logunit, LIS_verify
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
! !ARGUMENTS:
      integer,intent(in)          :: n
      integer,intent(in)          :: ftn
      integer,intent(in)          :: dimID(11)
      character(len=*),intent(in) :: wformat
! !DESCRIPTION:
!
! Writes dynamic deep soil temperature time lag information to an already
! opened restart file. Assumes LSM data has already been written out.
!EOP
      ! Local variables
      integer :: l,t
      integer :: tbot_lagdayid,tlagid,tdlyid,tyrid,tyraid,ndayid,nyearid
      real, allocatable :: tmptilen(:)

      if (LIS_rc%tbot_update_lag .eq. 0) return

      if ( LIS_masterproc ) then
         if ( wformat == "netcdf" ) then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
         call LIS_verify(nf90_redef(ftn), &
              'nf90_redef failed in LIS_writeTmnUpdateRestart')

         call LIS_verify(nf90_def_dim(ftn,'tbot_lagday',LIS_rc%tbot_lagday,&
              tbot_lagdayid),&
              'nf90_def_dim for tbot_lagday failed '//&
              'in LIS_writeTmnUpdateRestart')

         call LIS_writeHeader_restart(ftn,n,dimID,tlagid, "tlag", &
              "daily average skin temperature",&
              "K",vlevels=LIS_rc%tbot_lagday,valid_min=220.0, valid_max=330.0,&
               var_flag = "tbot_lagday")
         call LIS_writeHeader_restart(ftn,n,dimID,tdlyid, "tdly", &
              "daily mean skin temperature",&
              "K",vlevels=1,valid_min=220.0, valid_max=330.0)
         call LIS_writeHeader_restart(ftn,n,dimID,tyrid, "tyr", &
              "mean skin temperature of previous year",&
              "K",vlevels=1,valid_min=220.0, valid_max=330.0)
         call LIS_writeHeader_restart(ftn,n,dimID,tyraid, "tyra", &
              "skin temperature in current year",&
              "K",vlevels=1,valid_min=0.0, valid_max=120450.0)
         call LIS_writeHeader_restart(ftn,n,dimID,ndayid, "nday", &
              "Accumulated time (sec) during current UTC day",&
              "-",vlevels=1,valid_min=0.0, valid_max=31622400.0)
         call LIS_writeHeader_restart(ftn,n,dimID,nyearid, "nyear", &
              "Accumulated time (days) during current UTC year",&
              "-",vlevels=1,valid_min=0.0, valid_max=366.0)

         call LIS_verify(nf90_enddef(ftn), &
              'nf90_enddef failed in LIS_writeTmnUpdateRestart')
         endif
#endif
      endif

      if (LIS_masterproc) then
         if ( wformat == "netcdf" ) then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
            call LIS_verify(nf90_redef(ftn), &
                 'nf90_redef failed in LIS_writeTmnUpdateRestart')
            call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"tbot_update_lag", &
                 LIS_rc%tbot_update_lag),&
                 'nf90_put_att failed for tbot_update_lag')
            call LIS_verify(nf90_enddef(ftn), &
                 'nf90_enddef failed in LIS_writeTmnUpdateRestart')
#endif
         else
            write(ftn) LIS_rc%tbot_update_lag
         endif
      end if

      if (LIS_masterproc) then
         if ( wformat == "netcdf" ) then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
            call LIS_verify(nf90_redef(ftn), &
                 'nf90_redef failed in LIS_writeTmnUpdateRestart')
            call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"tbot_lagday", &
                 LIS_rc%tbot_lagday),&
                 'nf90_put_att failed for tbot_lagday')
            call LIS_verify(nf90_enddef(ftn), &
                 'nf90_enddef failed in LIS_writeTmnUpdateRestart')
#endif
         else
            write(ftn) LIS_rc%tbot_lagday
         endif
      end if         

      allocate(tmptilen(LIS_rc%npatch(n,LIS_rc%lsm_index)))
      tmptilen(:) = 0.0
      do l=1,LIS_rc%tbot_lagday
         do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
            tmptilen(t) = LIS_tmnUpdate(n)%tile(t)%tlag(l)
         end do
         call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,tmptilen, &
              varid=tlagid,dim=l,wformat=wformat)
      end do
      call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,&
           LIS_tmnUpdate(n)%tile%tdly,&
           varid=tdlyid,dim=1,wformat=wformat)
      call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,&
           LIS_tmnUpdate(n)%tile%tyr,&
           varid=tyrid,dim=1,wformat=wformat)
      call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,&
           LIS_tmnUpdate(n)%tile%tyra,&
           varid=tyraid,dim=1,wformat=wformat)
      call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,&
           LIS_tmnUpdate(n)%tile%nday,&
           varid=ndayid,dim=1,wformat=wformat)
      call LIS_writevar_restart(ftn,n,LIS_rc%lsm_index,&
           LIS_tmnUpdate(n)%tile%nyear,&
           varid=nyearid,dim=1,wformat=wformat)

      deallocate(tmptilen)
      return
   end subroutine LIS_writeTmnUpdateRestart

!------------------------------------------------------------------------------

!BOP
! !ROUTINE: LIS_readTmnUpdateRestart
! \label{LIS_readTmnUpdateRestart}
!
! !INTERFACE:
   subroutine LIS_readTmnUpdateRestart(n,ftn,wformat)
! !USES:
      use LIS_coreMod,    only : LIS_rc, LIS_masterproc
      use LIS_historyMod, only : LIS_readvar_restart
      use LIS_logMod,     only : LIS_logunit, LIS_endrun, LIS_verify
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
! !ARGUMENTS:
      integer,intent(in)          :: n ! Domain
      integer,intent(in)          :: ftn
      character(len=*),intent(in) :: wformat
! !DESCRIPTION:
!
! Reads in dynamic deep soil temperature time lag information from already
! opened restart file. Assumes LSM data has already been read in.
!EOP
      ! Local variables
      integer :: tbot_update_lag
      integer :: tbot_lagday
      integer :: l,t,rc
      real, allocatable :: tmptilen(:)

      if ( LIS_rc%tbot_update_lag == 0 ) then
         write(LIS_logunit,*) '[INFO] dynamic deep soil temperature updating '//&
                              'was disabled in lis.config file.'
         write(LIS_logunit,*) '[INFO] returning from LIS_readTmnUpdateRestart'
         return
      endif

      if ( wformat == "netcdf" ) then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
         rc = nf90_get_att(ftn,NF90_GLOBAL,"tbot_update_lag", &
              tbot_update_lag)
#endif
      else
         read(ftn, iostat=rc) tbot_update_lag
      endif

      ! tbot_update_lag will be missing from older restart files and
      ! from newer restart files where LIS_rc%tbot_update_lag was set to 0.
      ! In these cases, simply return.
      if ( rc /= 0 ) then
         write(LIS_logunit,*) '[WARN] restart file is missing the '       // &
                              'dynamic deep soil temperature time lag ' // &
                              'information'
         write(LIS_logunit,*) '[WARN] returning from LIS_readTmnUpdateRestart'
         return
      endif

!      if (tbot_update_lag .ne. LIS_rc%tbot_update_lag) then
!         write(LIS_logunit,*) 'ERROR, mismatch between tbot_update_lag settings!'
!         write(LIS_logunit,*) 'Expected ',LIS_rc%tbot_update_lag
!         write(LIS_logunit,*) 'Found ',tbot_update_lag
!         call LIS_endrun()
!      end if
!
!      if (tbot_update_lag .eq. 0) return

      if ( wformat == "netcdf" ) then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
         call LIS_verify(nf90_get_att(ftn,NF90_GLOBAL,"tbot_lagday", &
              tbot_lagday),&
              'nf90_get_att failed for tbot_lagday')
#endif
      else
         read(ftn) tbot_lagday
      endif

      if (tbot_lagday .ne. LIS_RC%tbot_lagday) then
         write(LIS_logunit,*) '[ERR] mismatch between tbot_lagday settings!'
         write(LIS_logunit,*) '[ERR] Expected ',LIS_rc%tbot_lagday
         write(LIS_logunit,*) '[ERR] Found ',tbot_lagday
         call LIS_endrun()
      end if

      allocate(tmptilen(LIS_rc%npatch(n,LIS_rc%lsm_index)))

      do l=1,LIS_rc%tbot_lagday
         call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,&
              tmptilen,varname="tlag",&
              dim=l,vlevels=LIS_rc%tbot_lagday,wformat=wformat)
         do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
            LIS_tmnUpdate(n)%tile(t)%tlag(l) = tmptilen(t)
         end do
      end do

      deallocate(tmptilen)

      call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,&
           LIS_tmnUpdate(n)%tile%tdly,&
           varname="tdly",wformat=wformat)
      call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,&
           LIS_tmnUpdate(n)%tile%tyr,&
           varname="tyr",wformat=wformat)
      call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,&
           LIS_tmnUpdate(n)%tile%tyra,&
           varname="tyra",wformat=wformat)
      call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,&
           LIS_tmnUpdate(n)%tile%nday,&
           varname="nday",wformat=wformat)
      call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,&
           LIS_tmnUpdate(n)%tile%nyear,&
           varname="nyear",wformat=wformat)
   
   end subroutine LIS_readTmnUpdateRestart

!------------------------------------------------------------------------------

#if 0   
!BOP
!
! !ROUTINE: LIS_tbotTerrainAdjustment
! \label{LIS_tbotTerrainAdjustment}
!
! !REVISION HISTORY:
!  30 May 2013: Eric Kemp; Initial Code
!
! !INTERFACE:
   subroutine LIS_tbotTerrainAdjustment(nx,ny,placetbot1,n)
! !USES:
      use LIS_coreMod,    only : LIS_rc, LIS_domain
      use LIS_logMod,     only : LIS_logunit, LIS_endrun
      use LIS_topoMod,    only : LIS_topo       
! !ARGUMENTS:
      integer,intent(in) :: nx ! Global x-dimension
      integer,intent(in) :: ny ! Global y-dimension
      real,intent(inout) :: placetbot1(nx,ny) ! 2-D climatological deep soil
                                              ! temperature (K)
      integer,intent(in) :: n ! Domain
! !DESCRIPTION:
!  This subroutine adjusts the input deep soil temperature for the Noah LSM.
!  This original source for this is a time-averaged skin temperature 
!  climatology from ECMWF, which was then extrapolated to mean sea level using
!  the standard atmosphere lapse rate (6.5 K / km). This process is reversed
!  to adjust the temperatures to the LIS terrain.  Based on WRF.
!EOP

! Local variables
      integer :: i
      real :: tmp_t, tmp_topo
      
! Sanity check dimensions
      if (nx .ne. LIS_rc%lnc(n) .or. &
          ny .ne. LIS_rc%lnr(n)) then
         write(LIS_logunit,*) '[ERR] dimension mismatch for LIS domain!'
         write(LIS_logunit,*) '[ERR] Expected ',nx,ny
         write(LIS_logunit,*) '[ERR] Found ',LIS_rc%lnc(n),LIS_rc%lnr(n)
         call LIS_endrun()
      end if

      ! Only adjust terrain if requested by the user
      if (LIS_rc%tbot_terrain_adj .eq. 1) then
         do i = 1, LIS_rc%npatch(n,LIS_rc%lsm_index)
            tmp_t = &
                 placetbot1(LIS_domain(n)%tile(i)%col, &
                            LIS_domain(n)%tile(i)%row)
            if (tmp_t .ne. -9999.00) then
               tmp_topo = &
                    LIS_topo(n)%elevation(LIS_domain(n)%tile(i)%col, &
                                          LIS_domain(n)%tile(i)%row)
               tmp_t = tmp_t - 0.0065 * tmp_topo
               
               placetbot1(LIS_domain(n)%tile(i)%col, &
                          LIS_domain(n)%tile(i)%row) = tmp_t
            end if
         end do
      end if
      
      return
   end subroutine LIS_tbotTerrainAdjustment
#endif

end module LIS_tbotAdjustMod
