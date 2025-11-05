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
! !ROUTINE: HYMAP3_routing_readrst
! \label{HYMAP3_routing_readrst}
!
! !REVISION HISTORY:
! 15 Nov 2011: Augusto Getirana;  Initial implementation
! 19 Jan 2016: Augusto Getirana;  Inclusion of four Local Inertia variables
! 10 Mar 2019: Sujay Kumar;       Added support for NetCDF and parallel
!                                 processing.
! 27 Apr 2020: Augusto Getirana;  Added support for urban drainage
!
! !INTERFACE:
subroutine HYMAP3_routing_readrst

! !USES:
  use ESMF
  use HYMAP3_routingMod, only : HYMAP3_routing_struc
  use LIS_constantsMod, only: LIS_CONST_PATH_LEN
  use LIS_coreMod
  use LIS_fileIOMod
  use LIS_logMod
  use LIS_timeMgrMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
!
! !DESCRIPTION:
!  This routine reads NetCDF formatted HYMAP3 restart files.
!
!EOP

  implicit none

  integer       :: n
  integer       :: ftn
  character(LIS_CONST_PATH_LEN) :: filename
  logical       :: read_restart
  integer           :: yr,mo,da,hr,mn,ss,doy
  real*8            :: time
  real              :: gmt
  integer :: status

  external :: HYMAP3_readvar_restart
  external :: HYMAP3_readbif_restart
  external :: HYMAP3_readvar_restart_ens
  external :: HYMAP3_readbif_restart_ens

  do n=1, LIS_rc%nnest

     read_restart = .false.

     if(HYMAP3_routing_struc(n)%startmode.eq."restart") then !cold start
        read_restart = .true.
     endif

     if(LIS_rc%runmode.eq."ensemble smoother") then
        if(LIS_rc%iterationId(n).gt.1) then
           read_restart = .true.

           if(HYMAP3_routing_struc(n)%rstInterval.eq.2592000) then
              !create the restart filename based on the timewindow start
              ! time
              call ESMF_TimeGet(LIS_twStartTime,yy=yr,mm=mo,&
                   dd=da,calendar=LIS_calendar,rc=status)
              hr = 0
              mn = 0
              ss = 0
              call LIS_tick(time,doy,gmt,yr,mo,da,hr,mn,ss,&
                   (-1)*(HYMAP3_routing_struc(n)%dt))
           else
              call ESMF_TimeGet(LIS_twStartTime,yy=yr,mm=mo,&
                   dd=da,calendar=LIS_calendar,rc=status)
              hr = 0
              mn = 0
              ss = 0
           endif

           call LIS_create_restart_filename(n,filename, &
                'ROUTING','HYMAP3_router', &
                yr, mo, da, hr, mn, ss, &
                wformat="netcdf")

           HYMAP3_routing_struc(n)%rstfile = filename
        endif
     endif

     if(read_restart) then
        write(LIS_logunit,*) '[INFO] HYMAP3 restart file used: ', &
             trim(HYMAP3_routing_struc(n)%rstfile)

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
        status = nf90_open(path=HYMAP3_routing_struc(n)%rstfile, &
             mode=NF90_NOWRITE, ncid=ftn)
        call LIS_verify(status, "Error opening file "//&
             trim(HYMAP3_routing_struc(n)%rstfile))
#endif
        if(HYMAP3_routing_struc(n)%useens.eq.0) then
           call HYMAP3_readvar_restart(ftn,n,&
                HYMAP3_routing_struc(n)%rivsto,&
                "RIVSTO")
           call HYMAP3_readvar_restart(ftn,n,&
                HYMAP3_routing_struc(n)%fldsto,&
                "FLDSTO")
           call HYMAP3_readvar_restart(ftn,n,&
                HYMAP3_routing_struc(n)%rnfsto,&
                "RNFSTO")
           call HYMAP3_readvar_restart(ftn,n,&
                HYMAP3_routing_struc(n)%bsfsto,&
                "BSFSTO")
           !ag(13Jan2021)
           !write out variables below if flowtype is not kinematic wave
           if(HYMAP3_routing_struc(n)%flowtype/=1)then
              call HYMAP3_readvar_restart(ftn,n,&
                   HYMAP3_routing_struc(n)%rivout_pre,&
                   "RIVOUT_PRE")
              call HYMAP3_readvar_restart(ftn,n,&
                   HYMAP3_routing_struc(n)%rivdph_pre,&
                   "RIVDPH_PRE")
              call HYMAP3_readvar_restart(ftn,n,&
                   HYMAP3_routing_struc(n)%fldout_pre,&
                   "FLDOUT_PRE")
              call HYMAP3_readvar_restart(ftn,n,&
                   HYMAP3_routing_struc(n)%flddph_pre,&
                   "FLDDPH_PRE")
           endif
           !ag (27Apr2020)
           !urban drainage storage
           if(HYMAP3_routing_struc(n)%flowtype==4)then
              call HYMAP3_readvar_restart(ftn,n,&
                   HYMAP3_routing_struc(n)%drsto,&
                   "DRSTO")
              call HYMAP3_readvar_restart(ftn,n,&
                   HYMAP3_routing_struc(n)%drout,&
                   "DROUT")
           endif
           !ag(8Sep2022)
           !bifurcation outflow previous time step
           if(HYMAP3_routing_struc(n)%bifflag==1)then
              call HYMAP3_readbif_restart(ftn,n,&
                   HYMAP3_routing_struc(n)%bifout_pre,&
                   "BIFOUT_PRE")
           endif
        else
           call HYMAP3_readvar_restart_ens(ftn,n,&
                HYMAP3_routing_struc(n)%rivsto,&
                "RIVSTO")
           call HYMAP3_readvar_restart_ens(ftn,n,&
                HYMAP3_routing_struc(n)%fldsto,&
                "FLDSTO")
           call HYMAP3_readvar_restart_ens(ftn,n,&
                HYMAP3_routing_struc(n)%rnfsto,&
                "RNFSTO")
           call HYMAP3_readvar_restart_ens(ftn,n,&
                HYMAP3_routing_struc(n)%bsfsto,&
                "BSFSTO")
           !ag(13Jan2021)
           !write out variables below if flowtype is not kinematic wave
           if(HYMAP3_routing_struc(n)%flowtype/=1)then
              call HYMAP3_readvar_restart_ens(ftn,n,&
                   HYMAP3_routing_struc(n)%rivout_pre,&
                   "RIVOUT_PRE")
              call HYMAP3_readvar_restart_ens(ftn,n,&
                   HYMAP3_routing_struc(n)%rivdph_pre,&
                   "RIVDPH_PRE")
              call HYMAP3_readvar_restart_ens(ftn,n,&
                   HYMAP3_routing_struc(n)%fldout_pre,&
                   "FLDOUT_PRE")
              call HYMAP3_readvar_restart_ens(ftn,n,&
                   HYMAP3_routing_struc(n)%flddph_pre,&
                   "FLDDPH_PRE")
           endif
           !ag (27Apr2020)
           !urban drainage storage
           if(HYMAP3_routing_struc(n)%flowtype==4)then
              call HYMAP3_readvar_restart_ens(ftn,n,&
                   HYMAP3_routing_struc(n)%drsto,&
                   "DRSTO")
              call HYMAP3_readvar_restart_ens(ftn,n,&
                   HYMAP3_routing_struc(n)%drout,&
                   "DROUT")
           endif
           !ag(8Sep2022)
           !bifurcation outflow previous time step
           if(HYMAP3_routing_struc(n)%bifflag==1)then
              call HYMAP3_readbif_restart_ens(ftn,n,&
                   HYMAP3_routing_struc(n)%bifout_pre,&
                   "BIFOUT_PRE")
           endif
        endif

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
        status = nf90_close(ftn)
        call LIS_verify(status, &
             "Error in nf90_close in HYMAP3_routing_readrst")
#endif
     endif
  enddo

end subroutine HYMAP3_routing_readrst

!BOP
! !ROUTINE: HYMAP3_readvar_restart
! \label{HYMAP3_readvar_restart}
!
! !INTERFACE:
subroutine HYMAP3_readvar_restart(ftn, n, var, varname)

! !USES:
  use HYMAP3_routingMod
  use LIS_coreMod
  use LIS_logMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none

! !ARGUMENTS:
  integer, intent(in)   :: ftn
  integer, intent(in)   :: n
  real, intent(inout)   :: var(LIS_rc%nroutinggrid(n))
  character(len=*), intent(in) :: varname

! !DESCRIPTION:
!  Reads a real variable from a NetCDF restart file.
!
!  The arguments are:
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
!EOP
  real, allocatable :: gtmp(:)
  integer           :: varid
  integer           :: i,ix,iy,ix1,iy1
  integer           :: status

  allocate(gtmp(LIS_rc%glbnroutinggrid(n)))
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  status = nf90_inq_varid(ftn,trim(varname),varid)
  call LIS_verify(status, &
       'Error in nf90_inq_varid in HYMAP3_readvar_restart')
  status = nf90_get_var(ftn,varid,gtmp)
  call LIS_verify(status, &
       'Error in nf90_get_var in HYMAP3_readvar_restart')
#endif

  do i=1,LIS_rc%nroutinggrid(n)
     ix = HYMAP3_routing_struc(n)%seqx(i)
     iy = HYMAP3_routing_struc(n)%seqy(i)
     ix1 = ix + LIS_ews_halo_ind(n,LIS_localPet+1) -1
     iy1 = iy + LIS_nss_halo_ind(n,LIS_localPet+1) -1
     var(i)  = gtmp(HYMAP3_routing_struc(n)%sindex(ix1,iy1))
  enddo

  deallocate(gtmp)
end subroutine HYMAP3_readvar_restart

!BOP
! !ROUTINE: HYMAP3_readvar_restart_ens
! \label{HYMAP3_readvar_restart_ens}
!
! !INTERFACE:
subroutine HYMAP3_readvar_restart_ens(ftn, n, var, varname)
! !USES:
  use HYMAP3_routingMod
  use LIS_coreMod
  use LIS_logMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none

! !ARGUMENTS:
  integer, intent(in)   :: ftn
  integer, intent(in)   :: n
  real, intent(inout)   :: var(LIS_rc%nroutinggrid(n), &
       LIS_rc%nensem(n))
  character(len=*), intent(in) :: varname

! !DESCRIPTION:
!  Reads a real variable from a NetCDF restart file.
!
!  The arguments are:
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
!EOP
  real, allocatable :: gtmp(:)
  integer           :: varid
  integer           :: i,m,ix,iy,ix1,iy1
  integer           :: status

  allocate(gtmp(LIS_rc%glbnroutinggrid(n)*LIS_rc%nensem(n)))
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  status = nf90_inq_varid(ftn,trim(varname),varid)
  call LIS_verify(status, &
       'Error in nf90_inq_varid in HYMAP3_readvar_restart_ens')
  status = nf90_get_var(ftn,varid,gtmp)
  call LIS_verify(status, &
       'Error in nf90_get_var in HYMAP3_readvar_restart_ens')
#endif
  do i=1,LIS_rc%nroutinggrid(n)
     do m=1,LIS_rc%nensem(n)
        ix = HYMAP3_routing_struc(n)%seqx(i)
        iy = HYMAP3_routing_struc(n)%seqy(i)
        ix1 = ix + LIS_ews_halo_ind(n,LIS_localPet+1) -1
        iy1 = iy + LIS_nss_halo_ind(n,LIS_localPet+1) -1

        var(i,m)  = gtmp(m + &
             (HYMAP3_routing_struc(n)%sindex(ix1,iy1)-1)*LIS_rc%nensem(n))

     enddo
  enddo
  deallocate(gtmp)
end subroutine HYMAP3_readvar_restart_ens

!BOP
! !ROUTINE: HYMAP3_readbif_restart
! \label{HYMAP3_readbif_restart}
!
! !INTERFACE:
subroutine HYMAP3_readbif_restart(ftn, n, var, varname)
! !USES:
  use HYMAP3_routingMod
  use LIS_coreMod
  use LIS_logMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none

! !ARGUMENTS:
  integer, intent(in)   :: ftn
  integer, intent(in)   :: n
  real, intent(inout)   :: var(HYMAP3_routing_struc(n)%nbif)
  character(len=*), intent(in) :: varname

! !DESCRIPTION:
!  Reads a real variable from a NetCDF restart file.
!
!  The arguments are:
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
!EOP
  integer           :: varid
  integer           :: status

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  status = nf90_inq_varid(ftn,trim(varname),varid)
  call LIS_verify(status, &
       'Error in nf90_inq_varid in HYMAP3_readbif_restart')
  status = nf90_get_var(ftn,varid,var)
  call LIS_verify(status, &
       'Error in nf90_get_var in HYMAP3_readbif_restart')
#endif

end subroutine HYMAP3_readbif_restart

!BOP
! !ROUTINE: HYMAP3_readbif_restart_ens
! \label{HYMAP3_readbif_restart_ens}
!
! !INTERFACE:
subroutine HYMAP3_readbif_restart_ens(ftn, n, var, varname)
! !USES:
  use HYMAP3_routingMod
  use LIS_coreMod
  use LIS_logMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none

! !ARGUMENTS:
  integer, intent(in)   :: ftn
  integer, intent(in)   :: n
  real, intent(inout)   :: var(HYMAP3_routing_struc(n)%nbif,&
       LIS_rc%nensem(n))
  character(len=*), intent(in) :: varname

! !DESCRIPTION:
!  Reads a real variable from a NetCDF restart file.
!
!  The arguments are:
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
!EOP

  write(LIS_logunit,*) &
       '[ERR] HYMAP3_readbif_restart_ens not implemented'
  call LIS_endrun

end subroutine HYMAP3_readbif_restart_ens

