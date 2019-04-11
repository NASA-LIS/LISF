!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.1
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"

!BOP
! !ROUTINE: HYMAP2_routing_readrst
! \label{HYMAP2_routing_readrst} 
!
! !REVISION HISTORY:
! 15 Nov 2011: Augusto Getirana;  Initial implementation
! 19 Jan 2016: Augusto Getirana;  Inclusion of four Local Inertia variables
! 10 Mar 2019: Sujay Kumar;       Added support for NetCDF and parallel 
!                                 processing. 
!  
! !INTERFACE: 
subroutine HYMAP2_routing_readrst
! !USES: 
  use ESMF
  use LIS_fileIOMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use HYMAP2_routingMod, only : HYMAP2_routing_struc

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
! 
! !DESCRIPTION: 
!  This routine reads NetCDF formatted HYMAP2 restart files. 
! 
!EOP

  implicit none

  integer       :: n 
  integer       :: ftn
  integer       :: i,j,k
  integer       :: ios,status
  character*100 :: filename
  logical       :: read_restart
  integer           :: yr,mo,da,hr,mn,ss,doy
  real*8            :: time
  real              :: gmt

  do n=1, LIS_rc%nnest
     
     read_restart = .false. 
     
     if(HYMAP2_routing_struc(n)%startmode.eq."restart") then !cold start
        read_restart = .true. 
     endif
     
     if(LIS_rc%runmode.eq."ensemble smoother") then 
        if(LIS_rc%iterationId(n).gt.1) then 
           read_restart = .true. 
           
           if(HYMAP2_routing_struc(n)%rstInterval.eq.2592000) then 
              !create the restart filename based on the timewindow start time
              call ESMF_TimeGet(LIS_twStartTime,yy=yr,mm=mo,&
                   dd=da,calendar=LIS_calendar,rc=status)
              hr = 0 
              mn = 0 
              ss = 0 
              call LIS_tick(time,doy,gmt,yr,mo,da,hr,mn,ss,&
                   (-1)*(HYMAP2_routing_struc(n)%dt))
           else
              call ESMF_TimeGet(LIS_twStartTime,yy=yr,mm=mo,&
                   dd=da,calendar=LIS_calendar,rc=status)
              hr = 0 
              mn = 0 
              ss = 0 
           endif
           
           call LIS_create_restart_filename(n,filename,&
                'ROUTING','HYMAP2_router',&
                yr, mo, da, hr, mn, ss, & 
                wformat="netcdf")
           
           HYMAP2_routing_struc(n)%rstfile = filename
        endif
     endif
     
     if(read_restart) then 
        write(LIS_logunit,*) 'HYMAP2 restart file used: ', &
             trim(HYMAP2_routing_struc(n)%rstfile)

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
        status = nf90_open(path=HYMAP2_routing_struc(n)%rstfile, &
             mode=NF90_NOWRITE, ncid=ftn)
        call LIS_verify(status, "Error opening file "//&
             trim(HYMAP2_routing_struc(n)%rstfile))
#endif
        call HYMAP2_readvar_restart(ftn,n,&
             HYMAP2_routing_struc(n)%rivsto,&
             "RIVSTO")
        call HYMAP2_readvar_restart(ftn,n,&
             HYMAP2_routing_struc(n)%fldsto,&
             "FLDSTO")
        call HYMAP2_readvar_restart(ftn,n,&
             HYMAP2_routing_struc(n)%rnfsto,&
             "RNFSTO")
        call HYMAP2_readvar_restart(ftn,n,&
             HYMAP2_routing_struc(n)%bsfsto,&
             "BSFSTO")

        if(HYMAP2_routing_struc(n)%flowtype==3)then
           call HYMAP2_readvar_restart(ftn,n,&
                HYMAP2_routing_struc(n)%rivout_pre,&
                "RIVOUT_PRE")
           call HYMAP2_readvar_restart(ftn,n,&
                HYMAP2_routing_struc(n)%rivdph_pre,&
                "RIVDPH_PRE")
           call HYMAP2_readvar_restart(ftn,n,&
                HYMAP2_routing_struc(n)%fldout_pre,&
                "FLDOUT_PRE")
           call HYMAP2_readvar_restart(ftn,n,&
                HYMAP2_routing_struc(n)%flddph_pre,&
                "FLDDPH_PRE")
        endif
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
        status = nf90_close(ftn)
        call LIS_verify(status, "Error in nf90_close in HYMAP2_routing_readrst")
#endif
     endif
  enddo

end subroutine HYMAP2_routing_readrst

!BOP
! !ROUTINE: HYMAP2_readvar_restart
! \label{HYMAP2_readvar_restart}
! 
! !INTERFACE:
  subroutine HYMAP2_readvar_restart(ftn, n, var, varname)
! !USES:
    use LIS_coreMod
    use LIS_logMod
    use HYMAP2_routingMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

    implicit none
! !ARGUMENTS: 
    integer, intent(in)   :: ftn
    integer, intent(in)   :: n
    real, intent(inout)   :: var(HYMAP2_routing_struc(n)%nseqall)
    character(len=*)      :: varname
    
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

    allocate(gtmp(HYMAP2_routing_struc(n)%nseqall_glb))
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    status = nf90_inq_varid(ftn,trim(varname),varid)
    call LIS_verify(status,'Error in nf90_inq_varid in HYMAP2_readvar_restart')
    status = nf90_get_var(ftn,varid,gtmp)
    call LIS_verify(status,'Error in nf90_get_var in HYMAP2_readvar_restart')
#endif       

    do i=1,HYMAP2_routing_struc(n)%nseqall
       ix = HYMAP2_routing_struc(n)%seqx(i)
       iy = HYMAP2_routing_struc(n)%seqy(i)
       ix1 = ix + LIS_ews_halo_ind(n,LIS_localPet+1) -1
       iy1 = iy + LIS_nss_halo_ind(n,LIS_localPet+1) -1
       var(i)  = gtmp(HYMAP2_routing_struc(n)%sindex(ix1,iy1))
    enddo
     
    deallocate(gtmp)   
  end subroutine HYMAP2_readvar_restart
