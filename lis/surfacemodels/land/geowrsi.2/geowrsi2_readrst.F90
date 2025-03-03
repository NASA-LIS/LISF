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
!
! !ROUTINE: geowrsi2_readrst
! \label{geowrsi2_readrst}
! 
! !REVISION HISTORY:
!  30 Jan 2014: KR Arsenault; Initial code set-up
!
! !INTERFACE:
subroutine geowrsi2_readrst()

! !USES:
  use LIS_coreMod,    only : LIS_rc, LIS_masterproc
  use LIS_historyMod, only : LIS_readvar_restart
  use LIS_logMod,     only : LIS_logunit, LIS_endrun, LIS_verify,&
          LIS_getNextUnitNumber, LIS_releaseUnitNumber
  use geowrsi2_lsmMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

! !ARGUMENTS: 
!
! !DESCRIPTION:
!  This program reads restart and start-of-season (SOS) files for 
!  GeoWRSI 2.0, including important tile information.
!
!  The following is the list of variables specified in the GeoWRSI2.0
!  restart file: 
!  \begin{verbatim}
!   nc,nr,ntiles - grid and tile space dimensions 
!   sos          - GeoWRSI2 SOS 
!   sosanom      - GeoWRSI2 SOSanom
!  \end{verbatim}
!
!  The routines invoked are: 
! \begin{description}
! \item[LIS\_readvar\_restart](\ref{LIS_readvar_restart}) \newline
!   Reads a variable from the restart file
! \end{description}
!
!EOP
  implicit none
  integer           :: n, t, l
  integer           :: nc,nr,npatch
  integer           :: ftn
  integer           :: file_status
  logical           :: file_exists
  character*20      :: wformat
  real, allocatable :: tmptilen(:)
  integer           :: year_temp

! __________________________________


  do n=1,LIS_rc%nnest   ! Loop over nests

  !- Read-in ensemble-based history files (temporary) of netcdf-format for SOS-input:
     if( geowrsi2_lsmRunMode == "WRSI" .and. LIS_rc%nensem(n) > 1  ) then

        wformat="netcdf"   ! Hard-coded to netcdf for now

        write(LIS_logunit,*) "... READING IN ENSEMBLE SOS/SOSa NETCDF (Restart) FILES "
        if( LIS_rc%startcode .ne. "restart" ) then
           write(LIS_logunit,*) "ERR: Need to select 'restart' if reading in Ensemble-based file ..."
           write(LIS_logunit,*) "Program stopping ..."
           call LIS_endrun()
        endif
       
        allocate(tmptilen(LIS_rc%npatch(n,LIS_rc%lsm_index)))
        inquire(file=trim(geowrsi2_struc(n)%rfile), exist=file_exists)

        if(.not.file_exists) then
           write(LIS_logunit,*) "GeoWRSI2.0 restart file ", &
                geowrsi2_struc(n)%rfile," does not exist."
           write(LIS_logunit,*) "Program stopping ..."
           call LIS_endrun()

        endif
        write(LIS_logunit,*)    &
            "GeoWRSI2.0 restart file used: ",geowrsi2_struc(n)%rfile

      ! Open netcdf file
        if( wformat == "netcdf" ) then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
          file_status=nf90_open(path=geowrsi2_struc(n)%rfile,&
              mode=NF90_NOWRITE,ncid=ftn)
          call LIS_verify(file_status,&
             "Error opening file "//geowrsi2_struc(n)%rfile)
#endif
        endif

      ! Read SOS-fields from netcdf file:
        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,geowrsi2_struc(n)%wrsi%SOS, &
             varname="SOS_inst",wformat=wformat)

        call LIS_readvar_restart(ftn,n,LIS_rc%lsm_index,geowrsi2_struc(n)%wrsi%SOSa, &
             varname="SOSa_inst",wformat=wformat)

!           geowrsi2_struc(n)%wrsi(t)%sos_write = geowrsi2_struc(n)%wrsi(t)%SOS
!           geowrsi2_struc(n)%wrsi(t)%sosa_write= geowrsi2_struc(n)%wrsi(t)%SOSa

!           do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
!              geowrsi2_struc(n)%wrsi(t)%SOS  = tmptilen(t)
!              geowrsi2_struc(n)%wrsi(t)%SOSa = tmptilen(t)
!           enddo

      ! Close netcdf file
        if( wformat == "netcdf" ) then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
           file_status = nf90_close(ftn)
           call LIS_verify(file_status,'Error in nf90_close in geowrsi2_readrst')
#endif     
        endif
        deallocate(tmptilen)

     endif  ! End ensemble restart read

  end do    ! End nest loop

end subroutine geowrsi2_readrst
