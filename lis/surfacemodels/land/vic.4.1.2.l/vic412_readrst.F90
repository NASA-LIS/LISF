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
! !ROUTINE: vic412_readrst
! \label{vic412_readrst}
!
! !REVISION HISTORY:
! 02 Aug 2011; James Geiger, Initial implementation of VIC 4.1.1 into LIS.
! 15 May 2014; Shugong Wang, Add NetCDF restart support to VIC 4.1.2.l for LIS 7
! 
! !INTERFACE:
subroutine vic412_readrst()

! !USES:
    use LIS_coreMod, only    : LIS_rc, LIS_masterproc, LIS_domain
    use LIS_historyMod, only : LIS_readvar_restart
    use LIS_logMod, only     : LIS_logunit, LIS_endrun, &
                               LIS_getNextUnitNumber,   &
                               LIS_releaseUnitNumber,   &
                               LIS_verify                
    use vic412_lsmMod
    use netcdf
  implicit none
! !ARGUMENTS: 
!
! !DESCRIPTION:
!  This program reads restart files for VIC.  This
!  includes all relevant water/energy storages and tile information. 
!EOP

    integer           :: t, l
    integer           :: nc, nr, npatch
    integer           :: n
    integer           :: ftn
    integer           :: status
    real, allocatable :: tmptilen(:)
    logical           :: file_exists
    character*20      :: wformat
    integer           :: ts, nest, vt_scheme, vegclass 
    integer           :: start_mode
    real              :: start_t, end_t

    do n=1, LIS_rc%nnest
        !wformat = "netcdf"
        !wformat = "binary"
        wformat = trim(vic412_struc(n)%rfile_format)
        if(LIS_rc%startcode .eq. "coldstart") then  
            start_mode = 0
        ! restart
        elseif(LIS_rc%startcode .eq. "restart") then
            start_mode = 1 
            allocate(tmptilen(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            ! check the existance of restart file
            inquire(file=vic412_struc(n)%rfile, exist=file_exists)
            If (.not. file_exists) then 
                write(LIS_logunit,*) "vic412 restart file ", vic412_struc(n)%rfile," does not exist "
                write(LIS_logunit,*) "Program stopping ..."
                call LIS_endrun
            endif
            write(LIS_logunit,*) "vic412 restart file used: ", vic412_struc(n)%rfile
  
            call cpu_time(start_t)
            if(trim(wformat) .eq. "binary") then 
              ftn = LIS_getNextUnitNumber()
              open(ftn,file=vic412_struc(n)%rfile,form='unformatted')        
              read(ftn) nc,nr,npatch  !time, veg class, no. tiles
           
              if(nc.ne.LIS_rc%gnc(n) .or. nr.ne.LIS_rc%gnr(n))then
                  write(LIS_logunit,*) vic412_struc(n)%rfile,                 &
                   'grid space mismatch - VIC 4.1.2 halted'
                call LIS_endrun
              endif
              
              if(npatch .ne. LIS_rc%glbnpatch_red(n,LIS_rc%lsm_index))then           
                write(LIS_logunit,*) 'restart tile space mismatch, halting..'
                call LIS_endrun
              endif
            elseif(wformat.eq."netcdf") then 
              ! open netcdf restart file 
              status = nf90_open(path=vic412_struc(n)%rfile, &
                                 mode=NF90_NOWRITE, ncid=ftn)
              call LIS_verify(status, "Error opening file "//vic412_struc(n)%rfile)
            endif

            do l=1, vic412_struc(n)%state_chunk_size ! 
                 call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index,                 &
                                          tmptilen, varname="state_chunk", dim=l,   &
                                          vlevels=vic412_struc(n)%state_chunk_size, &
                                          wformat=wformat)
                 do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                     vic412_struc(n)%vic(t)%state_chunk(l) = tmptilen(t)
                 enddo
            enddo
            
            if(trim(wformat) .eq. "binary") then 
              close(ftn)
            else
              ! close netcdf restart file
              status = nf90_close(ftn)
              call LIS_verify(status, "Error in nf90_close in vic412_readrst")
            endif
            call cpu_time(end_t)
            write(LIS_logunit,*) "LIC VIC restart reading time: ", end_t - start_t, "(s)" 
            deallocate(tmptilen)
        endif
        
        ! distribute state chunck into state variables 
        vt_scheme = vic412_struc(n)%veg_tiling_scheme
        ts = 1 
        nest = n
        do t = 1, LIS_rc%npatch(n,LIS_rc%lsm_index)
           if (vt_scheme == 1) then ! LIS-based tiling
               vegclass = LIS_domain(n)%tile(t)%vegt
           else                     ! VIC-based tiling
               vegclass = -1
           endif
           call vic412_initialize(t, ts, vt_scheme, vegclass, nest, &
                                  vic412_struc(n)%vic(t)%state_chunk, &
                                  start_mode)
        enddo
    enddo
end subroutine vic412_readrst
