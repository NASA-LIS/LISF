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
! !ROUTINE: NoahMP401_read_OPT_parameters
! \label{NoahMP401_read_OPT_parameters}
!
! !REVISION HISTORY:
!
!  This subroutine reads the optimized parameters generated typically 
!  from the LIS OPT/UE system. In general, this routine can be used
!  to overwrite the default lookup table values that the model uses. 
! 
!   5 May 2020: Sujay Kumar; Initial specification
!
! !INTERFACE:
!
subroutine NoahMP401_read_OPT_parameters()
! !USES:
  use LIS_coreMod
  use NoahMP401_lsmMod

  implicit none

  integer           :: mtype
  logical           :: var_found
  integer           :: t, k, n
  integer           :: col, row
  real, allocatable :: placeholder(:,:)

  mtype = LIS_rc%lsm_index

  do n=1,LIS_rc%nnest
     allocate(placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n)))

     call NOAHMP401_read_OPT_param(n, "ALBDRY1",  placeholder, var_found)

     if(var_found) then 
        do t = 1, LIS_rc%npatch(n, mtype)
           col = LIS_surface(n, mtype)%tile(t)%col
           row = LIS_surface(n, mtype)%tile(t)%row
           if(placeholder(col,row).ne.LIS_rc%udef) then 
              NOAHMP401_struc(n)%noahmp401(t)%param%albdry(1) =&
                   placeholder(col, row)
           endif
        enddo
     endif
     
     call NOAHMP401_read_OPT_param(n, "ALBDRY2",  placeholder, var_found)
     if(var_found) then 
        do t = 1, LIS_rc%npatch(n, mtype)
           col = LIS_surface(n, mtype)%tile(t)%col
           row = LIS_surface(n, mtype)%tile(t)%row
           if(placeholder(col,row).ne.LIS_rc%udef) then 
              NOAHMP401_struc(n)%noahmp401(t)%param%albdry(2) =&
                   placeholder(col, row)
           endif
        enddo
     endif

     call NOAHMP401_read_OPT_param(n, "ALBICE1", placeholder, var_found)
     if(var_found) then 
        do t = 1, LIS_rc%npatch(n, mtype)
           col = LIS_surface(n, mtype)%tile(t)%col
           row = LIS_surface(n, mtype)%tile(t)%row
           if(placeholder(col,row).ne.LIS_rc%udef) then 
              NOAHMP401_struc(n)%noahmp401(t)%param%albice(1) =&
                   placeholder(col, row)
           endif
        enddo
     endif

     call NOAHMP401_read_OPT_param(n, "ALBICE2", placeholder, var_found)
     if(var_found) then 
        do t = 1, LIS_rc%npatch(n, mtype)
           col = LIS_surface(n, mtype)%tile(t)%col
           row = LIS_surface(n, mtype)%tile(t)%row
           if(placeholder(col,row).ne.LIS_rc%udef) then 
              NOAHMP401_struc(n)%noahmp401(t)%param%albice(2) =&
                   placeholder(col, row)
           endif
        enddo
     endif

     call NOAHMP401_read_OPT_param(n, "ALBSAT1", placeholder, var_found)
     if(var_found) then 
        do t = 1, LIS_rc%npatch(n, mtype)
           col = LIS_surface(n, mtype)%tile(t)%col
           row = LIS_surface(n, mtype)%tile(t)%row
           if(placeholder(col,row).ne.LIS_rc%udef) then 
              NOAHMP401_struc(n)%noahmp401(t)%param%albsat(1) =&
                   placeholder(col, row)
           endif
        enddo
     endif

     call NOAHMP401_read_OPT_param(n, "ALBSAT2", placeholder, var_found)
     if(var_found) then 
        do t = 1, LIS_rc%npatch(n, mtype)
           col = LIS_surface(n, mtype)%tile(t)%col
           row = LIS_surface(n, mtype)%tile(t)%row
           if(placeholder(col,row).ne.LIS_rc%udef) then 
              NOAHMP401_struc(n)%noahmp401(t)%param%albsat(2) =&
                   placeholder(col, row)
           endif
        enddo
     endif

     call NOAHMP401_read_OPT_param(n, "BETADS", placeholder, var_found)
     if(var_found) then 
        do t = 1, LIS_rc%npatch(n, mtype)
           col = LIS_surface(n, mtype)%tile(t)%col
           row = LIS_surface(n, mtype)%tile(t)%row
           if(placeholder(col,row).ne.LIS_rc%udef) then 
              NOAHMP401_struc(n)%noahmp401(t)%param%betads =&
                   placeholder(col, row)
           endif
        enddo
     endif
     
     call NOAHMP401_read_OPT_param(n, "BETAIS", placeholder, var_found)
     if(var_found) then 
        do t = 1, LIS_rc%npatch(n, mtype)
           col = LIS_surface(n, mtype)%tile(t)%col
           row = LIS_surface(n, mtype)%tile(t)%row
           if(placeholder(col,row).ne.LIS_rc%udef) then 
              NOAHMP401_struc(n)%noahmp401(t)%param%betais =&
                   placeholder(col, row)
           endif
        enddo
     endif

     call NOAHMP401_read_OPT_param(n, "EG1", placeholder, var_found)
     if(var_found) then 
        do t = 1, LIS_rc%npatch(n, mtype)
           col = LIS_surface(n, mtype)%tile(t)%col
           row = LIS_surface(n, mtype)%tile(t)%row
           if(placeholder(col,row).ne.LIS_rc%udef) then 
              NOAHMP401_struc(n)%noahmp401(t)%param%eg(1) =&
                   placeholder(col, row)
           endif
        enddo
     endif

     call NOAHMP401_read_OPT_param(n, "EG2", placeholder, var_found)
     if(var_found) then 
        do t = 1, LIS_rc%npatch(n, mtype)
           col = LIS_surface(n, mtype)%tile(t)%col
           row = LIS_surface(n, mtype)%tile(t)%row
           if(placeholder(col,row).ne.LIS_rc%udef) then 
              NOAHMP401_struc(n)%noahmp401(t)%param%eg(2) =&
                   placeholder(col, row)
           endif
        enddo
     endif

     call NOAHMP401_read_OPT_param(n, "MFSNO", placeholder, var_found)
     if(var_found) then 
        do t = 1, LIS_rc%npatch(n, mtype)
           col = LIS_surface(n, mtype)%tile(t)%col
           row = LIS_surface(n, mtype)%tile(t)%row
           if(placeholder(col,row).ne.LIS_rc%udef) then 
              NOAHMP401_struc(n)%noahmp401(t)%param%mfsno =&
                   placeholder(col, row)
           endif
        enddo
     endif

     call NOAHMP401_read_OPT_param(n, "OMEGAS1", placeholder, var_found)
     if(var_found) then 
        do t = 1, LIS_rc%npatch(n, mtype)
           col = LIS_surface(n, mtype)%tile(t)%col
           row = LIS_surface(n, mtype)%tile(t)%row
           if(placeholder(col,row).ne.LIS_rc%udef) then 
              NOAHMP401_struc(n)%noahmp401(t)%param%omegas(1) =&
                   placeholder(col, row)
           endif
        enddo
     endif

     call NOAHMP401_read_OPT_param(n, "OMGEAS2", placeholder, var_found)
     if(var_found) then 
        do t = 1, LIS_rc%npatch(n, mtype)
           col = LIS_surface(n, mtype)%tile(t)%col
           row = LIS_surface(n, mtype)%tile(t)%row
           if(placeholder(col,row).ne.LIS_rc%udef) then 
              NOAHMP401_struc(n)%noahmp401(t)%param%omegas(2) =&
                   placeholder(col, row)
           endif
        enddo
     endif

     call NOAHMP401_read_OPT_param(n, "RSURF_SNOW", placeholder, var_found)
     if(var_found) then 
        do t = 1, LIS_rc%npatch(n, mtype)
           col = LIS_surface(n, mtype)%tile(t)%col
           row = LIS_surface(n, mtype)%tile(t)%row
           if(placeholder(col,row).ne.LIS_rc%udef) then 
              NOAHMP401_struc(n)%noahmp401(t)%param%rsurf_snow =&
                   placeholder(col, row)
           endif
        enddo
     endif

     call NOAHMP401_read_OPT_param(n, "SSI", placeholder, var_found)
     if(var_found) then 
        do t = 1, LIS_rc%npatch(n, mtype)
           col = LIS_surface(n, mtype)%tile(t)%col
           row = LIS_surface(n, mtype)%tile(t)%row
           if(placeholder(col,row).ne.LIS_rc%udef) then 
              NOAHMP401_struc(n)%noahmp401(t)%param%ssi =&
                   placeholder(col, row)
           endif
        enddo
     endif

     call NOAHMP401_read_OPT_param(n, "Z0SNO", placeholder, var_found)
     if(var_found) then 
        do t = 1, LIS_rc%npatch(n, mtype)
           col = LIS_surface(n, mtype)%tile(t)%col
           row = LIS_surface(n, mtype)%tile(t)%row
           if(placeholder(col,row).ne.LIS_rc%udef) then 
              NOAHMP401_struc(n)%noahmp401(t)%param%z0sno =&
                   placeholder(col, row)
           endif
        enddo
     endif

     call NOAHMP401_read_OPT_param(n, "MXSNALB", placeholder, var_found)
     if(var_found) then 
        do t = 1, LIS_rc%npatch(n, mtype)
           col = LIS_surface(n, mtype)%tile(t)%col
           row = LIS_surface(n, mtype)%tile(t)%row
           if(placeholder(col,row).ne.LIS_rc%udef) then 
              NOAHMP401_struc(n)%noahmp401(t)%param%mxsnalb =&
                   placeholder(col, row)
           endif
        enddo
     endif


     call NOAHMP401_read_OPT_param(n, "MNSNALB", placeholder, var_found)
     if(var_found) then 
        do t = 1, LIS_rc%npatch(n, mtype)
           col = LIS_surface(n, mtype)%tile(t)%col
           row = LIS_surface(n, mtype)%tile(t)%row
           if(placeholder(col,row).ne.LIS_rc%udef) then 
              NOAHMP401_struc(n)%noahmp401(t)%param%mnsnalb =&
                   placeholder(col, row)
           endif
        enddo
     endif


     call NOAHMP401_read_OPT_param(n, "SNDECAYEXP", placeholder, var_found)
     if(var_found) then 
        do t = 1, LIS_rc%npatch(n, mtype)
           col = LIS_surface(n, mtype)%tile(t)%col
           row = LIS_surface(n, mtype)%tile(t)%row
           if(placeholder(col,row).ne.LIS_rc%udef) then 
              NOAHMP401_struc(n)%noahmp401(t)%param%sndecayexp =&
                   placeholder(col, row)
           endif
        enddo
     endif

     call NOAHMP401_read_OPT_param(n, "T_ULIMIT", placeholder, var_found)
     if(var_found) then 
        do t = 1, LIS_rc%npatch(n, mtype)
           col = LIS_surface(n, mtype)%tile(t)%col
           row = LIS_surface(n, mtype)%tile(t)%row
           if(placeholder(col,row).ne.LIS_rc%udef) then 
              NOAHMP401_struc(n)%noahmp401(t)%param%t_ulimit =&
                   placeholder(col, row)
           endif
        enddo
     endif

     call NOAHMP401_read_OPT_param(n, "T_MLIMIT", placeholder, var_found)
     if(var_found) then 
        do t = 1, LIS_rc%npatch(n, mtype)
           col = LIS_surface(n, mtype)%tile(t)%col
           row = LIS_surface(n, mtype)%tile(t)%row
           if(placeholder(col,row).ne.LIS_rc%udef) then 
              NOAHMP401_struc(n)%noahmp401(t)%param%t_mlimit =&
                   placeholder(col, row)
           endif
        enddo
     endif


     call NOAHMP401_read_OPT_param(n, "T_LLIMIT", placeholder, var_found)
     if(var_found) then 
        do t = 1, LIS_rc%npatch(n, mtype)
           col = LIS_surface(n, mtype)%tile(t)%col
           row = LIS_surface(n, mtype)%tile(t)%row
           if(placeholder(col,row).ne.LIS_rc%udef) then 
              NOAHMP401_struc(n)%noahmp401(t)%param%t_llimit =&
                   placeholder(col, row)
           endif
        enddo
     endif

     call NOAHMP401_read_OPT_param(n, "SNOWF_SCALEF", placeholder, var_found)
     if(var_found) then 
        do t = 1, LIS_rc%npatch(n, mtype)
           col = LIS_surface(n, mtype)%tile(t)%col
           row = LIS_surface(n, mtype)%tile(t)%row
           if(placeholder(col,row).ne.LIS_rc%udef) then 
              NOAHMP401_struc(n)%noahmp401(t)%param%snowf_scalef =&
                   placeholder(col, row)
           endif
        enddo
     endif

     deallocate(placeholder)
    
  end do
end subroutine NoahMP401_read_OPT_parameters

!BOP
!
! !ROUTINE: NOAHMP401_read_OPT_param
!  \label{NOAHMP401_read_OPT_param}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification for read_laiclimo
!  30 Oct  2013: Shugong Wang; Generalization for reading OPT spatial parameter
!
! !INTERFACE:
subroutine NOAHMP401_read_OPT_param(n, ncvar_name, placeholder,var_found)
! !USES:
  use netcdf
  use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_localPet,   &   
       LIS_ews_halo_ind, LIS_ewe_halo_ind, &
       LIS_nss_halo_ind, LIS_nse_halo_ind   
  use LIS_logMod,  only : LIS_logunit, LIS_verify, LIS_endrun
  use LIS_fileIOMod, only: LIS_read_param
  implicit none
! !ARGUMENTS: 
  integer, intent(in)          :: n
  character(len=*), intent(in) :: ncvar_name 
  real, intent(out)            :: placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n))
  logical                      :: var_found
! !DESCRIPTION:
!  This subroutine reads OPT parameters from the LIS
!  NetCDF parameter data file
!  
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[array]
!    array containing returned values
!   \end{description}
!
!EOP      

  integer       :: ios1
  integer       :: ios, nid, param_ID, nc_ID, nr_ID, dimids(3)
  integer       :: nc, nr, t,  k
  real, pointer :: level_data(:, :)
  logical       :: file_exists
  
  placeholder = LIS_rc%udef
  var_found  = .false. 

  inquire(file=LIS_rc%paramfile(n), exist=file_exists)
  if(file_exists) then
     write(LIS_logunit, *) '[INFO] Reading '//trim(ncvar_name)//&
          ' map '
     
     ! open NetCDF parameter file
     ios = nf90_open(path=trim(LIS_rc%paramfile(n)), &
          mode=NF90_NOWRITE, ncid=nid)
     call LIS_verify(ios, 'Error in nf90_open in NOAHMP401_read_OPT_param')
     
     ! inquire the ID of east-west dimension
     ios = nf90_inq_dimid(nid, 'east_west', nc_ID)
     call LIS_verify(ios, &
          'Error in nf90_inq_dimid in NOAHMP401_read_OPT_param')
     
     ! inquire the ID of north-south dimension
     ios = nf90_inq_dimid(nid, 'north_south', nr_ID)
     call LIS_verify(ios, &
          'Error in nf90_inq_dimid in NOAHMP401_read_OPT_param')
     
     ! inquire the length of east-west dimension
     ios = nf90_inquire_dimension(nid, nc_ID, len=nc)
     call LIS_verify(ios, &
          'Error in nf90_inquire_dimension in NOAHMP401_read_OPT_param')
     
     ! inquire the length of north-south dimension
     ios = nf90_inquire_dimension(nid, nr_ID, len=nr)
     call LIS_verify(ios, &
          'Error in nf90_inquire_dimension in NOAHMP401_read_OPT_param')
     
     ! inquire the ID of parameter. 
     ios = nf90_inq_varid(nid, Trim(ncvar_name), param_ID)
     if(ios.eq.0) then 
        ! inquire the IDs of all dimensions. The third dimension is the level dimension
        ios = nf90_inquire_variable(nid, param_ID, dimids = dimids)
        call LIS_verify(ios, trim(ncvar_name)//' failed to inquire dimensions')
        
        ! allocate memory
        allocate(level_data (LIS_rc%gnc(n), LIS_rc%gnr(n)))
        
        ! inquire the variable ID of parameter 
        ios = nf90_inq_varid(nid, trim(ncvar_name), param_ID)
        call LIS_verify(ios, trim(ncvar_name)//&
             ' field not found in the LIS param file')
        
        ! read parameter 
        ios = nf90_get_var(nid, param_ID, level_data)
        call LIS_verify(ios, 'Error in nf90_get_var in NOAHMP401_read_OPT_param')
        
        ! grab parameter at specific level
        placeholder(:, :) = & 
             level_data(LIS_ews_halo_ind(n, LIS_localPet+1):LIS_ewe_halo_ind(n, LIS_localPet+1), &
             LIS_nss_halo_ind(n, LIS_localPet+1):LIS_nse_halo_ind(n, LIS_localPet+1))
        
        deallocate(level_data)

        var_found = .true. 
     else
        write(LIS_logunit,*) '[WARN] ', trim(ncvar_name)//&
             ' field not found in the LIS param file'
     endif

        ! close netcdf file 
     ios = nf90_close(nid)
     call LIS_verify(ios, 'Error in nf90_close in NOAHMP401_read_OPT_param')
     
  endif
end subroutine NOAHMP401_read_OPT_param
                                         
