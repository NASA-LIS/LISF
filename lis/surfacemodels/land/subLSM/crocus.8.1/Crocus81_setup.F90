!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: Crocus81_setup
! \label{Crocus81_setup}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   11/6/19: Mahdi Navari, Shugong Wang; initial implementation for LIS 7 and Crocus81
!
! !INTERFACE:
subroutine Crocus81_setup()
! !USES:
    use LIS_logMod,    only: LIS_logunit, LIS_verify, LIS_endrun
    use LIS_fileIOMod, only: LIS_read_param !, LIS_convertParamDataToLocalDomain
    use LIS_coreMod,   only: LIS_rc, LIS_surface
    use Crocus81_lsmMod

!
! !DESCRIPTION:
!
!  This routine is the entry point to set up the parameters
!  required for Crocus81.  These include: 
!    SLOPE        - angle between the normal to the surface and the vertical 
!    (MN:  replaced PDIRCOSZW with slope and computed the cosine in the driver) [Radians]
!    PERMSNOWFRAC - Fraction of permanet snow/ice [-]
!    SLOPE_DIR    - Typical slope aspect in the grid  (clockwise from N) [Radians]
!    SAND         - Soil sand fraction (-) [-]
!    SILT         - Soil silt fraction (-) [-]
!    CLAY         - Soil clay fraction (-) [-]
!    POROSITY     - Soil porosity (m3 m-3) [m3/m3]
! 
!  The routines invoked are:
!  \begin{description}
!  \item[LIS\_read\_param](\ref{LIS_read_param}) \\ 
!    retrieves LIS parameter data from NetCDF file
!  \end{description}
!EOP
    implicit none
    integer           :: mtype
    integer           :: t,k, n
    integer           :: col, row
    real, allocatable :: placeholder(:,:)
    
    mtype = LIS_rc%lsm_index
    
    do n=1, LIS_rc%nnest
        ! allocate memory for place holder for #n nest
        allocate(placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n)))
        
        !------------------------------------!
        ! reading spatial spatial parameters !
        !------------------------------------!


        ! read: SLOPE
        write(LIS_logunit,*) "Crocus81: reading parameter SLOPE from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(CROCUS81_struc(n)%LDT_ncvar_SLOPE), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            CROCUS81_struc(n)%crocus81(t)%slope = placeholder(col, row)
        enddo 

        ! read: PERMSNOWFRAC
        write(LIS_logunit,*) "Crocus81: reading parameter PERMSNOWFRAC from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(CROCUS81_struc(n)%LDT_ncvar_PERMSNOWFRAC), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            CROCUS81_struc(n)%crocus81(t)%permsnowfrac = placeholder(col, row)
        enddo

        ! read: SLOPE_DIR
        write(LIS_logunit,*) "Crocus81: reading parameter SLOPE_DIR from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(CROCUS81_struc(n)%LDT_ncvar_SLOPE_DIR), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            CROCUS81_struc(n)%crocus81(t)%slope_dir = placeholder(col, row)
        enddo 

        ! read: SAND
        write(LIS_logunit,*) "Crocus81: reading parameter SAND from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(CROCUS81_struc(n)%LDT_ncvar_SAND), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            CROCUS81_struc(n)%crocus81(t)%sand = placeholder(col, row)
        enddo 

        ! read: SILT
        write(LIS_logunit,*) "Crocus81: reading parameter SILT from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(CROCUS81_struc(n)%LDT_ncvar_SILT), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            CROCUS81_struc(n)%crocus81(t)%silt = placeholder(col, row)
        enddo 

        ! read: CLAY
        write(LIS_logunit,*) "Crocus81: reading parameter CLAY from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(CROCUS81_struc(n)%LDT_ncvar_CLAY), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            CROCUS81_struc(n)%crocus81(t)%clay = placeholder(col, row)
        enddo 

        ! read: POROSITY
        write(LIS_logunit,*) "Crocus81: reading parameter POROSITY from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(CROCUS81_struc(n)%LDT_ncvar_POROSITY), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            CROCUS81_struc(n)%crocus81(t)%porosity = placeholder(col, row)
        enddo 

        !----------------------------------------------!
        ! MULTILEVEL reading spatial spatial parameters !
        !----------------------------------------------!
          ! read: ALB
          write(LIS_logunit,*) "Crocus81: reading parameter ALB from ", trim(LIS_rc%paramfile(n))
          do k = 1, 12
              call CROCUS81_read_MULTILEVEL_param(n, CROCUS81_struc(n)%LDT_ncvar_ALB, k, placeholder)
              do t = 1, LIS_rc%npatch(n, mtype)
                  col = LIS_surface(n, mtype)%tile(t)%col
                  row = LIS_surface(n, mtype)%tile(t)%row
                  CROCUS81_struc(n)%crocus81(t)%alb(k) = placeholder(col, row)
              enddo 
          enddo 
          
        deallocate(placeholder)
    enddo


end subroutine Crocus81_setup

!BOP
  !
  ! !ROUTINE: CROCUS81_read_MULTILEVEL_param
  !  \label{read_MULTILEVEL_param}
  !
  ! !REVISION HISTORY:
  !  03 Sept 2004: Sujay Kumar; Initial Specification for read_laiclimo
  !  30 Oct  2013: Shugong Wang; Generalization for reading MULTILEVEL spatial parameter
  !
  ! !INTERFACE:
  subroutine CROCUS81_read_MULTILEVEL_param(n, ncvar_name, level, placeholder)
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
    integer, intent(in)          :: level
    character(len=*), intent(in) :: ncvar_name
    real, intent(out)            :: placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n))
! !DESCRIPTION:
!  This subroutine reads MULTILEVEL parameters from the LIS
!  NetCDF parameter data file
!  
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[level]
!    level index (month, quarter, soil layer, snow layer) of the data to be read
!   \item[array]
!    array containing returned values
!   \end{description}
!
!EOP      

    integer       :: ios1
    integer       :: ios, nid, param_ID, nc_ID, nr_ID, dimids(3)
    integer       :: nc, nr, t, nlevel, k
    real, pointer :: level_data(:, :, :)
    real, pointer :: level_data1(:, :, :)
    logical       :: file_exists

    inquire(file=LIS_rc%paramfile(n), exist=file_exists)
    if(file_exists) then
        write(LIS_logunit, *) 'Reading '//trim(ncvar_name)//' map for level ', level

        ! open NetCDF parameter file
        ios = nf90_open(path=trim(LIS_rc%paramfile(n)), mode=NF90_NOWRITE, ncid=nid)
        call LIS_verify(ios, 'Error in nf90_open in RUC37_read_MULTILEVEL_param')

        ! inquire the ID of east-west dimension
        ios = nf90_inq_dimid(nid, 'east_west', nc_ID)
        call LIS_verify(ios, 'Error in nf90_inq_dimid in RUC37_read_MULTILEVEL_param')

        ! inquire the ID of north-south dimension
        ios = nf90_inq_dimid(nid, 'north_south', nr_ID)
        call LIS_verify(ios, 'Error in nf90_inq_dimid in RUC37_read_MULTILEVEL_param')

        ! inquire the length of east-west dimension
        ios = nf90_inquire_dimension(nid, nc_ID, len=nc)
        call LIS_verify(ios, 'Error in nf90_inquire_dimension in RUC37_read_MULTILEVEL_param')

        ! inquire the length of north-south dimension
        ios = nf90_inquire_dimension(nid, nr_ID, len=nr)
        call LIS_verify(ios, 'Error in nf90_inquire_dimension in RUC37_read_MULTILEVEL_param')

        ! inquire the ID of parameter. 
        ios = nf90_inq_varid(nid, Trim(ncvar_name), param_ID)
        call LIS_verify(ios, trim(ncvar_name)//' field not found in the LIS param file')

        ! inquire the IDs of all dimensions. The third dimension is the level dimension
        ios = nf90_inquire_variable(nid, param_ID, dimids = dimids)
        call LIS_verify(ios, trim(ncvar_name)//' failed to inquire dimensions')

        ! inquire the length of the level dimension
        ios = nf90_inquire_dimension(nid, dimids(3), len=nlevel)
        call LIS_verify(ios, trim(ncvar_name)//' failed to inquire the length of the 3rd dimension')

        ! allocate memory
        allocate(level_data(LIS_rc%gnc(n), LIS_rc%gnr(n), nlevel))

        ! inquire the variable ID of parameter 
        ios = nf90_inq_varid(nid, trim(ncvar_name), param_ID)
        call LIS_verify(ios, trim(ncvar_name)//' field not found in the LIS param file')

        ! read parameter 
        ios = nf90_get_var(nid, param_ID, level_data)
        call LIS_verify(ios, 'Error in nf90_get_var in RUC37_read_MULTILEVEL_param')

        ! close netcdf file 
        ios = nf90_close(nid)
        call LIS_verify(ios, 'Error in nf90_close in RUC37_read_MULTILEVEL_param')

        ! grab parameter at specific level
        placeholder(:, :) = &
             level_data(LIS_ews_halo_ind(n, LIS_localPet+1):LIS_ewe_halo_ind(n, LIS_localPet+1), &
                        LIS_nss_halo_ind(n, LIS_localPet+1):LIS_nse_halo_ind(n, LIS_localPet+1), level)

        ! free memory 
        deallocate(level_data)

    else
        write(LIS_logunit, *) 'MULTILEVEL parameter data file: ', LIS_rc%paramfile(n), ' does not exist'
        write(LIS_logunit, *) 'program stopping ...'
        call LIS_endrun
    endif
  end subroutine CROCUS81_read_MULTILEVEL_param

