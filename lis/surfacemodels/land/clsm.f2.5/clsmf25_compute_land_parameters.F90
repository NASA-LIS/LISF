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
! !ROUTINE: clsmf25_read_land_parameters3
! \label{clsmf25_read_land_parameters3}
!
! !REVISION HISTORY:
! 12 May 2003 Rolf Reichle, Initial Specification
! 06 Jun 2005 Rolf Reichle, adapted to read "SiB2_V2" parameters  
! 10 Jul 2006 James Geiger, Implementation in LIS
! 23 Nov 2012: David Mocko, Added Catchment Fortuna-2.5
!
! !INTERFACE:
subroutine clsmf25_compute_land_parameters(nest) 
! !USES:
    !use tile_coord_types
    use clsmf25_constants, only : N_gt
    use clsmf25_lsmMod
    use LIS_coreMod, only : LIS_rc, LIS_surface
    use LIS_logMod,  only : LIS_logunit, LIS_endrun

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: nest
!
! !DESCRIPTION:
! Reads in the vegetation, soil properties and topographic 
! parameters from global files and outputs parameters for the domain.
!
! Additional parameters are derived from the ones that have been read 
! from files.
!
!  The arguments are: 
!  \begin{description}
!   \item [nest]
!      index of the nest
!  \end{description}
!
! -------------------------------------------------------------
! 
!EOP
    integer :: d2g
    integer :: dummy_int, k, month, m
    real :: z_in_m, term1, term2
    character(300) :: filename
    
    write(LIS_logunit,*)"[INFO] Computing derived CLSM land surface parameters ..."
    write(LIS_logunit,*)
    
    do k=1,LIS_rc%npatch(nest,LIS_rc%lsm_index)
       
       ! Three soil depths for soil moisture model: 
       !
       ! dzsf: surface layer
       ! dzrz: root zone        -> water capacity of the root zone
       ! dzpr: unsaturated zone -> approx depth-to-bedrock
       !
       ! NOTE: Units of dz** are [mm] while excess/deficits from catchment()
       !       are in SI units (ie kg/m^2) or loosely speaking, in mm of water.
       !       In other words, density of water (1000 kg/m^3) is built 
       !       into dz** (reichle, 5 Feb 04).

       clsmf25_struc(nest)%cat_param(k)%dzsf = &
            clsmf25_struc(nest)%dzsfcrd * 1000.
       clsmf25_struc(nest)%cat_param(k)%dzrz = 1000.
       
       ! changed re-setting of dzrz back to earlier value because
       ! Sarith parameters are in fact consistent that the earlier version
       ! reichle, 12 Sep 2007
       !
       ! cp(k)%dzpr = max(1500., cp(k)%dpth)
       !
       ! previously, root zone depth ranged from .75m to 1m, which
       ! is inconsistent with subroutine catchment(), where root
       ! zone depth is hard-wired to 1m, and with the time scale
       ! parameters, that have been derived for 1m root zone depth
       ! (THE LATTER IS IN FACT *NOT* TRUE - reichle, 12 Sep 2007)
       ! - reichle, 30 May 2003
       
       clsmf25_struc(nest)%cat_param(k)%dzpr = &
            max(1000., clsmf25_struc(nest)%cat_param(k)%dpth)

       if ( clsmf25_struc(nest)%cat_param(k)%dzrz > &
            0.75*clsmf25_struc(nest)%cat_param(k)%dzpr) &
            clsmf25_struc(nest)%cat_param(k)%dzrz = &
            0.75*clsmf25_struc(nest)%cat_param(k)%dzpr
       
       ! soil storages

       clsmf25_struc(nest)%cat_param(k)%vgwmax =                           &
            clsmf25_struc(nest)%cat_param(k)%poros*&
            clsmf25_struc(nest)%cat_param(k)%dzrz
       
       z_in_m = clsmf25_struc(nest)%cat_param(k)%dzpr/1000.
       
       term1 = -1.+((clsmf25_struc(nest)%cat_param(k)%psis-z_in_m)/        &
                     clsmf25_struc(nest)%cat_param(k)%psis)**              &
                   ((clsmf25_struc(nest)%cat_param(k)%bee-1.)/             &
                     clsmf25_struc(nest)%cat_param(k)%bee)
       
       term2 = clsmf25_struc(nest)%cat_param(k)%psis *                     &
               clsmf25_struc(nest)%cat_param(k)%bee /                      &
              (clsmf25_struc(nest)%cat_param(k)%bee-1)
       
       clsmf25_struc(nest)%cat_param(k)%cdcr1 = 1000.*&
            clsmf25_struc(nest)%cat_param(k)%poros*(z_in_m-(-term2*term1))
       
       clsmf25_struc(nest)%cat_param(k)%cdcr2 =        &
          (1.-clsmf25_struc(nest)%cat_param(k)%wpwet)* &
          clsmf25_struc(nest)%cat_param(k)%poros*&
          clsmf25_struc(nest)%cat_param(k)%dzpr

       ! Soil depths for ground temperature model
       
       if (N_gt/=6) then
          
          write(LIS_logunit,*) '[INFO] clsmf25_read_land_parameters:'
          write(LIS_logunit,*) '   using N_gt = ',N_gt,', but only'
          write(LIS_logunit,*) '   6 layer depths are specified.'
          call LIS_endrun()
                    
       end if
       
       clsmf25_struc(nest)%cat_param(k)%dzgt(1) =  0.0988  
       clsmf25_struc(nest)%cat_param(k)%dzgt(2) =  0.1952
       clsmf25_struc(nest)%cat_param(k)%dzgt(3) =  0.3859
       clsmf25_struc(nest)%cat_param(k)%dzgt(4) =  0.7626
       clsmf25_struc(nest)%cat_param(k)%dzgt(5) =  1.5071
       clsmf25_struc(nest)%cat_param(k)%dzgt(6) = 10.0000
       
    end do
    
  end subroutine clsmf25_compute_land_parameters
