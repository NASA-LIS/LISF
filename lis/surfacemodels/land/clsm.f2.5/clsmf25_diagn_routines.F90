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
! !MODULE: clsmf25_diagn_routines
! 
! !DESCRIPTION: 
! 
!  This module contains subroutinte for computing various diagnostics
!  in the Catchment model
! 
! !REVISION HISTORY: 
!  Feb 5,  2004, Rolf Reichle:  Initial Specification
!  Mar 19, 2004, Rolf Reichle:  Revised subroutine, calc_soil_moist
!  Aug 31, 2004, Rolf Reichle:  Added calc_tsurf
!  Jun 21, 2005, Rolf Reichle:  Added calc_arX
!
!EOP
module clsmf25_diagn_routines

  use clsmf25_model
  use clsmf25_types
  use clsmf25_constants  
  implicit none
  
contains

!BOP
! 
! !ROUTINE: calc_tsurf
! \label{calc_tsurf}
! 
! !INTERFACE:   
  subroutine calc_tsurf( N_cat, cat_param, cat_progn, tsurf )
!
! !DESCRIPTION:     
! Calculate diagnostic surface temperature "tsurf" from prognostics
!
!EOP
    
    implicit none
    
    integer,                  intent(in) :: N_cat
    
    type(cat_param_type), dimension(N_cat), intent(in)  :: cat_param
    type(cat_progn_type), dimension(N_cat), intent(in)  :: cat_progn
    
    real,                 dimension(N_cat), intent(out) :: tsurf
    
    ! ----------------------------
    !    
    ! local variables
    
    integer :: n
    
    real, dimension(N_cat) :: ar1, ar2, ar4
    
    real, dimension(N_cat) :: asnow 
    
    real, dimension(N_cat) :: tpsn1, real_dummy
    
    real, parameter :: tf = 273.16
    
    ! ------------------------------------------------------------------
    
    call calc_arX( N_cat, cat_param, cat_progn, ar1, ar2, ar4 )
    
    ! Compute snow covered area
    
    call calc_asnow( N_cat, cat_progn, asnow )
    
    ! Compute tsurf
    
    do n=1,N_cat
       
       tsurf(n) = ar1(n)*cat_progn(n)%tc1 + ar2(n)*cat_progn(n)%tc2 &
            + ar4(n)*cat_progn(n)%tc4
       
       if (asnow(n)>0.) then
          
          ! get_tf_nd returns snow temperature in deg Celsius
          
          call get_tf_nd( 1, cat_progn(n)%htsn(1), cat_progn(n)%wesn(1), &
               tpsn1(1), real_dummy(1) ) 
          
          tsurf(n) = (1. - asnow(n))*tsurf(n) + asnow(n)*(tpsn1(1) + tf)
          
       end if
       
    enddo
    
    return
    
  end subroutine calc_tsurf
  
!BOP
! 
! !ROUTINE: calc_arX
! \label{calc_arX}
! 
! !INTERFACE:   
  subroutine calc_arX( N_cat, cat_param, cat_progn, ar1, ar2, ar4 )
! 
! !DESCRIPTION: 
! 
!  ???
!EOP    
    implicit none
    
    integer, intent(in) :: N_cat
     
    type(cat_param_type), dimension(N_cat), intent(in) :: cat_param
    type(cat_progn_type), dimension(N_cat), intent(in) :: cat_progn
    
    real, dimension(N_cat), intent(out) :: ar1, ar2, ar4
    
    ! locals
    
    real, parameter :: dtstep_dummy = -9999.
    
    real, dimension(N_cat) :: rzeq, runsrf_dummy
    real, dimension(N_cat) :: catdef_dummy, rzexc_dummy, srfexc_dummy    
    real, dimension(N_cat) :: srfmx, srfmn, swsrf1, swsrf2, swsrf4, rzi

    ! ------------------------------------------------------------------
    !
    ! Call partition to get saturated/unsaturated/wilting 
    !  areas ar1, ar2, and ar4:
    !
    ! Need to calculate root zone equilibrium moisture for given 
    !  catchment deficit (needed for call to partition).
    
    call rzequil( &
         N_cat, cat_param%vegcls, cat_progn%catdef, cat_param%vgwmax,      &
         cat_param%cdcr1, cat_param%cdcr2, cat_param%wpwet,                &
         cat_param%ars1,  cat_param%ars2,  cat_param%ars3,                 &
         cat_param%ara1,  cat_param%ara2,  cat_param%ara3, cat_param%ara4, &
         cat_param%arw1,  cat_param%arw2,  cat_param%arw3, cat_param%arw4, &
         rzeq)
    
    ! Call partition with dtstep_dummy:
    !  In partition, dtstep is only used for a correction that
    !  puts water into runsrf (for which runsrf_dummy is used here).
    !  Also use catdef_dummy etc because partition() updates catdef
    !  whenever srfexc exceeds physical bounds, but this is not desired here.
    
    runsrf_dummy = 0.
    
    catdef_dummy = cat_progn%catdef          
    rzexc_dummy  = cat_progn%rzexc
    srfexc_dummy = cat_progn%srfexc
    
    call partition( &
         N_cat, dtstep_dummy, cat_param%vegcls, cat_param%dzsf, rzexc_dummy, &
         rzeq, cat_param%vgwmax, cat_param%cdcr1, cat_param%cdcr2,           &
         cat_param%psis,   cat_param%bee,  cat_param%poros, cat_param%wpwet, &
         cat_param%ars1,   cat_param%ars2, cat_param%ars3,                   &
         cat_param%ara1,   cat_param%ara2, cat_param%ara3,  cat_param%ara4,  &
         cat_param%arw1,   cat_param%arw2, cat_param%arw3,  cat_param%arw4,  &
         .false.,                                                            &
         srfexc_dummy, catdef_dummy, runsrf_dummy,                           &
         ar1, ar2, ar4, srfmx, srfmn, swsrf1, swsrf2, swsrf4,rzi )
    
  end subroutine calc_arX

!BOP
! 
! !ROUTINE: calc_asnow
! \label{calc_asnow}
! 
! !INTERFACE: 
  subroutine calc_asnow( N_cat, cat_progn, asnow )
!
! !DESCRIPTION:     
! Calculate diagnostic snow area from prognostic SWE
!
!EOP    
    implicit none
    
    integer,                  intent(in) :: N_cat
    
    type(cat_progn_type), dimension(N_cat), intent(in) :: cat_progn
    
    real, dimension(N_cat), intent(out) :: asnow
    
    ! local variables
    
    integer :: n

    ! "wemin" MUST BE CONSISTENT WITH "wemin" IN SUBROUTINE SNOWRT()!
    
!sm    real, parameter :: wemin = 13.    ! [kg/m^2]
    
    ! -----------------------------------------------------------

    do n=1,N_cat
       
       asnow(n) = min( sum(cat_progn(n)%wesn(1:N_snow))/wemin, 1.)
       
    end do
    
  end subroutine calc_asnow    

end module clsmf25_diagn_routines

! *************************** EOF ****************************************
