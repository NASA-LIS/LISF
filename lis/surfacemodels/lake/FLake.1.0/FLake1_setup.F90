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
! !ROUTINE: FLake1_setup
! \label{FLake1_setup}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   6/4/13: Shugong Wang; initial implementation for LIS 7 and FLake1
!
! !INTERFACE:
subroutine FLake1_setup()
! !USES:
    use LIS_logMod,    only: LIS_logunit
    use LIS_fileIOMod, only: LIS_read_param
    use LIS_coreMod,   only: LIS_rc, LIS_surface
    use FLake1_Mod

!
! !DESCRIPTION:
!
!  This routine is the entry point to set up the parameters
!  required for FLake1.  These include: 
!    lon          - longitude of lake center [-]
!    lat          - latitude of lake center [-]
!    depth_w      - lake depth [m]
!    fetch        - typical wind fetch [m]
!    depth_bs     - depth of the thermally active layer of the bottom sediments [m]
!    T_bs         - temperature at the outer edge of the thermally active layer of the bottom sediments [K]
! 
!  The routines invoked are:
!  \begin{description}
!  \item[LIS\_read\_param](\ref{LIS_read_param}) \newline
!    retrieves LIS parameter data fron NetCDF file
!  \end{description}
!EOP
    implicit none
    integer           :: mtype
    integer           :: i, n
    integer           :: col, row
    real, allocatable :: placeholder(:,:)
    
    mtype = LIS_rc%lake_index
    
    ! allocate memory for place holder
    do n=1, LIS_rc%nnest
       allocate(placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n)))
    
        ! read: lon
        write(LIS_logunit,*) "FLake1: reading parameter LON from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, "lon", placeholder)
        do i = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(i)%col
            row = LIS_surface(n, mtype)%tile(i)%row
            FLAKE1_struc(n)%flake1(i)%lon = placeholder(col, row)
        enddo 

        ! read: lat
        write(LIS_logunit,*) "FLake1: reading parameter LAT from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, "lat", placeholder)
        do i = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(i)%col
            row = LIS_surface(n, mtype)%tile(i)%row
            FLAKE1_struc(n)%flake1(i)%lat = placeholder(col, row)
        enddo 

        ! read: depth_w
        write(LIS_logunit,*) "FLake1: reading parameter DEPTH_W from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, "LAKEDEPTH", placeholder)
        do i = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(i)%col
            row = LIS_surface(n, mtype)%tile(i)%row
            FLAKE1_struc(n)%flake1(i)%depth_w = placeholder(col, row)
        enddo 
        ! read: fetch
        write(LIS_logunit,*) "FLake1: reading parameter FETCH from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, "LAKEWINDFETCH", placeholder)
        do i = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(i)%col
            row = LIS_surface(n, mtype)%tile(i)%row
            FLAKE1_struc(n)%flake1(i)%fetch = placeholder(col, row)
        enddo 

        ! read: depth_bs
        write(LIS_logunit,*) "FLake1: reading parameter DEPTH_BS from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, "LAKESEDIMDEPTH", placeholder)
        do i = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(i)%col
            row = LIS_surface(n, mtype)%tile(i)%row
            FLAKE1_struc(n)%flake1(i)%depth_bs = placeholder(col, row)
        enddo 

        ! read: T_bs
        write(LIS_logunit,*) "FLake1: reading parameter T_BS from ", trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, "LAKESEDIMTEMP", placeholder)
        do i = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(i)%col
            row = LIS_surface(n, mtype)%tile(i)%row
            FLAKE1_struc(n)%flake1(i)%t_bs = placeholder(col, row)
        enddo 
        deallocate(placeholder)


    enddo

end subroutine FLake1_setup

