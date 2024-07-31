!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: calc_SlopeAspect_module
! \label{calc_SlopeAspect_module}
!
! !REVISION HISTORY:
!  21 May  2014: KR Arsenault; Estimate slope and aspect from elevation
!
! !INTERFACE:
 module calc_SlopeAspect_module
!
! !DESCRIPTION:
!  This subroutine module contains routines which can estimate slope and
!   and aspect from elevation field input.
!
!EOP

contains 

 subroutine calc_Slope_fromElev( elev_input, slope_output, res, lon0, lat0, nx, ny )

! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_logunit, LDT_endrun

  implicit none
! !ARGUMENTS: 
  real,    intent(in) :: elev_input(nx,ny)
  real,    intent(out):: slope_output(nx,ny)
  real,    intent(in) :: res
  real,    intent(in) :: lat0, lon0
  integer, intent(in) :: nx, ny

! !DESCRIPTION:
!  This subroutine calculates slope from elevation input field.
!
!  The arguments are:
!  \begin{description}
!  \item[elev]
!    input elevation field
!  \item[slope]
!    output slope field
!  \item[res]
!    resolution of input and output field
!  \item[lat0]
!    lower-left latitude input gridcell value
!  \item[lon0]
!    lower-left longitude input gridcell value
!  \item[nx]
!    number of x-dir column points 
!  \item[ny]
!    number of y-dir column points 
!  \end{description}
!EOP

  integer :: ib, i, j, k, ii, jj, ie, iw
  real    :: lat, lon
  real    :: dx, dy, dzdx, dzdy, pi
  real    :: Az, Bz, Cz, Dz, Ez, Fz, Gz, Hz, Iz
  real    :: grid(-1:1, -1:1)
  real, allocatable :: elev(:, :), slp(:) !, asp(:)
  real, parameter   :: Re = 6378.1E3  ! Radius of Earth, m

  pi = acos(-1.0)

  write(LDT_logunit,*) "... Calculating Slope Field from Elevation field ... "

  allocate( elev(0:nx-1, -1:1) )    ! i starts from 0 for easier wrap-around calc.
  allocate( slp(0:nx-1) )

  slope_output = -9999.

! Dealing with south and north edges 
  Do j = 1, ny
     if (j .eq. 1) then   ! initial point
       elev(:, -1) = elev_input(:,1)
       elev(:,  0) = elev_input(:,1)
       elev(:,  1) = elev_input(:,2)
!       read(11, rec=1) elev(:, 0)
!       read(11, rec=2) elev(:, 1)
!       elev(:, -1) = elev(:, 0)
     else if (j .eq. ny) then  ! "wrap-around"
       elev(:, -1) = elev_input(:, ny-1)
       elev(:,  0) = elev_input(:, ny)
       elev(:,  1) = elev_input(:, ny)
!       elev(:, -1) = elev(:, 0)
!       elev(:, 0) = elev(:, 1)
!       elev(:, 1) = elev(:, 0)
     else     ! all other points
       elev(:, -1) = elev_input(:, j-1)
       elev(:,  0) = elev_input(:, j)
       elev(:,  1) = elev_input(:, j)
!       elev(:, -1) = elev(:, 0)
!       elev(:,  0) = elev(:, 1)
!       read(11, rec=j+1) elev(:, 1)
     end if

     Do i = 0, nx-1
        lat = lat0 + (j-1) * res
        lon = lon0 + i * res
        dx = Re * cos(lat*pi/180.0) * (res * pi /180.0)
        dy = Re * res * pi /180.0
        iw = mod(i - 1 + nx, nx)
        ie = mod(i + 1, nx)

        Do jj=-1, 1
           grid(-1, jj) = elev(iw, jj)
           grid(0, jj)  = elev(i, jj)
           grid(1, jj)  = elev(ie, jj)
       End Do

     ! Fill undefined values 
       if ( grid(0, 0) .eq. -9999.0) then
         slp(i) = -9999.0
!         slope_output(i+1,j) = -9999.0
       else
         do jj = -1, 1
           do ii = -1, 1
              if(grid(ii, jj) .eq. -9999.0) then
                grid(ii, jj) = grid(0, 0)
              end if
           end do
         end do
       end if

     ! 3x3 grid :   A  B  C 
     !              D  E  F 
     !              G  H  I 

       Az = grid(-1, 1)
       Bz = grid(0, 1)
       Cz = grid(1, 1)

       Dz = grid(-1, 0)
       Ez = grid(0, 0)
       Fz = grid(1, 0)

       Gz = grid(-1, -1)
       Hz = grid(0, -1)
       Iz = grid(1, -1)

       dzdx = - ( (Az - Cz) + 2.0 * (Dz - Fz) + (Gz - Iz) ) / (8.0 * dx )
       dzdy =   ( (Az - Gz) + 2.0 * (Bz - Hz) + (Cz - Iz) ) / (8.0 * dy )

       slp(i) = atan( Sqrt( dzdx*dzdx + dzdy*dzdy) )    ! radian
!       slope_output(i+1,j) = atan( Sqrt( dzdx*dzdx + dzdy*dzdy) )    ! radian

       slope_output(i+1,j) = slp(i)

    End do  ! i
  End Do    ! j

  deallocate(elev)
  deallocate(slp)

  write(LDT_logunit,*) "... Finished Calculating Slope Field ... "

 end subroutine calc_Slope_fromElev


 subroutine calc_Aspect_fromElev( elev_input, aspect_output, res, lon0, lat0, nx, ny )

! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_logunit, LDT_endrun

  implicit none
! !ARGUMENTS: 
  real,    intent(in) :: elev_input(nx,ny)
  real,    intent(out):: aspect_output(nx,ny)
  real,    intent(in) :: res
  real,    intent(in) :: lat0, lon0
  integer, intent(in) :: nx, ny

! !DESCRIPTION:
!  This subroutine calculates aspect from elevation input field.
!
!  The arguments are:
!  \begin{description}
!  \item[elev]
!    input elevation field
!  \item[aspect]
!    output aspect field
!  \item[res]
!    resolution of input and output field
!  \item[lat0]
!    lower-left latitude input gridcell value
!  \item[lon0]
!    lower-left longitude input gridcell value
!  \item[nx]
!    number of x-dir column points 
!  \item[ny]
!    number of y-dir column points 
!  \end{description}
!EOP

  integer :: ib, i, j, k, ii, jj, ie, iw
  real    :: lat, lon
  real    :: dx, dy, dzdx, dzdy, pi
  real    :: Az, Bz, Cz, Dz, Ez, Fz, Gz, Hz, Iz
  real    :: grid(-1:1, -1:1)
  real, allocatable :: elev(:, :),  asp(:)
  real, parameter   :: Re = 6378.1E3  ! Radius of Earth, m

  pi = acos(-1.0)

  write(LDT_logunit,*) "... Calculating Aspect Field from Elevation field ... "

  allocate( elev(0:nx-1, -1:1) )    ! i starts from 0 for easier wrap-around calc.
  allocate( asp(0:nx-1) )

  aspect_output = -9999.0

! Dealing with south and north edges 
  Do j = 1, ny
     if (j .eq. 1) then   ! initial point
       elev(:, -1) = elev_input(:,1)
       elev(:,  0) = elev_input(:,1)
       elev(:,  1) = elev_input(:,2)
!       read(11, rec=1) elev(:, 0)
!       read(11, rec=2) elev(:, 1)
!       elev(:, -1) = elev(:, 0)
     else if (j .eq. ny) then  ! "wrap-around"
       elev(:, -1) = elev_input(:, ny-1)
       elev(:,  0) = elev_input(:, ny)
       elev(:,  1) = elev_input(:, ny)
!       elev(:, -1) = elev(:, 0)
!       elev(:, 0) = elev(:, 1)
!       elev(:, 1) = elev(:, 0)
     else     ! all other points
       elev(:, -1) = elev_input(:, j-1)
       elev(:,  0) = elev_input(:, j)
       elev(:,  1) = elev_input(:, j)
!       elev(:, -1) = elev(:, 0)
!       elev(:,  0) = elev(:, 1)
!       read(11, rec=j+1) elev(:, 1)
     end if

     Do i = 0, nx-1
        lat = lat0 + (j-1) * res
        lon = lon0 + i * res
        dx = Re * cos(lat*pi/180.0) * (res * pi /180.0)
        dy = Re * res * pi /180.0
        iw = mod(i - 1 + nx, nx)
        ie = mod(i + 1, nx)

        Do jj=-1, 1
           grid(-1, jj) = elev(iw, jj)
           grid(0, jj)  = elev(i, jj)
           grid(1, jj)  = elev(ie, jj)
       End Do

     ! Fill undefined values 
       if ( grid(0, 0) .eq. -9999.0) then
         asp(i) = -9999.0
       else
         do jj = -1, 1
           do ii = -1, 1
              if(grid(ii, jj) .eq. -9999.0) then
                grid(ii, jj) = grid(0, 0)
              end if
           end do
         end do
       end if

     ! 3x3 grid :   A  B  C 
     !              D  E  F 
     !              G  H  I 

       Az = grid(-1, 1)
       Bz = grid(0, 1)
       Cz = grid(1, 1)

       Dz = grid(-1, 0)
       Ez = grid(0, 0)
       Fz = grid(1, 0)

       Gz = grid(-1, -1)
       Hz = grid(0, -1)
       Iz = grid(1, -1)

       dzdx = - ( (Az - Cz) + 2.0 * (Dz - Fz) + (Gz - Iz) ) / (8.0 * dx )
       dzdy =   ( (Az - Gz) + 2.0 * (Bz - Hz) + (Cz - Iz) ) / (8.0 * dy )

       if ( dzdx .ne. 0.0) then
          asp(i)  = pi - atan( dzdy / dzdx ) + pi * 0.5 * dzdx / abs(dzdx)
       else if (dzdy .lt. 0 )  then
          asp(i) = 0.0
       else
          asp(i) = pi
       end if

       aspect_output(i+1,j) = asp(i)

    End do  ! i
  End Do    ! j

  deallocate(elev)
  deallocate(asp)

  write(LDT_logunit,*) "... Finished Calculating Aspect Field ... "


 end subroutine calc_Aspect_fromElev

end module calc_SlopeAspect_module
