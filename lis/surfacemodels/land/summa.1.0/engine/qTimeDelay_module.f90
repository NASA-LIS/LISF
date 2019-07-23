! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module qTimeDelay_module
USE nrtype
implicit none
private
public::qOverland
contains


 ! *************************************************************************************************************
 ! public subroutine qOverland: compute the time delay in runoff in a basin (places runoff in future time steps)
 ! *************************************************************************************************************
 subroutine qOverland(&
                      ! input
                      ixRouting,             &  ! index for routing method
                      averageSurfaceRunoff,  &  ! surface runoff (m s-1)
                      averageSoilBaseflow,   &  ! baseflow from the soil profile (m s-1)
                      averageAquiferBaseflow,&  ! baseflow from the aquifer (m s-1)
                      fracFuture,            &  ! fraction of runoff in future time steps (m s-1)
                      qFuture,               &  ! runoff in future time steps (m s-1)
                      ! output
                      averageInstantRunoff,  &  ! instantaneous runoff (m s-1)
                      averageRoutedRunoff,   &  ! routed runoff (m s-1)
                      err,message)              ! error control
 ! look-up values for the sub-grid routing method
 USE mDecisions_module,only:      &
  timeDelay,&  ! time-delay histogram
  qInstant     ! instantaneous routing
 implicit none
 ! input
 integer(i4b),intent(in)    :: ixRouting              ! index for routing method
 real(dp),intent(in)        :: averageSurfaceRunoff   ! surface runoff (m s-1)
 real(dp),intent(in)        :: averageSoilBaseflow    ! baseflow from the soil profile (m s-1)
 real(dp),intent(in)        :: averageAquiferBaseflow ! baseflow from the aquifer (m s-1)
 real(dp),intent(in)        :: fracFuture(:)          ! fraction of runoff in future time steps (m s-1)
 real(dp),intent(inout)     :: qFuture(:)             ! runoff in future time steps (m s-1)
 ! output
 real(dp),intent(out)       :: averageInstantRunoff   ! instantaneous runoff (m s-1)
 real(dp),intent(out)       :: averageRoutedRunoff    ! routed runoff (m s-1)
 integer(i4b),intent(out)   :: err                    ! error code
 character(*),intent(out)   :: message                ! error message
 ! internal
 integer(i4b)               :: nTDH                   ! number of points in the time-delay histogram
 integer(i4b)               :: iFuture                ! index in time delay histogram
 ! initialize error control
 err=0; message='qOverland/'

 ! compute instantaneous runoff (m s-1)
 averageInstantRunoff = averageSurfaceRunoff + averageAquiferBaseflow + averageSoilBaseflow

 ! compute routed runoff (m s-1)
 select case(ixRouting)  ! (select option for sub-grid routing)
  ! ** instantaneous routing
  case(qInstant)
   averageRoutedRunoff = averageInstantRunoff

  ! ** time delay histogram
  case(timeDelay)
   ! identify number of points in the time-delay histogram
   nTDH = size(qFuture)
   ! place a fraction of runoff in future steps
   qFuture(1:nTDH) = qFuture(1:nTDH) + averageInstantRunoff*fracFuture(1:nTDH)
   ! save the routed runoff
   averageRoutedRunoff = qFuture(1)
   ! move array back
   do iFuture=2,nTDH
    qFuture(iFuture-1) = qFuture(iFuture)
   end do
   qFuture(nTDH) = 0._dp

   !print*, 'averageInstantRunoff, averageRoutedRunoff = ', averageInstantRunoff, averageRoutedRunoff
   !print*, 'qFuture(1:100) = ', qFuture(1:100)
   !pause

  ! ** error checking
  case default; err=20; message=trim(message)//'cannot find option for sub-grid routing'; return

 end select ! (select option for sub-grid routing)

 end subroutine qOverland


end module qTimeDelay_module
