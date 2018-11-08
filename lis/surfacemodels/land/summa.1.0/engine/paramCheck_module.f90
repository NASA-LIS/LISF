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

module paramCheck_module
! define numerical recipes data type
USE nrtype
! define look-up values for the choice of method to combine and sub-divide snow layers
USE mDecisions_module,only:&
 sameRulesAllLayers, & ! SNTHERM option: same combination/sub-dividion rules applied to all layers
 rulesDependLayerIndex ! CLM option: combination/sub-dividion rules depend on layer index
implicit none
private
public::paramCheck
contains


 ! ************************************************************************************************
 ! public subroutine paramCheck: check consistency of model parameters
 ! ************************************************************************************************
 subroutine paramCheck(mpar_data,err,message)
 ! model decisions
 USE globalData,only:model_decisions  ! model decision structure
 USE var_lookup,only:iLookDECISIONS   ! named variables for elements of the decision structure
 ! SUMMA look-up variables
 USE data_types,only:var_dlength      ! data vector with variable length dimension (dp): x%var(:)%dat(:)
 USE var_lookup,only:iLookPARAM       ! named variables for elements of the data structures
 implicit none
 ! define input
 type(var_dlength),intent(in)    :: mpar_data            ! model parameters
 ! define output
 integer(i4b),intent(out)        :: err                  ! error code
 character(*),intent(out)        :: message              ! error message
 ! local variables
 integer(i4b)                    :: iLayer               ! index of model layers
 real(dp),dimension(5)           :: zminLayer            ! minimum layer depth in each layer (m)
 real(dp),dimension(4)           :: zmaxLayer_lower      ! lower value of maximum layer depth
 real(dp),dimension(4)           :: zmaxLayer_upper      ! upper value of maximum layer depth
 ! Start procedure here
 err=0; message="paramCheck/"

 ! *****
 ! * check that the snow layer bounds are OK...
 ! ********************************************

 ! select option for combination/sub-division of snow layers
 select case(model_decisions(iLookDECISIONS%snowLayers)%iDecision)
  ! SNTHERM option
  case(sameRulesAllLayers)
   if(mpar_data%var(iLookPARAM%zmax)%dat(1)/mpar_data%var(iLookPARAM%zmin)%dat(1) < 2.5_dp)then
    message=trim(message)//'zmax must be at least 2.5 times larger than zmin: this avoids merging layers that have just been divided'
    err=20; return
   end if
  ! CLM option
  case(rulesDependLayerIndex)
   ! (build vectors of min/max)
   zminLayer       = (/mpar_data%var(iLookPARAM%zminLayer1)%dat(1),&
                       mpar_data%var(iLookPARAM%zminLayer2)%dat(1),&
                       mpar_data%var(iLookPARAM%zminLayer3)%dat(1),&
                       mpar_data%var(iLookPARAM%zminLayer4)%dat(1),&
                       mpar_data%var(iLookPARAM%zminLayer5)%dat(1)/)
   zmaxLayer_lower = (/mpar_data%var(iLookPARAM%zmaxLayer1_lower)%dat(1),&
                       mpar_data%var(iLookPARAM%zmaxLayer2_lower)%dat(1),&
                       mpar_data%var(iLookPARAM%zmaxLayer3_lower)%dat(1),&
                       mpar_data%var(iLookPARAM%zmaxLayer4_lower)%dat(1)/)
   zmaxLayer_upper = (/mpar_data%var(iLookPARAM%zmaxLayer1_upper)%dat(1),&
                       mpar_data%var(iLookPARAM%zmaxLayer2_upper)%dat(1),&
                       mpar_data%var(iLookPARAM%zmaxLayer3_upper)%dat(1),&
                       mpar_data%var(iLookPARAM%zmaxLayer4_upper)%dat(1)/)
   ! (check consistency)
   do iLayer=1,4  ! NOTE: the lower layer does not have a maximum value
    ! ensure that we have higher maximum thresholds for sub-division when fewer number of layers
    if(zmaxLayer_lower(iLayer) < zmaxLayer_upper(iLayer))then
     write(message,'(a,2(i0,a))') trim(message)//'expect the maximum threshold for sub-division in the case where there is only ', &
                                  iLayer,' layer(s) is greater than the maximum threshold for sub-division in the case where there are > ',&
                                  iLayer,' layer(s)'
     err=20; return
    end if
    ! ensure that the maximum thickness is 3 times greater than the minimum thickness
    if(zmaxLayer_upper(iLayer)/zminLayer(iLayer) < 2.5_dp .or. zmaxLayer_upper(iLayer)/zminLayer(iLayer+1) < 2.5_dp)then
     write(*,'(a,1x,3(f20.10,1x))') 'zmaxLayer_upper(iLayer), zminLayer(iLayer), zminLayer(iLayer+1) = ', &
                                     zmaxLayer_upper(iLayer), zminLayer(iLayer), zminLayer(iLayer+1)
     write(message,'(a,3(i0,a))') trim(message)//'zmaxLayer_upper for layer ',iLayer,' must be 2.5 times larger than zminLayer for layers ',&
                                  iLayer,' and ',iLayer+1,': this avoids merging layers that have just been divided'
     err=20; return
    end if
   end do  ! loop through layers
  case default; err=20; message=trim(message)//'unable to identify option to combine/sub-divide snow layers'; return
 end select ! (option to combine/sub-divide snow layers)

 ! -------------------------------------------------------------------------------------------------------------------------------------------

 ! *****
 ! * check parameter dependencies...
 ! *********************************

 ! associations
 associate(&
 ! canopy geometry
 heightCanopyTop        => mpar_data%var(iLookPARAM%heightCanopyTop)%dat(1),   & ! intent(in): [dp] height at the top of the vegetation canopy (m)
 heightCanopyBottom     => mpar_data%var(iLookPARAM%heightCanopyBottom)%dat(1),& ! intent(in): [dp] height at the bottom of the vegetation canopy (m)
 ! transpiration
 critSoilWilting        => mpar_data%var(iLookPARAM%critSoilWilting)%dat(1),   & ! intent(in): [dp] critical vol. liq. water content when plants are wilting (-)
 critSoilTranspire      => mpar_data%var(iLookPARAM%critSoilTranspire)%dat(1), & ! intent(in): [dp] critical vol. liq. water content when transpiration is limited (-)
 ! soil properties
 fieldCapacity          => mpar_data%var(iLookPARAM%fieldCapacity)%dat(1),     & ! intent(in): [dp]    field capacity (-)
 theta_sat              => mpar_data%var(iLookPARAM%theta_sat)%dat,            & ! intent(in): [dp(:)] soil porosity (-)
 theta_res              => mpar_data%var(iLookPARAM%theta_res)%dat             & ! intent(in): [dp(:)] soil residual volumetric water content (-)
 ) ! associations to parameters

 ! check canopy geometry
 if(heightCanopyTop < heightCanopyBottom)then
  write(message,'(a,i0,a)') trim(message)//'height of canopy top is less than the height of the canopy bottom'
  err=20; return
 endif

 ! check that the maximum transpiration limit is within bounds
 if( any(critSoilTranspire > theta_sat) .or. any(critSoilTranspire < theta_res) )then
  print*, 'theta_res         = ', theta_res
  print*, 'theta_sat         = ', theta_sat
  print*, 'critSoilTranspire = ', critSoilTranspire
  message=trim(message)//'critSoilTranspire parameter is out of range '// &
                         '[NOTE: if overwriting Noah-MP soil table values in paramTrial, must overwrite all soil parameters]'
  err=20; return
 end if

 ! check that the soil wilting point is within bounds
 if( any(critSoilWilting > theta_sat) .or. any(critSoilWilting < theta_res) )then
  print*, 'theta_res       = ', theta_res
  print*, 'theta_sat       = ', theta_sat
  print*, 'critSoilWilting = ', critSoilWilting
  message=trim(message)//'critSoilWilting parameter is out of range '// &
                         '[NOTE: if overwriting Noah-MP soil table values in paramTrial, must overwrite all soil parameters]'
  err=20; return
 end if

 ! check that the field capacity is within bounds
 if( any(fieldCapacity > theta_sat) .or. any(fieldCapacity < theta_res) )then
  print*, 'theta_res     = ', theta_res
  print*, 'theta_sat     = ', theta_sat
  print*, 'fieldCapacity = ', fieldCapacity
  message=trim(message)//'fieldCapacity parameter is out of range '// &
                         '[NOTE: if overwriting Noah-MP soil table values in paramTrial, must overwrite all soil parameters]'
  err=20; return
 end if

 ! check transpiration
 if(critSoilTranspire < critSoilWilting)then
  write(message,'(a,i0,a)') trim(message)//'critical point for transpiration is less than the wilting point'
  err=20; return
 endif

 ! check porosity
 if( any(theta_sat < theta_res) )then
  write(message,'(a,i0,a)') trim(message)//'porosity is less than the residual liquid water content'
  err=20; return
 endif

 ! end associations to parameter dependencies
 end associate

 end subroutine paramCheck


end module paramCheck_module
