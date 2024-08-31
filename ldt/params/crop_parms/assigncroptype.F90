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
! !ROUTINE: assigncroptype
!  \label{assigncroptype}
!
! !REVISION HISTORY:
!  23  Jul 2014: KR Arsenault; Implementing crop layers
!  13  Mar 2020: HK Beaudoing; Added MIRCA dataset handling
!
! !INTERFACE:
subroutine assigncroptype( nest, crop_classification, &
                           num_types, croptype, crop_index)

! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_logunit, LDT_getNextUnitNumber, &
          LDT_releaseUnitNumber, LDT_verify, LDT_endrun
  use LDT_LSMCropModifier_Mod

  implicit none

! !ARGUMENTS: 
  integer,       intent(in)  :: nest
  integer,       intent(in)  :: num_types
  character(20), intent(in)  :: crop_classification
  character(20), intent(in)  :: croptype
  integer,       intent(out) :: crop_index
!
! !DESCRIPTION:
!  This subroutine reads a CROP classification and single-crop
!   entry and returns the classification's index value for that 
!   crop type.
!
!  The arguments are:
!  \begin{description}
!   \item[num_types]
!     number of crop types
!   \item[crop_classification]
!     Crop classification
!   \item[croptype]
!     the crop parameter name entry
!   \item[crop_index]
!     the returned index value for the entered crop type
!   \end{description}
!EOP      
!- Local:
   integer        :: ftn, ios1
   integer        :: i, k
   character(100) :: cropinv_file
   character(100) :: header1
   character(20)  :: read_cropname
   character(40)  :: read_fullname
! _____________________________________

  write(LDT_logunit,*)"[INFO] Obtaining the crop type for given classification: ",&
                       trim(crop_classification)

!- Create filename for the crop classification inventory:
!  e.g., FAOSTAT05_Crop.Inventory 
   cropinv_file = trim(LDT_LSMCrop_struc(nest)%croplib_dir)//&
                  trim(crop_classification)//"_Crop.Inventory"

   write(LDT_logunit,*) "- Reading Crop Library File: ",&
                         trim(cropinv_file)

!- Open file:
   ftn = LDT_getNextUnitNumber()
   open(ftn, file=cropinv_file, status='old', form='formatted',&
        iostat=ios1)

!- Read crop inventory file:
   read(ftn,fmt=*) header1
   do i = 1, num_types
      read(ftn,fmt=*) k, read_cropname, read_fullname 
!-- special case handling for MIRCA
      if ( trim(crop_classification) .eq. "MIRCA" .or. &
           trim(crop_classification) .eq. "MIRCA52" ) then
           call filetocropname (read_cropname)
      endif
      if( croptype == read_cropname ) then
        crop_index = k    ! Crop_index assignment
        exit
      endif
   end do
   
   call LDT_releaseUnitNumber(ftn)

end subroutine assigncroptype

 subroutine filetocropname (cropname)
! MIRCA filename contains number instead of cropname, so it needs to be 
! mapped out to those found in CROPMAP and FAO
 implicit none
 character(20),intent(inout)  :: cropname
 character(20) :: temp
 
 temp = cropname
 select case (trim(temp))
  case ("01")
    cropname = "wheat"
  case ("02")
    cropname = "maize"
  case ("03")
    cropname = "rice"
  case ("04")
    cropname = "barley"
  case ("05")
    cropname = "rye"
  case ("06")
    cropname = "millet"
  case ("07")
    cropname = "sorghum"
  case ("08")
    cropname = "soybean"
  case ("09")
    cropname = "sunflower"
  case ("10")
    cropname = "potato"
  case ("11")
    cropname = "cassava"
  case ("12")
    cropname = "sugarcane"
  case ("13")
    cropname = "sugarbeet"
  case ("14")
    cropname = "oilpalm"
  case ("15")
    cropname = "rapeseed"
  case ("16")
    cropname = "groundnut"
  case ("17")
    cropname = "pulsenes"
  case ("18")
    cropname = "citrusnes"
  case ("19")
    cropname = "datepalm"
  case ("20")
    cropname = "grape"
  case ("21")
    cropname = "cotton"
  case ("22")
    cropname = "cocoa"
  case ("23")
    cropname = "coffee"
  case ("24")
    cropname = "othersper"
  case ("25")
    cropname = "grassnes"
  case ("26")
    cropname = "othersann"
 end select
 end subroutine filetocropname
