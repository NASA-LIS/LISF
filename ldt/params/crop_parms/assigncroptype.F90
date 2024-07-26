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
      if( croptype == read_cropname ) then
        crop_index = k    ! Crop_index assignment
        exit
      endif
   end do
   
   call LDT_releaseUnitNumber(ftn)

end subroutine assigncroptype
