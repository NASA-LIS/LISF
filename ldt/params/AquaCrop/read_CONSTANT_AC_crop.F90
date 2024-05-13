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
! !ROUTINE: read_CONSTANT_AC_crop
! \label{read_CONSTANT_AC_crop}
!
! !REVISION HISTORY:
!  13 May 2024; Louise Busschaert: Initial implementation
!
! !INTERFACE:
subroutine read_CONSTANT_AC_crop(n, array)

! !USES:
  use ESMF
  use LDT_coreMod,  only : LDT_rc, LDT_config
  use LDT_logMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  real, intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n))

! !LOCALS:
  integer        :: ftn, ios1
  integer        :: i, k, rc, num_types, crop_index
  character(100) :: croptype
  character(200) :: dir
  character(100) :: cropinv_file
  character(100) :: header1
  character(100) :: read_cropname
  character(100) :: read_fullname

! !DESCRIPTION:
!  This subroutine sets a constant crop type defined in the
!  configuration file.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[array]
!   output field with the crop type
!  \end{description}
!EOP
! ____________________________

!- Finds chosen crop type
  call ESMF_ConfigFindLabel(LDT_config,"AC crop type:",rc=rc)
  call ESMF_ConfigGetAttribute(LDT_Config, croptype, rc=rc)
  call LDT_verify(rc,"AC crop type: not defined")

!- Finds crop inventory
  call ESMF_ConfigFindLabel(LDT_config,"AC crop library directory:",rc=rc)
  call ESMF_ConfigGetAttribute(LDT_Config, dir,rc=rc)
  call LDT_verify(rc,"AC crop library directory: not defined")

   cropinv_file = trim(dir)//"AC_Crop.Inventory"

   write(LDT_logunit,*) "- Reading Crop Library File: ",&
                         trim(cropinv_file)

!- Open file:
   ftn = LDT_getNextUnitNumber()
   open(ftn, file=cropinv_file, status='old', form='formatted',&
        iostat=ios1)

!- Read crop inventory file:
   read(ftn,fmt=*) num_types
   read(ftn,fmt=*) header1
   do i = 1, num_types
      read(ftn,fmt=*) k, read_cropname, read_fullname 
      if( trim(croptype) == trim(read_cropname) ) then
        crop_index = k    ! Crop_index assignment
        exit
      endif
   end do
   
   call LDT_releaseUnitNumber(ftn)
  array = crop_index    ! Set all values to chosen crop type

  write(LDT_logunit, *) "[INFO] Setting a CONSTANT crop type value ..."
  

end subroutine read_CONSTANT_AC_crop
