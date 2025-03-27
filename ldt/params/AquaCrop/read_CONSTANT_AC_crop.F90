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
  use LDT_constantsMod, only: LDT_CONST_PATH_LEN
  use LDT_coreMod,      only: LDT_rc, LDT_config
  use LDT_paramDataMod, only: LDT_LSMparam_struc
  use LDT_logMod, only: LDT_verify, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_logunit, LDT_endrun

  implicit none
  ! !ARGUMENTS: 
  integer, intent(in) :: n
  real, intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n))

  ! !LOCALS:
  integer        :: ftn, ios1
  integer        :: i, k, rc, num_types, crop_index
  integer        :: r,c
  character(100) :: croptype
  character(len=LDT_CONST_PATH_LEN) :: dir
  character(len=LDT_CONST_PATH_LEN) :: cropinv_file
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
  call ESMF_ConfigFindLabel(LDT_config,"AquaCrop crop type:",rc=rc)
  call ESMF_ConfigGetAttribute(LDT_Config, croptype, rc=rc)
  call LDT_verify(rc,"AquaCrop crop type: not defined")

  !- Finds crop inventory
  call ESMF_ConfigFindLabel(LDT_config,"AquaCrop crop library directory:",rc=rc)
  call ESMF_ConfigGetAttribute(LDT_Config, dir,rc=rc)
  call LDT_verify(rc,"AquaCrop crop library directory: not defined")

  cropinv_file = trim(dir)//"AC_Crop.Inventory"

  write(LDT_logunit,*) "- Reading Crop Library File: ",&
       trim(cropinv_file)

  !- Open file:
  ftn = LDT_getNextUnitNumber()
  open(ftn, file=cropinv_file, status='old', form='formatted',&
       iostat=ios1)
  call LDT_verify(ios1,"AquaCrop crop inventory does not exist")

  !- Read crop inventory file:
  read(ftn,fmt=*) num_types
  read(ftn,fmt=*) header1
  do i = 1, num_types
     read(ftn,fmt=*) k, read_cropname, read_fullname 
     if( trim(croptype) == trim(read_cropname) ) then
        crop_index = k    ! Crop_index assignment
        exit
     endif
     if ( (i.eq.num_types).and.(crop_index.ne.k) ) then
        write(LDT_logunit,*) "[ERR] AC72: ", trim(croptype), &
             " not in AC_Crop.Inventory"
        call LDT_endrun
     endif
  end do

  call LDT_releaseUnitNumber(ftn)

  ! fill the array with the crop type in case not masked
  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        if (LDT_LSMparam_struc(n)%landmask%value(c,r,1).gt.0) then
           array(c,r) = crop_index
        else
           array(c,r) = LDT_rc%udef
        endif
     enddo
  enddo

  write(LDT_logunit, *) "[INFO] Setting a CONSTANT crop type value ..."


end subroutine read_CONSTANT_AC_crop
