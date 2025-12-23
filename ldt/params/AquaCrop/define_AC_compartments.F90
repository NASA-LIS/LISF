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
! !ROUTINE: define_AC_compartments
! \label{define_AC_compartments}
!
! !REVISION HISTORY:
!  29 May 2024; Louise Busschaert: Initial implementation
!
! !INTERFACE:
subroutine define_AC_compartments(n, array)

  ! !USES:
  use ac_utils
  use AquaCrop_parmsMod
  use ESMF
  use LDT_constantsMod, only: LDT_CONST_PATH_LEN
  use LDT_coreMod,      only: LDT_rc, LDT_config
  use LDT_logMod, only: LDT_verify, LDT_logunit, &
       LDT_getNextUnitNumber, LDT_releaseUnitNumber
  use LDT_paramDataMod, only: LDT_LSMparam_struc

  implicit none

  ! !ARGUMENTS:
  integer, intent(in) :: n
  real, intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),12)

  ! !LOCALS:
  integer           :: i,j,k,num_types,ct,rc,nrcomp
  integer           :: ftn, ftn1, r,c, ios1
  real, allocatable :: depths(:)
  real, allocatable :: rootingd(:)
  real, allocatable :: csizes(:,:)
  real              :: TotalDepthL, TotDepthC, DeltaZ, fAdd
  character(len=LDT_CONST_PATH_LEN)    :: dir
  character(len=LDT_CONST_PATH_LEN)    :: cropinv_file, crop_file
  character(100)    :: header1, str1
  character(100)    :: read_cropname
  character(100)    :: read_fullname

  ! !DESCRIPTION:
  !  This subroutine sets a constant crop type defined in the
  !  configuration file.
  !
  !  The arguments are:
  !  \begin{description}
  !  \item[n]
  !   index of the nest
  !  \item[array]
  !   output field with the compartment size for the 12 levels
  !  \end{description}
  !EOP
  ! ____________________________

  !- Get number of layers
  call ESMF_ConfigFindLabel(LDT_config,"AquaCrop number of soil layers:",rc=rc)
  call ESMF_ConfigGetAttribute(LDT_config,Aquacrop_struc(n)%nlayers,rc=rc)
  call LDT_verify(rc,"AquaCrop number of soil layers: not defined")

  !- Get thickness
  allocate(depths(Aquacrop_struc(n)%nlayers))
  call ESMF_ConfigFindLabel(LDT_config,"AquaCrop soil layer thickness:",rc=rc)
  call ESMF_ConfigGetAttribute(LDT_config,depths,rc=rc)
  call LDT_verify(rc,"AquaCrop soil layer thickness: not defined")

  do i = 1, AquaCrop_struc(n)%nlayers
     AquaCrop_struc(n)%lthickness(i) = depths(i)
  enddo

  write(LDT_logunit,*) "[INFO] AquaCrop number of soil layers: ", Aquacrop_struc(n)%nlayers
  do i=1,Aquacrop_struc(n)%nlayers
     write (str1, *) i
     write(LDT_logunit,*) "Depth layer "//trim(str1)//" =", AquaCrop_struc(n)%lthickness(i)
  enddo

  !- Get max rooting depth
  !- Finds crop inventory
  call ESMF_ConfigFindLabel(LDT_config,"AquaCrop crop library directory:",rc=rc)
  call ESMF_ConfigGetAttribute(LDT_Config, dir,rc=rc)
  call LDT_verify(rc,"AquaCrop crop library directory: not defined")

  cropinv_file = trim(dir)//"AC_Crop.Inventory"
  !- Open file:
  ftn = LDT_getNextUnitNumber()
  open(ftn, file=cropinv_file, status='old', form='formatted',&
       iostat=ios1)
  call LDT_verify(ios1,"AquaCrop crop inventory does not exist")

  !- Read crop inventory file:
  read(ftn,fmt=*) num_types
  read(ftn,fmt=*) header1
  allocate(rootingd(num_types))
  do i = 1, num_types
     read(ftn,fmt=*) k, read_cropname, read_fullname 
     ! Open corresponding crop file
     crop_file = trim(dir)//"AC_Crop.Files/"//trim(read_cropname)//".CRO"
     ftn1 = LDT_getNextUnitNumber()
     open(ftn1, file=crop_file, status='old', form='formatted',&
          iostat=ios1)
     call LDT_verify(ios1,"AquaCrop crop file does not exist")
     do j=1,37
        read(ftn1,fmt=*)
     enddo
     read(ftn1,fmt=*) rootingd(i)
     ! Check validity
     if((rootingd(i).le.0).or.(rootingd(i).gt.3))then
        write(LDT_logunit,*) "[WARN] AquaCrop rooting depth not in valid range"
        write(LDT_logunit,*) "Program stopping ..."
     endif
     call LDT_releaseUnitNumber(ftn1)
  end do
  call LDT_releaseUnitNumber(ftn)

  ! calculate the compartment sizes for each crop type
  ! Max 12 compartments and by default 10 cm each
  allocate(csizes(num_types,12))

  do ct = 1, num_types
     ! Adjust number of compartments based on the soil profile
     ! depth
     ! (based on DetermineNrandThicknessCompartments from global.f90)
     ! Does not use rootingd
     TotalDepthL = 0.
     do i = 1, AquaCrop_struc(n)%nlayers
        TotalDepthL = TotalDepthL + AquaCrop_struc(n)%lthickness(i)
     end do
     TotDepthC = 0.
     nrcomp = 0
     loop: do
        DeltaZ = (TotalDepthL - TotDepthC)
        nrcomp = nrcomp + 1
        if (DeltaZ > 0.10) then
           csizes(ct,nrcomp) = 0.10
        else
           csizes(ct,nrcomp) = DeltaZ
        end if
        TotDepthC = TotDepthC + csizes(ct,nrcomp)
        if ((nrcomp == 12) &
             .or. (abs(TotDepthC - TotalDepthL) < 0.0001)) exit loop
     end do loop

     ! AquaCrop procedure to adjust the size of the compartments 
     ! (based on AdjustSizeCompartments in global.f90)
     ! Requires rootingd

     ! Calculate new depth of compartments
     TotDepthC = 0.
     do i = 1, nrcomp
        TotDepthC = TotDepthC + csizes(ct,i)
     enddo

     if (nrcomp < 12) then
        loop_bis: do
           nrcomp = nrcomp + 1
           if ((rootingd(ct) - TotDepthC) > 0.10) then
              csizes(ct,nrcomp) = 0.10
           else
              csizes(ct,nrcomp) = rootingd(ct) - TotDepthC
           end if
           TotDepthC = TotDepthC + csizes(ct,nrcomp)
           if ((nrcomp == 12) &
                .or. ((TotDepthC + 0.00001) >= rootingd(ct))) exit loop_bis
        end do loop_bis
     end if

     if ((TotDepthC + 0.00001) < rootingd(ct)) then
        fAdd = (rootingd(ct)/0.1 - 12.)/78.
        do i = 1, 12
           csizes(ct,i) = 0.1 * (1. + i*fAdd)
           csizes(ct,i) = 0.05 * real(roundc(csizes(ct,i)*20., mold=1))
        end do
        TotDepthC = 0.
        do i = 1, 12
           TotDepthC = TotDepthC + csizes(ct,i)
        end do
        if (rootingd(ct) - TotDepthC > ac_zero_threshold) then
           loop2: do
              csizes(ct,12) = csizes(ct,12) + 0.05
              TotDepthC = TotDepthC + 0.05
              if (TotDepthC >= rootingd(ct)) exit loop2
           end do loop2
        else
           do while ((TotDepthC - 0.04999999) >= rootingd(ct))
              csizes(ct,12) = csizes(ct,12) - 0.05
              TotDepthC = TotDepthC - 0.05
           end do
        end if
     end if
  enddo
  ! End AquaCrop subroutines

  ! fill the array with the compartment size
  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        if (LDT_LSMparam_struc(n)%landmask%value(c,r,1).gt.0)then
           array(c,r,:) = csizes(int(AquaCrop_struc(n)%cropt%value(c,r,1)),:)
        else
           array(c,r,:) = LDT_rc%udef
        endif
     enddo
  enddo
  deallocate(csizes)
  deallocate(rootingd)

end subroutine define_AC_compartments


