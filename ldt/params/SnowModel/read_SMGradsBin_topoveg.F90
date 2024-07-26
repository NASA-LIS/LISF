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
! !ROUTINE: read_SMGradsBin_topoveg
! \label{read_SMGradsBin_topoveg}
!
! !REVISION HISTORY:
!  16Jul2020: Kristi Arsenault; Added SnowModel topo-veg reader
!  16Dec2021: Kristi Arsenault; Added parallel read of parameters
!
! !INTERFACE:
subroutine read_SMGradsBin_topoveg(n, array1, array2)

! !USES:
  use LDT_coreMod
  use LDT_paramDataMod
  use LDT_logMod,   only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use SnowModel_parmsMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n       ! Input nest index
  real, intent(inout) :: array1(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
  real, intent(inout) :: array2(LDT_rc%lnc(n),LDT_rc%lnr(n),1)

! !DESCRIPTION:
!  This subroutine retrieves SnowModel's topo-veg data and reprojects
!  it to the latlon projection. 
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[array1]
!   output field with the retrieved topo data
!  \item[array2]
!   output field with the retrieved vege data
!  \end{description}
!EOP

  integer :: ftn
  integer :: c, r
  logical :: file_exists
  real, allocatable :: read_topomap(:,:)   ! Read input parameters
  real, allocatable :: read_vegmap(:,:)

! ____________________________________________________________________________________

  array1 = LDT_rc%udef
  array2 = LDT_rc%udef

  inquire(file=trim(SnowModel_struc(n)%topoveg_file), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) "SnowModel topo-veg map, ",&
           trim(SnowModel_struc(n)%topoveg_file),", not found."
     call LDT_endrun
  endif
  select case ( SnowModel_struc(n)%topoveg_gridtransform )
    case( "none" ) 
      write(LDT_logunit,*) "[INFO] Reading Grads_binary topoveg file: ",&
            trim(SnowModel_struc(n)%topoveg_file)
  case default
     write(LDT_logunit,*) "[ERR] Since the Topo-Veg field involves discrete data values,"
     write(LDT_logunit,*) "  Please check your entries for this parameter."
     call LDT_endrun
  end select

  ftn = LDT_getNextUnitNumber()
  open( ftn, file=SnowModel_struc(n)%topoveg_file, &
        form="unformatted", access='direct',status='old', &
        recl=4*LDT_rc%gnc(n)*LDT_rc%gnr(n) )

!     Original code from Snowmodel's preprocess.f
!       open (unit=37,file=topoveg_fname, &
!             form='unformatted',access='direct',recl=4*nx*ny)
!       read (37,rec=1) ((topo_land(i,j),i=1,nx),j=1,ny)

  ! Read in topographic map
  allocate(read_topomap(LDT_rc%gnc(n),LDT_rc%gnr(n)))
  read(ftn,rec=1) ((read_topomap(c,r),c=1,LDT_rc%gnc(n)),r=1,LDT_rc%gnr(n))

  array1(1:LDT_rc%lnc(n),1:LDT_rc%lnr(n),1) = &
        read_topomap((LDT_ews_halo_ind(n,LDT_localPet+1)):(LDT_ewe_halo_ind(n,LDT_localPet+1)),&
                     (LDT_nss_halo_ind(n,LDT_localPet+1)):(LDT_nse_halo_ind(n,LDT_localPet+1)))
  deallocate(read_topomap)

  ! Read in vegetation type map
  allocate(read_vegmap(LDT_rc%gnc(n),LDT_rc%gnr(n)))
  read(ftn,rec=2) ((read_vegmap(c,r),c=1,LDT_rc%gnc(n)),r=1,LDT_rc%gnr(n))

  array2(1:LDT_rc%lnc(n),1:LDT_rc%lnr(n),1) = &
        read_vegmap((LDT_ews_halo_ind(n,LDT_localPet+1)):(LDT_ewe_halo_ind(n,LDT_localPet+1)),&
                    (LDT_nss_halo_ind(n,LDT_localPet+1)):(LDT_nse_halo_ind(n,LDT_localPet+1)))
  deallocate(read_vegmap)
 
  call LDT_releaseUnitNumber(ftn)

  ! Set landmask all to "=1" since SnowModel expects to have a land/water type
  !  at each point in the domain:
  LDT_LSMparam_struc(n)%landmask%value = 1
  LDT_LSMparam_struc(n)%dommask%value = 1

  write(LDT_logunit, *) "[INFO] Done reading Grads_binary topo-veg file"

end subroutine read_SMGradsBin_topoveg
