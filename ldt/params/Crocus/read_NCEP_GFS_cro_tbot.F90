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
! !ROUTINE: read_NCEP_cro_tbot
! \label{read_NCEP_cro_tbot}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  23 Sep 2019: Mahdi Navari, modiffied for Crocus 
!
! !INTERFACE:
subroutine read_NCEP_GFS_cro_tbot(n, array)

! !USES:
  use LDT_coreMod,     only : LDT_rc, LDT_domain
  use LDT_logMod,      only : LDT_logunit, LDT_getNextUnitNumber, &
          LDT_releaseUnitNumber, LDT_endrun
  use LDT_gridmappingMod
  use Crocus_parmsMod

  implicit none

! !ARGUMENTS: 
  integer,   intent(in) :: n
  real,   intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),1)

! !DESCRIPTION:
!  This subroutine retrieves the bottom temperature climatology for the 
!  specified month and returns the values in the latlon projection
!  
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[array]
!   output field with the retrieved greenness fraction
!  \end{description}
!
!EOP      
  integer :: ftn
  integer :: c,r
  logical :: file_exists
  integer :: subpnc, subpnr, glpnc, glpnr
  real    :: subparam_gridDesc(20)       ! Input parameter grid desc array
  integer, allocatable :: lat_line(:,:), lon_line(:,:)
  real,    allocatable :: read_tbot(:,:), dummy_lat(:,:), dummy_lon(:,:)
! _______________________________________________

  inquire(file=trim(Crocus_struc(n)%tbotFile), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) "[ERR] GFS TBOT map ",trim(Crocus_struc(n)%tbotFile)," not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif

  write(LDT_logunit,*)"[INFO] Reading NCEP-GFS Bottom Temperature file: ",&
        trim(Crocus_struc(n)%tbotfile)

  ftn = LDT_getNextUnitNumber()
  open(ftn, file=trim(Crocus_struc(n)%tbotFile), access='sequential',status='old', &
       form="unformatted", recl=4)
  
  subparam_gridDesc = 0.
  call LDT_RunDomainPts(n, Crocus_struc(n)%tbot_proj, Crocus_struc(n)%tbot_gridDesc(:), &
                         glpnc,glpnr,subpnc,subpnr,  &
                         subparam_gridDesc,lat_line,lon_line )

  allocate(dummy_lat(glpnc,glpnr))
  allocate(dummy_lon(glpnc,glpnr))
  allocate(read_tbot(glpnc,glpnr))

  read(ftn) dummy_lat
  read(ftn) dummy_lon
  read(ftn) read_tbot

  close(ftn)
  deallocate(dummy_lat)
  deallocate(dummy_lon)

  select case (Crocus_struc(n)%tbot_gridtransform)
    case("none")
      do r = 1,subpnr
        do c = 1,subpnc
           array(c,r,1) = read_tbot(lon_line(c,r),lat_line(c,r))
        enddo ! columns
      enddo ! rows
      deallocate(read_tbot)
    case default
       write(LDT_logunit,*) "[ERR] The spatial transform, ",&
                            trim(Crocus_struc(n)%tbot_gridtransform),&
                            ", for GFS bottom temperature is not available at this time ... "
       write(LDT_logunit,*) "Program stopping ..."
       call LDT_endrun
   end select

   do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        if(array(c,r,1).lt.0) then
           array(c,r,1) = LDT_rc%udef
        endif
     enddo
   enddo

  call LDT_releaseUnitNumber(ftn)
end subroutine read_NCEP_GFS_cro_tbot
