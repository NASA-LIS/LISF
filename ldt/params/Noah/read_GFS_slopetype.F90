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
! !ROUTINE: read_GFS_slopetype
! \label{read_GFS_slopetype}
!
! !REVISION HISTORY:
!  05 2014: Grey Nearing; Initial Specification
!
! !INTERFACE:
subroutine read_GFS_slopetype(n, array)

! !USES:
  use LDT_coreMod,        only : LDT_rc
  use LDT_logMod,         only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod,      only : readLISdata 
  use LDT_gridmappingmod
  use Noah_parmsMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  real, intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),1)

! !DESCRIPTION:
!  This subroutine retrieves static, slope type data and reprojects
!  it to the gaussian projection. 
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[array]
!   output field with the retrieved slope type data
!  \end{description}
!EOP

  integer :: ftn
  integer :: c, r
  logical :: file_exists
  integer :: glpnc, glpnr             ! Parameter (global) total columns and rows
  integer :: subpnc, subpnr           ! Parameter subsetted columns and rows
  real    :: subparam_gridDesc(20)    ! Input parameter grid desc array
  integer, allocatable  :: lat_line(:,:), lon_line(:,:)
  real,    allocatable  :: dummy_lat(:,:), dummy_lon(:,:), read_slope(:,:)
! ____________________________________________________

  array = LDT_rc%udef

  inquire(file=trim(Noah_struc(n)%slopetypefile), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) "[ERR] SLOPETYPE map ",trim(Noah_struc(n)%slopetypefile)," not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif
  select case ( Noah_struc(n)%slopetype_gridtransform )
    case("none")  
      write(LDT_logunit,*) "[INFO] Reading GFS slopetype file: ",&
            trim(Noah_struc(n)%slopetypefile)
  case default
     write(LDT_logunit,*) "[ERR] only the 'none' spatial transform type is currently allowed."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  end select

! - find subparam grid
   subparam_gridDesc = 0.
   call LDT_RunDomainPts(n,Noah_struc(n)%slopetype_proj,&
             Noah_struc(n)%slopetype_gridDesc(:), &
             glpnc,glpnr,subpnc,subpnr,subparam_gridDesc,lat_line,lon_line)

! - open file  
  ftn = LDT_getNextUnitNumber()
  open(ftn, file=Noah_struc(n)%slopetypefile, access='sequential',status='old', &
       form="unformatted", recl=4)

! - read file
  allocate(dummy_lat(glpnc,glpnr))
  allocate(dummy_lon(glpnc,glpnr))
  allocate(read_slope(glpnc,glpnr))
  read(ftn) dummy_lat
  read(ftn) dummy_lon
  read(ftn) read_slope
  close(ftn)
  deallocate(dummy_lat)
  deallocate(dummy_lon)

!- Transform parameter grid to LIS run domain:
  select case( Noah_struc(n)%slopetype_gridtransform )
    case( "none" )
      do r = 1, subpnr
        do c = 1, subpnc
          array(c,r,1) = read_slope(lon_line(c,r),lat_line(c,r))
        enddo ! columns
      enddo ! rows
      deallocate(read_slope)
   case default
     write(LDT_logunit,*)" This spatial transform, ",&
          trim(Noah_struc(n)%slopetype_gridtransform),&
          ", for GFS slope type is not available at this time ... "
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  end select
  close(ftn)

  do r = 1,LDT_rc%lnr(n)
     do c = 1,LDT_rc%lnc(n)
        if(array(c,r,1).gt.9.) array(c,r,1) = 9.  ! Set slopetype upper limit to 9  
     enddo
  enddo

  call LDT_releaseUnitNumber(ftn)
  write(LDT_logunit, *) "[INFO] Done reading GFS slope type file"

end subroutine read_GFS_slopetype
