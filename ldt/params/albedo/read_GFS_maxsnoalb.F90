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
! !ROUTINE: read_GFS_mxsnoalb
!  \label{read_GFS_mxsnoalb}

! !REVISION HISTORY:
!  05 2014; Grey Nearing; Initial Specification
!
! !INTERFACE:
subroutine read_GFS_mxsnoalb(n,array) 

! !USES:
  use LDT_coreMod,       only : LDT_rc, LDT_domain
  use LDT_logMod,        only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod,     only : readLISdata 
  use LDT_albedoMod
  use LDT_gridmappingmod

!EOP      
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  real, intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
!
! !DESCRIPTION:
!  This subroutine retrieves the maximum snow expected over 
!  deep snow and returns the values in a latlon projection
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[array]
!   output field with the retrieved max snow albedo
!  \end{description}
!
!EOP      
  integer :: ftn
  integer :: c, r
  logical :: file_exists
  integer :: glpnc, glpnr             ! Parameter (global) total columns and rows
  integer :: subpnc, subpnr           ! Parameter subsetted columns and rows
  real    :: subparam_gridDesc(20)    ! Input parameter grid desc array
  integer, allocatable  :: lat_line(:,:), lon_line(:,:)
  real,    allocatable  :: dummy_lat(:,:), dummy_lon(:,:), read_albedo(:,:)
! ____________________________

  array = LDT_rc%udef

  inquire(file=trim(LDT_albedo_struc(n)%mxsnoalbfile), exist=file_exists)
  if(.not.file_exists) then
     write(LDT_logunit,*) "MaxSnowAlb map ",trim(LDT_albedo_struc(n)%mxsnoalbfile)," not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif
  select case ( LDT_albedo_struc(n)%mxsnoalb_gridtransform)
  case("none")
     write(LDT_logunit,*) "[INFO] Reading GFS max snow alb file: ",&
          trim(LDT_albedo_struc(n)%mxsnoalbfile)
  case default
     write(LDT_logunit,*) "[ERR] only the none spatial transform type is currently allowed."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  end select

! - find subparam grid   subparam_gridDesc = 0.
  call LDT_RunDomainPts(n,LDT_albedo_struc(n)%mxsnoalb_proj,&
       LDT_albedo_struc(n)%mxsnoalb_gridDesc, &
       glpnc,glpnr,subpnc,subpnr,subparam_gridDesc,lat_line,lon_line)
  
  ! - open file  
  ftn = LDT_getNextUnitNumber()
  open(ftn, file=LDT_albedo_struc(n)%mxsnoalbfile, &
       access='sequential',status='old', &
       form="unformatted", recl=4)

! - read file
  allocate(dummy_lat(glpnc,glpnr))
  allocate(dummy_lon(glpnc,glpnr))
  allocate(read_albedo(glpnc,glpnr))
  read(ftn) dummy_lat
  read(ftn) dummy_lon
  read(ftn) read_albedo
  close(ftn)
  deallocate(dummy_lat)
  deallocate(dummy_lon)

!- Transform parameter grid to LIS run domain:
  select case (LDT_albedo_struc(n)%mxsnoalb_gridtransform)
    case( "none" )
      do r = 1,subpnr
        do c = 1,subpnc
          array(c,r,1) = read_albedo(lon_line(c,r),lat_line(c,r))
        enddo ! columns
      enddo ! rows
      deallocate(read_albedo)
   case default
     write(LDT_logunit,*) "[ERR] This spatial transform, ",&
          trim(LDT_albedo_struc(n)%mxsnoalb_gridtransform),&
          ", for GFS max snow albedo is not available at this time ... "
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  end select
  close(ftn)

  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        if(array(c,r,1).lt.0) then
           array(c,r,1) = LDT_rc%udef
        endif
     enddo
  enddo
  
  call LDT_releaseUnitNumber(ftn)
  write(LDT_logunit, *) "[INFO] Done reading GFS max snow albedo file"

end subroutine read_GFS_mxsnoalb
