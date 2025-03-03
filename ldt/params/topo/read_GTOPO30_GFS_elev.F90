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
! !ROUTINE: read_GTOPO30_GFS_elev
! \label{read_GTOPO30_GFS_elev}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  01 Aug  2012: KR Arsenault; Expanded for elevation tiling
!
! !INTERFACE:
subroutine read_GTOPO30_GFS_elev( n, num_bins, fgrd, elevave )

! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_gridmappingMod
  use LDT_fileIOMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n  
  integer, intent(in) :: num_bins
  real,    intent(out):: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real,    intent(out):: elevave(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)

! !DESCRIPTION:
!  This subroutine retrieves static elevation data from the GFS 
!  GTOPO30 source and reprojects it to the latlon projection. 
!
!  Source information:
!   http://webgis.wr.usgs.gov/globalgis/gtopo30/gtopo30.htm
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[num_bins]
!   number of bins (or bands)
!  \item[fgrd]
!   output grid fractions for elevation bins (or bands)
!  \item[elevave]
!   elevation average for each bin (or band)
!  \end{description}
!EOP
   integer   :: ftn
   logical   :: file_exists
   integer   :: c, r
   integer   :: subpnc, subpnr, glpnc, glpnr
   real      :: subparam_gridDesc(20)       ! Input parameter grid desc array
   integer, allocatable :: lat_line(:,:), lon_line(:,:)
   real,    allocatable :: read_elev(:,:), dummy_lat(:,:), dummy_lon(:,:) 
   real      :: elev(LDT_rc%lnc(n),LDT_rc%lnr(n))

!________________________________________________________________________

  elev = 0.
  fgrd = 0.
  elevave = 0.
  
  inquire(file=trim(LDT_rc%elevfile(n)), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) "[ERR] GFS Elevation map ",trim(LDT_rc%elevfile(n))," not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif
  write(LDT_logunit,*) "[INFO] Reading GFS GTOPO30 elevation file: ",&
                       trim(LDT_rc%elevfile(n))
  if (num_bins.ne.1) then 
     write(LDT_logunit,*) "[ERR] Elevation map for GFS must be 2-D"
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif

  ! open file 
  ftn = LDT_getNextUnitNumber()
  open(ftn,file=LDT_rc%elevfile(n),form='unformatted', &
       access='sequential',recl=1,status='old',convert='big_endian')

  subparam_gridDesc = 0.
  call LDT_RunDomainPts(n,LDT_rc%topo_proj,LDT_rc%topo_gridDesc(n,:), &
                         glpnc,glpnr,subpnc,subpnr,  &
                         subparam_gridDesc,lat_line,lon_line )

  allocate(dummy_lat(glpnc,glpnr))
  allocate(dummy_lon(glpnc,glpnr))
  allocate(read_elev(glpnc,glpnr))
  
  read(ftn) dummy_lat
  read(ftn) dummy_lon
  read(ftn) read_elev
      
  close(ftn) 
  deallocate(dummy_lat)
  deallocate(dummy_lon)
! -------------------------------------------------------------------
!     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
! -------------------------------------------------------------------
  select case (LDT_rc%topo_gridtransform(n))
    case("none")
      do r = 1,subpnr
        do c = 1,subpnc
           elev(c,r) = read_elev(lon_line(c,r),lat_line(c,r))
        enddo ! columns
      enddo ! rows
      deallocate(read_elev)
  
       do r = 1,LDT_rc%lnr(n)
         do c = 1,LDT_rc%lnc(n)
           fgrd(c,r,1)    = 1.0 
           elevave(c,r,1) = elev(c,r)
         enddo
       enddo

    case default
       write(LDT_logunit,*) "[ERR] The spatial transform, ",&
                            trim(LDT_rc%topo_gridtransform(n)),&
                            ", for GFS elevation is not available at this time ... "
       write(LDT_logunit,*) "Program stopping ..."
       call LDT_endrun
   end select

   call LDT_releaseUnitNumber(ftn)
   write(LDT_logunit, *) "[INFO] Done reading GFS GTOPO30 elevation file."

end subroutine read_GTOPO30_GFS_elev

