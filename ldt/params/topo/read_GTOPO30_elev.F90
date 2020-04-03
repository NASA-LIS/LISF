!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: read_GTOPO30_elev
! \label{read_GTOPO30_elev}
!
! !REVISION HISTORY:
!  03 Sept 2004: Sujay Kumar; Initial Specification
!  01 Aug  2012: KR Arsenault; Expanded for elevation tiling
!
! !INTERFACE:
subroutine read_GTOPO30_elev( n, num_bins, fgrd, elevave )

! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_gridmappingMod
  use LDT_fileIOMod
  use LDT_paramTileInputMod, only: param_1dbin_areacalc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n  
  integer, intent(in) :: num_bins
  real,    intent(out):: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real,    intent(out):: elevave(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)

! !DESCRIPTION:
!  This subroutine retrieves static elevation data from the GTOPO30 source
!  and reprojects it to the latlon projection. 
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
   integer   :: c, r, t, i, line
   integer   :: subpnc, subpnr, glpnc, glpnr
   integer   :: mi                          ! Total number of input param grid array points
   integer   :: mo                          ! Total number of output LIS grid array points
   real      :: subparam_gridDesc(20)       ! Input parameter grid desc array
   integer, allocatable :: lat_line(:,:), lon_line(:,:)
   integer, allocatable :: n11(:)           ! Maps each input grid point to output grid.
   real,    allocatable :: gi1(:)           ! input parameter 1d grid
   logical*1,allocatable:: li1(:)           ! input logical mask (to match gi)
   real,    allocatable :: read_elev(:,:)   ! Read input parameter
   real      :: elev(LDT_rc%lnc(n),LDT_rc%lnr(n),1)

!________________________________________________________________________

  elev = 0.
  fgrd = 0.
  elevave = 0.
  
  inquire(file=trim(LDT_rc%elevfile(n)), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) "Elevation map ",trim(LDT_rc%elevfile(n))," not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif
  write(LDT_logunit,*) "Reading GTOPO30 (LIS-based) elevation file: ",&
                       trim(LDT_rc%elevfile(n))
  
  ftn = LDT_getNextUnitNumber()
  open(ftn,file=LDT_rc%elevfile(n),form='unformatted', &
       access='direct',convert='big_endian',recl=4,status='old')

! -------------------------------------------------------------------
!     AGGREGATING FINE-SCALE GRIDS TO COARSER LIS OUTPUT GRID
! -------------------------------------------------------------------

!- Transform parameter grid to LIS run domain:
   select case ( LDT_rc%topo_gridtransform(n) )

 !- Transforming 2-D elevation field: 
    case( "none", "neighbor", "average", "bilinear", "budget-bilinear" )

       call readLISdata(n, ftn, LDT_rc%topo_proj, LDT_rc%topo_gridtransform(n), &
                        LDT_rc%topo_gridDesc(n,:), 1, elev )  ! 1 indicates 2D layer

       do r=1,LDT_rc%lnr(n)
          do c=1,LDT_rc%lnc(n)
             do t = 1, num_bins
             !- Single elevation layer, write to first bin:
                if( t == 1 ) then
                   fgrd(c,r,t)    = 1.0 
                   elevave(c,r,t) = elev(c,r,1)
                endif
             end do
          enddo
       enddo


 !- Transforming 3-D elevation field: 
    case( "tile" )

      subparam_gridDesc = 0.
      call LDT_RunDomainPts( n, LDT_rc%topo_proj, LDT_rc%topo_gridDesc(n,:), &
               glpnc, glpnr, subpnc, subpnr, subparam_gridDesc, lat_line, lon_line )

   !- Initialize parameter read-in array:
      allocate( read_elev(subpnc,subpnr) )
      read_elev = LDT_rc%udef
      line = 0
      do r = 1, subpnr
         do c = 1, subpnc
            line = (lat_line(c,r)-1)*glpnc + lon_line(c,r)
            read(ftn,rec=line) read_elev(c,r)
         enddo
      enddo
      mi = subpnc*subpnr
      mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
      allocate( gi1(mi), li1(mi), n11(mi) )
      gi1 = LDT_rc%udef; li1 = .false.

   !- Assign 2-D array to 1-D for aggregation routines:
      i = 0
      do r = 1, subpnr
         do c = 1, subpnc;  i = i + 1
            gi1(i) = read_elev(c,r)
            if( gi1(i) .ne. LDT_rc%udef )  li1(i) = .true.
         enddo
      enddo

   !- Create mapping between parameter domain and LIS grid domain:
      call upscaleByAveraging_input( subparam_gridDesc, LDT_rc%gridDesc(n,:), &
                                     mi, mo, n11 )

      call param_1dbin_areacalc( n, num_bins, mi, mo, n11, &
                                 LDT_rc%udef, gi1, fgrd, elevave )

      deallocate( gi1, li1, n11 )

    case default
       write(LDT_logunit,*) "[ERR] The spatial transform, ",trim(LDT_rc%topo_gridtransform(n)),&
                            ", for elevation is not available at this time ... "
       write(LDT_logunit,*) "Program stopping ..."
       call LDT_endrun
   end select

   call LDT_releaseUnitNumber(ftn)

   write(LDT_logunit, *) "[INFO] Done reading GTOPO30 elevation file."

end subroutine read_GTOPO30_elev

