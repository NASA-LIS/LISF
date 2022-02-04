#define VERIFY_(A)   IF(A/=0)THEN;PRINT *,'ERROR AT LINE ', __LINE__;STOP;ENDIF
#define ASSERT_(A)   if(.not.A)then;print *,'Error:',__FILE__,__LINE__;stop;endif

! -----------------------------------------------------------------
#include "LIS_misc.h"
module getCropIrrigTypes

  use ESMF
  use LIS_coreMod
  use LIS_numerRecipesMod, ONLY : LIS_indexArrayReal, LIS_reverse, LIS_quicksort
  use LIS_mpiMod
  use LIS_logMod,ONLY : LIS_verify
  
  implicit none

  PRIVATE
  
  type, public :: itype_params

     REAL                 :: ITYPE_MIN_FRAC = 0.1
     REAL                 :: CTYPE_AREA_TOL = 2.
     INTEGER              :: N_CROPTYPES  = 26
     INTEGER              :: N_IRRIGTYPES = 3
     REAL, allocatable    :: CROP_IRRIGTYPE (:,:,:)
     INTEGER, pointer     :: COUNTRY_CROP_IRRAVAIL (:)
     integer, allocatable :: crop_rank(:,:)

  end type itype_params
     
  type, public, extends (itype_params) :: matt_algorithm

   contains
     
     procedure, public  :: init_thres
     procedure, public  :: git => get_irrig_type
     procedure, private :: fpt => find_preferred_type
     procedure, private :: dti => dissolve_tiny_itypes
     procedure, private :: lct => least_common_type     
     procedure, private :: hcr => highest_crop_rank
     procedure, private :: ait => assign_irrig_type
     procedure, private :: fpt_two => find_preferred_types_fortwo
     
  end type matt_algorithm

  ! ---------------------------------------------------------  

contains

  
  ! ---------------------------------------------------------    

  SUBROUTINE init_thres (ma, ncrops, nitypes, global)

    use LIS_fileIOMod  
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    use LIS_fileIOMod
    
    implicit none

    class(matt_algorithm), intent (inout) :: ma
    type (itype_params)                   :: ip
    integer                               :: rc, i, j
    integer, intent (in), optional        :: ncrops, nitypes
    logical, intent (in), optional        :: global
    character*300                         :: crop_irrig_file
    integer                               :: ios, nid, vid, n_country, g_ncrops, g_nitypes
    allocate (ma%crop_rank(1:ncrops, 1:nitypes))
    call init_crop_rank(ma%crop_rank)
    
    if (present (global)) then
       call ESMF_ConfigGetAttribute(LIS_config, LABEL = 'ITYPE_MIN_FRAC_GLOBAL:', VALUE = ma%ITYPE_MIN_FRAC, DEFAULT = ip%ITYPE_MIN_FRAC, rc = rc)
       call ESMF_ConfigGetAttribute(LIS_config, LABEL = 'CTYPE_AREA_TOL_GLOBAL:', VALUE = ma%CTYPE_AREA_TOL, DEFAULT = ip%CTYPE_AREA_TOL, rc = rc)
    else
       call ESMF_ConfigGetAttribute(LIS_config, LABEL = 'ITYPE_MIN_FRAC_USA:', VALUE = ma%ITYPE_MIN_FRAC, DEFAULT = ip%ITYPE_MIN_FRAC, rc = rc)
       call ESMF_ConfigGetAttribute(LIS_config, LABEL = 'CTYPE_AREA_TOL_USA:', VALUE = ma%CTYPE_AREA_TOL, DEFAULT = ip%CTYPE_AREA_TOL, rc = rc)
    endif
    
    if(present (ncrops )) ma%N_CROPTYPES  = ncrops
    if(present (nitypes)) ma%N_IRRIGTYPES = nitypes

    if (present (global)) then
       
       call ESMF_ConfigGetAttribute(LIS_config, LABEL = 'Country level crop irrigtype file:', VALUE = crop_irrig_file, rc = rc)
       call LIS_verify(rc,'ERROR crop irrigtype file not defined')
       
#if (defined USE_NETCDF3 || defined USE_NETCDF4)

       ios = nf90_open(path=trim(crop_irrig_file), mode=NF90_NOWRITE,ncid=nid) ; call LIS_verify(ios,'ERROR opening '//trim(crop_irrig_file))
       ios = nf90_inq_dimid        (nid, "country", vid)       ; call LIS_verify(ios,'Problem with country dimension id.')
       ios = nf90_inquire_dimension(nid, vid, len = n_country) ; call LIS_verify(ios,'Problem with country dimension.')
       ios = nf90_inq_dimid        (nid, "crops", vid)         ; call LIS_verify(ios,'Problem with crops dimension id.')
       ios = nf90_inquire_dimension(nid, vid, len = g_ncrops)  ; call LIS_verify(ios,'Problem with crops dimension.')
       ios = nf90_inq_dimid        (nid, "irrigtype", vid)     ; call LIS_verify(ios,'Problem with irrigtypes dimension id.')
       ios = nf90_inquire_dimension(nid, vid, len = g_nitypes) ; call LIS_verify(ios,'Problem with irrigtypes dimension.')

       call LIS_verify(g_ncrops - ma%N_CROPTYPES ,'Number of croptypes in crop_irrig_file does not match.' )
       call LIS_verify(g_nitypes -ma%N_IRRIGTYPES,'Number of irrigtypes in crop_irrig_file does not match.')

       allocate (ma%CROP_IRRIGTYPE (1: g_ncrops, g_nitypes,1: n_country))
       
       ios = nf90_inq_varid(nid,'CROP_IRRIGTYPE',vid) ; call LIS_verify(ios,'nf90_inq_varid failed for CROP_IRRIGTYPE')
       ios = nf90_get_var(nid,vid, ma%CROP_IRRIGTYPE) ; call LIS_verify(ios,'nf90_get_var failed for CROP_IRRIGTYPE'  )

       ios = nf90_close(nid)

#endif           

    endif

  contains
    SUBROUTINE init_crop_rank  (crop_rank)

    integer,dimension(:,:),  intent(inout) :: crop_rank
    
     ! From Matt's Excel sheet
    crop_rank (:,1) = (/ 7, 5,26,11,12,13, 6,15,16, 2, 4,17, 3,22,10, 1, 8,25,21,19,14,24,23,20,18, 9/)
    crop_rank (:,2) = (/17,18,26,23,21,22,16,24, 9,13,15,10,14, 4,20,25,11, 7, 3, 1,12, 6, 5, 2,19, 8/)
    crop_rank (:,3) = (/12,14, 1, 8, 9, 7,13, 5, 4,26,24, 3,25,18,10,23,11,15,19,21, 6,16,17,20, 2,22/)
    
  end subroutine init_crop_rank
    
  END SUBROUTINE init_thres

  ! ---------------------------------------------------------

  SUBROUTINE get_irrig_type (this, GNC, GNR, CELLAREA, ADMU_GRID, CROPTYPE_IN, IRRIGFRAC, PREFTYPE, usa)

    !A.	For each model grid pixel: 
    !1.	Determine the vegetation/crop types of the tiles (LDT)
    !2.	Assign each pixel to a county, based on whatever county makes up the highest fractional area of the pixel (COUNTY) (LDT).
    !B.	For each county 
    !1.	Calculate the fraction of each irrigation type (IRRITYPE).
    !2.	Start with the least common irrigation type (**  least_common_type **)
    !3.	If that irrigation type covers less than 10% of the area of the county divided by the number of pixels in the county, then ignore it (for example, county area = 6000 km2, drip fraction = 0.008, number of pixels in county = 10; then threshold = 6000 km2 * 10% / 10 pixels = 60 km2; drip area = 0.008 * 6000 km2 = 48 km2; therefore ignore drip) (LDT). (** dissolve_tiny_itypes **)
    !4.	For the least common irrigation type that passes the test (#2):
    !5.	Use the data in the attached spreadsheet to determine the crop with the highest Crop Rank corresponding to that irrigation type that also exists in the county (after tiles have been assigned).  
    !6.	Does the area of the largest tile of that crop type in the county exceed the area of that irrigation type in the county by more than 2% of the county area?  If no, then assign that irrigation type to the tile and, if not within 2% (county area) of the irrigation type area target, repeat with the next largest tile of that type.  If yes, then go to the second largest tile of that crop type in the county and repeat the test.  Continue until the assigned area of that irrigation t1ype is within +- 2% of the desired fractional area.  If that target cannot be achieved with the tile of this crop type, then go to the next highest ranked (Crop Rank) tile for the irrigation type and repeat this step (5).
    !7.	If all tiles in the county are tested before the +-2% goal is achieved for the irrigation type, then move to step 7, keeping any tile assignments that were made in step 5.
    !8.	Move to the second most common irrigation type in the county and repeat step 5, being sure not to re-use any tiles already assigned.
    !9.	All remaining tiles in the county are assigned to the most common irrigation type for that crop type (Irr Rank in spreadsheet).
    !10.Return to step 1 for the next county.

    implicit none
    
    class(matt_algorithm), intent(inout) :: this
    integer, intent (in)                 :: GNC, GNR       ! number of columns/rows in LIS global domain
    real, dimension (:,:), intent (in)   :: ADMU_GRID      ! COUNTRY/COUNTY grid in the local domain
    real, dimension (:,:), intent (in)   :: CELLAREA      ! LIS grid cell area
    real, dimension (:,:,:), intent (in) :: CROPTYPE_IN    ! land cover fractions in crop tiles
    real, dimension (:,:,:), intent (in) :: IRRIGFRAC      ! irrigtype fractions  1=Sprinkler, 2=Drip, 3=Floodg 
    real, dimension (:,:,:),intent(inout):: PREFTYPE       ! preferred Irrig type for each CROPTYPE_IN
    logical, intent (in), optional       :: USA

    real, allocatable, dimension (:,:,:) :: CROPTYPE
    integer, allocatable, dimension(:)   :: NCELLS_PROC, displ, buf_cnt
    logical            :: isusa = .false.
    integer            :: DX, DY, NX, NY, C1, C2, R1,R2, n, i_cnt, i,j,NCNTY, NPLUS, &
         NCNTY_MAX,n_AdminUnit, NT_CNTY, nt, t1, t2
    integer            :: comm_rank, comm_size, status, mpistatus(MPI_STATUS_SIZE),INFOS, req,index_of_country,n_lctypes
    logical            :: root_proc = .true., run_alg
    real, allocatable, dimension (:,:):: RDUMMY, GDUMMY, G_COUNTY, CNTY_CROPTYPE, CNTY_IRRIGTYPE, CNTY_IRRIGFRAC
    real, allocatable, dimension (:)  :: CNTY_CELLAREA
    real, ALLOCATABLE, DIMENSION (:)  :: loc_int, AdminUnit, cnty_var, cnty_var1,cnty_var2, AIRRIGFRAC
    logical, dimension(:),allocatable :: unq_mask
    real                              :: thisAdmUnit, itype_thres, county_area, local_area, total_itype_area
    integer, allocatable, dimension (:,:), target :: COUNTRY_CROP_IRRAVAIL
    
    DX = SIZE (ADMU_GRID,1)
    DY = SIZE (ADMU_GRID,2)
     
    call MPI_COMM_Size(LIS_mpi_comm,comm_size,status)  ; VERIFY_(STATUS)
    call MPI_COMM_Rank(LIS_mpi_comm,comm_rank,status)  ; VERIFY_(STATUS)
    if (comm_rank /= 0) root_proc = .false.
    
    NCNTY_MAX = 0
    isusa = .false.
    
    if(present (usa)) then
       ! check whether global domain contains any US county
       isusa = .true.
       NCNTY     = 0
       NCNTY_MAX = 0
       if(isusa) then
          NCNTY = COUNT (ADMU_GRID > 0.)
       else
          NCNTY = COUNT (ADMU_GRID >= 0.)
       endif

       ! Check whether US counties actually exist within the global domain
       call MPI_AllReduce(NCNTY,NCNTY_MAX,1,MPI_INTEGER,MPI_MAX,LIS_mpi_comm,STATUS)
       VERIFY_(STATUS)
       
       if (NCNTY_MAX == 0) then
          ! Leave SUBROUTINE get_irrig_type because there is nothing to do in the US.
          RETURN
       endif
              
    endif

    NX = GNC / DX
    NY = GNR / DY
    
    allocate (RDUMMY   (1: DX, 1: DY))
    ! temperoroy array that gets continually adjusted during irrigation type assigning 
    allocate (CROPTYPE (1: DX, 1: DY, 1: this%N_CROPTYPES))
    CROPTYPE = CROPTYPE_IN

    N_LCTYPES = size (PREFTYPE, 3) - this%N_CROPTYPES
           
    ! STEP 1: root_proc gathers available countries (or US counties) from the global domain.
    ! --------------------------------------------------------------------------------------

    if(root_proc) then
       allocate (GDUMMY      (1: GNC, 1: GNR))
       GDUMMY = -9999.
       if(.not.isusa) allocate (COUNTRY_CROP_IRRAVAIL (1:this%N_CROPTYPES, 1:this%N_IRRIGTYPES))
    endif
       
    PROC_LOOP_1: do j = 1, comm_size
          
       N = CEILING(REAL(j)/REAL(NX))
       i = j - (N-1)*NX
       
       ! Below C1, C2 and R1, R2 are minimum and maximum global i,j indices
       ! for the rectangular region in jth processor respectively.
       c1 = (i-1)*DX + 1
       c2 = i*DX
       R1 = (N-1)*DY + 1
       R2 = N*DY
       
       if (comm_rank == 0) GDUMMY (1:DX, 1:DY) = ADMU_GRID
       if (j > 1) then
          if(j-1 == comm_rank) then
             ! Sends ADMU_GRID(DX,DY) to root
             call MPI_ISend(ADMU_GRID, DX*DY,MPI_REAL,0,990,LIS_mpi_comm,req,STATUS); VERIFY_(STATUS) 
             call MPI_WAIT (req,MPI_STATUS_IGNORE,STATUS) ; VERIFY_(STATUS)
          else if (comm_rank == 0) then
                ! root receives ADMU_GRID from the jth processor
             call MPI_RECV(GDUMMY (c1:c2,R1:R2),DX*DY , MPI_REAL, j-1,990,LIS_mpi_comm,MPI_STATUS_IGNORE,STATUS); VERIFY_(STATUS)
          endif
       endif
    end do PROC_LOOP_1
    
    if(root_proc) then
       if(isusa) then
          nplus = count(GDUMMY > 0.)
       else
          nplus = count(GDUMMY >= 0.)
       endif
       allocate (loc_int (1:NPLUS))
       allocate (unq_mask(1:NPLUS))
       if(isusa) then
          loc_int = pack(GDUMMY,mask = (gdummy > 0)) ! loc_int contains available COUNTY IDs
       else
          loc_int = pack(GDUMMY,mask = (gdummy >= 0)) ! loc_int contains available COUNTRY IDs
       endif
       
       call LIS_quicksort (loc_int)

       unq_mask = .true.

       do n = 2,NPLUS 
          unq_mask(n) = .not.(loc_int(n) == loc_int(n-1)) ! to count number of unique numbers in loc_int 
       end do
       n_AdminUnit = count(unq_mask)
       allocate(AdminUnit (1:n_AdminUnit))
       AdminUnit = pack(loc_int,mask =unq_mask)
    endif

    call MPI_BCAST (n_AdminUnit, 1, MPI_INTEGER, 0,LIS_mpi_comm, STATUS); VERIFY_(STATUS)
    
    allocate (NCELLS_PROC (1: comm_size))
    allocate (DISPL       (1: comm_size))
    allocate (BUF_CNT     (1: comm_size))
    
    do n = 1, comm_size
       DISPL (n) = n-1
    end do
    
    BUF_CNT = 1
    N = comm_rank + 1
    j = CEILING(REAL(N)/REAL(NX))
    i = N - (j-1)*NX
    
    ! Below c1,c2, and r1,r2 are column,row indices of the global array
    ! -----------------------------------------------------------------
    
    c1 = (i-1)*DX + 1
    c2 = i*DX
    R1 = (j-1)*DY + 1
    R2 = j*DY
    
  COUNTY_LOOP: do i_cnt = 1, n_AdminUnit

     run_alg     = .false.
     local_area  = 0.
     nplus       = 0
     county_area = 0.
     NCELLS_PROC = 0

     if(root_proc) thisAdmUnit = AdminUnit(i_cnt)
     call MPI_BCAST (thisAdmUnit, 1, MPI_REAL, 0,LIS_mpi_comm, STATUS); VERIFY_(STATUS)
         
     ! compute sum of cell areas that make up the county fraction in the local processor

     nplus     = count (ADMU_GRID == thisAdmUnit)

     if(NPLUS > 0)  then
        allocate (cnty_var (1:NPLUS))
        allocate (cnty_var1(1:NPLUS))
        allocate (cnty_var2(1:NPLUS))
        cnty_var = pack (cellarea, mask = (ADMU_GRID == thisAdmUnit))
        
        do n = 1, nplus
           local_area = local_area + cnty_var(n)
        end do
     endif

     ! total county area and number of grid cells provided by each processor to form the county

     call MPI_Reduce(local_area,county_area,1,MPI_REAL,MPI_SUM,0,LIS_mpi_comm,STATUS)
     call MPI_ALLGATHERV(nplus, BUF_CNT, MPI_integer, NCELLS_PROC, BUF_CNT,displ, MPI_integer, &
          LIS_mpi_comm, status ) ; VERIFY_(STATUS)

     NT_CNTY = SUM (NCELLS_PROC)
     if(NT_CNTY > 0) run_alg = .true.
     
     call MPI_Barrier(LIS_mpi_comm, STATUS)
     
     ! (4) Assign PREFTYPE for CROPTYPES
      
     RUN_MATALG: if (run_alg) then
        
        if(root_proc) then
           
           allocate (CNTY_CROPTYPE  (1: NT_CNTY, 1 : this%N_CROPTYPES))
           allocate (CNTY_IRRIGTYPE (1: NT_CNTY, 1 : this%N_CROPTYPES))
           allocate (CNTY_IRRIGFRAC (1: NT_CNTY, 1 : this%N_IRRIGTYPES))
           allocate (CNTY_CELLAREA  (1: NT_CNTY))
           allocate (AIRRIGFRAC (1 : this%N_IRRIGTYPES))
           CNTY_IRRIGTYPE = -1.
           CNTY_IRRIGFRAC = 0.
           CNTY_CELLAREA  = 0.
           AIRRIGFRAC     = 0.
                                 
        endif
        
        GATHER_FROM_PROCS: do j = 1, comm_size
           CHECK_ACTIVE_PROCS: if(NCELLS_PROC(j) > 0) then
              NUM_CROPS1: do nt =  1, this%N_CROPTYPES

                 ! 1) Gather CROPTYPE
                 ! ------------------
                 RDUMMY = CROPTYPE (:,:,Nt)
                 cnty_var = pack (RDUMMY, mask = (ADMU_GRID == thisAdmUnit))

                 if (nt <= this%N_IRRIGTYPES) then

                    ! 2) Gather IRRIGFRAC
                    ! -------------------
                    RDUMMY = IRRIGFRAC (:,:,Nt)
                    cnty_var1 = pack (RDUMMY, mask = (ADMU_GRID == thisAdmUnit))

                 endif
                 
                 if (nt == 1) then

                    ! 3) Gather CELLAREA
                    ! ------------------
                    RDUMMY = CELLAREA (:,:)
                    cnty_var2 = pack (RDUMMY, mask = (ADMU_GRID == thisAdmUnit))

                 endif
 
                 ! if root is active itself
                 if((comm_rank == 0).and.(nplus > 0)) then
                    CNTY_CROPTYPE (1:NPLUS,NT) = cnty_var
                    if (nt <= this%N_IRRIGTYPES) CNTY_IRRIGFRAC (1:NPLUS,NT) = cnty_var1
                    if (nt == 1)                 CNTY_CELLAREA  (1:NPLUS)    = cnty_var2 
                 endif
                 
                 
                 if (j > 1) then                    
                    ! inter processor communication                    
                    if(j-1 == comm_rank) then
                       ! send to root
                       call MPI_ISend(cnty_var, nplus,MPI_REAL,0,991,LIS_mpi_comm,req,STATUS); VERIFY_(STATUS) 
                       call MPI_WAIT (req,MPI_STATUS_IGNORE,STATUS) ; VERIFY_(STATUS)
                       if (nt <= this%N_IRRIGTYPES) then
                          call MPI_ISend(cnty_var1, nplus,MPI_REAL,0,994,LIS_mpi_comm,req,STATUS); VERIFY_(STATUS) 
                          call MPI_WAIT (req,MPI_STATUS_IGNORE,STATUS) ; VERIFY_(STATUS)
                       endif
                       if (nt == 1) then
                          call MPI_ISend(cnty_var2, nplus,MPI_REAL,0,995,LIS_mpi_comm,req,STATUS); VERIFY_(STATUS) 
                          call MPI_WAIT (req,MPI_STATUS_IGNORE,STATUS) ; VERIFY_(STATUS)
                       endif
                    else if (comm_rank == 0) then
                       ! root receives
                       t1 = SUM (NCELLS_PROC(1:j-1)) + 1
                       t2 = SUM (NCELLS_PROC(1:j))
                       call MPI_RECV(CNTY_CROPTYPE (t1:t2, nt), NCELLS_PROC(j), MPI_REAL, j-1,991,LIS_mpi_comm,MPI_STATUS_IGNORE,STATUS); VERIFY_(STATUS)
                       if (nt <= this%N_IRRIGTYPES) then
                          call MPI_RECV(CNTY_IRRIGFRAC (t1:t2, nt), NCELLS_PROC(j), MPI_REAL, j-1,994,LIS_mpi_comm,MPI_STATUS_IGNORE,STATUS); VERIFY_(STATUS)
                       endif
                       if (nt == 1) then
                          call MPI_RECV(CNTY_CELLAREA (t1:t2),      NCELLS_PROC(j), MPI_REAL, j-1,995,LIS_mpi_comm,MPI_STATUS_IGNORE,STATUS); VERIFY_(STATUS)
                       endif
                    endif
                 endif
                 
              end do NUM_CROPS1
           endif CHECK_ACTIVE_PROCS
        end do GATHER_FROM_PROCS
        
        ! Run Matt's algorithm at ROOT
        if(root_proc) then

           if(.not.isusa) then
              index_of_country = thisAdmUnit + 1 ! FINDLOC (this%COUNTRY_CODE, NINT(thisAdmUnit))
           endif
           
           ! County IRRIGFRAC
           do j = 1, this%N_IRRIGTYPES
              AIRRIGFRAC(j) = SUM (CNTY_IRRIGFRAC (:,J))/NT_CNTY
!              if (AIRRIGFRAC(j) /= CNTY_IRRIGFRAC (1,j)) print *, 'IRRIGFRAC should have been uniform within the county, please check', AIRRIGFRAC(j), CNTY_IRRIGFRAC (1,j), NT_CNTY
           end do

           !(3). that irrigation type covers less than 10%  (ITYPE_MIN_FRAC) of the area of the county divided by the number of pixels in the county, then ignore it
           ! (for example, county area = 6000 km2, drip fraction = 0.008, number of pixels in county = 10; then threshold = 6000 km2 * 10% / 10 pixels = 60 km2;
           ! drip area = 0.008 * 6000 km2 = 48 km2; therefore ignore drip). (** dissolve_tiny_itypes **)
           
           itype_thres = this%ITYPE_MIN_FRAC * county_area/NT_CNTY /county_area
           
           ! Use itype_thres to dissolve tiny irrigtype fractions 
           call this%dti (itype_thres, AIRRIGFRAC)
           ! CHANGE CHANGE CHANGE
           if(count (AIRRIGFRAC > 0) == 1.) then
              ! Single irrig type in the county
              where (CNTY_CROPTYPE > 0.)
                 CNTY_IRRIGTYPE = maxloc (AIRRIGFRAC,1)
              endwhere
           else
              ! more than 1 IRRIGFRACS (irrigation types) passed threshold test
              ! Inputs: CNTY_IRRIGFRAC(NOF_LIS_CELLS_IN_COUNTY,3), CNTY_CROPTYPE (NOF_LIS_CELLS_IN_COUNTY, 26)
              ! OutPut: CNTY_IRRIGTYPE(NOF_LIS_CELLS_IN_COUNTY, 26) - then that will be scattered back to proessors
              !         to populate PREFTYPE (DX, DY, LIS_rc%nsurfacetypes) DX and DY are cols and rows in locally in processors

              ! CNTY_IRRIGFRAC contains area [km2] of irrigtypes within the tile/grid cell whose sum is equal to the sum of croptype in that tile/grid cell.
              ! For a given tile 'j', scale to ensure SUM (CNTY_IRRIGTYPE)  = SUM (CNTY_CROPTYPE)
              do j = 1, NT_CNTY
                 CNTY_IRRIGFRAC (j,:) = AIRRIGFRAC * SUM (CNTY_CROPTYPE (j,:))*CNTY_CELLAREA(j) ! change  km^2
                 CNTY_CROPTYPE  (j,:) = CNTY_CROPTYPE (j,:)*CNTY_CELLAREA(j)                    ! change  km^2
              end do
              
              if(isusa) then
                 NT = count (AIRRIGFRAC > 0)
                 ! A special case for single-tile counties with only 2 reported irrigtation types at 50-50 split
                 if ((NT_CNTY==1).AND.(NT == 2).AND.(SUM(AIRRIGFRAC)/2. == MAXVAL (AIRRIGFRAC))) THEN
                    call this%fpt_two (CNTY_CROPTYPE(1,:), CNTY_IRRIGFRAC(1,:), CNTY_IRRIGTYPE(1,:))
                    AIRRIGFRAC = 0.
                 endif
              endif
              
              NT = count (AIRRIGFRAC > 0)
              do j = 1, NT

                 ! 4.For the least common irrigation type that passes the test (#2):
                 ! We call Matt's algorytithm here.
                 i = this%lct (AIRRIGFRAC)

                 !8. (loop through available irrigation types)
                 !   Move to the second most common irrigation type in the county and
                 !   repeat step 5, being sure not to re-use any tiles already assigned.
                 
                 if (j < NT) then
                    !6.	Does the area of the largest tile of that crop type in the county
                    !  exceed the area of that irrigation type in the county by more than
                    ! CTYPE_AREA_TOL% of the county area? 
                    total_itype_area = SUM(CNTY_IRRIGFRAC(:,i)) ! [km^2]  
                    if(isusa) then
                       call this%fpt (i, NT_CNTY, CNTY_CROPTYPE, CNTY_IRRIGTYPE, total_itype_area, first = .true.)                       
                    else
                       if(SUM (this%CROP_IRRIGTYPE (:,i,index_of_country)) >  0) &
                            call this%fpt (i, NT_CNTY, CNTY_CROPTYPE, CNTY_IRRIGTYPE, total_itype_area, first = .true., IRRAVAIL=  NINT(this%CROP_IRRIGTYPE (:,i,index_of_country)))
                    endif
                    AIRRIGFRAC (i) = 0.
                 else
                    !9.	All remaining tiles in the county are assigned to the
                    !   most common irrigation type for that crop type (Irr Rank in spreadsheet).
                    where (CNTY_CROPTYPE > 0.)
                       CNTY_IRRIGTYPE = maxloc (AIRRIGFRAC,1)
                    endwhere
                 endif
              end do              
           endif
        endif

        ! Unpack CNTY_IRRIGTYPE, populate 2D array RDUMMY and then global PREFTYPE

        UNPACK_TO_PROCS: do j = 1, comm_size
           CHECK_ACTIVE_PROCS2: if(NCELLS_PROC(j) > 0) then

              N = CEILING(REAL(j)/REAL(NX))
              i = j - (N-1)*NX
              ! Below C1, C2 and R1, R2 are minimum and maximum global i,j indices
              ! for the rectangular region in jth processor respectively.
              c1 = (i-1)*DX + 1
              c2 = i*DX
              R1 = (N-1)*DY + 1
              R2 = N*DY
              t1 = SUM (NCELLS_PROC(1:j-1)) + 1
              t2 = SUM (NCELLS_PROC(1:j))

              NUM_CROPS2: do nt =  1, this%N_CROPTYPES
                 RDUMMY = 0 
                 ! if root is active itself 
                 if ((comm_rank == 0).AND. (NPLUS > 0)) then
                    RDUMMY = UNPACK(CNTY_IRRIGTYPE(1:NCELLS_PROC(1),nt),mask = (ADMU_GRID == thisAdmUnit), field=-1.)
                    where (RDUMMY > 0)
                       PREFTYPE (:, :, N_LCTYPES + NT) = RDUMMY
                    end where
                 endif

                 if (j > 1) then
                    ! inter processor communication
                    if(j-1 == comm_rank) then
                       
                       !   1st receives the subset of CNTY_IRRIGTYPE  from root
                       call MPI_RECV(cnty_var,NPLUS, MPI_REAL,0,992,LIS_mpi_comm,MPI_STATUS_IGNORE,STATUS); VERIFY_(STATUS) 
                       
                       ! unpack and send back the 2D array to root
                       RDUMMY = UNPACK(cnty_var,mask = (ADMU_GRID == thisAdmUnit), field=-1.)
                       where (RDUMMY > 0)
                          PREFTYPE (:, :, N_LCTYPES + NT) = RDUMMY
                       end where

                    else if (comm_rank == 0) then
                       !  root sends jth portion of CNT_IRRIGTYPE to the jth processor
                       call MPI_ISend(CNTY_IRRIGTYPE(t1:t2, nt), NCELLS_PROC(j),MPI_REAL,j-1,992,LIS_mpi_comm,req,STATUS); VERIFY_(STATUS) 
                       call MPI_WAIT (req,MPI_STATUS_IGNORE,STATUS) ; VERIFY_(STATUS)    

                    endif
                 endif                 
              end do NUM_CROPS2
           endif CHECK_ACTIVE_PROCS2
           call MPI_Barrier(LIS_mpi_comm, STATUS)
        end do UNPACK_TO_PROCS
                
        if(root_proc) then
           if (allocated (CNTY_CROPTYPE)) then 
              deallocate (CNTY_CROPTYPE, CNTY_IRRIGTYPE, CNTY_IRRIGFRAC)
              deallocate (CNTY_CELLAREA, AIRRIGFRAC)
           endif
        endif
     endif RUN_MATALG

     if(allocated(cnty_var )) deallocate (cnty_var)
     if(allocated(cnty_var1)) deallocate (cnty_var1)
     if(allocated(cnty_var2)) deallocate (cnty_var2)
     if(allocated(COUNTRY_CROP_IRRAVAIL)) deallocate (COUNTRY_CROP_IRRAVAIL)
     
  end do COUNTY_LOOP

  call MPI_Barrier(LIS_mpi_comm, STATUS)
  deallocate (RDUMMY, NCELLS_PROC, DISPL, BUF_CNT)   

END SUBROUTINE get_irrig_type
  
  ! ---------------------------------------------------------  

  RECURSIVE SUBROUTINE find_preferred_type( this, this_lct, NT_CNTY, &
       CROPTYPE, PREFTYPE, total_itype_area, first, IRRAVAIL)

    implicit none

    ! INPUTS
    ! ------
    ! this_lct : IRRIGATION TYPE that is being processed.
    ! NT_CNTY  : Number of tiles in the county
    ! CROPTYPE (NT_CNTY, N_CROPTYPES) : CROP TYPES fractions
    ! IRRIGTYPE(NT_CNTY)              : this_lct fraction in each tile

    ! OUTPUTS
    ! -------
    ! PREFTYPE (NT_CNTY, N_CROPTYPES) : preferred irrigation type 
    
    class(matt_algorithm),  intent(inout) :: this
    real,    dimension(:,:),intent(inout) :: PREFTYPE
    real,    dimension(:,:),intent(inout) :: CROPTYPE
    integer, intent (in)                  :: this_lct, NT_CNTY
    real,    intent (inout)               :: total_itype_area
    logical, optional, intent (in)        :: first
    integer, optional, dimension(:), target, intent(in) :: IRRAVAIL
    real,    dimension (NT_CNTY    , this%N_CROPTYPES),      target :: CROPTYPE_TMP, PREFTYPE_TMP
    integer, allocatable                  :: cr_lct (:)
    real, save                            :: original_area
    integer    :: hcr, hcr_nxt, i
    logical                               :: look_nearest
    real, dimension (:), allocatable      :: tile_size
    
    if (present (first)) then
       original_area = total_itype_area
    else
       if (abs(total_itype_area)/original_area < this%CTYPE_AREA_TOL/100.)  RETURN
    endif
    
    !5.	Use the data in the attached spreadsheet to determine the crop
    !   with the highest Crop Rank (CROP_RANK) corresponding to that irrigation.
    allocate (cr_lct (this%N_CROPTYPES))
    cr_lct = this%crop_rank (:,this_lct)
    if(present (IRRAVAIL)) then
       this%COUNTRY_CROP_IRRAVAIL => IRRAVAIL
    endif
    hcr = this%hcr (cr_lct,CROPTYPE)    ! highest crop rank for this_lct
    
    if(hcr > 0) then

       if(.not.present (IRRAVAIL)) then
          if(present (first)) then
             if((NT_CNTY > 1).AND.(NT_CNTY < 5)) then
                if(count(CROPTYPE > 0.)/NT_CNTY == 1) then                
                   ! A single crop type is grown in the entire county. Instead the processing in the ascending order of tile area, select the county 
                   ! whose irrigaated crop area is closest to the area of this_lct.
                   allocate (tile_size (1: NT_CNTY))
                   tile_size = ABS(CROPTYPE(:,hcr) - total_itype_area)
                   PREFTYPE (minloc (tile_size,1),hcr) = this_lct
                   total_itype_area = total_itype_area - CROPTYPE (minloc (tile_size,1),hcr)
                   CROPTYPE (minloc (tile_size,1),hcr) = 0.
                   deallocate (tile_size)
                   return
                endif
             endif
          endif
       endif
       
       !6. Does the area of the largest tile of that crop type in the county exceed the area of that irrigation type in the county by
    !   more than 2% (CTYPE_AREA_TOL) of the county area?
       call this%ait(this_lct, total_itype_area, CROPTYPE(:,hcr), PREFTYPE(:,hcr))
    else
       RETURN
    endif
    
    if(.not.present (IRRAVAIL)) then
       if(present (first)) then
          if(total_itype_area == original_area) then
             call this%ait(this_lct, total_itype_area, CROPTYPE(:,hcr), PREFTYPE(:,hcr), second_attempt = .true.)
             return
          endif
       endif
    endif
    
    !   If that target cannot be achieved with the tile of this crop type, then go to
    !   the next highest ranked (Crop Rank) tile for the irrigation type and repeat this step (5).
    if (abs(total_itype_area)/original_area > this%CTYPE_AREA_TOL/100.) then
       CROPTYPE_TMP = CROPTYPE
       CROPTYPE_TMP(:,hcr) = 0.
       PREFTYPE_TMP = PREFTYPE
       hcr_nxt = this%hcr (cr_lct,CROPTYPE)
       if(hcr_nxt > 0) then
          if (present (IRRAVAIL))then 
             call this%fpt (this_lct, NT_CNTY, CROPTYPE_TMP, PREFTYPE_TMP, total_itype_area, IRRAVAIL=IRRAVAIL)
          else
             call this%fpt (this_lct, NT_CNTY, CROPTYPE_TMP, PREFTYPE_TMP, total_itype_area)
          endif
          CROPTYPE_TMP(:,hcr) = CROPTYPE (:,hcr)
          PREFTYPE_TMP(:,hcr) = PREFTYPE (:,hcr)
          CROPTYPE = CROPTYPE_TMP
          PREFTYPE = PREFTYPE_TMP
       endif
    else
       RETURN
    endif
    deallocate (cr_lct)
    
  END SUBROUTINE find_preferred_type

  ! -------------------------------------------------------------------------------------
  
  SUBROUTINE find_preferred_types_fortwo (this, CROPTYPE, IRRIGFRAC, PREFTYPE)

    implicit none

    class(matt_algorithm),  intent(inout) :: this
    real,    dimension(:),intent(inout)   :: PREFTYPE
    real,    dimension(:),intent(inout)   :: CROPTYPE, IRRIGFRAC

    integer                               :: n,i, this_itype, sel_itype, this_rank
    real, dimension (:), allocatable      :: tile_size
    integer, dimension (:), allocatable   :: tile_index
    real                                  :: tile_area, ifrac(3)
        
    n = size (CROPTYPE)
    allocate (tile_size (1:n))
    allocate (tile_index(1:n))
    
    call LIS_indexArrayReal (n, CROPTYPE, tile_index)
    call LIS_reverse (tile_index)
    tile_size = CROPTYPE (tile_index)
    
    ifrac = IRRIGFRAC(:)/SUM(IRRIGFRAC)
    tile_area = SUM(IRRIGFRAC)
    
    do i = 1,n
       if(tile_size(i) > 0.) then
          ! assign the irrigation types from the largest crop fractions
          sel_itype = 0
          this_rank = 100
          do this_itype = 1,3
             if(IRRIGFRAC(this_itype) > 0) then
                if(this%crop_rank(tile_index(i),this_itype) < this_rank) then
                   sel_itype = this_itype
                   this_rank = this%crop_rank(tile_index(i),this_itype)
                endif
             endif
          end do

          if(sel_itype > 0) then
             PREFTYPE (tile_index(i)) = sel_itype           
          else
             PREFTYPE (tile_index(i)) = 1
          endif

          IRRIGFRAC (NINT(PREFTYPE (tile_index(i)))) = IRRIGFRAC (NINT(PREFTYPE (tile_index(i)))) - CROPTYPE (tile_index(i))
          CROPTYPE (tile_index(i)) = 0.
          if(IRRIGFRAC (NINT(PREFTYPE (tile_index(i)))) < 0.) IRRIGFRAC (NINT(PREFTYPE (tile_index(i)))) = 0.
          
       endif
       
    end do
  
  END SUBROUTINE find_preferred_types_fortwo

  ! ---------------------------------------------------------  

  SUBROUTINE assign_irrig_type (this, this_lct, total_itype_area,CROPTYPE, PREFTYPE, second_attempt)
    ! (6) ... If no, then assign that irrigation type to the tile and, if not within
    !      CTYPE_AREA_TOL % (county area) of the irrigation type area target, repeat with the
    !     next largest tile of that type.  If yes, then go to the second largest tile of that
    !     crop type in the county and repeat the test.  Continue until the assigned area of
    !     that irrigation t1ype is within +- 2% of the desired fractional area. 
    implicit none

    class(matt_algorithm),intent(inout) :: this
    integer, intent(in)                 :: this_lct
    real,    intent(inout)              :: total_itype_area
    real,    dimension(:),intent(inout) :: PREFTYPE
    real,    dimension(:),intent(inout) :: CROPTYPE
    integer                             :: n, i
    real, dimension (:), allocatable    :: tile_size
    integer, dimension (:), allocatable :: tile_index
    logical, intent (in), optional      :: second_attempt
    
    n = size (CROPTYPE)
    allocate (tile_size (1:n))
    allocate (tile_index(1:n))

    call LIS_indexArrayReal (n, CROPTYPE, tile_index)
    if(present (second_attempt)) then
       tile_size = CROPTYPE (tile_index)
       do i = 1,n
          if(tile_size(i) > 0.) then
             if((tile_size(i) - total_itype_area)/total_itype_area <= 0.3) then
             PREFTYPE (tile_index(i)) = this_lct
             CROPTYPE (tile_index(i)) = 0.
             total_itype_area = 0. !total_itype_area - tile_size(i)
             return
          endif
          endif
       end do
       return
    endif    
    call LIS_reverse (tile_index)
    tile_size = CROPTYPE (tile_index)

    do i = 1,n
       if(tile_size(i) > 0.) then
          if(tile_size(i) <= total_itype_area) then
             CROPTYPE (tile_index(i)) = 0.
             PREFTYPE (tile_index(i)) = this_lct
             total_itype_area = total_itype_area - tile_size(i)
          endif
       endif
    end do

    deallocate (tile_size, tile_index)
     
  END SUBROUTINE assign_irrig_type
  
  ! ---------------------------------------------------------
  
  SUBROUTINE dissolve_tiny_itypes (this, itype_thres, IRRIGTYPE)

    !(3). that irrigation type covers less than itype_thres of the area of the county divided by the number of pixels
    !     in the county, then ignore it (dissolve it and distribute the fraction evenly)
    
    implicit none
    class(matt_algorithm),intent(inout) :: this
    real, intent (inout), dimension (:) :: IRRIGTYPE
    real, intent (in)                   :: itype_thres
    real, allocatable                   :: ITYPE(:)
    integer :: n,i
    real    :: dis_fr

    n = size(IRRIGTYPE)
    allocate (ITYPE (1:n))

    ITYPE = IRRIGTYPE
    where (IRRIGTYPE < itype_thres)
       ITYPE = 0.
    end where
    n = count (IRRIGTYPE > 0.)
    i = count (ITYPE > 0.    )

    if (n /= i) then
       dis_fr = SUM (IRRIGTYPE) - SUM (ITYPE)
       if (i == 2) then
          where (itype > 0.)
             itype = itype + dis_fr / 2.
          endwhere
       endif
       if (i == 1) then
          where (itype > 0.)
             itype = itype + dis_fr 
          endwhere
       endif       
    endif
    
    IRRIGTYPE = ITYPE
    deallocate (itype)
    
  END SUBROUTINE dissolve_tiny_itypes

  ! ----------------------------------------------------------------
  
  integer function highest_crop_rank (this, cr_lct,CROPTYPE)

    !5.	Use the data in the attached spreadsheet to determine the crop with the highest Crop Rank
    ! corresponding to that irrigation type that also exists in the county
    ! (after tiles have been assigned).
    
    implicit none
    
    class(matt_algorithm),intent(inout) :: this
    real,    dimension(:,:),intent(in)  :: CROPTYPE
    integer, dimension(:),  intent(in)  :: cr_lct 
    integer :: cr (this%N_CROPTYPES), i

    cr = cr_lct
    
    ! void all ranks if  the croptype fraction is zero.
    do i = 1, this%N_CROPTYPES
       if (SUM (CROPTYPE (:,i)) == 0) cr (i) = 100
       if(associated (this%COUNTRY_CROP_IRRAVAIL))then
          if(this%COUNTRY_CROP_IRRAVAIL(i) == 0) cr (i) = 100
       endif
    end do
    highest_crop_rank = 0
    if(minval (cr) < 100) highest_crop_rank = minloc (cr,1)
    
  end function highest_crop_rank

  ! ----------------------------------------------------------------
 
  integer function least_common_type (this, itype)
    
    ! 4.For the least common irrigation type that passes the test (#2)
    ! irrig type: 1: sprinkler; 2: drip; 3: flood
    
    implicit none
    class(matt_algorithm),intent(inout) :: this
    real, intent (in),dimension(:)      :: itype
    integer                             :: i
    real                                :: min_frac
       
    min_frac = 1.
    do i = 1,size(itype)
       if ((itype(i) > 0.).and.(min_frac > itype(i))) then
          least_common_type = i
          min_frac          = itype(i)
       endif
    end do
    
  end function least_common_type

end module getCropIrrigTypes

