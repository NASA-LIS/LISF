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
! !ROUTINE: createtiles_catch
!  \label{createtiles_catch}
! 
! !REVISION HISTORY:
!  10 Jul 2006: James Geiger; Initial version
!
! !INTERFACE:
subroutine createtiles_catch
! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_localPet, LIS_masterproc,  &
       LIS_npes, LIS_ews_ind, LIS_ewe_ind, LIS_nss_ind, LIS_nse_ind, &
       LIS_deltas, LIS_offsets
  use LIS_mpiMod
  use LIS_logMod, only :LIS_logunit, LIS_verify, LIS_endrun
  use LIS_domainMod, only : LIS_domain_setup
  use tile_coord_types
  use tile_coord_routines, only : read_tile_coord, read_tile_coord_header

   implicit none
! !ARGUMENTS: 
!
! !DESCRIPTION:
!  This routine computes the grid and tile spaces based on the 
!  input distribution of the catchment-based tiles.
!  The routine 
!  also computes the following structures necessary for domain 
!  decomposition and other parallel routines. \newline
!  {\tt  ntiles\_pergrid} : Number of tiles per each grid \newline
!  {\tt  s\_tid} : Index of the first tile in each grid cell  \newline
! 
!EOP
   type(grid_def_type)                          :: atm_grid
   type(tile_coord_type), allocatable, dimension(:) :: tile_coord


   integer :: N_tiles, n, t, count1, tilec, tiler, c, r, m, &
              index, ierr, max_tpg, nt

   integer, allocatable, dimension(:,:)    :: tiles_per_grid, acount
   integer, allocatable, dimension(:,:,:)  :: gtn, sindex
   integer, allocatable, dimension(:)      :: tindex_g2l
   real, allocatable, dimension(:,:,:)     :: fgrd
   real, allocatable, dimension(:,:)       :: sum
   real, allocatable, dimension(:)         :: temp_sort
   integer, allocatable, dimension(:)      :: temp_index
   integer                             :: ftn
   integer                             :: N_tile

   real :: locallat, locallon

   integer :: cat_row_offset, row_offset, col_offset
   real    :: atm_lat, atm_lon

   ! THIS BELONGS SOMEWHERE ELSE!
   integer, allocatable, dimension(:) :: tmptileid
   integer                            :: dummy_int

   do n = 1, LIS_rc%nnest

      ftn = 10 
      open (ftn, file=trim(LIS_rc%tile_coord_file(n)), form='formatted', action='read') 

      call read_tile_coord_header(ftn, &
           trim(LIS_rc%tile_coord_file(n)), atm_grid, N_tile)
      allocate(tile_coord(N_tile))
      call read_tile_coord(ftn, N_tile, tile_coord)
      close(ftn, status='keep')


      cat_row_offset = atm_grid%N_lat - LIS_rc%pnr(n)
      row_offset = nint((LIS_rc%gridDesc(n,4)-LIS_rc%gridDesc(n,34))/&
                        LIS_rc%gridDesc(n,10))
      col_offset = nint((LIS_rc%gridDesc(n,5)-LIS_rc%gridDesc(n,35))/&
                        LIS_rc%gridDesc(n,9))

      !print*,'cat_row_offset = ',cat_row_offset
      !print*,'row_offset = ',row_offset
      !print*,'col_offset = ',col_offset

      N_tiles = size(tile_coord)

      allocate( tindex_g2l(N_tiles) )
      do t = 1, N_tiles
         tindex_g2l(tile_coord(t)%tile_id) = t
      enddo

      ! Impose surface_minp
      do t = 1, N_tiles
         if ( tile_coord(t)%frac_atm < LIS_rc%surface_minp ) then
            tile_coord(t)%frac_atm = 0.0
         endif
      enddo

      ! Create gindex -- this maps (column, row) to grid-cell index.
      allocate( LIS_domain(n)%gindex(LIS_rc%lnc(n), LIS_rc%lnr(n)) )
      allocate( tiles_per_grid(LIS_rc%lnc(n), LIS_rc%lnr(n)) )

      ! Count all land tiles
      count1 = 0
      do t = 1, N_tiles
         if ( tile_coord(t)%typ == tile_typ_land ) then
            count1 = count1 + 1
         endif
      enddo
      LIS_rc%ncatg(n) = count1

      LIS_domain(n)%gindex = -1
      tiles_per_grid = 0
      do t = 1, N_tiles
         if ( tile_coord(t)%typ == tile_typ_land .and. &
              tile_coord(t)%frac_atm > 0.0) then
            tilec = tile_coord(t)%i_atm
            tiler = tile_coord(t)%j_atm-cat_row_offset
            atm_lat = LIS_rc%gridDesc(n,34) + (tiler-1)*LIS_rc%gridDesc(n,39)
            atm_lon = LIS_rc%gridDesc(n,35) + (tilec-1)*LIS_rc%gridDesc(n,40)
            if ( atm_lat >= LIS_rc%gridDesc(n,4) .and. &
                 atm_lat <= LIS_rc%gridDesc(n,7) .and. &
                 atm_lon >= LIS_rc%gridDesc(n,5) .and. &
                 atm_lon <= LIS_rc%gridDesc(n,8) ) then
               c = tilec - col_offset 
               r = tiler - row_offset 
               LIS_domain(n)%gindex(c, r) = 1
               tiles_per_grid(c, r) = tiles_per_grid(c, r) + 1
            endif
         endif
      enddo

      ! Count all grid boxes containing land
      count1 = 1
      do r = 1, LIS_rc%lnr(n)
         do c = 1, LIS_rc%lnc(n)
            if ( LIS_domain(n)%gindex(c,r) == 1 ) then
               LIS_domain(n)%gindex(c,r) = count1
               count1 = count1 + 1
            endif
         enddo
      enddo
      LIS_rc%ngrid(n) = count1 - 1

      allocate( LIS_domain(n)%grid(LIS_rc%ngrid(n)) )

      LIS_domain(n)%minLat = 200.0
      LIS_domain(n)%minLon = 200.00
      LIS_domain(n)%maxLat = -200.0
      LIS_domain(n)%maxLon = -200.00

      do r = 1, LIS_rc%lnr(n)
         do c = 1, LIS_rc%lnc(n)
            t = LIS_domain(n)%gindex(c,r)
            locallat = LIS_rc%gridDesc(n,4)+(r-1)*LIS_rc%gridDesc(n,10)
            locallon = LIS_rc%gridDesc(n,5)+(c-1)*LIS_rc%gridDesc(n,9)

            if(localLat.lt.LIS_domain(n)%minLat) LIS_domain(n)%minLat = localLat
            if(localLat.gt.LIS_domain(n)%maxLat) LIS_domain(n)%maxLat = localLat
            if(localLon.lt.LIS_domain(n)%minLon) LIS_domain(n)%minLon = localLon
            if(localLon.gt.LIS_domain(n)%maxLon) LIS_domain(n)%maxLon = localLon

            if ( t /= -1 ) then
               LIS_domain(n)%grid(t)%lat = locallat
               LIS_domain(n)%grid(t)%lon = locallon
               LIS_domain(n)%grid(t)%col = c
               LIS_domain(n)%grid(t)%row = r
            endif
         enddo
      enddo

      max_tpg = maxval(tiles_per_grid)
      allocate( fgrd(LIS_rc%lnc(n),LIS_rc%lnr(n),max_tpg), stat=ierr)
      call LIS_verify(ierr)
      allocate( gtn(LIS_rc%lnc(n),LIS_rc%lnr(n),max_tpg) , stat=ierr)
      call LIS_verify(ierr)
      allocate( sindex(LIS_rc%lnc(n),LIS_rc%lnr(n),max_tpg), stat=ierr)
      call LIS_verify(ierr)
      allocate( acount(LIS_rc%lnc(n),LIS_rc%lnr(n)), stat=ierr)
      call LIS_verify(ierr)
      fgrd = 0.0
      sindex = -1
      gtn = -1
      acount = 1
      do t = 1, N_tiles
         if ( tile_coord(t)%typ == tile_typ_land .and. &
              tile_coord(t)%frac_atm > 0.0) then
            tilec = tile_coord(t)%i_atm
            tiler = tile_coord(t)%j_atm-cat_row_offset
            c = tilec - col_offset
            r = tiler - row_offset
            if ( r >= 1 .and. r <= LIS_rc%lnr(n) .and. &
                 c >= 1 .and. c <= LIS_rc%lnc(n) ) then
               if ( LIS_domain(n)%gindex(c,r) /= -1 ) then
                  fgrd(c,r,acount(c,r)) = tile_coord(t)%frac_atm
                  gtn(c,r,acount(c,r))  = tile_coord(t)%tile_id
                  acount(c,r) = acount(c,r) + 1
               endif
            endif
         endif
      enddo

      ! For each grid-cell (c,r) sort the fraction-of-grid values.
      ! Note that the sorting routine used sorts in ascending order.
      ! We want to use the values in descending order, hence the -1.0*fgrd.
      ! Note that this routine does not actually sort the given array, rather
      ! it returns a list of indices.
      allocate(temp_sort(max_tpg))
      allocate(temp_index(max_tpg))
      do r = 1, LIS_rc%lnr(n)
         do c = 1, LIS_rc%lnc(n)
            temp_sort = -1.0*fgrd(c,r,:)
            call sortrx(max_tpg, temp_sort, temp_index)
            sindex(c,r,:) = temp_index
         enddo
      enddo
      deallocate(temp_sort)
      deallocate(temp_index)

      ! Impose SURFACE_MAXT
      ! Recall that the fgrd array has been ``sorted'', and the tiles we
      ! want are at the beginning of the array.
      ! If there are more tiles per grid than we need, then reset the smaller
      ! ones to 0.
      if ( max_tpg > LIS_rc%surface_maxt ) then
         nt = max_tpg - LIS_rc%surface_maxt
         do r = 1, LIS_rc%lnr(n)
            do c = 1, LIS_rc%lnc(n)
               do t = LIS_rc%surface_maxt+1, max_tpg
                  fgrd(c,r,sindex(c,r,t)) = 0.0
               enddo
            enddo
         enddo
      endif

      LIS_rc%ntiles(n) = 0
      do m=1,LIS_rc%nensem(n)
         do r = 1, LIS_rc%lnr(n)
            do c = 1, LIS_rc%lnc(n)
               do t = 1, max_tpg
                  if ( fgrd(c,r,sindex(c,r,t)) > 0 ) then
                     LIS_rc%ntiles(n) = LIS_rc%ntiles(n) + 1
                  endif
               enddo
            enddo
         enddo
      enddo

      LIS_rc%glbntiles(n) = LIS_rc%ntiles(n)
      allocate( LIS_domain(n)%tile(LIS_rc%ntiles(n)) )
      count1 = 1

      do r = 1, LIS_rc%lnr(n)
         do c = 1, LIS_rc%lnc(n)
            do t = 1, max_tpg
               do m=1,LIS_rc%nensem(n)
                  if ( fgrd(c,r,sindex(c,r,t)) > 0.0 ) then
                     index = tindex_g2l(gtn(c,r,sindex(c,r,t)))
                     LIS_domain(n)%tile(count1)%ensem   = m
                     LIS_domain(n)%tile(count1)%row     = r
                     LIS_domain(n)%tile(count1)%col     = c
                     LIS_domain(n)%tile(count1)%index   = LIS_domain(n)%gindex(c,r)
! INCORRECT USE OF VEGT
                     LIS_domain(n)%tile(count1)%vegt    = tile_coord(index)%typ 
                     LIS_domain(n)%tile(count1)%fgrd    = tile_coord(index)%frac_atm
                     LIS_domain(n)%tile(count1)%tile_id = tile_coord(index)%tile_id
                     LIS_domain(n)%tile(count1)%com_lon = tile_coord(index)%com_lon
                     LIS_domain(n)%tile(count1)%com_lat = tile_coord(index)%com_lat
                     count1 = count1 + 1
                  endif
               enddo
            enddo
         enddo
      enddo
      ! renormalize fgrd
      allocate( sum(LIS_rc%lnc(n),LIS_rc%lnr(n)) )
      sum = 0.0
      do count1 = 1, LIS_rc%ntiles(n)
         sum(LIS_domain(n)%tile(count1)%col,LIS_domain(n)%tile(count1)%row)    = &
            sum(LIS_domain(n)%tile(count1)%col,LIS_domain(n)%tile(count1)%row) + &
            LIS_domain(n)%tile(count1)%fgrd
      enddo
      do count1 = 1, LIS_rc%ntiles(n)
         LIS_domain(n)%tile(count1)%fgrd = LIS_domain(n)%tile(count1)%fgrd *&
              LIS_rc%nensem(n)/  &
              sum(LIS_domain(n)%tile(count1)%col,LIS_domain(n)%tile(count1)%row)
      enddo

      ! create d2g
      ! THIS BELONGS SOMEWHERE ELSE!
      allocate(tmptileid(LIS_rc%ncatg(n)))
      open (10,file=trim(LIS_rc%tile_veg_file(n)),form='formatted',status='old')
      tmptileid = 0
      do t=1,LIS_rc%ncatg(n)
         read (10,*) tmptileid(t), dummy_int, dummy_int
      end do
      close(10)
      do t = 1, LIS_rc%ntiles(n)
         index = LIS_domain(n)%tile(t)%tile_id
         do count1 = 1, LIS_rc%ncatg(n)
            if ( tmptileid(count1) == index ) then
               exit
            endif
         enddo
         if ( count1 > LIS_rc%ncatg(n) ) then
            print*,'ERR creating d2g'
            call LIS_endrun
         else
            LIS_domain(n)%tile(t)%d2g = count1
         endif
      enddo

      deallocate(sum)
      deallocate(tiles_per_grid)
      deallocate(fgrd)
      deallocate(sindex)
      deallocate(gtn)
      deallocate(acount)
      deallocate(tindex_g2l)
      deallocate(tmptileid)


      
      write(unit=LIS_logunit,fmt=*)'global domain',LIS_localPet,':(',LIS_rc%gnc(n),LIS_rc%gnr(n),')'
      write(unit=LIS_logunit,fmt=*)'local domain',LIS_localPet,':(',LIS_rc%lnc(n),LIS_rc%lnr(n),')'
      write(unit=LIS_logunit,fmt=*)'num grids',LIS_localPet,':(',LIS_rc%ngrid(n),')'
      write(unit=LIS_logunit,fmt=*)'num tiles',LIS_localPet,':(',LIS_rc%ntiles(n),')'
!open(unit=10,file='mask.dat',form='unformatted')
!write(10) LIS_domain(n)%gindex
!close(10)
      !stop 666
      call LIS_domain_setup(n)

   enddo

end subroutine createtiles_catch

