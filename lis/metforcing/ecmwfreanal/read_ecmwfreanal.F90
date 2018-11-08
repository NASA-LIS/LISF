!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: read_ecmwfreanal
! \label{read_ecmwfreanal}
!
! !REVISION HISTORY:
!  09 Apr 2002: Urszula Jambor; Original code, based on readgeos.f
!  20 Dec 2002: Matt Rodell; Don't ipolate if running 1/2 deg, also fix
!               error in setting v-wind to zero
!  25 Mar 2003: Urszula Jambor; Modified argument list passed to GEOGFILL2_ECMWFREANAL.
!  24 Nov 2003; Sujay Kumar; Included in LIS
!  15 Mar 2004; Matt Rodell; Fix test for beginning of restart run
!
! !INTERFACE:
subroutine read_ecmwfreanal(order, n, findex, yr, mon, da, hr, ferror)
! USES:
  use LIS_coreMod,        only : LIS_rc,LIS_domain, LIS_localPet, LIS_masterproc
  use LIS_metforcingMod, only : LIS_forc
  use LIS_logMod,         only : LIS_logunit, LIS_endrun
  use ecmwfreanal_forcingMod, only : ecmwfreanal_struc


  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: order
  integer, intent(in)    :: n 
  integer, intent(in)    :: findex
  integer, intent(in)    :: yr
  integer, intent(in)    :: mon
  integer, intent(in)    :: da
  integer, intent(in)    :: hr
  integer, intent(inout) :: ferror
!
! !DESCRIPTION:
!  For the given time, reads 8 parameters from 1/2 degree
!  Reanalysis ECMWF data, transforms into 9 LIS forcing 
!  parameters and interpolates to the LIS domain.
!
! Reanal. ECMWF FORCING VARIABLES available every 6 hours: \newline
!  1. 2T      2 metre air temperature [K], instantaneous \newline
!  2. 2D      2 metre dew point temperature [K], instantaneous \newline
!  3. SSRD    Downward shortwave flux, surface [W/m2s], accumulation \newline
!  4. STRD    Downward longwave radiation, surface [W/m2s], accumulation \newline
!  5. WIND    Calculated absolute 10 metre wind speed [m/s], instantaneous \newline
!  6. P       Calculate surface pressure [Pa], instantaneous \newline
!  7. LSP+CP  Total precipitation [m], accumulation \newline
!  8. CP      Convective precipatation [m], accumulation \newline
!
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    6 hourly instance, order=2, read the next 6 hourly instance)
!  \item[n]
!    index of the nest
!  \item[yr]
!    current year
!  \item[mon]
!    current month
!  \item[da]
!    current day of the year
!  \item[hr]
!    current hour of day
!  \item[ferror]
!    flag to indicate success of the call (=0 indicates success)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[ecmwfreanalgrid\_2\_lisgrid](\ref{ecmwfreanalgrid_2_lisgrid}) \newline
!    transform the ECMWF REANALYSIS data to the LIS grid 
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
!  \item[geogfill2\_ecmwfreanal](\ref{geogfill2_ecmwfreanal}) \newline
!    fill in any grid points due to mismatched grids. 
!  \end{description}
!EOP
  integer :: mo
  ! Specify Reanalysis ECMWF REANALYSIS forcing parameters & file paths
  real, parameter :: no_data_value = -9999.0
  integer, parameter :: N_EF=8, NF=9  ! # ECMWF & GLDAS forcing variables
  real, parameter :: epsln = 0.622        !Constants' values taken from 
  real, parameter :: A = 2.53E8, B=5.42E3 !Roger&Yau, A Short Course
                                          !in Cloud Physics, pp.12-17
  character(25), dimension(N_EF), parameter :: ecmwf_fv = (/  &
       'TEMP-AIR   ',    &
       'TEMP-DEW   ',    & 
       'RAD-SW-DWN ',    &
       'RAD-LW-DWN ',    &
       'WIND       ',    &
       'PRES-SRF   ',    &
       'PREC-TOTL  ',    &
       'PREC-CONV  '     /)

  character(100), dimension(N_EF), parameter :: v_inpath = (/ &
       './FORCING/GRID/     ',  &
       './FORCING/GRID/     ',  &
       './FORCING/GRID/     ',  &
       './FORCING/GRID/     ',  &
       './FORCING/GRID/     ',  &
       './FORCING/GRID/     ',  &
       './FORCING/GRID/     ',  &
       './FORCING/GRID/     '   /)

  integer, dimension(N_EF), parameter :: fnum = (/ &
       31, 33, 34, 35, 36, 37, 38, 39 /) 

  character*200 :: infile

  real :: tempecmwf(ecmwfreanal_struc(n)%ncold,ecmwfreanal_struc(n)%nrold,N_EF)
  real :: temp2ecmwf(ecmwfreanal_struc(n)%ncold,ecmwfreanal_struc(n)%nrold,NF)
  real,allocatable :: templdas(:,:, :)
  real :: p_kPa(ecmwfreanal_struc(n)%ncold,ecmwfreanal_struc(n)%nrold)
  character :: cyr*4, cmo*2 
  integer :: timestep
  integer :: i, j, t, v, ii, r, kk, eindex, x, y
  integer :: ios     ! set to non-zero if there's an error
  integer :: istat

  integer :: gldas,necmwf     ! Size of I/O 1D fields
  real    :: f(ecmwfreanal_struc(n)%ncold*ecmwfreanal_struc(n)%nrold),go(LIS_rc%lnc(n)*LIS_rc%lnr(n))  ! 1D I/O fields
  real    :: tg(LIS_rc%lnc(n),LIS_rc%lnr(n))           ! Interpolated 2D data field
  integer :: iret,km,ip,c
  real :: gridDesco(50)       ! Input,output grid info arrays
  logical*1 :: lb(ecmwfreanal_struc(n)%ncold*ecmwfreanal_struc(n)%nrold)
  logical*1 :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))      ! Input and output bitmaps
  logical*1 :: geogmask(LIS_rc%lnc(n),LIS_rc%lnr(n))! 2D output bitmap
  logical*1 :: tmask(LIS_rc%lnc(n),LIS_rc%lnr(n))   ! 2d valid temperature bitmap
  logical :: file_exists

  ! if a problem, ferror is set to zero
  ferror = 1
  mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)
  !----------------------------------------------------------------
  ! Establish which file timestep the date & hour correspond to
  ! and which file to open.
  !----------------------------------------------------------------

  timestep = 4*(da - 1) + (1 + hr/6)
  if(LIS_masterproc) &
       write(LIS_logunit,*) 'order and month-timestep', order, timestep
  write(cyr, '(i4.4)') yr
  write(cmo, '(i2.2)') mon
  allocate(templdas(LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%nf), stat=ios)
  if(ios .ne.0) then 
     write(LIS_logunit,*) 'Error allocating templdas.',LIS_localPet
     stop 344
  endif
  !=== Open ECMWF forcing files

  if ((ecmwfreanal_struc(n)%findtime1==1) .and. (order==1)) then  !beginning of run
     do v = 1, 7 !N_EF

        ! File name for data variable(v)/year/yearmo
        infile=trim(ecmwfreanal_struc(n)%ecmwfreanaldir)//trim(ecmwf_fv(v))//'/'//cyr//'/'//cyr//cmo

        open(fnum(v), file=trim(infile), form='unformatted', iostat=istat)
        inquire(unit=fnum(v), exist=file_exists)
        
        if (.not.file_exists) then
          write(*,*) 'READ_ECMWFREANAL(1): unit ', fnum(v), '= ', infile,file_exists
          stop
        end if

        if (istat/=0) then
           if(LIS_masterproc) then 
              write(LIS_logunit,*) 'Problem opening file: ', trim(infile),istat
              write(LIS_logunit,*) 'Stopping...'
           endif
           call LIS_endrun
        else
           if(LIS_masterproc) write(LIS_logunit,*) 'Opened file: ', trim(infile)
        end if
        ! Fast-forward to desired location within data file
        if (timestep>1) then
           do t = 1, (timestep-1)
              do i = 1, ecmwfreanal_struc(n)%nrold
                 read(fnum(v)) (tempecmwf(j,i,v), j=1, ecmwfreanal_struc(n)%ncold)
              end do 
           end do
        end if
     end do !v
     
  else if ( timestep==1 ) then !close previous month files & open current
     do v = 1, 7!N_EF

        close(fnum(v))
        ! New file name for data variable(v)/year/yearmo
        infile=trim(ecmwfreanal_struc(n)%ecmwfreanaldir)//trim(ecmwf_fv(v))//'/'//cyr//'/'//cyr//cmo
        open(fnum(v), file=trim(infile), form='unformatted', iostat=istat)
        inquire(unit=fnum(v), exist=file_exists)
        if (.not.file_exists) then
           write(*,*) 'READ_ECMWFREANAL(2): unit ', fnum(v), '= ', infile
           stop
        end if
        if (istat/=0) then
           if(LIS_masterproc) then 
              write(LIS_logunit,*) 'Problem opening file: ', trim(infile)
              write(LIS_logunit,*) 'Stopping...'
              stop
           endif
        else
           if(LIS_masterproc) write(LIS_logunit,*) 'Opened file: ', trim(infile)
        end if

     end do !v

  end if
  do v = 1, 7!N_EF
     !----------------------------------------------------------------
     ! Extract grid array for chosen time
     !----------------------------------------------------------------
     do i = 1, ecmwfreanal_struc(n)%nrold
        read(fnum(v)) (tempecmwf(j,i,v), j=1, ecmwfreanal_struc(n)%ncold)
     end do

     !----------------------------------------------------------------
     ! Change data from ECMWF grid convention to GLDAS one 
     ! Swap latitude band and shift longitudes by 180deg. 
     !----------------------------------------------------------------
!     if(v==2) then 
!        do j=1,berg_struc(n)%nrold
!           do i=1,berg_struc(n)%ncold
!              write(LIS_logunit,*) 'bef ',i,j,tempecmwf(i,j,2)
!           enddo
!        enddo
!     endif
     call ecmwfreanalgrid_2_lisgrid(ecmwfreanal_struc(n)%ncold,ecmwfreanal_struc(n)%nrold,tempecmwf(:,:,v))
!     if(v==2) then 
!        do j=1,ecmwfreanal_struc(n)%nrold
!           do i=1,ecmwfreanal_struc(n)%ncold
!              write(LIS_logunit,*) 'aft ',i,j,tempecmwf(i,j,2)
!           enddo
!        enddo
!     endif
  end do !v
  
  !-----------------------------------------------------------------
  ! Calculate specific humidity from Dew point Temp. & Sfc Pressure.
  ! Approximate: q~epsilon*e/p, p~p_sfc, e=e_s(Td), e_s=A*exp(-B/T)
  !-----------------------------------------------------------------
  do j=1,ecmwfreanal_struc(n)%nrold
     do i=1,ecmwfreanal_struc(n)%ncold
        if(tempecmwf(i,j,6)>=0.0) &
             p_kPa(i,j) = tempecmwf(i,j,6) / 1000.0
        if (tempecmwf(i,j,6)< 0.0) then 
           p_kPa(i,j) = no_data_value
        endif
        if (tempecmwf(i,j,2).gt.no_data_value) then 
           tempecmwf(i,j,2) = epsln*(A * EXP(-B/tempecmwf(i,j,2))) /&
                p_kPa(i,j) 
        endif
  !-----------------------------------------------------------------
  ! Filter out any unrealistic forcing values.
  !-----------------------------------------------------------------

  ! Shortwave
        if (tempecmwf(i,j,3) < 0.0001) &
             tempecmwf(i,j,3) = 0.0001
  ! Wind
        if (tempecmwf(i,j,5) < 0.0001) &
             tempecmwf(i,j,5) = 0.0001
  ! Total precipitation
        if (tempecmwf(i,j,7) < 0.0) &
             tempecmwf(i,j,7) = 0.0
  ! Convective precipitation
        if (tempecmwf(i,j,8) < 0.0) &
             tempecmwf(i,j,8) = 0.0
     
  !----------------------------------------------------------------
  ! Transfer ECMWF forcing fields to GLDAS format
  !----------------------------------------------------------------

        temp2ecmwf(i,j,1) = tempecmwf(i,j,1)
        temp2ecmwf(i,j,2) = tempecmwf(i,j,2)
        temp2ecmwf(i,j,3) = tempecmwf(i,j,3)
        temp2ecmwf(i,j,4) = tempecmwf(i,j,4)
        temp2ecmwf(i,j,5) = tempecmwf(i,j,5)    !Since absolute wind speed 
        temp2ecmwf(i,j,6) = 0.0                 ! already provided as WIND,
        temp2ecmwf(i,j,7) = tempecmwf(i,j,6)    ! let U=WIND and V=0.0
        temp2ecmwf(i,j,8) = tempecmwf(i,j,7)/(6*3600)
        temp2ecmwf(i,j,9) = tempecmwf(i,j,8)/(6*3600)
     end do
  enddo
!  open(112,file="forcing.bin",form='unformatted')
!  write(112) temp2ecmwf(:,:,3)
!  write(LIS_logunit,*) temp2ecmwf(:,:,3)
  !-----------------------------------------------------------------
  ! Interpolate each forcing variable to GLDAS domain
  !-----------------------------------------------------------------

  !=== Initialize input & output grid arrays
  gridDesco = 0
           
  !=== Set input & output grid array values (reanlECMWF to GLDAS)
        
  gridDesco = LIS_rc%gridDesc(n,:)

  !=== Define input & output data bitmaps
  necmwf = ecmwfreanal_struc(n)%ncold*ecmwfreanal_struc(n)%nrold
  gldas  = LIS_rc%lnc(n)*LIS_rc%lnr(n)
  do i=1,necmwf
     if ( ecmwfreanal_struc(n)%remask1d(i) > 0 ) then
        lb(i)=.true.
     else
        lb(i)=.false.
     end if
  enddo
!<debug>
!  tempecmwf(:,:,1) = -9999.0
!  v = 0
!  do r=1,ecmwfreanal_struc(n)%nrold
!     do c=1,ecmwfreanal_struc(n)%ncold
!        if(lb(c+v)) then 
!           tempecmwf(c,r,1) = 1
!        endif
!     enddo
!     v= v + ecmwfreanal_struc(n)%ncold
!  enddo
!  open(112,file="forcing.bin",form='unformatted')
!  write(112) tempecmwf(:,:,1)
!  stop
!</debug>
  lo = .false.
  tmask = .false.

  do v=1,NF
     if (v .ne. 6) then ! not the v-wind component, which is set to zero.
     
      !=== Transferring current data to 1-D array for interpolation
      c=0
      do i=1,ecmwfreanal_struc(n)%nrold
          do j=1,ecmwfreanal_struc(n)%ncold
             c = c + 1
             f(c) = temp2ecmwf(j,i,v)
!             if(v.eq.1 .and. f(c).ne.-9999.0) write(LIS_logunit,*) c,f(c)
          enddo
      enddo

      !=== Interpolate if forcing and model grids are not both half deg.
!      if (LIS_rc%domain .ne. 3) then
      if (LIS_rc%gridDesc(n,9).ne.0.5) then
          ip = 0
          km = 1
          if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 
             call bilinear_interp(gridDesco,lb,f,lo,go,ecmwfreanal_struc(n)%mi,mo, & 
                  LIS_domain(n)%lat, LIS_domain(n)%lon,&
                  ecmwfreanal_struc(n)%w111,ecmwfreanal_struc(n)%w121,&
                  ecmwfreanal_struc(n)%w211,ecmwfreanal_struc(n)%w221, & 
                  ecmwfreanal_struc(n)%n111,ecmwfreanal_struc(n)%n121,&
                  ecmwfreanal_struc(n)%n211,ecmwfreanal_struc(n)%n221,LIS_rc%udef, iret)
          elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 
             if(v.eq.7) then 
                call conserv_interp(gridDesco,lb,f,lo,go,ecmwfreanal_struc(n)%mi,mo, & 
                     LIS_domain(n)%lat, LIS_domain(n)%lon,&
                     ecmwfreanal_struc(n)%w112,ecmwfreanal_struc(n)%w122,&
                     ecmwfreanal_struc(n)%w212,ecmwfreanal_struc(n)%w222, & 
                     ecmwfreanal_struc(n)%n112,ecmwfreanal_struc(n)%n122,&
                     ecmwfreanal_struc(n)%n212,ecmwfreanal_struc(n)%n222,LIS_rc%udef,iret)
             else
                call bilinear_interp(gridDesco,lb,f,lo,go,ecmwfreanal_struc(n)%mi,mo, & 
                     LIS_domain(n)%lat, LIS_domain(n)%lon,&
                     ecmwfreanal_struc(n)%w111,ecmwfreanal_struc(n)%w121,&
                     ecmwfreanal_struc(n)%w211,ecmwfreanal_struc(n)%w221, & 
                     ecmwfreanal_struc(n)%n111,ecmwfreanal_struc(n)%n121,&
                     ecmwfreanal_struc(n)%n211,ecmwfreanal_struc(n)%n221,LIS_rc%udef, iret)
             endif
          endif
      else ! forcing and model grids both half degree
        kk = 0
        do r=1,LIS_rc%lnr(n)
          y = r + (LIS_rc%gridDesc(n,4) + 89.750) / 0.500
          do c=1,LIS_rc%lnc(n)
            x = c + (LIS_rc%gridDesc(n,5) + 179.750) / 0.500
            kk = kk + 1
            eindex = ((y - 1) * 720) + x
            go(kk) = f(eindex)
            lo(kk) = lb(eindex)
          end do
        end do

      end if ! LIS_rc%domain 
    
       !=== Convert data to original 3D array & a 2D array to 
       !=== fill in of missing points due to geography difference  
!      tg = -9999.0
      tg = 0.0 
      c = 0
      geogmask = .false.
      do j = 1, LIS_rc%lnr(n)
         do i = 1, LIS_rc%lnc(n)
!            write(LIS_logunit,*) i,j,gindex(i,j),lo(i+c)
            if(LIS_domain(n)%gindex(i,j).ne.-1) then 
               geogmask(i,j) = lo(i+c)
               tg(i,j) = go(i+c)
            endif
         enddo
         c = c + LIS_rc%lnc(n)
      enddo
!      if(v==3) then 
!         open(112,file="forcing.bin",form='unformatted')
!         write(112) tg
!         stop
!      endif
      if(v==2) then 
!         do j=1,LIS_rc%lnr(n)
!            do i=1, LIS_rc%lnc(n)
!               write(LIS_logunit,*) 'before ',tg(14,1)
!            enddo
!         enddo
      endif

      call geogfill2_ecmwfreanal(n, LIS_rc%lnc(n),LIS_rc%lnr(n),geogmask,tg,v,tmask)
      !       write(LIS_logunit,*) gindex(21,3),tg(21,3),v
#if 0 
      if(v==3) then 
         do j=1,LIS_rc%lnr(n)
            do i=1, LIS_rc%lnc(n)
               write(LIS_logunit,*) 'after ',i,j, tg(i,j)
            enddo
         enddo
      endif
#endif
      do j = 1, LIS_rc%lnr(n)
         do i = 1, LIS_rc%lnc(n)
            if(LIS_domain(n)%gindex(i,j).ne.-1) then 
               if (tg(i,j) .eq. -1.0) then
                  write(LIS_logunit,*) 'No nearest neighbours, v, i, j',v,i,j
                  stop
               endif
               templdas(i,j,v) = tg(i,j)
            endif
         end do !c
       enddo ! r

    else ! v==6, v-wind component, always zero
       templdas(:,:,6) = 0.0
    endif
    
 enddo !v

  do v= 1,NF
    !=== Fill in undefined and ocean points
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1)then
                 ecmwfreanal_struc(n)%metdata1(v,LIS_domain(n)%gindex(c,r))=templdas(c,r,v)
              else
                 ecmwfreanal_struc(n)%metdata2(v,LIS_domain(n)%gindex(c,r))=templdas(c,r,v)
              endif
           endif
        enddo
     enddo !r
  enddo !v
  deallocate(templdas)

end subroutine read_ecmwfreanal


!BOP
! 
! !ROUTINE: ecmwfreanalgrid_2_lisgrid
! \label{ecmwfreanalgrid_2_lisgrid}
!
! !REVISION HISTORY:
!  10 Apr 2002: Urszula Jambor;  Code adapted from 
!               ecmwfgrid_2_grid2catgrid, by R. Reichle
! !INTERFACE:
subroutine ecmwfreanalgrid_2_lisgrid( nx, ny, grid_data )

  implicit none
! !ARGUMENTS:   
  integer, intent(in)                   :: nx, ny
  real, intent(inout), dimension(nx,ny) :: grid_data
!
! !DESCRIPTION:
! Changes grid data from ECMWF data convention to LIS convention
!
! ECMWF: North-to-South around Greenwich Meridian
! Global grid. Data are written as flat binary from ``upper left to 
! lower right'' starting at 0.5-degree grid point center coordinates: 
! 0.25E,89.75N and going to 0.25W,89.75S. Here is the write statement:
! 
! \begin{verbatim}
! do i = 1,360   
!   write(14) (val(j,i),j=1,720)
! end do
! \end{verbatim}
!
! LIS: South-to-North around Date Line
! Full global grid.  Starts at the southernmost latitude and date line, 
! going east and then north.
!EOP
  
  integer :: i, j, m, n
  real :: tmp, tmp_data1(nx), tmp_data2(nx)
  
  ! ------------------------------------------------------------------
  ! some checks
  
  if ((nx /= 720) .or. (ny /= 360)) then
     write (*,*) 'Recmwfgrid_2_gldasgrid(): This routine has only been'
     write (*,*) 'checked for nx=720 and ny=360. Make sure you know'
     write (*,*) 'what you are doing. STOPPING.'
     stop
  end if  
  if ((mod(nx,2) /= 0) .or. (mod(ny,2) /= 0)) then
     write (*,*) 'Recmwfgrid_2_gldasgrid(): This routine can only work'
     write (*,*) 'for even nx and ny. Make sure you know'
     write (*,*) 'what you are doing. STOPPING.'
     stop
  end if
  
  !-------------------------------------------------------------------
  
  do j=1,ny/2
     
     ! swap latitude bands (North-to-South becomes South-to-North)
     n = ny-j+1
     tmp_data1      = grid_data(:,j)
     tmp_data2      = grid_data(:,n)
     
     do i=1,nx/2

        ! shift longitudes (wrapping around Greenwhich Meridian becomes
        !  wrapping around Date Line)
        m = i + nx/2
        tmp          = tmp_data1(i)
        tmp_data1(i) = tmp_data1(m)
        tmp_data1(m) = tmp
        
        tmp          = tmp_data2(i)
        tmp_data2(i) = tmp_data2(m)
        tmp_data2(m) = tmp
        
     end do
     
     grid_data(:,j) = tmp_data2
     grid_data(:,n) = tmp_data1
     
  end do

end subroutine Ecmwfreanalgrid_2_lisgrid

!BOP
! 
! !ROUTINE: fillgaps
! 
! !DESCRIPTION:
! Fills in values for NSIPP tilespace land points
! where no ECMWF reanalysis data is available via GEOGFILL 
! by assigning most appropriate land-point value along
! the latitudinal circle of original tilespace point.
! Developed manually with 17 points in mind.
!
! !REVISION HISTORY:
!  23 Jul 2002: Urszula Jambor
!  03 Sep 2002: Urszula Jambor, revised points to reflect 
!               NSIPP land mask correction.
! !INTERFACE:
subroutine fillgaps( nc, nr, v, arr )
!EOP
  implicit none

  integer :: nc, nr, v
  real :: arr(nc,nr)

  arr( 55, 1) = arr( 59, 4)
  arr( 62, 2) = arr( 59, 4)
  arr( 88, 8) = arr( 46, 8)
  arr( 94, 8) = arr( 46, 8)
  arr(  2, 9) = arr(142, 9)
  arr(  3,20) = arr(133,20)
  arr(139,20) = arr(133,20)
  arr(140,20) = arr(133,20)
  arr( 95,29) = arr( 89,29)
  arr( 37,30) = arr( 41,30)
  arr( 36,31) = arr( 41,31)
  arr( 37,31) = arr( 41,31)
  arr(136,34) = arr(123,34)
  arr( 62,50) = arr( 69,50)
  arr( 63,50) = arr( 69,50)
  arr(142,57) = arr(144,57)
  arr( 70,67) = arr( 83,67)

  if (v == 3) then
    arr( 88, 8) = arr( 84,14)
    arr( 94, 8) = arr( 84,14)
  endif

end subroutine fillgaps



