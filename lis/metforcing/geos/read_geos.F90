!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: read_geos
! \label{read_geos}
! 
! !REVISION HISTORY:
!  11 Apr 2000: Brian Cosgrove; Added read statements for forcing interpolation
!  17 Apr 2001: Jon Gottschalck; Added code to perform initialization of
!                                Mosaic with GEOS forcing and new intp. scheme
!  14 Aug 2001: Urszula Jambor; Added ferror flag as a routine argument
!  07 Dec 2001: Urszula Jambor; Began used of LDAS$%$LDAS_GRIDDESC array
!  10 Feb 2003: Sujay Kumar: Initial Specification in LIS.    
!
! !INTERFACE:      
subroutine read_geos(order,n, findex, &
     name,tscount,ferror)
! !USES:
  use LIS_coreMod, only : LIS_rc,LIS_domain, LIS_masterproc
  use LIS_logMod,only : LIS_logunit, LIS_endrun, LIS_verify, &
                        LIS_getNextUnitNumber, LIS_releaseUnitNumber
  use LIS_metforcingMod, only: LIS_forc
  use geos_forcingMod, only : geos_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in)          :: n
  integer, intent(in)          :: findex
  character(len=*), intent(in) :: name
  integer, intent(in)          :: order, tscount
  integer, intent(out)         :: ferror          

!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  GEOS data, transforms into 9 LIS forcing 
!  parameters and interpolates to the LIS domain. \newline
!
! GEOS FORCING VARIABLES (unless noted, fields are 3-hr upstream averaged): \newline
!  1. T 2m    Temperature interpolated to 2 metres [$K$] \newline
!  2. q 2m    Instantaneous specific humidity interpolated to 2 metres[$kg/kg$] \newline
!  3. radswg  Downward shortwave flux at the ground [$W/m^2$] \newline
!  4. lwgdwn  Downward longwave radiation at the ground [$W/m^2$] \newline
!  5. u 10m   Instantaneous zonal wind interpolated to 10 metres [$m/s$] \newline
!  6. v 10m   Instantaneous meridional wind interpolated to 10 metres[$m/s$] \newline
!  7. ps      Instantaneous Surface Pressure [$Pa$] \newline
!  8. preacc  Total precipitation [$mm/s$] \newline
!  9. precon  Convective precipatation [$mm/s$] \newline
! 10. albedo  Surface albedo (0-1) \newline
!
!  The arguments are: 
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    3 hourly instance, order=2, read the next 3 hourly instance)
!  \item[n]
!    index of the nest
!  \item[name]
!    name of the 3 hour GDAS forecast file
!  \item[tscount]
!    time step count
!  \item[ferror]
!    return error code (0 indicates success)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
!  \end{description}
!EOP

  integer :: mo
  integer :: nforce,i,j,v
  integer :: ios,ioerror      ! set to non-zero if there's an error
  real, allocatable :: tempgeos(:,:,:) 
                                  ! lf%nmif = Max.# parameters to retrieve
  real, allocatable :: tempvar(:,:,:)
  integer :: glis,ngeos                ! Size of I/O 1D fields
  real, allocatable :: f(:),go(:) ! 1D I/O fields
  real, allocatable :: gmask(:)  ! GEOS forcing mask
  integer :: count1
  integer :: iret,c,r
  real :: gridDesco(50)       ! Input,output grid info arrays
  logical*1, allocatable :: lb(:)
  logical*1, allocatable :: lo(:)      ! Input and output bitmaps
 
  integer :: ftn 
  integer :: nr_index, nc_index
  nr_index = geos_struc(n)%nrold
  nc_index = geos_struc(n)%ncold

  allocate(tempvar(nc_index,nr_index,geos_struc(n)%nmif), stat=ios)
  call LIS_verify(ios,'Error allocating tempvar.')
  
  tempvar = 0.0

  allocate(tempgeos(LIS_rc%lnc(n),LIS_rc%lnr(n),geos_struc(n)%nmif), stat=ios)
  call LIS_verify(ios,'Error allocating tempgeos.')
  
  allocate(f(nc_index*nr_index), stat=ios)
  call LIS_verify(ios,'Error allocating f.')
  
  allocate(go(LIS_rc%lnc(n)*LIS_rc%lnr(n)), stat=ios)
  call LIS_verify(ios,'Error allocating go.')
  
  allocate(gmask(nc_index*nr_index), stat=ios)
  call LIS_verify(ios,'Error allocating gmask.')
  
  allocate(lb(nc_index*nr_index), stat=ios)
  call LIS_verify(ios,'Error allocating lb.')
  
  allocate(lo(LIS_rc%lnc(n)*LIS_rc%lnr(n)), stat=ios)
  call LIS_verify(ios,'Error allocating lo.')
  
  gmask = 0.0
  ngeos = nc_index*nr_index
  glis = LIS_rc%lnc(n)*LIS_rc%lnr(n)
  ferror = 1                
!------------------------------------------------------------------------
! Open GEOS forcing file
!------------------------------------------------------------------------
  if(LIS_masterproc) then 
     write(LIS_logunit,*) 'MSG: read_geos -- Reading GEOS forcing file -', & 
          trim(name)
  endif
  
  ftn = LIS_getNextUnitNumber()
  open(ftn,file=name,form='unformatted',iostat=ios)
  !         open(40,file=name,form='unformatted', &
  !              access="direct",recl=geos_struc(n)%ncold*geos_struc(n)%nrold*4,iostat=ios)
  
  if ( ios /= 0 ) then
     write(LIS_logunit,*)'ERR: read_geos -- Error opening ', & 
          trim(name),'. Stopping.'
     call LIS_endrun
  endif

  read(ftn,iostat=ioerror)tempvar
  
  if ( ioerror /= 0 ) then
     write(LIS_logunit,*)'ERR: read_geos -- Error reading ', &
          trim(name),'. ioerror = ',ioerror,'. Stopping.'
     call LIS_endrun
  else
     if(LIS_masterproc) then 
        write(LIS_logunit,*) 'MSG: read_geos -- Read GEOS forcing file -', & 
             trim(name)
     endif
  endif
!------------------------------------------------------------------------
! Finding number of forcing variables 
! (13 if time step is 0, otherwise the normal 10)
!------------------------------------------------------------------------
  if (tscount .eq. 0) then
     nforce = geos_struc(n)%nmif
  else
     nforce = 10 
  endif
  do v=1,nforce
!------------------------------------------------------------------------
! Transferring current data to 1-D array for interpolation
!------------------------------------------------------------------------

     c=0

     do j=1,nr_index    !i=1,lf%nrold
        do i=1,nc_index !j=1,lf%ncold
           c = c + 1
           f(c) = tempvar(i,j,v)
           if (tscount .eq. 0 .and. order .eq. 1 & 
                .and. v .eq. 11) then
              gmask(c) = f(c) ! Storing geos land mask for later use
           endif
        enddo
     enddo
!------------------------------------------------------------------------     
! Initializing input and output grid arrays
!------------------------------------------------------------------------
     gridDesco = 0
     gridDesco = LIS_rc%gridDesc(n,:)
     mo = LIS_rc%lnc(n) * LIS_rc%lnr(n)
!------------------------------------------------------------------------
! Defining input data bitmap
!------------------------------------------------------------------------
     do i=1,ngeos
        lb(i)=.true.
     enddo
!------------------------------------------------------------------------     
! Alter default bitmap prescribed above for 
! surface parameters (soil wetness, snow)
!------------------------------------------------------------------------
     if (v .eq. 12 .or. v .eq. 13) then
        do i=1,ngeos
           if (gmask(i)==100.0 .or. gmask(i)==101.0) then
              lb(i)=.false.
           else
              lb(i)=.true.
           endif
        enddo
     endif
!------------------------------------------------------------------------
! Defining output data bitmap
!------------------------------------------------------------------------
     do i=1,glis
        lo(i)=.true.
     enddo
!------------------------------------------------------------------------
! Interpolate data from GEOS grid to GLDAS grid
!------------------------------------------------------------------------
     if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 
        call bilinear_interp(gridDesco,lb,f,lo,go,geos_struc(n)%mi,mo, & 
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             geos_struc(n)%w111,geos_struc(n)%w121,&
             geos_struc(n)%w211,geos_struc(n)%w221,&
             geos_struc(n)%n111,geos_struc(n)%n121,&
             geos_struc(n)%n211,geos_struc(n)%n221,LIS_rc%udef, iret)
     elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then 
        if(v.eq.8 .or. v.eq. 9) then            
           call conserv_interp(gridDesco,lb,f,lo,go,geos_struc(n)%mi,mo,& 
                LIS_domain(n)%lat, LIS_domain(n)%lon,&
                geos_struc(n)%w112,geos_struc(n)%w122,&
                geos_struc(n)%w212,geos_struc(n)%w222,&
                geos_struc(n)%n112,geos_struc(n)%n122,&
                geos_struc(n)%n212,geos_struc(n)%n222,LIS_rc%udef, iret)
        else 
           call bilinear_interp(gridDesco,lb,f,lo,go,geos_struc(n)%mi,mo, & 
                LIS_domain(n)%lat, LIS_domain(n)%lon,&
                geos_struc(n)%w111,geos_struc(n)%w121,&
                geos_struc(n)%w211,geos_struc(n)%w221,&
                geos_struc(n)%n111,geos_struc(n)%n121,&
                geos_struc(n)%n211,geos_struc(n)%n221,LIS_rc%udef, iret)
        endif
     endif

!------------------------------------------------------------------------
! Convert data to original 3D array & a 2D array to 
! fill in of missing points due to geography difference  
!------------------------------------------------------------------------
     count1 = 0
     do j = 1, LIS_rc%lnr(n) !j = 1, LIS_rc%nr
        do i = 1, LIS_rc%lnc(n)                !i = 1, LIS_rc%nc
           if(LIS_domain(n)%gindex(i,j) .ne. -1) then
              tempgeos(i,j,v) = go(i+count1)
           endif
        enddo
        count1 = count1 + LIS_rc%lnc(n)
     enddo

!------------------------------------------------------------------------
! Fill in undefined and ocean points
!------------------------------------------------------------------------
     if(v.le.LIS_rc%met_nf(findex)) then 
        do r=1,LIS_rc%lnr(n)
           do c=1,LIS_rc%lnc(n)
              if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                 if(order.eq.1)then
                    geos_struc(n)%metdata1(v,LIS_domain(n)%gindex(c,r))=tempgeos(c,r,v)
                 else
                    geos_struc(n)%metdata2(v,LIS_domain(n)%gindex(c,r))=tempgeos(c,r,v)
                 endif
              endif
           enddo
        enddo
     endif
  enddo                  !v
  
  deallocate(tempvar, stat=ios)
  call LIS_verify(ios,'Error deallocating tempvar.')
  
  deallocate(tempgeos, stat=ios)
  call LIS_verify(ios,'Error deallocating tempgeos.')
  
  deallocate(f, stat=ios)
  call LIS_verify(ios,'Error deallocating f.')
  
  deallocate(go, stat=ios)
  call LIS_verify(ios,'Error deallocating go.')
  
  deallocate(gmask, stat=ios)
  call LIS_verify(ios,'Error deallocating gmask.')
  
  deallocate(lb, stat=ios)
  call LIS_verify(ios,'Error deallocating lb.')
  
  deallocate(lo, stat=ios)
  call LIS_verify(ios,'Error deallocating lo.')
  if(LIS_masterproc) then 
     write(LIS_logunit,*) 'MSG: read_geos -- Closing GEOS forcing file -', & 
          trim(name)
  endif
  close(ftn, iostat=ios)
  call LIS_releaseUnitNumber(ftn)
  if ( ios /= 0 ) then
     write(LIS_logunit,*)'ERR: read_geos -- Error closing ', trim(name), & 
          '. Stopping.'
     call LIS_endrun
  endif
  
  
  return

end subroutine read_geos
