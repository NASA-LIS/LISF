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
! !ROUTINE: read_gswp1
! \label{read_gswp1}
!
! !REVISION HISTORY:
! 11 Dec 2003: Sujay Kumar, Initial Specification
!
! !INTERFACE:
subroutine read_gswp1(order,n,findex, name,tscount,ferror)
! !USES:
  use LIS_coreMod,        only : LIS_rc, LIS_domain, LIS_localPet
  use LIS_metforcingMod, only : LIS_forc
  use LIS_logMod, only         : LIS_logunit, LIS_endrun, LIS_verify
  use gswp1_forcingMod,   only : gswp1_struc

!
! !DESCRIPTION:
!  Reads in GSWP-1 data and performs interpolation to the LIS domain.
!
! \subsection{Core Functions of read\_gswp1}
!  \begin{description}
!  \item[bilinear\_interp]
!      Interpolates GSWP-1 data to LIS grid using bilinear interpolation
!  \end{description}
!
! GSWP-1 FORCING VARIABLES (unless noted, fields are 3-hr upstream averaged): \newline
!  1. T 2m    Temperature interpolated to 2 metres [$K$] \newline
!  2. q 2m    Instantaneous specific humidity interpolated to 2 metres[$kg/kg$] \newline
!  3. radswg  Downward shortwave flux at the ground [$W/m^2$] \newline
!  4. lwgdwn  Downward longwave radiation at the ground [$W/m^2$] \newline
!  5. wind speed  wind speed \newline
!  6. wind east   dummy field, always 0 \newline
!  7. ps      Instantaneous Surface Pressure [$Pa$] \newline
!  8. preacc  Total precipitation [$mm/s$] \newline
!  9. precon  Convective precipatation [$mm/s$] \newline
!

!EOP
  implicit none
  integer,      intent(in) :: n 
  integer,      intent(in) :: findex
  character(len=*), intent(in) :: name
  integer, intent(in)      :: order, tscount
  integer, intent(out)     :: ferror
  integer :: mo
  integer :: nforce,i,j,v
  integer :: ios,ioerror    ! set to non-zero if there's an error
  real, allocatable :: tempgswp1(:,:) ! lf%nmif = Max.# parameters to retrieve
  real, allocatable :: tempvar(:,:,:)
  integer :: glis,ngswp1     ! Size of I/O 1D fields
  real, allocatable :: f(:),go(:) ! 1D I/O fields
  real, allocatable :: gmask(:) ! GSWP-1 forcing mask
  integer :: count1
  integer :: iret,c
  real :: gridDesco(50)     ! Input,output grid info arrays
  logical*1, allocatable :: lb(:)
  logical*1, allocatable :: lo(:) ! Input and output bitmaps


  integer :: nr_index, nc_index
  nr_index = gswp1_struc(n)%nrold
  nc_index = gswp1_struc(n)%ncold

!-----------------------------------------------------------------------
! Finding number of forcing variables
!-----------------------------------------------------------------------
  nforce = LIS_rc%met_nf(findex)

  allocate(tempvar(nc_index,nr_index,nforce),stat=ios)
  call LIS_verify(ios,'Error allocating tempvar.')
  
  allocate(tempgswp1(LIS_rc%ngrid(n),nforce),stat=ios)
  call LIS_verify(ios,'Error allocating tempgswp1.')
  
  allocate(f(nc_index*nr_index),stat=ios)
  call LIS_verify(ios,'Error allocating f.')
  
  allocate(go(LIS_rc%lnc(n)*LIS_rc%lnr(n)),stat=ios)
  call LIS_verify(ios,'Error allocating go.')
  
  allocate(gmask(nc_index*nr_index),stat=ios)
  call LIS_verify(ios,'Error allocating gmask.')
  
  allocate(lb(nc_index*nr_index),stat=ios)
  call LIS_verify(ios,'Error allocating lb.')
  
  allocate(lo(LIS_rc%lnc(n)*LIS_rc%lnr(n)),stat=ios)
  call LIS_verify(ios,'Error allocating lo.')
  
  gmask = 0.0
  ngswp1 = nc_index*nr_index
  glis = LIS_rc%lnc(n)*LIS_rc%lnr(n)
  ferror = 1
!-----------------------------------------------------------------------
! Open GSWP-1 forcing file
!-----------------------------------------------------------------------

  write(LIS_logunit,*)'MSG: read_gswp1 -- Reading GSWP-1 forcing file: ', &
       trim(name),' (',LIS_localPet,')'
  
  open(40,file=name,form='unformatted',iostat=ios)
  if (ios.ne.0) then
     write(LIS_logunit,*)'ERR: read_gswp1 -- Error opening: ',trim(name),'.'
     write(LIS_logunit,*)'Stopping.',' (', LIS_localPet,')'
     call LIS_endrun
  endif
  

  do i = 1,nforce
     read(40,iostat=ioerror) tempvar(:,:,i)
  enddo

  if (ioerror.ne.0) then
     write(LIS_logunit,*)'ERR: read_gswp1 -- Error reading: ',trim(name),'.'
     write(LIS_logunit,*)'ioerror = ',ioerror,'.  Stopping.',' (', LIS_localPet,')'
     call LIS_endrun
  endif

!-----------------------------------------------------------------------
! Transferring current data to 1-D array for interpolation
!-----------------------------------------------------------------------
  do v = 1,nforce
     c=0
     do i=1,nr_index        !i=1,lf%nrold
        do j=1,nc_index     !j=1,lf%ncold
           c = c + 1
           f(c) = tempvar(j,i,v)
           if ((tscount.eq.0).and.(order.eq.1).and.(v.eq.9)) then
              gmask(c) = f(c) ! Storing gswp1 land mask for later use
           endif
        enddo
     enddo
     
!-----------------------------------------------------------------------
! Initializing input and output grid arrays
!-----------------------------------------------------------------------
     gridDesco = 0
     gridDesco = LIS_rc%gridDesc(n,:)
     mo = LIS_rc%lnc(n) * LIS_rc%lnr(n)

!-----------------------------------------------------------------------
! Defining input data bitmap
!-----------------------------------------------------------------------
     do i = 1,ngswp1
        lb(i) = .true.
     enddo

!-----------------------------------------------------------------------
! Alter default bitmap prescribed above for
! surface parameters (soil wetness, snow)
!-----------------------------------------------------------------------
     if ((v.eq.12).or.(v.eq.13)) then
        do i = 1,ngswp1
           if (gmask(i).le.-998.0) then
              lb(i) = .false.
           else
              lb(i) = .true.
           endif
        enddo
     endif

!-----------------------------------------------------------------------
! Defining output data bitmap
!-----------------------------------------------------------------------
     do i = 1,glis
        lo(i) = .true.
     enddo
     
!-----------------------------------------------------------------------
! Interpolate data from GSWP-1 grid to LIS grid
!-----------------------------------------------------------------------
     if (trim(LIS_rc%met_interp(findex)).eq."bilinear") then
        call bilinear_interp(gridDesco,lb,f,lo,go,gswp1_struc(n)%mi,mo,   &
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             gswp1_struc(n)%w111,gswp1_struc(n)%w121,&
             gswp1_struc(n)%w211,gswp1_struc(n)%w221,&
             gswp1_struc(n)%n111,gswp1_struc(n)%n121,&
             gswp1_struc(n)%n211,gswp1_struc(n)%n221,LIS_rc%udef, iret)
     elseif (trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then
        if ((v.eq.8).or.(v.eq.9)) then
           call conserv_interp(gridDesco,lb,f,lo,go,gswp1_struc(n)%mi,mo, &
                LIS_domain(n)%lat, LIS_domain(n)%lon,&
                gswp1_struc(n)%w112,gswp1_struc(n)%w122,&
                gswp1_struc(n)%w212,gswp1_struc(n)%w222,&
                gswp1_struc(n)%n112,gswp1_struc(n)%n122,&
                gswp1_struc(n)%n212,gswp1_struc(n)%n222,LIS_rc%udef, iret)
        else
           call bilinear_interp(gridDesco,lb,f,lo,go,gswp1_struc(n)%mi,mo,   &
                LIS_domain(n)%lat, LIS_domain(n)%lon,&
                gswp1_struc(n)%w111,gswp1_struc(n)%w121,&
                gswp1_struc(n)%w211,gswp1_struc(n)%w221,&
                gswp1_struc(n)%n111,gswp1_struc(n)%n121,&
                gswp1_struc(n)%n211,gswp1_struc(n)%n221,LIS_rc%udef, iret)
        endif
     endif

!-----------------------------------------------------------------------
! Convert data to original 3D array & a 2D array to
! fill in of missing points due to geography difference
!-----------------------------------------------------------------------
     count1 = 0
     do j = 1,LIS_rc%lnr(n)
        do i = 1,LIS_rc%lnc(n)
           if (LIS_domain(n)%gindex(i,j).ne.-1) then
              tempgswp1(LIS_domain(n)%gindex(i,j),v) = go(i+count1)
           endif
        enddo
        count1 = count1 + LIS_rc%lnc(n)
     enddo

!-----------------------------------------------------------------------
! Fill in undefined and ocean points
!-----------------------------------------------------------------------
     do i = 1,LIS_rc%ngrid(n)
        if (tempgswp1(i,v).ge.9.9e+14) then
           tempgswp1(i,v) = LIS_rc%udef
        endif
        if (order.eq.1) then
           gswp1_struc(n)%metdata1(v,i) = tempgswp1(i,v)
        else
           gswp1_struc(n)%metdata2(v,i) = tempgswp1(i,v)
        endif
     enddo                  !i
  enddo                     !v

!      write(LIS_logunit,*)'DBG: read_gswp1 -- Deallocating arrays'
  
  deallocate(tempvar,stat=ios)
  call LIS_verify(ios,'Error deallocating tempvar.')
  
  deallocate(tempgswp1,stat=ios)
  call LIS_verify(ios,'Error deallocating tempgswp1.')
  
  deallocate(f,stat=ios)
  call LIS_verify(ios,'Error deallocating f.')
  
  deallocate(go,stat=ios)
  call LIS_verify(ios,'Error deallocating go.')
  
  deallocate(gmask,stat=ios)
  call LIS_verify(ios,'Error deallocating gmask.')
  
  deallocate(lb,stat=ios)
  call LIS_verify(ios,'Error deallocating lb.')
  
  deallocate(lo,stat=ios)
  call LIS_verify(ios,'Error deallocating lo.')
  
  write(LIS_logunit,*)'MSG: read_gswp1 -- Closing GSWP-1 forcing file: ', &
       trim(name),' (',LIS_localPet,')'
  close(40,iostat=ios)
  if (ios.ne.0) then
     write(LIS_logunit,*)'ERR: read_gswp1 -- Error closing: ',trim(name),'.'
     write(LIS_logunit,*)'Stopping.',' (', LIS_localPet, ')'
     call LIS_endrun
  endif
  
  !      write(LIS_logunit,*)'DBG: read_gswp1 -- leaving',' (', LIS_localPet, ')'
  
  return
end subroutine read_gswp1


