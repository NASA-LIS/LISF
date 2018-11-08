!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readrhoneAGG
! \label{readrhoneAGG}
!
! !REVISION HISTORY:
!  5 Nov 2003: Dave Mocko, Initial Specification
!
! !INTERFACE:
subroutine readrhoneAGG(n, order,findex, name,tscount,ferror)
! !USES:
  use LIS_coreMod,        only : LIS_rc,LIS_domain
  use LIS_metforcingMod, only : LIS_forc
  use LIS_logMod,         only : LIS_logunit, LIS_verify, LIS_endrun
  use rhoneAGG_forcingMod,   only : rhoneAGG_struc

!
! !DESCRIPTION:
!  Reads in RHONEAGG data and performs interpolation to the LIS domain.
!
! \subsection{Core Functions of readrhoneAGG}
!  \begin{description}
!  \item[bilinear\_interp]
!      Interpolates RHONEAGG data to LIS grid using bilinear interpolation
!  \end{description}
!
! RHONEAGG FORCING VARIABLES (unless noted, fields are 3-hr upstream averaged): \newline
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

!EOP
  implicit none
  integer,      intent(in) :: n
  integer,      intent(in) :: findex
  character*80, intent(in) :: name
  integer, intent(in)      :: order, tscount
  integer, intent(out)     :: ferror
  integer :: mo
  integer :: nforce,i,j,v
  integer :: ios,ioerror    ! set to non-zero if there's an error
  real, allocatable :: temprhoneAGG(:,:) ! lf%nmif = Max.# parameters to retrieve
  real, allocatable :: tempvar(:,:,:)
  integer :: glis,nrhoneAGG     ! Size of I/O 1D fields
  real, allocatable :: f(:),go(:) ! 1D I/O fields
  real, allocatable :: gmask(:) ! RHONEAGG forcing mask
  integer :: count1
  integer :: c
  real :: gridDesco(50)     ! Input,output grid info arrays
  logical*1, allocatable :: lb(:)
  logical*1, allocatable :: lo(:) ! Input and output bitmaps
  integer :: nr_index, nc_index
  
  nr_index = rhoneAGG_struc(n)%nrold
  nc_index = rhoneAGG_struc(n)%ncold
  
  
  allocate(tempvar(nc_index,nr_index,rhoneAGG_struc(n)%nmif),stat=ios)
  call LIS_verify(ios,'Error allocating tempvar.')
  
  allocate(temprhoneAGG(LIS_rc%ngrid(n),rhoneAGG_struc(n)%nmif),stat=ios)
  call LIS_verify(ios,'Error allocating temprhoneAGG.')
  
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
  nrhoneAGG = nc_index*nr_index
  glis = LIS_rc%lnc(n)*LIS_rc%lnr(n)
  ferror = 1
!-----------------------------------------------------------------------
! Open RHONEAGG forcing file
!-----------------------------------------------------------------------

  write(LIS_logunit,*)'MSG: readrhoneAGG -- Reading RHONEAGG forcing file: ', &
       trim(name)
  
  open(40,file=trim(name),form='unformatted',iostat=ios)
  if (ios.ne.0) then
     write(LIS_logunit,*)'ERR: readrhoneAGG -- Error opening: ',trim(name)
     write(LIS_logunit,*)'Stopping.'
     call LIS_endrun
  endif
  
  do i = 1,rhoneAGG_struc(n)%nmif
     read(40,iostat=ioerror) tempvar(:,:,i)
  enddo
  
  if (ioerror.ne.0) then
     write(LIS_logunit,*)'ERR: readrhoneAGG -- Error reading: ',trim(name),'.'
     write(LIS_logunit,*)'ioerror = ',ioerror,'.  Stopping.'
     call LIS_endrun
  else
  endif
!-----------------------------------------------------------------------
! Finding number of forcing variables
!-----------------------------------------------------------------------
  nforce = rhoneAGG_struc(n)%nmif

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
              gmask(c) = f(c) ! Storing rhoneAGG land mask for later use
           endif
        enddo
     enddo
     
!-----------------------------------------------------------------------
! Initializing input and output grid arrays
!-----------------------------------------------------------------------
     gridDesco = 0
!     gridDesco = LIS_rc%gridDesc
     mo = LIS_rc%lnc(n) * LIS_rc%lnr(n)
     
!-----------------------------------------------------------------------
! Defining input data bitmap
!-----------------------------------------------------------------------
     do i = 1,nrhoneAGG
        lb(i) = .true.
     enddo
     
!-----------------------------------------------------------------------
! Alter default bitmap prescribed above for
! surface parameters (soil wetness, snow)
!-----------------------------------------------------------------------
     if ((v.eq.12).or.(v.eq.13)) then
        do i = 1,nrhoneAGG
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
! Interpolate data from RHONEAGG grid to LIS grid
!-----------------------------------------------------------------------
!         if (LIS_rc%f%interp.eq.1) then
!            call bilinear_interp(gridDesco,ibi,lb,f,ibo,lo,go,mi,mo,   &
!               LIS_domain(n)%lat, LIS_domain(n)%lon,&
!               w110,w120,w210,w220,n110,n120,n210,n220,iret)
!         elseif (LIS_rc%f%interp.eq.2) then
!            if ((v.eq.8).or.(v.eq.9)) then
!               call conserv_interp(gridDesco,ibi,lb,f,ibo,lo,go,mi,mo, &
!               LIS_domain(n)%lat, LIS_domain(n)%lon,&
!               w113,w123,w213,w223,n113,n123,n213,n223,iret)
!            else
!               call bilinear_interp(gridDesco,ibi,lb,f,ibo,lo,go,mi,mo,&
!               LIS_domain(n)%lat, LIS_domain(n)%lon,&
!               w110,w120,w210,w220,n110,n120,n210,n220,iret)
!            endif
!         endif

!-----------------------------------------------------------------------
! Convert data to original 3D array & a 2D array to
! fill in of missing points due to geography difference
!-----------------------------------------------------------------------
     count1 = 1
     do j = 1,LIS_rc%lnr(n)
        do i = 1,LIS_rc%lnc(n)
           if (LIS_domain(n)%gindex(i,j).ne.-1) then
              temprhoneAGG(LIS_domain(n)%gindex(i,j),v) = f(count1)
           endif
           count1 = count1 + 1
        enddo
     enddo

!-----------------------------------------------------------------------
! Fill in undefined and ocean points
!-----------------------------------------------------------------------
     do i = 1,LIS_rc%ngrid(n)
        if (temprhoneAGG(i,v).ge.9.9e+14) then
           temprhoneAGG(i,v) = LIS_rc%udef
        endif
        if (order.eq.1) then
           rhoneAGG_struc(n)%metdata1(v,i) = temprhoneAGG(i,v)
        else
           rhoneAGG_struc(n)%metdata2(v,i) = temprhoneAGG(i,v)
        endif
     enddo                  !i
  enddo                     !v

!      write(LIS_logunit,*)'DBG: readrhoneAGG -- Deallocating arrays'

  deallocate(tempvar,stat=ios)
  call LIS_verify(ios,'Error deallocating tempvar.')
  
  deallocate(temprhoneAGG,stat=ios)
  call LIS_verify(ios,'Error deallocating temprhoneAGG.')
  
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

  write(LIS_logunit,*)'MSG: readrhoneAGG -- Closing RHONEAGG forcing file: ', &
       trim(name)
  close(40,iostat=ios)
  if (ios.ne.0) then
     write(LIS_logunit,*)'ERR: readrhoneAGG -- Error closing: ',trim(name),'.'
     write(LIS_logunit,*)'Stopping.'
     call LIS_endrun
  endif
  
  return
end subroutine readrhoneAGG



