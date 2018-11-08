!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: rdhm356_read_precip
! \label{rdhm356_read_precip}
!
! !REVISION HISTORY:
!  25 May 2006: Kristi Arsenault;  Data and code implementation
!  18 Dec 2013: Shugong Wang; implementation for RDHM 356
!  
! !INTERFACE:
subroutine rdhm356_read_precip(n, fname, findex, order, ferror_rdhm356)

! !USES:
  use LIS_coreMod,        only : LIS_rc, LIS_domain
  use LIS_logMod,         only : LIS_logunit, LIS_endrun
  use LIS_metforcingMod,  only : LIS_forc
  use rdhm356_forcingMod, only : rdhm356_struc_precip
  use LIS_XMRG_READER

 implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  character(len=*)   :: fname     
  integer, intent(in) :: findex
  integer, intent(in) :: order
  integer             :: ferror_rdhm356

! !DESCRIPTION:
!  For the given time, reads RDHM 356 precipitation 
!  and interpolates to a designated user-domain.
!  NOTE:: These subroutines use the LIS XMRG Reader 
!         for opening and reading the xmrg files.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[fname]
!    name of the hourly xmrg file
!  \item[ferror\_rdhm356]
!    flag to indicate success of the call (=0 indicates success)
!  \item[filehr]
!    current file hour
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_rdhm356\_precip](\ref{interp_rdhm356_precip}) \newline
!    spatially interpolates the rdhm356 data
!  \end{description}
!EOP

  integer            :: i, j,ierr ! Loop indicies and error flags
  integer            :: ndata, index1
  real, allocatable      :: precip_regrid(:,:) 
  real, allocatable  :: rdhm356in(:)
  logical*1          :: lb(rdhm356_struc_precip(n)%ncol*rdhm356_struc_precip(n)%nrow)
  integer            :: ksec1(100),c
  logical            :: file_exists  
  real, allocatable  :: xmrg_data(:,:)
  real    :: x_or, y_or, nodata_value, cell_size 
  integer :: xmrg_nrow, xmrg_ncol, len_rec2, status
  integer :: row, col, kk, xmrg_row, xmrg_col, row_offset, col_offset  

  if(order.eq.1) then 
     rdhm356_struc_precip(n)%metdata1(:) = LIS_rc%udef
  elseif(order.eq.2) then 
     rdhm356_struc_precip(n)%metdata2(:) = LIS_rc%udef
  endif

  ! allocate memory 
  allocate ( precip_regrid(LIS_rc%lnc(n), LIS_rc%lnr(n)) )
  precip_regrid = -1.0

  ! check initially if file exists:
  inquire (file=fname, exist=file_exists )
  if (.not. file_exists)  then
    write(LIS_logunit,*)"** Missing precipitation file: ", fname
    return      
  endif

  ! initialize and allocate arrays
  ndata = rdhm356_struc_precip(n)%ncol * rdhm356_struc_precip(n)%nrow
  allocate ( rdhm356in(ndata) )
  rdhm356in = rdhm356_struc_precip(n)%undef_value
  ! call read_xmrg2(trim(fname)//char(0),rdhm356in)
  
  call xmrg_read_header(trim(fname)//char(0), x_origin=x_or, y_origin=y_or,   &
                        num_row=xmrg_nrow, num_col=xmrg_ncol,        &
                        nodata_value=nodata_value, len_rec2=len_rec2,& 
                        cell_size=cell_size, status=status)
  call xmrg_read_data(trim(fname)//char(0), xmrg_data, xmrg_nrow=xmrg_nrow, xmrg_ncol=xmrg_ncol, &
                      status=status)
  
!  open(unit=1001, file='rdhm356_in', status='unknown');
  kk = 1
  do row=1, rdhm356_struc_precip(n)%nrow
    do col=1, rdhm356_struc_precip(n)%ncol
      ! calculate the offsets of row and column 
      row_offset = nint((rdhm356_struc_precip(n)%lower_left_hrapy - y_or)/ &
                         rdhm356_struc_precip(n)%hrap_resolution)
      col_offset = nint((rdhm356_struc_precip(n)%lower_left_hrapx - x_or)/ &
                         rdhm356_struc_precip(n)%hrap_resolution) 
      
      xmrg_row = row_offset + row
      xmrg_col = col_offset + col 
                 
      if((xmrg_row .gt. xmrg_nrow) .or. (xmrg_col .gt. xmrg_ncol)) then
        rdhm356in(kk) = LIS_rc%udef
      else
        if((xmrg_data(xmrg_row, xmrg_col) .eq. rdhm356_struc_precip(n)%undef_value) &
           .or. (xmrg_data(xmrg_row, xmrg_col)<0)) then
          rdhm356in(kk) = LIS_rc%udef
        else
          rdhm356in(kk) = xmrg_data(xmrg_row, xmrg_col)/3600.0  ! mm/hr -> kg/m2s
        endif
      endif
!      write(1001, *) row, col, rdhm356in(kk), xmrg_row, xmrg_col, xmrg_data(row, col)
      kk = kk + 1
    enddo
  enddo
!  close(1001)


  ksec1 = 0

  call interp_rdhm356_precip(n, findex, ksec1, ndata, rdhm356in, lb,     &
                             LIS_rc%gridDesc,                            &
                             LIS_rc%lnc(n), LIS_rc%lnr(n), precip_regrid )

!  open(unit=1001, file='precip_regrid', status='unknown');
  do j = 1, LIS_rc%lnr(n)
     do i = 1, LIS_rc%lnc(n)
        if ( precip_regrid(i,j) .ge. LIS_rc%udef ) then
           index1 = LIS_domain(n)%gindex(i,j)
           if(index1 .ne. -1) then
              if(order.eq.1) then
                 rdhm356_struc_precip(n)%metdata1(index1) = precip_regrid(i,j) 
              elseif(order.eq.2) then 
                 rdhm356_struc_precip(n)%metdata2(index1) = precip_regrid(i,j) 
              endif
              xmrg_row = row_offset + j
              xmrg_col = col_offset + i 
!              write(1001, '(3I9,2F16.6,3I6)') j, i, index1, precip_regrid(i,j), xmrg_data(xmrg_row, xmrg_col), xmrg_row, xmrg_col, order
           endif
        endif
     enddo
  enddo
!  close(unit=1001)
!  open(unit=1001, file='precip_regrid2', status='unknown');
!  write(unit=1001, fmt='(1043F14.6)') ((precip_regrid(i, j), i=1,1043), j=1, 774)
!  close(unit=1001)
  !call LIS_endrun();  
  ferror_rdhm356 = 0 
  write(LIS_logunit,*) "- Obtained precipitation file: ", fname
  write(LIS_logunit,*) " ------------------------------------ "

  deallocate(precip_regrid)
  deallocate(rdhm356in)
  deallocate(xmrg_data) 
end subroutine rdhm356_read_precip

