!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: rdhm356_read_temper
! \label{rdhm356_read_temper}
!
! !REVISION HISTORY:
!  25 May 2006: Kristi Arsenault;  Data and code implementation
!  19 Dec 2013: Shugong Wang; Impmementation for RDHM 356 
!
! !INTERFACE:
subroutine rdhm356_read_temper( n, findex, order, fname, ferror_rdhm356 )

! !USES:
  use LDT_coreMod,        only : LDT_rc, LDT_domain
  use LDT_logMod,         only : LDT_logunit
  use LDT_metforcingMod,  only : LDT_forc
  use rdhm356_forcingMod, only : rdhm356_struc_temper
  use LDT_xmrg_reader

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  character(len=*)   :: fname     
  integer, intent(in) :: findex
  integer, intent(in) :: order
  integer             :: ferror_rdhm356

! !DESCRIPTION:
!  For the given time, reads RDHM 356 temperature
!  and interpolates to a designated user-domain.
!  NOTE:: These subroutines use the LDT XMRG Reader 
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
!  \item[interp\_rdhm356\_temper](\ref{interp_rdhm356_temper}) \newline
!    spatially interpolates the rdhm356 data
!  \end{description}
!EOP

  integer            :: i, j, ierr ! Loop indicies and error flags
  integer            :: ndata, index1
  real, allocatable      :: temp_regrid(:,:)     
  real, allocatable  :: rdhm356in(:)
  logical*1          :: lb(rdhm356_struc_temper(n)%ncol*rdhm356_struc_temper(n)%nrow)
  integer            :: ksec1(100),c
  logical            :: file_exists  
  real, allocatable  :: xmrg_data(:,:)
  real    :: x_or, y_or
  integer :: xmrg_nrow, xmrg_ncol
  integer :: row, col, kk, xmrg_row, xmrg_col, row_offset, col_offset  


  if(order.eq.1) then 
     LDT_forc(n,findex)%metdata1(2,:) = LDT_rc%udef
  elseif(order.eq.2) then 
     LDT_forc(n,findex)%metdata2(2,:) = LDT_rc%udef
  endif

  ! allocate memory
  allocate ( temp_regrid(LDT_rc%lnc(n), LDT_rc%lnr(n)) )
  temp_regrid = -1.0

  ! check initially if file exists:
  inquire (file=fname, exist=file_exists )
  if (.not. file_exists)  then
    write(LDT_logunit,*)"** Missing temperature file: ", fname
    return
  endif

  ! initialize and allocate arrays
  ndata = rdhm356_struc_temper(n)%ncol * rdhm356_struc_temper(n)%nrow
  allocate ( rdhm356in(ndata) )
  rdhm356in = rdhm356_struc_temper(n)%undef_value 
  !call read_xmrg2temp(trim(fname)//char(0),rdhm356in,nx,ny)
  
  call xmrg_read_header(trim(fname)//char(0), x_origin=x_or, y_origin=y_or, &
                        num_row=xmrg_nrow, num_col=xmrg_ncol)
  call xmrg_read_data(trim(fname)//char(0), xmrg_data)

  kk = 1
  do row=1, rdhm356_struc_temper(n)%nrow
    do col=1, rdhm356_struc_temper(n)%ncol
      ! calculate the offsets of row and column 
      row_offset = nint((rdhm356_struc_temper(n)%lower_left_hrapy - y_or)/ &
                         rdhm356_struc_temper(n)%hrap_resolution)
      col_offset = nint((rdhm356_struc_temper(n)%lower_left_hrapx - x_or)/ &
                         rdhm356_struc_temper(n)%hrap_resolution) 
      xmrg_row = row_offset + row
      xmrg_col = col_offset + col 
      if((xmrg_row .gt. xmrg_nrow) .or. (xmrg_col .gt. xmrg_ncol)) then
        rdhm356in(kk) = LDT_rc%udef
      else
        if(xmrg_data(xmrg_row, xmrg_col) .eq. rdhm356_struc_temper(n)%undef_value) then
          rdhm356in(kk) = LDT_rc%udef
        else
          rdhm356in(kk) = (xmrg_data(xmrg_row, xmrg_col)-32.0)*5.0/9.0 + 273.16 ! F -> K 
        endif
      endif
      kk = kk + 1
    enddo
  enddo

  ksec1 = 0

  call interp_rdhm356_temper(n, findex, ksec1, ndata, rdhm356in, lb,     &
                             LDT_rc%gridDesc,                            &
                             LDT_rc%lnc(n), LDT_rc%lnr(n), temp_regrid )

  do j = 1, LDT_rc%lnr(n)
     do i = 1, LDT_rc%lnc(n)
        if ( temp_regrid(i,j) .ne. LDT_rc%udef ) then 
           index1 = LDT_domain(n)%gindex(i,j)
           if(index1 .ne. -1) then
              if(order.eq.1) then
                 LDT_forc(n,findex)%metdata1(2,index1) = temp_regrid(i,j) 
              elseif(order.eq.2) then 
                 LDT_forc(n,findex)%metdata2(2,index1) = temp_regrid(i,j) 
              endif
           endif
        endif
     enddo
  enddo
  ferror_rdhm356 = 0 
  write(LDT_logunit,*) "- Obtained RDHM 356 temperation file: ", fname
  write(LDT_logunit,*) " ------------------------------------ "

  deallocate(temp_regrid)
  deallocate(rdhm356in)
  deallocate(xmrg_data) 
end subroutine rdhm356_read_temper

