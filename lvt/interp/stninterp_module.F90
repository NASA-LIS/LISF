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
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
! by Matthew Garcia, GEST Research Associate
!    Hydrological Sciences Branch (Code 614.3)
!    NASA-GSFC 
!
! Calculation of measurement station weights from locations
!   for interpolation processing of measurement time series and
!   -- rasterization (visualization)
!   -- distributed modeling
!
! Primarily for use in conjunction with the
!   NASA-GSFC Land surface Verification Toolkit (LVT) V1.0
!
! Version History
!  xxNov03  Matthew Garcia  Program hpgrid2.f by Matthew Garcia, M.S.
!                                                Dept. of Civil Engineering
!                                                Colorado State University 
!  12Aug04  Matthew Garcia  Program idwgrid specified for ARMS project
!  12Aug04  Matthew Garcia  Delivered to USACE TEC for ARMS preprocessor use  
!  01Mar05  Matthew Garcia  Revised for integration into NASA/GSFC LISv4.1a
!  21Apr05  Matthew Garcia  Added MQ-B interpolation functionality and routines
!
! Notes:
!   Array declarations use the following convention:  Name(stnnum:cols:rows)
!
module stninterp_module
!
  use idw_module
  use mqb_module
!
  implicit none
!
  type interpgridsdec
    integer, allocatable :: stngrid(:,:,:)
    real, allocatable :: wtsgrid(:,:,:)
    real, allocatable :: mqbstns(:,:)
    real, allocatable :: mqbgrid(:,:)
 end type interpgridsdec
 !
 type(interpgridsdec) :: interpgrids
 !
  contains
!
!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
subroutine calcgrids(dirname,stnfile,stns,nnn,method,order,cols,rows,&
                     llx,lly,urx,ury,incr,npts,locarr)
!
  implicit none
!
! argument variables
  character*40, intent(IN) :: dirname,stnfile
  integer, intent(INOUT) :: stns,nnn,method,order
  integer, intent(INOUT) :: cols,rows
  real*4, intent(INOUT) :: llx,lly,urx,ury,incr
  integer, intent(INOUT), optional :: npts
  real*4, intent(INOUT), optional :: locarr(:,:)
!
! local variables
  integer :: s
  logical :: lerr = .FALSE.
  character*80 :: infile
  character*80 :: outfile
  integer :: funit = 10
  real(4), allocatable :: stndata(:,:)
  integer(4), allocatable :: W1grid(:,:,:) ! station numbers for each grid cell
  real(4), allocatable :: W2grid(:,:,:) ! station weights for each grid cell 
  real(4), allocatable :: W3grid(:,:) ! MQ-B station matrix for simulation domain 
!
  call checkvalues(lerr,stns,nnn,method,order)
  if (lerr) then
    print *,'ERR: calcgrids -- error returned from checkvalues'
    stop
  else
    allocate(stndata(3,stns))
    stndata = 0.0
    infile = trim(adjustl(dirname))//'/'//trim(adjustl(stnfile))
    open(unit=funit,file=infile,status='old')
    do s = 1,stns
      read(funit,'(F3.0,1X,2(F10.2,1X))') &
           stndata(1,s),stndata(2,s),stndata(3,s)
!       gauge ID     x-coord      y-coord
    end do
    close(funit)
    do s = 1,stns
      if (stndata(2,s).lt.llx .or. stndata(2,s).gt.urx) then
        print *,'WRN: calcgrids -- station',stndata(1,s),' occurs outside the specified domain.'
!        print *,' x =',stndata(2,s),' llx =',llx,' urx =',urx
!        stop
      endif
      if (stndata(3,s).lt.lly .or. stndata(3,s).gt.ury) then
        print *,'WRN: calcgrids -- station',stndata(1,s),' occurs outside the specified domain.'
!        print *,' y =',stndata(3,s),' lly =',lly,' ury =',ury
!        stop
      endif
    end do
    if (cols.ne.(urx - llx)/incr) then
      print *,'ERR: interp_module -- inconsistent domain, cell size and # of columns.'
      stop  
    endif
    if (rows.ne.(ury - lly)/incr) then
      print *,'ERR: interp_module -- inconsistent domain, cell size and # of rows.'
    stop
    endif
!
    if (method.eq.1 .or. method.eq.2) then
      allocate(W1grid(nnn+1,cols,rows))
      allocate(W2grid(nnn+1,cols,rows))
      W1grid = 0
      W2grid = 0.0
      call interpmethod(method,stns,nnn,cols,rows,llx,ury,incr,order,stndata,&
                        W2grid,W1grid)
    else if (method.eq.3) then
      allocate(W2grid(nnn+1,cols,rows))
      W2grid = 0.0
      call interpmethod(method,stns,nnn,cols,rows,llx,ury,incr,order,stndata,&
                        W2grid)
    else if (method.eq.4) then
      allocate(W1grid(1,1,1))
      allocate(W2grid(1,npts,nnn+1))
      allocate(W3grid(nnn+1,nnn+1))
      W1grid = 0
      W2grid = 0.0
      W3grid = 0.0
      call interpmethod(method,stns,nnn,cols,rows,llx,ury,incr,order,stndata,&
                        W2grid,W1grid,W3grid,npts,locarr)
    end if
!
    if (method.eq.1 .or. method.eq.2) then
      print *,'MSG: interp_module -- writing gridded station assignments to data file.'
      outfile = trim(adjustl(dirname))//'/stngrid.bin' 
      open(unit=30,file=outfile,status='replace',form='unformatted')
      write(30) W1grid
      close(30)
      deallocate(W1grid)
    end if
    print *,'MSG: interp_module -- writing gridded station weights to data file.'
    outfile = trim(adjustl(dirname))//'/stnwgts.bin' 
    open(unit=30,file=outfile,status='replace',form='unformatted')
    write(30) W2grid
    close(30)
    deallocate(W2grid)
    if (method.eq.4) then
      print *,'MSG: interp_module -- writing MQ-B station matrix to data file.'
      outfile = trim(adjustl(dirname))//'/stnmtx.bin' 
      open(unit=30,file=outfile,status='replace',form='unformatted')
      write(30) W3grid
      close(30)
      deallocate(W3grid)
    end if
  end if
!
  return
end subroutine calcgrids
!
!
!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
subroutine checkvalues(error,ns,nn,mthd,ord)
!
! argument variables
  logical, intent(INOUT) :: error
  integer, intent(IN) :: ns,mthd,ord
  integer, intent(INOUT) :: nn
!
! local variables
!
  if (nn.le.2.or.nn.gt.ns) then
    print *,'ERR: checkvalues -- # of nearest neighbors must satisfy 2 < Nnn < stns.'
!    print *,'   Reconfigure station data file and try again.'
    close(10)
    stop
  endif
!  print *,ns,' stations'
  select case (mthd)
    case (1)
      nn = 4
!     print *,nn,' weighted quadrant neighbors selected.'
    case (2) 
!     print *,nn,' weighted nearest neighbors selected.'
    case (3)
!     print *,'Full-domain weights (stns = nnn+1) selected.'
      if (ns.ne.nn+1) then
        print *,'ERR: checkvalues -- check specification of stns and nnn'
        close(10)
        stop
      endif
    case (4) 
!     print *,'Full-domain MQ-B method selected.'
    case default
      print *,'ERR: checkvalues -- check specification of interpolation method'
      close(10)
      stop
  end select
  select case (ord)
    case (1)
!     print *,'1st order (linear) IDW.'
    case (2)
!     print *,'2nd order (quadratic) IDW.'
    case (3)
!      print *,'3rd order (cubic) IDW.'
    case default
      print *,'ERR: checkvalues -- check specification of IDW order'
!     print *,' -- up to 3rd order only in current program version.'
      close(10)
      stop
  end select
!
  return
end subroutine checkvalues
!
!
!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
subroutine interpmethod(mthd,ns,nn,nc,nr,lx,ry,inc,ord,sdata,wgtarr,stnarr,stnmtx,npts,locarr)
!
! argument variables
  integer, intent(IN) :: mthd
  integer, intent(INOUT) :: ns,nn,nc,nr,ord
  real(4), intent(INOUT) :: lx,ry,inc
  real(4), intent(INOUT) :: sdata(:,:)
  real(4), intent(INOUT) :: wgtarr(:,:,:) ! station weights for each grid cell 
  integer(4), intent(INOUT), optional :: stnarr(:,:,:) ! station numbers for each grid cell
  real(4), intent(INOUT), optional :: stnmtx(:,:) ! MQ-B station matrix 
  integer, intent(INOUT), optional :: npts
  real(4), intent(INOUT), optional :: locarr(:,:) ! domain location indices
!
  select case (mthd)
    case (1)
      call idwquads(ns,nn,nc,nr,lx,ry,inc,ord,sdata,stnarr,wgtarr)
      call normalize(nn,nc,nr,wgtarr)
    case (2)
      call idwranked(ns,nn,nc,nr,lx,ry,inc,ord,sdata,stnarr,wgtarr)
      call normalize(nn,nc,nr,wgtarr)
    case (3)
      call idwalln(ns,nc,nr,lx,ry,inc,ord,sdata,wgtarr)
      call normalize(nn,nc,nr,wgtarr)
    case (4)
      call mqbmatrix(nn,npts,lx,ry,inc,sdata,stnmtx,wgtarr,locarr)
  end select
!
  return
end subroutine interpmethod
!
!
end module stninterp_module
