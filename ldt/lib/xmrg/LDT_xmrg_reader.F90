!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_xmrg_reader
!BOP
!
! !MODULE: LDT_xmrg_reader
!
! This module defines routines to read NOAA XMRG and XMRG-like data. 
! It automatically recognize big endian and little endian, compressed (gzip)
! and uncompressed (plain binary), XMRG and XMRG-like formats. All meta data 
! XMRG and XMRG-like data formats are defined as module variables, which are
! described as below:
!
!
! \begin{description}
!   \item[xor]  hrap-x coordinate of southwest corner of grid         
!   \item[yor]  hrap-x coordinate of southwest corner of grid
!   \item[maxx] number of hrap grid boxes in x direction 
!   \item[maxy] number of hrap grid boxes in y direction
!   \item[scale\_factor] scale factor to convert integer into real number (XMRG-like)
!   \item[cell\_size] scale factor of cell size (XMRG-like)
!   \item[nodata\_value] default value for missing data (XMRG-like)
!   \item[num\_byte] number of bytes of data (2 for short integer and 4 for real, XMRG-like)
!   \item[reverse\_byte] flag for little endian (0) and big endian (1)
!   \item[rec2\_len] length of the second record for format recognition 
! \end{description}
! 
! !REVISION HISTORY:
! 28 Nov 2013: Shugong Wang, Initial Code
! 02 Dec 2013: Shugong Wang, Add support to Sub-HRAP coordinate
!
    implicit none
    private

! !PUBLIC MEMBER FUNCTIONS:
    public :: xmrg_read_header
    public :: xmrg_read_data
    public :: xmrg_subset_data
    public :: xmrg_hrap_xy
    public :: xmrg_grid_latlon 
    public :: latlon_to_hrap
    public :: hrap_to_latlon
    public :: LDT_transform_xmrgparam
!EOP
    real    :: xor          
    real    :: yor          
    integer :: maxx        
    integer :: maxy        
    real    :: scale_factor
    real    :: size_cell 
    real    :: nodata_default
    integer :: num_byte
    integer :: reverse_byte 
    integer :: rec2_len     

contains

!BOP
! 
! !ROUTINE: xmrg_read_header
! \label{xmrg_read_header}
!
! !INTERFACE: 
  subroutine xmrg_read_header(xmrg_name,     &
                              x_origin,      &
                              y_origin,      &
                              num_row,       &
                              num_col,       &
                              nodata_value,  &
                              len_rec2,      &
                              cell_size,     &
                              status)
!
! !DESCRIPTION:
!
!  This routine reads the meta information of a XMRG or XMRG-like file
!  from its headers. It retrieves the HRAP origin, numbers of colum and 
!  rows, and nodata value of XMRG or XMRG-lile data. 
!
!  Description of arguments:
!  \begin{description}
!    \item[xmrg\_name] file name of XMRG/XMRG-like file
!    \item[x\_origin]  HRAP X coordiante of the lower left corner of data domain, optional output
!    \item[y\_origin]  HRAP Y coordinate of the lower left corner of data domain, optional output
!    \item[num\_row]   row number of data domain, optional output
!    \item[num\_col]   col number of data domain, optional output
!    \item[nodata\_value] no-data value, optional output
!    \item[len\_rec2] length of the second record, optional output
!    \item[cell\_size] HRAP unit of cell size, optional output 
!    \item[status] return status: 0 - success, 1 - fail
!  \end{description}
! 
!  The routine invoked is:
!  \begin{description}
!    \item{xmrg\_read\_header\_c} C function reading XMRG/XMRG-like
!      header records (1st and 2nd records). (\ref{xmrg_read_header_c}
!  \end{description}
!EOP

   implicit none 

   character(*), intent(in)       :: xmrg_name
   real, intent(out),    optional :: x_origin
   real, intent(out),    optional :: y_origin
   integer, intent(out), optional :: num_row
   integer, intent(out), optional :: num_col
   real, intent(out),    optional :: nodata_value
   integer, intent(out), optional :: len_rec2
   real, intent(out),    optional :: cell_size
   integer, intent(out), optional :: status

   integer :: xmrg_read_header_c
   integer :: rv

   rv = xmrg_read_header_c(trim(xmrg_name)//char(0),    &
                           xor, yor, maxx, maxy,        &
                           size_cell, scale_factor,     &
                           nodata_default, num_byte,    &
                           reverse_byte, rec2_len)
   if(present(status)) status = rv
   if(present(x_origin)) x_origin = xor
   if(present(y_origin)) y_origin = yor
   if(present(num_row)) num_row = maxy
   if(present(num_col)) num_col = maxx
   if(present(nodata_value)) nodata_value = nodata_default
   if(present(len_rec2)) len_rec2 = rec2_len
   if(present(cell_size)) cell_size = size_cell 

 end subroutine xmrg_read_header

!BOP
! 
! !ROUTINE: xmrg_read_data
! \label{xmrg_read_data}
!
! !INTERFACE: 
  subroutine xmrg_read_data(xmrg_name,   &
                            xmrg_data,   &
                            xmrg_nrow,   &
                            xmrg_ncol,   &
                            nodata_value,&
                            status)
!
! !DESCRIPTION:
!
!  This routine reads the data of a XMRG/XMRG-like file. It retrieves
!  data in a 2D array. The row number and the column number of XMRG/XMRG-lie
!  data are optional output. No-data value is an optional input to exclude
!  no-data values from scaling operation. Default no-data value is -1 if 
!  nodata\_value is not presented. 
!
!  Description of arguments:
!  \begin{description}
!    \item[xmrg\_name] name of XMRG/XMRG-like file
!    \item[xmrg\_data] 2D array of XMRG/XMRG-like data starting from the 
!                      lower left corner of the data domain
!    \item[xmrg\_nrow] row number of XMRG/XMRG-like data, optional output 
!    \item[xmrg\_ncol] column number of XMRG/XMRG-like data, optional output
!    \item[nodata\_value] no-data value of XMRG/XMRG-like data, optional input. 
!                         default value is -1 when not explicitly specified. 
!    \item[status] return status: 0 - success, 1 - fail
!  \end{description}
! 
!  The routine invoked is:
!  \begin{description}
!    \item[xmrg\_read\_data\_c](\ref{xmrg_read_data_c}) C function reading data from XMRG/XMRG-like
!     file. 
!    \item[xmrg\_read\_header] (\ref{xmrg_read_header}) read row and column numbers of XMRG/XMRG-lie data 
!  \end{description}
!EOP

   implicit none 
   character(*), intent(in)       :: xmrg_name
   real, allocatable, intent(out) :: xmrg_data(:,:)
   integer, optional, intent(out) :: xmrg_nrow
   integer, optional, intent(out) :: xmrg_ncol
   integer, optional, intent(out) :: status 
   real, optional, intent(in)     :: nodata_value 

   real, allocatable :: temp_data(:)
   integer :: xmrg_read_data_c
   integer :: row, col, k
   integer :: return_status
   real    :: nodata_default
! ________________________________________________________

   call xmrg_read_header(xmrg_name, nodata_value=nodata_default, status=return_status)

   if(return_status .eq. 1) then
      if(present(status)) status = return_status
      return
   else
      if(rec2_len .eq. 16) then
         if(present(nodata_value) .and. abs(nodata_value-nodata_default)>1E-5) then
            write(*,*) "xmrg_read_data fatal error: input no-data value is different XMRG-like default no-data value"
            write(*,*) "Input no-data value: ", nodata_value
            write(*,*) "Default no-data value: ", nodata_default
         endif
      else
         write(*,*) "xmrg_read_data: default no-data value (-1) is used when reading ", &
                    trim(xmrg_name)
         nodata_default = -1.0
      endif

      if(present(xmrg_nrow)) xmrg_nrow = maxy
      if(present(xmrg_ncol)) xmrg_ncol = maxx

      if(allocated(xmrg_data)) then
         deallocate(xmrg_data)
      end if

      allocate(xmrg_data(maxy , maxx))
      allocate(temp_data(maxx * maxy))

      return_status = xmrg_read_data_c(trim(xmrg_name)//char(0), temp_data, nodata_default)
      if(present(status)) status = return_status
      if(return_status .eq. 1) return
      k = 1
      do row=1, maxy
         do col=1, maxx
            xmrg_data(row, col) = temp_data(k)
            k = k + 1
         end do
      end do
   endif

   deallocate(temp_data)

  end subroutine xmrg_read_data

!BOP
! 
! !ROUTINE: xmrg_subset_data
! \label{xmrg_subset_data}
! 
! INTERFACE:
    subroutine xmrg_subset_data(xmrg_name,   &
                                subset_data, &
                                ll_hrap_x,   &
                                ll_hrap_y,   &
                                sub_nrow,    &
                                sub_ncol,    &
                                cell_size,   &
                                nodata_value,&
                                status)
! !DESCRIPTION:
! 
! This routine read a subset of XMRG/XMRG-like data. Caller need to 
! specify the HRAP coordinate of the low left coorner of the subset
! and the row and column numbers of the subset.
! 
! \begin{description}
!  \item[xmrg_name] name of XMRG/XMRG-like file
!  \item[xmrg_data] subset data, allocable 2D real array
!  \item[ll_hrap_x] HRAP X coordiante of the lower left corner of subset domain, real number
!  \item[ll_hrap_y] HRAP Y coordinate of the lower left corner of subset domain, real number
!  \item[sub_nrow] row number of subset domain, integer
!  \item[sub_ncol] col number of subset domain, integer
!  \item[resolution] cell size of XMRG/XMRG-like data, real number, optional 
!  \item[nodata_value] no-data value, real number, optional
!  \item[status] status of subsetting XMRG/XMRG-like data, optional, integer, 0: success, 1: fail
! \end{description}
!
! This routine invokes xmrg_read_data and xmrg_read_header  
!EOP

    character(len=*), intent(in)  :: xmrg_name
    real, allocatable, intent(out):: subset_data(:,:)
    real, intent(in)              :: ll_hrap_x
    real, intent(in)              :: ll_hrap_y
    integer, intent(in)           :: sub_nrow
    integer, intent(in)           :: sub_ncol
    real, optional, intent(in)    :: cell_size
    real, optional, intent(in)    :: nodata_value
    integer, intent(out), optional:: status

    real, allocatable :: xmrg_data(:,:)
    real              :: nodata, xor, yor, cell_size0, res
    integer           :: xmrg_nrow, xmrg_ncol, return_status
    integer           :: start_row, start_col, end_row, end_col
! __________________________________________________________________

    if(present(nodata_value)) then
        nodata = nodata_value
    else
        nodata = -1
    endif

    if(present(cell_size)) then
        res = cell_size
    else
        res = 1.0
    endif

    call xmrg_read_header(trim(xmrg_name), x_origin=xor, y_origin=yor, &
                          cell_size=cell_size0, nodata_value=nodata,  status=return_status)

    if(return_status .eq. 0) then
        if(present(status)) status = 0
        if(abs(res-cell_size0)>1E-6) then
            write(*,*) "xmrg_subset_data fatal error: input cell_size is different from XMRG cell size"
            status = 1
            return
        endif
    else
        if(present(status)) status = 1
        write(*,*) "xmrg_subset_data fatal error: call xmrg_read_header fail"
        return
    endif

    call xmrg_read_data(trim(xmrg_name),            &
                        xmrg_data,                  &
                        xmrg_nrow=xmrg_nrow,        &
                        xmrg_ncol=xmrg_ncol,        &
                        nodata_value=nodata_value,  &
                        status=return_status)

    if(return_status .eq. 0) then
        if(present(status)) status = 0
    else
        if(present(status)) status = 1
        write(*,*) "xmrg_subset_data fatal error: call xmrg_read_data fail"
        return
    endif

    start_row = int((ll_hrap_y-yor)/res) + 1
    start_col = int((ll_hrap_x-xor)/res) + 1
    end_row = start_row + sub_nrow - 1
    end_col = start_col + sub_ncol - 1

    if(allocated(subset_data)) deallocate(subset_data)
    allocate(subset_data(sub_nrow, sub_ncol))

    subset_data(:,:) = xmrg_data(start_row:end_row, start_col:end_col)

    deallocate(xmrg_data)

end subroutine xmrg_subset_data


!BOP
! 
! !ROUTINE: xmrg_hrap_xy
! \label{xmrg_hrap_xy}
!
! !INTERFACE: 
  subroutine xmrg_hrap_xy(xmrg_name, hrapx, hrapy, resolution, status)
!
! !DESCRIPTION:
!
!  This routine retrieves HRAP coordinates (HRAPX, HRAPY) of 
!  XMRG/XMRG-like data boxes. HRAP resolution is an optional input. The default
!  resolution is 1.0 HRAP unit. 
!  
!  Description of arguments: 
!  \begin{description}
!   \item[xmrg\_name] name of XMRG/XMRG-like file
!   \item[hrapx] HRAP X coordinates, allocatable 2D array 
!   \item[hrapy] HRAP Y coordinates, allocatable 2D array 
!   \item[resolution] HRAP resolution of XMRG/XMRG-like data, optional input
!   \item[status] return status: 0 - success, 1 - fail, optional output 
!  \end{description}
!  
!  The routine invoked is:
!  \begin{description}
!    \item[xmrg\_read\_header] (\ref{xmrg_read_header}) Read header information of
!    XMRG/XMRG-like file.
!  \end{description}
!EOP

    implicit none 

    character(*), intent(in)       :: xmrg_name
    real, allocatable, intent(out) :: hrapx(:,:)
    real, allocatable, intent(out) :: hrapy(:,:)
    real, optional, intent(in)     :: resolution
    integer, optional, intent(out) :: status

    integer :: row, col
    integer :: return_status
    real    :: hrap_resolution

  ! retrieve header information 
    call xmrg_read_header(xmrg_name, status=return_status)

    if(return_status .eq. 1) then
       if(present(status)) status = return_status
       return
    else
     ! cell size is specified when rec2 length = 16
       if(rec2_len .eq. 16) then
         if(present(resolution)) then
           if(abs(resolution-size_cell)<1E-5) then
              hrap_resolution = resolution
           else
              write(*, *) "xmrg_hrap_xy fatal error: HRAP resolution is different from XMRG-like cell size"
              write(*, *) "HRAP resolution: ", resolution
              write(*, *) "XMRG-like cell size: ", size_cell
              if(present(status)) status = 1 
              return
           endif
         else
            hrap_resolution = size_cell 
         endif
       else
          write(*,*) "xmrg_hrap_xy: default HRAP resolution (1.0) is used when reading ", &
                     trim(xmrg_name)
          hrap_resolution = 1.0
       endif

       if(allocated(hrapx)) deallocate(hrapx)
       if(allocated(hrapy)) deallocate(hrapy)
       allocate(hrapx(maxy, maxx))
       allocate(hrapy(maxy, maxx))
       do row=1, maxy
          do col=1, maxx
             hrapx(row, col) = xor + (col - 1) * hrap_resolution
             hrapy(row, col) = yor + (row - 1) * hrap_resolution
          end do
       end do
    endif

  end subroutine xmrg_hrap_xy

!BOP
! 
! !ROUTINE: xmrg_grid_latlon
! \label{xmrg_grid_latlon}
!
! !INTERFACE: 
  subroutine xmrg_grid_latlon(xmrg_name, lat, lon, resolution, status)

! !DESCRIPTION:
!
!  This routine retrieves lower left latitude-longitude coordinates of 
!  XMRG/XMRG-like data boxes. HRAP resolution is an optional input. The 
!  default resolution is 1.0 HRAP unit. 
!  
!  Description of arguments: 
!  \begin{description}
!   \item[xmrg\_name] name of XMRG/XMRG-like file
!   \item[lat] latitude of data domain, allocatable 2D array 
!   \item[lon] longitude of data domain, allocatable 2D array 
!   \item[resolution] HRAP resolution of XMRG/XMRG-like data, optional input
!   \item[status] return status: 0 - success, 1 - fail, optional output 
!  \end{description}
!  
!  The routines invoked are:
!  \begin{description}
!    \item{xmrg\_read\_header} (\ref{xmrg_read_header})
!    \item{hrap\_to\_latlon} (\ref{hrap_to_latlon})
!  \end{description}
!EOP

    implicit none 
    character(*), intent(in)       :: xmrg_name
    real, allocatable, intent(out) :: lat(:,:)
    real, allocatable, intent(out) :: lon(:,:)
    real, optional, intent(in)     :: resolution
    integer, optional, intent(out) :: status

    integer :: row, col, return_status, len_rec2
    real    :: rlon, rlat, hrap_resolution
    real    :: hrapx, hrapy, cell_size

  ! retrieve header information 
    call xmrg_read_header(xmrg_name,               &
                          cell_size=cell_size,     &
                          len_rec2=len_rec2,       &
                          status=return_status)
        
    if(return_status .eq. 1) then
       if(present(status)) status = return_status
       return
    else
       ! cell size is specified when rec2 length = 16
       if(len_rec2 .eq. 16) then
          if(present(resolution)) then
             if(abs(resolution-cell_size)<1E-6) then
                hrap_resolution = resolution
             else
                write(*, *) "xmrg_grid_latlon fatal error: HRAP resolution is different from XMRG/XMRG-like cell size"
                write(*, *) "HRAP resolution: ", resolution
                write(*, *) "XMRG/XMRG-like cell size: ", cell_size 
                if(present(status)) status = 1 
                return
             endif
          else
             hrap_resolution = cell_size
          endif
       else
          write(*,*) "xmrg_hrap_xy: default HRAP resolution (1.0) is used when reading ", trim(xmrg_name) 
          hrap_resolution = 1.0
       endif

       if(allocated(lat)) deallocate(lat)
       if(allocated(lon)) deallocate(lon)
       allocate(lat(maxy, maxx))
       allocate(lon(maxy, maxx))
       do row=1, maxy
          do col=1, maxx
             hrapx = xor + (col - 1) * hrap_resolution
             hrapy = yor + (row - 1) * hrap_resolution
             call hrap_to_latlon(hrapx, hrapy, rlon, rlat)
             lat(row, col) = rlat
             lon(row, col) = rlon
          end do
       end do
    endif

  end subroutine xmrg_grid_latlon

    
!BOP
! 
! !ROUTINE: latlon_to_hrap
! \label{latlon_to_hrap}
! !INTERFACE:
!
    subroutine latlon_to_hrap(rlon, rlat, hrap_x, hrap_y)
!
! !DESCRIPTION:
!   This routine convert latitude and longitude into HRAP coordinates 
!   (HRAPX and HRAPY). The code is adapted from
!   \url{http://www.nws.noaa.gov/oh/hrl/dmip/lat_lon.txt}
!EOP
    implicit none 
    real, intent(in)     :: rlon
    real, intent(in)     :: rlat
    real, intent(out)    :: hrap_x
    real, intent(out)    :: hrap_y
    real :: pi, d2rad, earthr, ref_lat, ref_lon
    real :: rmesh, tlat, re, flat, flon, r, x, y, rlon2

  ! west -> negative 
    rlon2 = -1.0 * rlon

    pi = 3.141592654
    d2rad = pi/180.
    earthr = 6371.2
    ref_lat = 60.
    ref_lon = 105.
    rmesh = 4.7625
    tlat = ref_lat*d2rad
    re = (earthr*(1.+sin(tlat)))/rmesh
    flat = rlat*d2rad
    flon = (rlon2+180.-ref_lon)*d2rad
    r = re*cos(flat)/(1.+sin(flat))
    x = r*sin(flon)
    y = r*cos(flon)
    hrap_x = x+401.
    hrap_y = y+1601.

  ! try to preserve a integer form, accurate to the 3rd decimal point
    hrap_x = 0.001*floor(1000*hrap_x+0.5)
    hrap_y = 0.001*floor(1000*hrap_y+0.5)

 end subroutine latlon_to_hrap
    
    
!BOP
! 
! !ROUTINE: hrap_to_latlon
! \label{hrap_to_latlon}
!
! !INTERFACE:
  subroutine hrap_to_latlon(hrap_x, hrap_y, rlon, rlat)
!
! !DESCRIPTION:
!   This routine converts HRAP coordinates into  latitude and longitude
!    The code is adapted from \url{http://www.nws.noaa.gov/oh/hrl/dmip/lat_lon.txt}
!EOP
    implicit none
    real, intent(in)    :: hrap_x
    real, intent(in)    :: hrap_y
    real, intent(out)   :: rlon
    real, intent(out)   :: rlat

    real :: earthr, stlon, pi, raddeg, xmesh 
    real :: tlat, x, y, rr, gi, ang
 
        earthr = 6371.2
        stlon = 105.
        pi = 3.141592654
        raddeg = 180./pi
        xmesh = 4.7625
        tlat = 60./raddeg
        x = hrap_x-401.
        y = hrap_y-1601.
        rr = x*x+y*y
        gi = ((earthr*(1.+sin(tlat)))/xmesh)
        gi = gi*gi
        rlat = asin((gi-rr)/(gi+rr))*raddeg
        ang = atan2(y,x)*raddeg
        if(ang.lt.0.) ang = ang+360.
        rlon = 270.+stlon-ang
        if(rlon.lt.0.) rlon = rlon+360.
        if(rlon.gt.360.) rlon = rlon-360.
       
        ! west -> negative 
        rlon = -1.0 * rlon

  end subroutine hrap_to_latlon


!BOP
! 
! !ROUTINE: LDT_transform_xmrgparam
! \label{LDT_transform_xmrgparam}
!
! !INTERFACE: 
  subroutine LDT_transform_xmrgparam( nest, target_nc, target_nr, &
                 target_gridDesc, xmrg_name, fill_value, final_data )

   implicit none

   integer,      intent(in)  :: nest
   integer,      intent(in)  :: target_nc
   integer,      intent(in)  :: target_nr
   real,         intent(in)  :: target_gridDesc(20)
   character(*), intent(in)  :: xmrg_name
   real,         intent(in)  :: fill_value
   real,        intent(out)  :: final_data(target_nc,target_nr)
!
! !DESCRIPTION:
!
!  This routine reads in an XMRG or XMRG-like file with HRAP domain and
!   transforms the read-in file to a new output target grid.
!
!  Description of arguments:
!  \begin{description}
!    \item[nest]        LIS/LDT domain nest index  
!    \item[target\_nc]  Target grid number of columns
!    \item[target\_nr]  Target grid number of rows
!    \item[target\_gridDesc]  Target grid description
!    \item[xmrg\_name]  File name of XMRG/XMRG-like file
!    \item[fill\_value] Value to fill missing/undefined grid points
!    \item[final\_data] Final data on the target grid 
!  \end{description}
! 
!  The routine invoked is:
!  \begin{description}
!    \item{xmrg\_read\_header\_c} C function reading XMRG/XMRG-like
!      header records (1st and 2nd records). (\ref{xmrg_read_header_c}
!  \end{description}
!EOP

   integer  :: rv
   integer  :: pi, pj, ti, tj
   integer  :: r, c
   integer  :: param_nc, param_nr
   real     :: target_xll, target_yll
   real     :: param_xll, param_yll
   real     :: xoffset, yoffset
   real     :: rlon, rlat
   real     :: cell_size, res
   real, allocatable :: data2d(:,:)
   real, allocatable :: data2d_transp(:,:)
   integer  :: file_status
   logical  :: file_exists

! ____________________________________________

   inquire(file=trim(xmrg_name), exist=file_exists)
   if( file_exists ) then
      final_data = fill_value

      call xmrg_read_header( xmrg_name, x_origin=param_xll, &
                             y_origin=param_yll, cell_size=cell_size, &
                             STATUS=file_status)

      call xmrg_read_data( xmrg_name, data2d, &
                           xmrg_nrow=param_nr, xmrg_ncol=param_nc, &
                           STATUS=file_status)

      rlat   = target_gridDesc(4)
      rlon   = target_gridDesc(5)

      call latlon_to_hrap(rlon, rlat, target_xll, target_yll)

      allocate( data2d_transp(param_nc, param_nr) )
      data2d_transp = transpose(data2d)

!      print *, "ll_x, ll_y, num_cols, num_rows"
!      print *, "Parameter grid: ", param_xll, param_yll, param_nc, param_nr
!      print *, "   Target grid: ", target_xll, target_yll, target_nc, target_nr
!      print *, "param c:  ", param_xll,  (param_xll+param_nc-1)
!      print *, "param r:  ", param_yll,  (param_yll+param_nr-1)
!      print *, "target c: ", target_xll, (target_xll+target_nc-1)
!      print *, "target r: ", target_yll, (target_yll+target_nr-1)

   !- Estimate x- and y-offsets for parameter index-values:
   ! Note!:  If resolution < 1.0, then need to figure out how to remove nint()
      if( nint(target_xll) > nint(param_xll) ) then
!        xoffset = (target_xll-param_xll)
        xoffset = param_xll
      elseif( nint(target_xll) <= nint(param_xll) ) then
        xoffset = param_xll
      endif
!      print *, "xoffset: ",xoffset

      if( nint(target_yll) > nint(param_yll) ) then
!        yoffset = (target_yll-param_yll)
        yoffset = param_yll
      elseif( nint(target_yll) <= nint(param_yll) ) then
        yoffset = param_yll
      endif
!      print *, "yoffset:", yoffset

      tj = 0
      do r = nint(target_yll), nint(target_yll+target_nr-1)
         tj = tj + 1
         pj = ( r - yoffset ) + 1
         ti = 0; pi = 0
         do c = nint(target_xll), nint(target_xll+target_nc-1)
            ti = ti + 1
            pi = ( c - xoffset ) + 1
            if( pj <= 0 .or. pi <= 0 ) then
               final_data(ti,tj) = -1.
!               final_data(ti,tj) = fill_value
            elseif( r > nint(param_yll+param_nr-1) .or. &
                    c > nint(param_xll+param_nc-1) ) then
               final_data(ti,tj) = -1.
!               final_data(ti,tj) = fill_value
            else 
               final_data(ti,tj) = data2d_transp(pi,pj)
            endif
         end do
      end do

      deallocate( data2d_transp )

   end if

 end subroutine LDT_transform_xmrgparam

end module LDT_xmrg_reader
