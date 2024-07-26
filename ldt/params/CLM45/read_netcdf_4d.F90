!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! modified version of Yudong's code translate_to_LDT.F90
! read from CLM4.5 native parameter file 
! Usage:
! read_netcdf_4d nest inputfile varname #dim dim1_name dim2_name dim3_name dim4_name
! Three data_type: int, float, double 
! fixed to four dimensions 


#include "LDT_misc.h"
subroutine read_netcdf_4d (n,infile,varname,ndims,dname1,dname2,dname3,dname4,nx,ny,nv,nt,subset_param) 

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LDT_coreMod,  only : LDT_rc
  use LDT_gridmappingMod
  use CLM45_parmsMod

  implicit none
!  include "netcdf.inc"

  integer, intent(in) :: n
  character (len=*),intent(in) :: infile
  character (len=*),intent(in) :: varname
  character (len=*),intent(in) :: dname1,dname2,dname3,dname4
  integer, intent(in) :: ndims
  integer, intent(in) :: nx,ny,nv,nt
  real*8, dimension(nx,ny,nv,nt),intent(inout) :: subset_param
  character (len=8),dimension(4) ::  dim_name
  integer :: start(4), count(4), dimid(4) ! max 4-dimensional  
  real*8, allocatable  :: ddata4d(:, :, :, : )

  integer :: ncid, varid, iargc, i, j, k, l, len, ic 
  integer :: ncid2, varid2

  integer :: r,c
  integer :: glpnc, glpnr           ! Parameter (global) total columns and rows
  integer :: subpnc, subpnr         ! Parameter subsetted columns and rows
  real    :: subparam_gridDesc(20)  ! Input parameter grid desc array
  integer, allocatable  :: lat_line(:,:), lon_line(:,:)
  real    :: shifted_gridDesc(20)  ! target parameter grid desc array

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  call check( nf90_open(infile, NF90_NOWRITE, ncid) )

  !Note: if variable in ncdump output is MONTHLY_SAI(time, lsmpft, lsmlat, lsmlon)
  !  the fortran data array should be data(lsmlon, lsmlat, lsmpft, time),
  !  and so should count array, ie., count=(/lsmlon, lsmlat, lsmpft, time/). 
  ! get all the dimension sizes 
  start(:) = 1
  count(:) = 1
  len = 1 
  dim_name(1) = dname1
  dim_name(2) = dname2
  dim_name(3) = dname3
  dim_name(4) = dname4

  do i=ndims, 1, -1  ! swap the dim order for fortran for count() array
    ic = ndims - i + 1  ! reverse order
    call check( nf90_inq_dimid(ncid, trim(dim_name(i)), dimid(i) ) ) 
    call check( nf90_inquire_dimension(ncid, dimid(i), len=count(ic) ) ) 
!    write(*, *) trim(dim_name(i)), " dimid= ", dimid(i), "  dim size=", count(ic) 
    len = len * count(ic) 
  end do 
!  write(*, *) "count1-4",(count(i), i=1, ndims) 
!  write(*, *) "nx ny nv nt",nx,ny,nv,nt 

  call check( nf90_inq_varid(ncid, trim(varname), varid) )

! forget data type matching. Convert everything to double when reading whatever types
! see https://www.unidata.ucar.edu/software/netcdf/docs/netcdf/Type-Conversion.html#Type-Conversion

  allocate(ddata4d(count(1), count(2), count(3), count(4)))
  call check( nf90_get_var(ncid, varid, ddata4d, start=start, count=count ) )
  ! Close the file, freeing all resources.
  call check( nf90_close(ncid) )

  ! make sure array size fits into odata4d structure
  if ( count(1) .eq. CLM45_struc(n)%clm45parms_gridDesc(2) ) then

  ! Shift longitude from 0-360 to -180~+180. Assume Global input
    if ( CLM45_struc(n)%clm45parms_gridDesc(5) .gt. 0.0 .or.&
          CLM45_struc(n)%clm45parms_gridDesc(8) .gt. 180.0 ) then
     call shift_lon( count(1), count(2), count(3), count(4),  &
                     CLM45_struc(n)%clm45parms_gridDesc, shifted_gridDesc, ddata4d )
    endif

! -------------------------------------------------------------------
!    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
!  - Map Parameter Grid Info to LIS Target Grid/Projection Info --
    subparam_gridDesc = 0.
    call LDT_RunDomainPts( n, CLM45_struc(n)%clm45parms_proj, shifted_gridDesc(:), &
             glpnc, glpnr, subpnc, subpnr, subparam_gridDesc, lat_line, lon_line )
!    write(*, *) "subpnc, subpnr:",subpnc, subpnr
!    write(*, *) "lon_line:",lon_line(1,1),lon_line(subpnc,subpnr)
!    write(*, *) "lat_line:",lat_line(1,1),lat_line(subpnc,subpnr)
!    write(*, *) "subparam_gridDesc:",subparam_gridDesc

!- Subset parameter read-in array:
!  array size: subset_param(subpnc, subpnr, count(3), count(4)) 
    do r = 1, subpnr
       do c = 1, subpnc
          subset_param(c,r,:,:) = ddata4d(lon_line(c,r),lat_line(c, r),:,:)
       enddo
    enddo

    deallocate(lon_line)
    deallocate(lat_line)
! -------------------------------------------------------------------
  else
! 1D data
!   write(*,*) "1D data: ",trim(varname)
   subset_param = ddata4d
!   write(*,*) "[ERR]: need to restructure odata4d ",trim(varname)
!   stop
  endif
      

! save to output file 
!     call check( nf90_create("outfile.nc", NF90_CLOBBER, ncid2) )
!     call check( nf90_inq_varid(ncid2, trim(varname), varid2) )
!     call check( nf90_put_var(ncid2, varid2, odata4d ) )
!     call check( nf90_close(ncid2) )

!  print*,trim(varname),maxval(ddata4d),minval(ddata4d)
!  print*,'ddata4d:',ddata4d(1,1,1,1)
!  print*,trim(varname),maxval(subset_param),minval(subset_param)
!  print*,'subset_param:',subset_param(1,1,1,1)
  deallocate(ddata4d)  

#endif
end subroutine read_netcdf_4d 

subroutine read_netcdf_4d_xvyv (n,infile,varname,ndims,dname1,dname2,dname3,dname4,nx,ny,nv,nt,subset_param) 

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LDT_coreMod,  only : LDT_rc
  use LDT_gridmappingMod
  use CLM45_parmsMod

  implicit none
!  include "netcdf.inc"

  integer, intent(in) :: n
  character (len=*),intent(in) :: infile
  character (len=*),intent(in) :: varname
  character (len=*),intent(in) :: dname1,dname2,dname3,dname4
  integer, intent(in) :: ndims
  integer, intent(in) :: nx,ny,nv,nt
  real*8, dimension(nx,ny,nv,nt),intent(inout) :: subset_param
  character (len=8),dimension(4) ::  dim_name
  integer :: start(4), count(4), dimid(4) ! max 4-dimensional  
  real*8, allocatable  :: ddata4d(:, :, :, : ),odata4d(:, :, :, : )

  integer :: ncid, varid, iargc, i, j, k, l, len, ic 
  integer :: ncid2, varid2

  integer :: r,c
  integer :: glpnc, glpnr           ! Parameter (global) total columns and rows
  integer :: subpnc, subpnr         ! Parameter subsetted columns and rows
  real    :: subparam_gridDesc(20)  ! Input parameter grid desc array
  integer, allocatable  :: lat_line(:,:), lon_line(:,:)
  real    :: shifted_gridDesc(20)  ! target parameter grid desc array

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  call check( nf90_open(infile, NF90_NOWRITE, ncid) )

  !Note: if variable in ncdump output is MONTHLY_SAI(time, lsmpft, lsmlat, lsmlon)
  !  the fortran data array should be data(lsmlon, lsmlat, lsmpft, time),
  !  and so should count array, ie., count=(/lsmlon, lsmlat, lsmpft, time/). 
  ! get all the dimension sizes 
  start(:) = 1
  count(:) = 1
  len = 1 
  dim_name(1) = dname1   !nj
  dim_name(2) = dname2   !ni
  dim_name(3) = dname3   !nv
  dim_name(4) = dname4   !nt

  call check( nf90_inq_dimid(ncid, trim(dim_name(3)), dimid(3) ) ) 
  call check( nf90_inquire_dimension(ncid, dimid(3), len=count(1) ) ) 
  len = len * count(1) 
!  write(*, *) trim(dim_name(3)), " dimid= ", dimid(3), "  dim size=", count(1) 
  call check( nf90_inq_dimid(ncid, trim(dim_name(2)), dimid(2) ) ) 
  call check( nf90_inquire_dimension(ncid, dimid(2), len=count(2) ) ) 
  len = len * count(2) 
!  write(*, *) trim(dim_name(2)), " dimid= ", dimid(2), "  dim size=", count(2) 
  call check( nf90_inq_dimid(ncid, trim(dim_name(1)), dimid(1) ) ) 
  call check( nf90_inquire_dimension(ncid, dimid(1), len=count(3) ) ) 
  len = len * count(3) 
!  write(*, *) trim(dim_name(1)), " dimid= ", dimid(1), "  dim size=", count(3) 
!  write(*, *) "xvyv: count1-4",count
!  write(*, *) "xvyv: nx ny nv nt",nx,ny,nv,nt 

  call check( nf90_inq_varid(ncid, trim(varname), varid) )

! forget data type matching. Convert everything to double when reading whatever types
! see https://www.unidata.ucar.edu/software/netcdf/docs/netcdf/Type-Conversion.html#Type-Conversion

  allocate(ddata4d(count(1), count(2), count(3), count(4)))
  allocate(odata4d(count(2),count(3),count(1),count(4))) ! global LonxLatxLevxT
  call check( nf90_get_var(ncid, varid, ddata4d, start=start, count=count ) )
  ! Close the file, freeing all resources.
  call check( nf90_close(ncid) )

  ! make sure array size fits into odata4d structure
  if ( count(1) .eq. CLM45_struc(n)%clm45parms_gridDesc(2) ) then
   odata4d = ddata4d
  else
!   write(*,*) "[WARNING]: xvyv need to restructure odata4d ",trim(varname)
   do j = 1, count(3)
    do i = 1, count(2)
      do k = 1, count(1)
         odata4d(i,j,k,1) = ddata4d(k,i,j,1)
      end do
    end do
   end do
  endif 
  ! Shift longitude from 0-360 to -180~+180. Assume Global input
  if ( CLM45_struc(n)%clm45parms_gridDesc(5) .gt. 0.0 .or.&
        CLM45_struc(n)%clm45parms_gridDesc(8) .gt. 180.0 ) then
   call shift_lon( count(2), count(3), count(1), count(4),  &
                   CLM45_struc(n)%clm45parms_gridDesc, shifted_gridDesc ,odata4d )
  endif
      
! -------------------------------------------------------------------
!    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
!  - Map Parameter Grid Info to LIS Target Grid/Projection Info --
    subparam_gridDesc = 0.
    call LDT_RunDomainPts( n, CLM45_struc(n)%clm45parms_proj, shifted_gridDesc(:), &
             glpnc, glpnr, subpnc, subpnr, subparam_gridDesc, lat_line, lon_line )
!    write(*, *) "subpnc, subpnr:",subpnc, subpnr
!    write(*, *) "lon_line:",lon_line(1,1),lon_line(subpnc,subpnr)
!    write(*, *) "lat_line:",lat_line(1,1),lat_line(subpnc,subpnr)
!    write(*, *) "subparam_gridDesc:",subparam_gridDesc

!- Subset parameter read-in array:
!  array size: subset_param(subpnc, subpnr, count(3), count(4)) 
    do r = 1, subpnr
       do c = 1, subpnc
          subset_param(c,r,:,:) = odata4d(lon_line(c,r),lat_line(c, r),:,:)
       enddo
    enddo

    deallocate(lon_line)
    deallocate(lat_line)

  deallocate(ddata4d)  
  deallocate(odata4d)  

#endif
end subroutine read_netcdf_4d_xvyv

subroutine read_netcdf_4d_global (n,infile,varname,ndims,dname1,dname2,dname3,dname4,nx,ny,nv,nt,odata4d,shifted_gridDesc)
! This routine read the CLM-4.5 NetCDF input file at global scale

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LDT_coreMod,  only : LDT_rc

  implicit none
!  include "netcdf.inc"

  integer, intent(in) :: n
  character (len=*),intent(in) :: infile
  character (len=*),intent(in) :: varname
  character (len=*),intent(in) :: dname1,dname2,dname3,dname4
  integer, intent(in) :: ndims
  integer, intent(in) :: nx,ny,nv,nt
  real*8, dimension(nx,ny,nv,nt),intent(inout) :: odata4d
  real, intent(out)     :: shifted_gridDesc(20) 

  character (len=8),dimension(4) ::  dim_name
  integer :: start(4), count(4), dimid(4) ! max 4-dimensional
  real*8, allocatable  :: ddata4d(:, :, :, : )

  integer :: ncid, varid, iargc, i, j, k, l, len, ic
  integer :: ncid2, varid2

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  call check( nf90_open(infile, NF90_NOWRITE, ncid) )

  !Note: if variable in ncdump output is MONTHLY_SAI(time, lsmpft, lsmlat, lsmlon)
  !  the fortran data array should be data(lsmlon, lsmlat, lsmpft, time),
  !  and so should count array, ie., count=(/lsmlon, lsmlat, lsmpft, time/).
  ! get all the dimension sizes
  start(:) = 1
  count(:) = 1
  len = 1
  dim_name(1) = dname1
  dim_name(2) = dname2
  dim_name(3) = dname3
  dim_name(4) = dname4

  do i=ndims, 1, -1  ! swap the dim order for fortran for count() array
    ic = ndims - i + 1  ! reverse order
    call check( nf90_inq_dimid(ncid, trim(dim_name(i)), dimid(i) ) )
    call check( nf90_inquire_dimension(ncid, dimid(i), len=count(ic) ) )
!    write(*, *) trim(dim_name(i)), " dimid= ", dimid(i), "  dim size=", count(ic)
    len = len * count(ic)
  end do
!  write(*, *) "count1-4",(count(i), i=1, ndims)
!  write(*, *) "nx ny nv nt",nx,ny,nv,nt

  call check( nf90_inq_varid(ncid, trim(varname), varid) )

! forget data type matching. Convert everything to double when reading whatever types
! see https://www.unidata.ucar.edu/software/netcdf/docs/netcdf/Type-Conversion.html#Type-Conversion

  allocate(ddata4d(count(1), count(2), count(3), count(4)))
  call check( nf90_get_var(ncid, varid, ddata4d, start=start, count=count ) )
  ! Close the file, freeing all resources.
  call check( nf90_close(ncid) )

  ! make sure array size fits into odata4d structure
  if ( count(1) .eq. LDT_rc%mask_gridDesc(n,2) ) then
  ! Shift longitude from 0-360 to -180~+180. Assume Global input
    if ( LDT_rc%mask_gridDesc(n,5) .gt. 0.0 .or.&
          LDT_rc%mask_gridDesc(n,8) .gt. 180.0 ) then
     call shift_lon( count(1), count(2), count(3), count(4),  &
                     LDT_rc%mask_gridDesc(n,:), shifted_gridDesc, ddata4d )
    endif
   odata4d = ddata4d
  else
   write(*,*) "[ERR]: CLM-4.5 array size mismatch ",trim(varname)
   stop
  endif

  deallocate(ddata4d)

#endif
end subroutine read_netcdf_4d_global

  subroutine check(status)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
    integer, intent ( in) :: status

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
#endif
  end subroutine check  

  subroutine shift_lon( nx, ny, nv, nt, clm45parms_gridDesc, shifted_gridDesc, sdata )
! this routine shifts longitude ranging 0-360 to -180~ +180.
  implicit none
  integer, intent (in) :: nx, ny, nv, nt
  real, intent(in)     :: clm45parms_gridDesc(20) 
  real, intent(out)     :: shifted_gridDesc(20) 
  real*8, intent(inout):: sdata(nx,ny,nv,nt)  
  integer              :: c, r, shiftx, t, v
  real                 :: lon1, lon, dx
  real*8, dimension(nx):: temp, shifted
  logical              :: found

  found = .false.
  ! copy grid coord info
  shifted_gridDesc = clm45parms_gridDesc
  dx = clm45parms_gridDesc(9)
  lon1 = clm45parms_gridDesc(5) - dx/2.0   ! corner of native (0-360)
  do c = 1, nx
   lon = lon1 + (c-1)*dx
   if ( lon .ge. 180.0 ) then
    shiftx = c - 1   
    shifted_gridDesc(5) = (lon + dx/2.0 ) - 360.     ! -179.375
    shifted_gridDesc(8) = lon1 + dx/2.0 + (nx-1)*dx  !  179.375
    if ( shifted_gridDesc(8) .gt. 180.0 ) then
     shifted_gridDesc(8) = shifted_gridDesc(8) - 180.0
    endif
    found = .true.
    exit
   endif
   if ( found ) exit
  end do
!  print*,'shift lon by ',shiftx,' new lon range: ',shifted_gridDesc(5),shifted_gridDesc(8)
  do t = 1, nt
   do v = 1, nv
    do r = 1, ny
       temp = sdata(1:nx,r,v,t)
       shifted = cshift(temp,shift=-shiftx)
       sdata(1:nx,r,v,t) = shifted
    enddo
   enddo
  enddo

  return

  end subroutine shift_lon

