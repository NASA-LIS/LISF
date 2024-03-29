!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
#include "LDT_NetCDF_inc.h"
!BOP
! !ROUTINE: read_Precip_climo
! \label{read_Precip_climo}

! !REVISION HISTORY:
! 19Jan2022: Mahdi Navari ; Initial Specification
!
! !INTERFACE:
subroutine read_Precip_climo(ncol, nrow, filename, precip)

  use LDT_coreMod
  use LDT_logMod
  use LDT_historyMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
  ! !ARGUMENTS:
  integer,   intent(in)    :: ncol, nrow
  character(len=*), intent(in)  :: filename
  real                     :: precip(ncol,nrow,12)
  logical                  :: file_exists
!
! !DESCRIPTION:
!  This routine reads the input CDF file (generated by LDT in NETCDF format)
!  The xrange values and the corresponding CDFs are read for each grid point.
!  Both these fields are expected to be in the 1-d grid vector dimension.
!
!  The arguments are:
!  \begin{description}
!  \item[n]             index of the nest
!  \item[nbins]         number of bins used to compute the model and obs CDFs
!  \item[filename]      name of the CDF file
!  \item[varname]       name of the variable being extracted.
!  \item[xrange]        x-axis values corresponding to the CDF
!  \item[cdf]           y-axis (CDF) values corresponding to the CDF
! \end{description}
!EOP
  real, allocatable        :: precip_climo(:,:,:)
  integer                  :: ios,nid,ncId,nrId,varId
  integer                  :: nc,nr,i,c,r
  character*3              :: month_name(12)

  month_name = (/"JAN","FEB","MAR","APR","MAY","JUN",&
       "JUL","AUG","SEP","OCT","NOV","DEC"/)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  inquire(file=trim(filename), exist=file_exists)

  if (.not. file_exists) then
     write(LDT_logunit,*) '[ERR] (',filename, &
          ') does not exist. Please provide a monthly precipitation climatology in CDF format '
     call LDT_endrun()
  end if


  write(LDT_logunit,*) '[INFO] Reading Precipitation climatology from CDF file ',trim(filename)
  call LDT_verify(nf90_open(path=trim(filename),mode=NF90_NOWRITE,&
       ncid=nid),'failed to open file '//trim(filename))

  ios = nf90_inq_dimid(nid,"east_west",ncId)
  call LDT_verify(ios,'Error in nf90_inq_dimid in read_precip_climo')

  ios = nf90_inq_dimid(nid,"north_south",nrId)
  call LDT_verify(ios,'Error in nf90_inq_dimid in read_precip_climo')

  ios = nf90_inquire_dimension(nid,ncId, len=nc)
  call LDT_verify(ios,'Error in nf90_inquire_dimension in read_precip_climo')

  ios = nf90_inquire_dimension(nid,nrId, len=nr)
  call LDT_verify(ios,'Error in nf90_inquire_dimension in read_precip_climo')

  if (ncol .ne. nc .or. nrow .ne.nr)  then
     write(LDT_logunit,*) &
          '[ERR] The number of columns or rows specified in the file '//&
          trim(filename)
     write(LDT_logunit,*) '[ERR] (', nc, nr,  &
          ') is different from the number of columns and rows specified'
     write(LDT_logunit,*) '[ERR] in the ldt.config file (', ncol, nrow, ')'
     call LDT_endrun()
  endif

  allocate(precip_climo(nc,nr,12))

  do i=1,12
     ios = nf90_inq_varid(nid,'TotalPrecip_'//trim(month_name(i)),varId)
     call LDT_verify(ios,'Precipitation climo field not found in the file')

     ios = nf90_get_var(nid,varId,precip_climo(:,:,i))
     call LDT_verify(ios,'Error in nf90_get_var in LDT_procDataForREAL')
  enddo

  ios = nf90_close(nid)
  call LDT_verify(ios,'Error in nf90_close in readldtparam_real_2d')
  write(LDT_logunit,*) '[INFO] Successfully read Precipitation climo file ',trim(filename)

  precip = precip_climo
  deallocate(precip_climo)

! end USE_NETCDF4
#endif

end subroutine read_Precip_climo

