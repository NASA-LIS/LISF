!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
#include "LDT_NetCDF_inc.h"
!BOP
! !ROUTINE: read_Drange
! \label{read_Drange}

! !REVISION HISTORY: 
! 2Dec2021: Mahdi Navari ; Initial Specification
!
! !INTERFACE: 
subroutine read_Drange(ngrid, filename, varname, xrange)

  use LDT_coreMod
  use LDT_logMod
  use LDT_historyMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

     implicit none
! !ARGUMENTS:      
!     integer,   intent(in)    :: n
!     integer,   intent(in)    :: k
!     integer,   intent(in)    :: nbins
!     integer,   intent(in)    :: ntimes
     integer,   intent(in)    :: ngrid
     character(len=*)         :: filename
     character(len=*)         :: varname
     real                     :: xrange(ngrid,2)
     !real        :: cdf(ngrid,ntimes, nbins)

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
     integer                  :: j,kk
     integer                  :: ngridId, nbinsId, nlevsId,ntimesId
     integer                  :: ngrid_file, nbins_file, nlevs_file, ntimes_file
     integer                  :: xid, cdfid
     real, allocatable        :: xrange_file(:,:,:)
     !real, allocatable        :: cdf_file(:,:,:)
     integer                  :: nid

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
     write(LDT_logunit,*) '[INFO] Reading Drange form CDF file ',trim(filename)
     !if(ngrid.gt.0) then
        call LDT_verify(nf90_open(path=trim(filename),mode=NF90_NOWRITE,&
             ncid=nid),'failed to open file '//trim(filename))

        call LDT_verify(nf90_inq_dimid(nid,"ngrid",ngridId), &
             'nf90_inq_dimid failed for ngrid')
        call LDT_verify(nf90_inq_dimid(nid,"nbins",nbinsId), &
             'nf90_inq_dimid failed for nbins')
        call LDT_verify(nf90_inq_dimid(nid,trim(varname)//"_levels",nlevsId), &
             'nf90_inq_dimid failed for '//trim(varname)//"_levels")
        call LDT_verify(nf90_inq_dimid(nid,"ntimes",ntimesId), &
             'nf90_inq_dimid failed for ntimes')


        call LDT_verify(nf90_inquire_dimension(nid,ngridId, len=ngrid_file),&
             'nf90_inquire_dimension failed for ngrid')
        call LDT_verify(nf90_inquire_dimension(nid,nbinsId, len=nbins_file),&
             'nf90_inquire_dimension failed for nbins')
        call LDT_verify(nf90_inquire_dimension(nid,nlevsId, len=nlevs_file),&
             'nf90_inquire_dimension failed for nbins')
        call LDT_verify(nf90_inquire_dimension(nid,ntimesId, len=ntimes_file),&
             'nf90_inquire_dimension failed for ntimes')

        !if(nbins.ne.nbins_file) then
        !   write(LDT_logunit,*) '[ERR] The number of bins specified in the file '//&
        !        trim(filename)
        !   write(LDT_logunit,*) '[ERR] (',nbins_file, &
        !        ') is different from the number of bins specified'
        !   write(LDT_logunit,*) '[ERR] in the lis.config file (',nbins,')'
        !   call LDT_endrun()
        !endif
        
        if (ntimes_file .gt. 1) then 
            write(LDT_logunit,*) '[ERR] The number of times specified in the file '//&
                trim(filename)
           write(LDT_logunit,*) '[ERR] (',ntimes_file, &
                ') should be 1 set the Temporal resolution of precipitation CDFs to "yearly" '
           call LDT_endrun()
        endif

        !allocate(xrange(ngrid_file,2))
        allocate(xrange_file(ngrid_file,nlevs_file, nbins_file))
        !allocate(cdf_file(ngrid_file,nlevs_file, nbins))

        do j=1,ntimes_file
           call LDT_verify(nf90_inq_varid(nid,trim(varname)//'_xrange',xid),&
                'nf90_inq_varid failed for for '//trim(varname)//'_xrange')
           !call LDT_verify(nf90_inq_varid(nid,trim(varname)//'_CDF',cdfid),&
           !     'nf90_inq_varid failed for '//trim(varname)//'_CDF')

           call LDT_verify(nf90_get_var(nid,xid,xrange_file, &
                start=(/1,j,1,1/), count=(/ngrid_file,1,nlevs_file,nbins_file/)),&
                'nf90_get_var failed for '//trim(varname)//'_xrange')
           !call LDT_verify(nf90_get_var(nid,cdfid,cdf_file,&
           !     start=(/1,j,1,1/), count=(/ngrid_file,1,nlevs_file,nbins/)),&
           !     'nf90_get_var failed for '//trim(varname)//'_CDF')

! Note: CDF is generated by LDT --> therefore data is in the LDT domain
#if 0 
           if(ngrid.gt.0) then
              do kk=1,nbins
                 call LIS_convertObsVarToLocalSpace(n,k,xrange_file(:,1,kk), &
                      xrange(:,j,kk))
                 call LIS_convertObsVarToLocalSpace(n,k,cdf_file(:,1,kk), &
                      cdf(:,j,kk))
              enddo
           endif
#endif
        enddo

        xrange(:,1) = xrange_file(:,1,1) ! min value for each gridcell
        xrange(:,2) = xrange_file(:,1,nbins_file) ! max value for each gridcell
!        xrange = xrange_file
!        cdf = cdf_file

        deallocate(xrange_file)
        !deallocate(cdf_file)

        call LDT_verify(nf90_close(nid),&
             'failed to close file '//trim(filename))
        write(LDT_logunit,*) '[INFO] Successfully read CDF file ',trim(filename)
     !endif
#endif

end subroutine read_Drange
