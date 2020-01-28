!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_ASO_SWE
! \label{read_ASO_SWE}
!
! !REVISION HISTORY:
!  01 Jul 2010: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_ASO_SWE(n, k, OBS_State, OBS_Pert_State)
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_dataAssimMod
  use LIS_DAobservationsMod
  use map_utils
  use UTM_utils
  use ASO_SWE_Mod, only : ASO_SWE_struc
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the Level 3 AMSR-E snow water equivalent observations 
!  from the HDF-EOS files and packages into an ESMF state object. 
!  The routine reads the snow data at 0z, performs spatial 
!  interpolation to the LIS grid and keeps it in memory. 
!  At 10.30 AM localtime for each grid point, the code then 
!  packages the interpolated observations into an ESMF state object. 
!  This routine handles the AMSR-E SWE retrievals available in the 
!  HDF-EOS format. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  type(ESMF_Field)    :: sweField

  integer             :: iret
  integer             :: i,j
  integer             :: ftn,nx,ny,dim1id,dim2id
  integer             :: xid,yid,sweid
  integer             :: stn_col, stn_row
  real                :: col, row
  real,  allocatable  :: var(:,:)
  real,  allocatable  :: xval(:)
  real,  allocatable  :: yval(:)
  real                :: swe_ip(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  integer             :: nswe_ip(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real,    pointer    :: obsl(:)
  integer             :: gid(LIS_rc%obs_ngrid(k))
  integer             :: assimflag(LIS_rc%obs_ngrid(k))

  character*100       :: sweobsdir
  logical             :: data_update
  logical             :: file_exists
  character*200       :: name
  real                :: latdeg, londeg
  logical             :: alarmCheck

  logical             :: readflag
  integer             :: status
  logical             :: dataflag(LIS_npes)
  logical             :: dataflag_local
  logical             :: data_upd
  integer             :: t,c,r,p
  real                :: lon, lhour
  real                :: timenow
  integer             :: zone
  integer             :: ierr


  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       sweobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  data_upd = .false. 

  swe_ip = LIS_rc%udef
!-------------------------------------------------------------------------
!   Read the data at 0z
!-------------------------------------------------------------------------

  timenow = float(LIS_rc%hr)*3600 + 60.0*LIS_rc%mn + LIS_rc%ss
  alarmCheck = (mod(timenow, 86400.0).eq.0)

  if(alarmCheck) then 

     call ASO_SWE_filename(name,sweobsdir,&
          LIS_rc%yr,LIS_rc%mo,LIS_rc%da)   
     
     inquire(file=name,exist=file_exists)
     if(file_exists) then 
        readflag = .true. 
     else 
        readflag = .false.
     endif
     
     if (readflag) then 
        write(LIS_logunit,*) 'Reading NASA ASO file ',name

#if(defined USE_NETCDF3 || defined USE_NETCDF4)        
        call LIS_verify(nf90_open(path=name,mode=NF90_NOWRITE, ncid=ftn),&
             'Error opening file '//trim(name))
        call LIS_verify(nf90_inq_dimid(ftn,'x',dim1Id),&
             'Error in nf90_inq_dimid: x')
        call LIS_verify(nf90_inquire_dimension(ftn,dim1Id,len=nx),&
             'Error in nf90_inquire_dimension: x')

        call LIS_verify(nf90_inq_dimid(ftn,'y',dim2Id),&
             'Error in nf90_inq_dimid: y')
        call LIS_verify(nf90_inquire_dimension(ftn,dim2Id,len=ny),&
             'Error in nf90_inquire_dimension: y')
        allocate(var(nx,ny))
        allocate(xval(nx))
        allocate(yval(ny))

        call LIS_verify(nf90_inq_varid(ftn,'x',xid),&
             'Error in nf90_inq_varid: x')
        call LIS_verify(nf90_get_var(ftn,xid,xval),&
             'Error in nf90_get_var: x')
        call LIS_verify(nf90_inq_varid(ftn,'y',yid),&
             'Error in nf90_inq_varid: y')
        call LIS_verify(nf90_get_var(ftn,yid,yval),&
             'Error in nf90_get_var: y')
        call LIS_verify(nf90_inq_varid(ftn,'Band1',sweid),&
             'Error in nf90_inq_varid: Band1')
        call LIS_verify(nf90_get_var(ftn,sweid,var),&
             'Error in nf90_get_var: Band1')
        zone = 11
        swe_ip = 0
        nswe_ip = 0 
        
        do i=1,nx
           do j=1,ny
              call UTM2Geo(zone, yval(j), xval(i), latdeg, londeg)
              
              call latlon_to_ij(LIS_domain(n)%lisproj, latdeg, londeg,&
                   col, row)
              stn_col = nint(col)
              stn_row = nint(row)
              
              if(.not.isNaN(var(i,j)).and.var(i,j).ge.0.and.&
                   stn_col.gt.0.and.stn_col.le.LIS_rc%obs_lnc(k).and.&
                   stn_row.gt.0.and.stn_row.le.LIS_rc%obs_lnr(k)) then 
                 swe_ip(stn_col,stn_row) = swe_ip(stn_col,stn_row) + & 
                      var(i,j)*1000.0 !to mm
                 nswe_ip(stn_col,stn_row) = nswe_ip(stn_col,stn_row) + 1
              endif
           enddo
        enddo

        call LIS_verify(nf90_close(ftn), 'Error nf90_close in read_ASO_SWE')
        
        deallocate(var)
        deallocate(xval)
        deallocate(yval)
        
        do r=1,LIS_rc%obs_lnr(k)
           do c=1,LIS_rc%obs_lnc(k)
              if(nswe_ip(c,r).ne.0) then
                 swe_ip(c,r) = swe_ip(c,r)/nswe_ip(c,r)
! Because of the boundary of the observed image, we screen out 
! zero snow values
                 if(swe_ip(c,r).eq.0) then
                    swe_ip(c,r) = LIS_rc%udef
                 endif

              else
                 swe_ip(c,r) = LIS_rc%udef
              endif
           enddo
        enddo

#endif

     endif
  endif

  call ESMF_StateGet(OBS_State,"Observation01",sweField,&
       rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status)
  
  obsl = LIS_rc%udef 
  dataflag_local = .false. 
  
  do r =1,LIS_rc%obs_lnr(k)
     do c =1,LIS_rc%obs_lnc(k)
        if (LIS_obs_domain(n,k)%gindex(c,r) .ne. -1)then
           
           if(LIS_rc%hr.eq.0) then   !assimilate at 0z         
              obsl(LIS_obs_domain(n,k)%gindex(c,r))=&
                   swe_ip(c,r)
              if(obsl(LIS_obs_domain(n,k)%gindex(c,r)).ne.-9999.0) then 
                 dataflag_local = .true.
              endif
           endif
           
        end if
     end do
  end do

#if (defined SPMD)
  call MPI_ALLGATHER(dataflag_local,1, MPI_LOGICAL, dataflag(:),&
       1, MPI_LOGICAL, LIS_mpi_comm, ierr)
#endif
  data_upd = .false.
  
  do p=1,LIS_npes
     data_upd = data_upd.or.dataflag(p)
  enddo
  
  if(data_upd) then 
     do t=1,LIS_rc%obs_ngrid(k)
        gid(t) = t
        if(obsl(t).ne.-9999.0) then 
           assimflag(t) = 1
        else
           assimflag(t) = 0
        endif
     enddo
     
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .true., rc=status)
     call LIS_verify(status, 'Error: AttributeSet in Data Update Status')
     
     call ESMF_AttributeSet(sweField,"Grid Number",&
          gid,itemCount=LIS_rc%ngrid(n),rc=status)
     call LIS_verify(status,'Error: AttributeSet in Grid Number')
     
     call ESMF_AttributeSet(sweField,"Assimilation Flag",&
          assimflag,itemCount=LIS_rc%ngrid(n),rc=status)
     call LIS_verify(status, 'Error: AttributeSet in Assimilation Flag')
           
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status, "Error: AttributeSet Data Update Status")
     return
  end if
  
end subroutine read_ASO_SWE



!BOP
! !ROUTINE: ASO_SWE_filename
! \label{ASO_SWE_filename}
! 
! !INTERFACE: 
subroutine ASO_SWE_filename(name, ndir, yr, mo,da)
! !USES:   
  use LIS_coreMod,only : LIS_rc

  implicit none
! !ARGUMENTS: 
  character*200      :: name
  integer           :: yr, mo, da, hr,mn
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the Level 3 ASO SWE filename based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[name] name of the ASO SWE filename
!  \item[ndir] name of the ASO SWE root directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') LIS_rc%yr
  write(unit=fmo, fmt='(i2.2)') LIS_rc%mo
  write(unit=fda, fmt='(i2.2)') LIS_rc%da
  
  name = trim(ndir)//'/ASO_swe_'//&
       trim(fyr)//trim(fmo)//trim(fda)//'.nc'

end subroutine ASO_SWE_filename




