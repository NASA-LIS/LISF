!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_IMSsca
! \label{read_IMSsca}
!
! !REVISION HISTORY:
!  1 Jun 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_IMSsca(n, OBS_State,OBS_Pert_State) 
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_pluginIndices, only : LIS_IMSscaobsId
  use LIS_constantsMod,  only : LIS_CONST_PATH_LEN
  use IMSsca_Mod, only : IMSsca_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads ANSA snow observations
!  and packages it into an ESMF State with certain predefined 
!  attributes. The routine reads the snow data at 0z, performs 
!  spatial interpolation to the LIS grid and keeps it in memory. 
!  At 14:00 localtime for each grid point, the code then 
!  packages the interpolated observations into the ESMF state. 
!  Note that the data files are in HDF5 format and the routine
!  employs the HDF5-subsetting tools to extract and subset data. 
!
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  integer,         parameter    :: IMSnc = 1500,IMSnr = 375
  type(ESMF_Field)              :: snowField,pertfield
  logical                       :: alarmCheck
  logical                       :: data_upd, file_exists
  logical                       :: dataflag(LIS_npes)
  logical                       :: dataflag_local
  integer                       :: c,r, p, t
  integer                       :: ftn
  character(len=LIS_CONST_PATH_LEN):: obsdir, ansa_filename, imsfile,MODISfile
  integer                       :: file_id, snwd_field_id,snwd_flag_field_id
  integer                       :: memrank = 2
  real                          :: tsnow(IMSsca_struc(n)%nc*&
       IMSsca_struc(n)%nr)
  real                          :: tsnow_flag(IMSsca_struc(n)%nc*&
       IMSsca_struc(n)%nr)
  logical*1                     :: li(IMSsca_struc(n)%nc*IMSsca_struc(n)%nr)
  real                          :: lon, lhour
  integer                       :: zone
  real                          :: ssdev(LIS_rc%ngrid(n))
  logical*1                     :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real,             pointer     :: obsl(:)
  integer                       :: gid(LIS_rc%ngrid(n))
  integer                       :: assimflag(LIS_rc%ngrid(n))
  integer                       :: status, iret, ierr
  real                          :: IMSdata(IMSnc*IMSnr)
  logical*1                     :: IMS_li(IMSnc*IMSnr)
  real                          :: IMSdata_ip(IMSnc*IMSnr)


  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       obsdir, rc=status)
  call LIS_verify(status,'Error in AttributeGet: Data Directory')

!-------------------------------------------------------------------------
!   Read the data at 0z daily. 
!-------------------------------------------------------------------------
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "IMS read alarm")

  if(alarmCheck.or.IMSsca_struc(n)%startMode) then 
     IMSsca_struc(n)%startMode = .false.
     
     call create_IMSsca_filename(imsfile,obsdir,&
          LIS_rc%yr,LIS_rc%doy)  
     inquire(file=imsfile,exist=file_exists)
     
     if(file_exists) then 
        
        ftn = LIS_getNextUnitNumber()
        
        write(LIS_logunit,*) 'Reading ',trim(imsfile)
        open(ftn,file=trim(imsfile),form='unformatted')
        read(ftn) IMSdata
        close(ftn)
        
        call LIS_releaseUnitNumber(ftn)
        
        ims_li  = .false.
        do c=1,IMSnc*IMSnr
           if(IMSdata(c).ge.0) then 
              ims_li(c) = .true. 
           endif
        enddo
        
        call bilinear_interp(LIS_rc%gridDesc(n,:),ims_li,IMSdata,&
             lo,IMSdata_ip,&
             IMSnc*IMSnr,LIS_rc%lnc(n)*LIS_rc%lnr(n),&
             LIS_domain(n)%lat,LIS_domain(n)%lon,&
             IMSsca_struc(n)%ims_w11,IMSsca_struc(n)%ims_w12, &
             IMSsca_struc(n)%ims_w21,IMSsca_struc(n)%ims_w22, &
             IMSsca_struc(n)%ims_n11,IMSsca_struc(n)%ims_n12, &
             IMSsca_struc(n)%ims_n21,IMSsca_struc(n)%ims_n22, &
             LIS_rc%udef,iret)           
        
        do r=1,LIS_rc%lnr(n)
           do c=1,LIS_rc%lnc(n)
              if(IMSdata_ip(c+LIS_rc%lnc(n)*(r-1)).lt.0.0) then 
                 IMSsca_struc(n)%snwd(c+LIS_rc%lnc(n)*(r-1)) = &
                      LIS_rc%udef
              elseif(IMSdata_ip(c+LIS_rc%lnc(n)*(r-1)).eq.0.0) then 
                 IMSsca_struc(n)%snwd(c+LIS_rc%lnc(n)*(r-1)) = &
                      0.0
              endif
           enddo
        enddo
     endif
  endif

!-------------------------------------------------------------------------
!  Update the OBS_State 
!-------------------------------------------------------------------------     

  call ESMF_StateGet(OBS_State,"Observation01",snowfield,&
       rc=status)
  call LIS_verify(status, 'Error: StateGet Observation01')
  
  call ESMF_FieldGet(snowfield,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status, 'Error: FieldGet')
  
  obsl = LIS_rc%udef 
  
  do r =1,LIS_rc%lnr(n)
     do c =1,LIS_rc%lnc(n)
        if (LIS_domain(n)%gindex(c,r) .ne. -1)then
           lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(c,r))%lon
           call LIS_localtime(LIS_rc%gmt,lon,lhour,zone)
           
           if(lhour.ge.(14-nint(LIS_rc%ts/3600.0)).and.lhour.le.14) then  
              obsl(LIS_domain(n)%gindex(c,r))=&
                   IMSsca_struc(n)%snwd(c+LIS_rc%lnc(n)*(r-1))
           endif

        end if
     end do
  end do

  dataflag_local = .false. 

! LSM-based QC
  call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
       //trim(LIS_IMSscaobsId)//char(0), & 
       n, OBS_state)

  call ESMF_StateGet(OBS_State,"Observation01",snowField,&
       rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(snowField,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status)

  do r =1,LIS_rc%lnr(n)
     do c =1,LIS_rc%lnc(n)
        if (LIS_domain(n)%gindex(c,r) .ne. -1)then
           if(obsl(LIS_domain(n)%gindex(c,r)).ne.-9999.0) then 
              dataflag_local = .true. 
           endif
        endif
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
     do t=1,LIS_rc%ngrid(n)
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
     

     call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
          rc=status)
     call LIS_verify(status, 'ESMF_StateGet for Observation01 for OBS_Pert_State failed in read_IMSsca')
     
     if(LIS_rc%ngrid(n).gt.0) then 

        call ESMF_AttributeSet(snowfield,"Grid Number",&
             gid,itemCount=LIS_rc%ngrid(n),rc=status)
        call LIS_verify(status,'Error: AttributeSet in Grid Number')
        
        call ESMF_AttributeSet(snowfield,"Assimilation Flag",&
             assimflag,itemCount=LIS_rc%ngrid(n),rc=status)
        call LIS_verify(status, 'Error: AttributeSet in Assimilation Flag')
        
     endif
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status, "Error: AttributeSet Data Update Status")
     return
  end if
  
end subroutine read_IMSsca


!BOP
!
! !ROUTINE: create_IMSsca_filename
! \label{create_IMSsca_filename}
! 
! !INTERFACE: 
subroutine create_IMSsca_filename(name, ndir, yr, doy)
  
  implicit none
! !ARGUMENTS: 
  character(len=*)      :: name
  integer           :: yr, doy
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates a timestamped IMS filename
!  
!  The arguments are: 
!  \begin{description}
!  \item[name] name of the NESDIS AMSRE soil moisture filename
!  \item[ndir] name of the NESDIS AMSRE soil moisture directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=3) :: fdoy
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fdoy, fmt='(i3.3)') doy

  name = trim(ndir)//'/'//trim(fyr)//'/ims'//trim(fyr)//trim(fdoy)//'.bin'
    
end subroutine create_IMSsca_filename

