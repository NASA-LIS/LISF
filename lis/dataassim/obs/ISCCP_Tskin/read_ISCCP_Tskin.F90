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
! !ROUTINE: read_ISCCP_Tskin
! \label{read_ISCCP_Tskin}
!
! !REVISION HISTORY:
!  21Jun2006: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_ISCCP_Tskin(n, OBS_State, OBS_Pert_State) 
! !USES: 
  use ESMF
  use LIS_coreMod,  only : LIS_rc, LIS_domain, LIS_masterproc, LIS_localPet
  use LIS_logMod,     only : LIS_logunit, LIS_endrun, & 
       LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify
  use ISCCP_Tskin_module, only : isccp_tskin_struc
  use LIS_historyMod,     only : LIS_readvar_gridded
  use LIS_pluginIndices 
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the ISCCP Tskin observations
!  and packages it into an ESMF State with certain predefined 
!  attributes. The implementation is adopted from Rolf Reichle's
!  routines in the GMAO driver
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  real, parameter     :: tskin_min = 200.0
  real, parameter     :: tskin_max = 400.0
  integer, parameter  :: N_gswp2_compressed = 15238
  integer, parameter  :: N_gswp2            = 64800

  ! land_i_gswp2 and land_j_gswp2 as stored in 
  ! ISCCP_Tskin_GSWP2_grid_V1 files (by Sarith) follow the GSWP2 convention
  ! for grid orientation, that is counting from north-to-south
  ! and from west-to-east

  real, parameter     :: minlon_gswp2 = -180.5
  real, parameter     :: maxlat_gswp2 =   90.5
  
  real, parameter     :: dx_gswp2 = 1.
  real, parameter     :: dy_gswp2 = 1.
  
  integer             :: land_i_gswp2, land_j_gswp2
  real                :: tsclr
!  real, dimension(N_gswp2_compressed) :: tmp_obs, tmp_lat, tmp_lon
  real                     :: tmp_lat, tmp_lon
  real, dimension(N_gswp2) :: tmp_obs

  type(ESMF_Field)    :: tskinfield
  
  integer             :: i, j, istat
  real,    pointer    :: obsl(:)
  integer             :: gid(LIS_rc%ngrid(n))
  integer             :: assimflag(LIS_rc%ngrid(n))

  character(len=LIS_CONST_PATH_LEN) :: tskinobsdir, name
  logical             :: data_update
  logical             :: file_exists

  logical             :: readflag
  integer             :: status

  integer             :: t
  integer             :: ftn
  integer             :: index1, iret
  logical*1           :: li(N_gswp2)
  logical*1           :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real                :: go(LIS_rc%lnc(n)*LIS_rc%lnr(n))

  integer             :: count1
  integer             :: c,r, cnt
  character(len=LIS_CONST_PATH_LEN) :: filename1, filename2
  real                :: mean_v1(LIS_rc%ngrid(n),8)
  real                :: std_v1(LIS_rc%ngrid(n),8)
  real                :: mean_v2(LIS_rc%ngrid(n),8)
  real                :: std_v2(LIS_rc%ngrid(n),8)
  
  character (len=4) :: fyr
  character (len=2) :: fmo,fhr
  
  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       tskinobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  call ISCCP_Tskin_filename(name,tskinobsdir,&
       LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr)
  
  inquire(file=name,exist=file_exists)
  
  if(file_exists.and.LIS_rc%mn.eq.0) then 
     readflag = .true. 
  else 
     readflag = .false.
  endif
  
  if (readflag) then 
     write(LIS_logunit,*)  'Reading ISCCP Tskin data ',trim(name)
     
     call ESMF_StateGet(OBS_State,"Observation01",tskinfield,&
          rc=status)
     call LIS_verify(status)

     call ESMF_FieldGet(tskinfield,localDE=0,farrayPtr=obsl,rc=status)
     call LIS_verify(status)
     
!-------------------------------------------------------------------------
!   Reading and mapping ISCCP data 
!-------------------------------------------------------------------------
     li  = .false.
     tmp_obs = -9999.0
     
     ftn = LIS_getNextUnitNumber()
     open(ftn,file = trim(name), form='unformatted', status='old',&
          iostat = istat)
     if(istat.eq.0) then 
        j = 0 
        do i=1,N_gswp2_compressed
           read(ftn) land_i_gswp2, land_j_gswp2, tsclr

           if ( (tsclr > tskin_min) .and. &
                (tsclr < tskin_max) ) then 
              j = j+1
              index1 = land_i_gswp2+(180-land_j_gswp2)*360
              tmp_obs(index1) = tsclr
              li(index1) = .true. 
              tmp_lon = minlon_gswp2 + land_i_gswp2*dx_gswp2
              tmp_lat = maxlat_gswp2 - land_j_gswp2*dy_gswp2
!              if(i.eq.5920) write(LIS_logunit,*) i, tmp_lat, tmp_lon, tsclr          
           endif
        enddo
!        open(100,file='isccp.bin',form='unformatted')
!        write(100) tmp_obs
!        close(100)
!        stop
        call LIS_releaseUnitNumber(ftn)
        
        call neighbor_interp(LIS_rc%gridDesc(n,:),li,tmp_obs,&
             lo,go,&
             isccp_tskin_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),&
             LIS_domain(n)%lat,LIS_domain(n)%lon,&
             isccp_tskin_struc(n)%n11,LIS_rc%udef,iret)       
!---------------------------------------------------------------------------------
! Impose LSM-based QC on the observations (flag data for rain, snow on the ground
!---------------------------------------------------------------------------------
        call ESMF_StateGet(OBS_State,"Observation01",tskinfield,&
             rc=status)
        call LIS_verify(status, 'Error: StateGet Observation01')
        
        call ESMF_FieldGet(tskinfield,localDE=0,farrayPtr=obsl,rc=status)
        call LIS_verify(status, 'Error: FieldGet')
        
        obsl = LIS_rc%udef 
        
        do r =1,LIS_rc%lnr(n)
           do c =1,LIS_rc%lnc(n)
              if (LIS_domain(n)%gindex(c,r) .ne. -1)then
                 obsl(LIS_domain(n)%gindex(c,r))=go(c+LIS_rc%lnc(n)*(r-1))
              end if
           end do
        end do

        call lsmdaqcobsstate(LIS_rc%lsm, LIS_isccpTskinId, &
             n, OBS_State)
        if(isccp_tskin_struc(n)%scal.eq.1) then 
           t = (LIS_rc%hr/3)+1
           write(fmo,fmt='(i2.2)') LIS_rc%mo
           write(fhr,fmt='(i2.2)') t
           filename1=trim(isccp_tskin_struc(n)%modelmean)//'_'//fmo//fhr//'.dat'
           inquire(file=filename1,exist=file_exists)
           if(file_exists) then
              
              ftn = LIS_getNextUnitNumber()
              open(ftn,file=filename1,form='unformatted')
              call LIS_readvar_gridded(ftn,n,mean_v1(:,t), status)              
              call LIS_releaseUnitNumber(ftn)

           else
              write(LIS_logunit,*) 'File not found: ',trim(filename1)
              call LIS_endrun
           endif

           filename2=trim(isccp_tskin_struc(n)%modelstd)//'_'//fmo//fhr//'.dat'
           inquire(file=filename2,exist=file_exists)
           if(file_exists) then
              
              ftn = LIS_getNextUnitNumber()
              open(ftn,file=filename2,form='unformatted')
              call LIS_readvar_gridded(ftn,n,std_v1(:,t), status)
              call LIS_releaseUnitNumber(ftn)          

           else
              write(LIS_logunit,*) 'File not found: ',trim(filename2)
              call LIS_endrun
           endif
           
           filename1=trim(isccp_tskin_struc(n)%obsmean)//'_'//fmo//fhr//'.dat'
           inquire(file=filename1,exist=file_exists)
           if(file_exists) then
              
              ftn = LIS_getNextUnitNumber()
              open(ftn,file=filename1,form='unformatted')
              call LIS_readvar_gridded(ftn,n,mean_v2(:,t), status)
              call LIS_releaseUnitNumber(ftn)              

           else
              write(LIS_logunit,*) 'File not found: ',trim(filename1)
              call LIS_endrun
           endif
           
           filename2=trim(isccp_tskin_struc(n)%obsstd)//'_'//fmo//fhr//'.dat'
           inquire(file=filename2,exist=file_exists)
           if(file_exists) then
              
              ftn = LIS_getNextUnitNumber()
              open(ftn,file=filename2,form='unformatted')              
              call LIS_readvar_gridded(ftn,n,std_v2(:,t), status)
              call LIS_releaseUnitNumber(ftn)

           else
              write(LIS_logunit,*) 'File not found: ',trim(filename2)
              call LIS_endrun
           endif
           
!           if(LIS_localPet.eq.1) then 
!              open(100,file='test.bin',form='unformatted')
!              write(100) go
!              close(100)
!              stop
!           endif


!           write(LIS_logunit,*) 'bef scal',go
           do c=1,LIS_rc%ngrid(n)
              if(mean_v1(c,t).ne.-9999.0.and.&
                   mean_v2(c,t).ne.-9999.0.and.&
                   std_v1(c,t).ne.-9999.0.and.&
                   std_v2(c,t).ne.0.0.and.&
                   obsl(c).ne.-9999.0) then 
                 obsl(c) = mean_v1(c,t) + ((obsl(c)-mean_v2(c,t))*&
                      std_v1(c,t))/std_v2(c,t)
              else
                 obsl(c) = -9999.0
              endif
           enddo
        endif
     endif

!     write(unit=fyr, fmt='(i4.4)') LIS_rc%yr
!     write(unit=fmo, fmt='(i2.2)') LIS_rc%mo
!     write(unit=fda, fmt='(i2.2)') LIS_rc%da
!     write(unit=fhr, fmt='(i2.2)') LIS_rc%hr
     
!     open(100,file='isccp_interp.bin',form='unformatted')
!     open(100,file='../../ISCCP_TSKIN/isccp_'//trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)//'00.bin',form='unformatted')
!     write(100) go
!     close(100)
!     stop

!     open(100,file='test.bin',form='unformatted')
!     write(100) go
!     close(100)
!     stop
!-------------------------------------------------------------------------
!  Done reading ISCCP data
!-------------------------------------------------------------------------
     readflag = .false.
     
  !     write(LIS_logunit,*) 'obs ',obsl
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
     call LIS_verify(status)
     
     call ESMF_AttributeSet(tskinfield,"Grid Number",&
          gid,itemCount=LIS_rc%ngrid(n),rc=status)
     call LIS_verify(status)
     
     call ESMF_AttributeSet(tskinfield,"Assimilation Flag",&
          assimflag,itemCount=LIS_rc%ngrid(n),rc=status)
     call LIS_verify(status)
     
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)
     return
  end if

end subroutine read_ISCCP_Tskin

subroutine ISCCP_Tskin_filename(name, ndir, yr, mo,da,hr)
  
  implicit none
  character(len=*)      :: name
  integer           :: yr, mo, da, hr
  character (len=*) :: ndir
  character (len=4) :: fyr
  character (len=2) :: fmo,fda,fhr
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr
  
  name = trim(ndir)//'/Y'//trim(fyr)//'/M'//trim(fmo)//'/isccpdx_tskin.'//trim(fyr)//trim(fmo)//trim(fda)//'_'//trim(fhr)//'z.bin'

end subroutine ISCCP_Tskin_filename



