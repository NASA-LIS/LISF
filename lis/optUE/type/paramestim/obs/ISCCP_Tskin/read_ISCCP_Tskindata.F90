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
! !ROUTINE: read_ISCCP_Tskindata
! \label{read_ISCCP_Tskindata}
!
! !REVISION HISTORY:
!  21Jun2006: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_ISCCP_Tskindata(Obj_Space)
! !USES: 
  use ESMF
  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod,     only : LIS_logunit, LIS_verify, &
       LIS_getNextUnitNumber, LIS_releaseUnitNumber
  use ISCCP_Tskinobs_module, only : isccp_tskin_struc
  use LIS_constantsMod,      only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)    :: Obj_Space
!
! !DESCRIPTION:
!  
!  reads the synthetic soil moisture observations 
!  produced from a LIS control run and packages it 
!  into an ESMF State with certain predefined 
!  attributes
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
  real, dimension(N_gswp2) :: tmp_obs

  type(ESMF_Field)    :: tskinField
  real,    pointer    :: obsl(:)
  character(len=LIS_CONST_PATH_LEN) :: tskinobsdir
  logical             :: data_update
  logical             :: file_exists
  character(len=LIS_CONST_PATH_LEN) :: name
  integer             :: i, j, istat
  integer             :: index1, iret
  logical*1           :: li(N_gswp2)
  logical*1,allocatable   :: lo(:)
  real,     allocatable   :: go(:)
  logical             :: readflag
  integer             :: status
  integer             :: ftn
  integer             :: n 

  n = 1

  call ESMF_AttributeGet(Obj_Space,"Data Directory",&
       tskinobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(Obj_Space,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  call ISCCP_Tskin_filename1(name, tskinobsdir, &
       LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr)

  inquire(file=name,exist=file_exists)

  if(file_exists.and.LIS_rc%mn.eq.0) then 
     readflag = .true. 
  else 
     readflag = .false.
  endif
  
  if (readflag) then 
     allocate(lo(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
     allocate(go(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
     write(LIS_logunit,*)  'Reading ISCCP Tskin data ',trim(name)
     
     call ESMF_StateGet(Obj_Space,"Skin Temperature",tskinField,&
          rc=status)
     call LIS_verify(status)

     call ESMF_FieldGet(tskinField,localDE=0,farrayPtr=obsl,rc=status)
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
!              tmp_lon(j) = minlon_gswp2 + land_i_gswp2*dx_gswp2
!              tmp_lat(j) = maxlat_gswp2 - land_j_gswp2*dy_gswp2
              li(index1) = .true. 
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
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             isccp_tskin_struc(n)%n11,LIS_rc%udef,iret)
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
     
!-------------------------------------------------------------------------
!  Done reading ISCCP data
!-------------------------------------------------------------------------
     readflag = .false.
     obsl = go
     deallocate(lo)
     deallocate(go)

     call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
          .true., rc=status)
     call LIS_verify(status)

  else
     call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)
     return
  end if

end subroutine read_ISCCP_Tskindata

subroutine ISCCP_Tskin_filename1(name, ndir, yr, mo,da,hr)
  
  implicit none
  character(len=*)  :: name
  integer           :: yr, mo, da, hr
  character (len=*) :: ndir
  character (len=4) :: fyr
  character (len=2) :: fmo,fda,fhr
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr
  
  name = trim(ndir)//'/Y'//trim(fyr)//'/M'//trim(fmo)//'/isccpdx_tskin.'//trim(fyr)//trim(fmo)//trim(fda)//'_'//trim(fhr)//'z.bin'

end subroutine ISCCP_Tskin_filename1



