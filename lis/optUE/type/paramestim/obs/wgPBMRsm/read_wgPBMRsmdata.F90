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
! !ROUTINE: read_wgPBMRsmdata
! \label{read_wgPBMRsmdata}
!
! !REVISION HISTORY:
!  09 Jul 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_wgPBMRsmdata(Obj_Space)
! !USES: 
  use ESMF
  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod,     only : LIS_logunit, LIS_verify, &
       LIS_getNextUnitNumber, LIS_releaseUnitNumber
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN
  use LIS_fileIOMod,      only : LIS_readData
  use wgPBMRsmobs_module, only : wgPBMRsm_struc

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)    :: Obj_Space
!
! !DESCRIPTION:
!  
!  reads the Walnut Gulch PBMR soil moisture data
!  and packages it into an ESMF State with certain predefined 
!  attributes
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[Obj\_State] observations state
!  \end{description}
!
!EOP
  real,    allocatable    :: sm(:,:)
  real,    pointer    :: obsl(:)
  type(ESMF_Field)    :: smField
  character(len=LIS_CONST_PATH_LEN) :: smobsdir, name
  logical             :: data_update
  logical             :: file_exists
  logical             :: readflag
  integer             :: status
  integer             :: ftn
  integer             :: c,r
  integer             :: istat
  real                :: gridDesc(6)
  integer             :: n 

  n = 1

  call ESMF_AttributeGet(Obj_Space,"Data Directory",&
       smobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(Obj_Space,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  call wgPBMRsm_filename(name, smobsdir, wgPBMRsm_struc(n)%site,&
       LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn, LIS_rc%ss)

  inquire(file=name,exist=file_exists)

!assume that the obs is available at 12Z (strategy used in Santanello et al. 2007)
  if(file_exists.and.LIS_rc%hr.eq.12.and.LIS_rc%mn.eq.0) then 
     readflag = .true. 
  else 
     readflag = .false.
  endif

  if (readflag) then 
     write(LIS_logunit,*)  'Reading WG PBMR soil moisture data ',trim(name)
     
     call ESMF_StateGet(Obj_Space,"PBMR soil moisture",smField,&
          rc=status)
     call LIS_verify(status)

     call ESMF_FieldGet(smField,localDE=0,farrayPtr=obsl,rc=status)
     call LIS_verify(status)

!-------------------------------------------------------------------------
!   Reading and mapping PBMR data (The LIS domain is assumed to be the 
!   WG domain in UTM projection. No interpolation is done. 
!-------------------------------------------------------------------------
     ftn = LIS_getNextUnitNumber()
     open(ftn,file = trim(name), form='unformatted', status='old',&
          access='direct',recl=4, iostat = istat)

     gridDesc(1) = 12
     gridDesc(2) = 3507393.0
     gridDesc(3) = 586018.0
     gridDesc(4) = 660
     gridDesc(5) = 333
     gridDesc(6) = 40.0

     allocate(sm(LIS_rc%lnc(n),LIS_rc%lnr(n)))

     if(istat.eq.0) then 

        call LIS_readData(n,ftn,gridDesc,sm)

        call LIS_releaseUnitNumber(ftn)

     endif

!-------------------------------------------------------------------------
!  Done reading PBMR data
!-------------------------------------------------------------------------
     readflag = .false.
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              obsl(LIS_domain(n)%gindex(c,r)) = sm(c,r)
           endif
        enddo
     enddo

     call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
          .true., rc=status)
     call LIS_verify(status)

  else
     call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)
     return
  end if

end subroutine read_wgPBMRsmdata

!BOP
! 
! !ROUTINE: wgPBMRsm_filename
! \label{wgPBMRsm_filename}
! 
! !INTERFACE: 
subroutine wgPBMRsm_filename(name, ndir, site, yr, mo,da,hr,mn,ss)
! !USES:   
  use LIS_timeMgrMod, only : LIS_date2time
!
! !DESCRIPTION:
!  This method generates the filename for the PBMR 
!  observational data. 
!EOP

  implicit none
  character(len=*)  :: name
  integer           :: yr, mo, da, hr,mn,ss
  integer           :: site
  character (len=*) :: ndir

  character (len=4) :: fdoy 
  character (len=1) :: fsite
  real*8            :: time
  integer           :: doy
  real              :: gmt

  call LIS_date2time(time, doy, gmt, yr, mo, da, hr,mn,ss)

  write(unit=fdoy, fmt='(i3.3)') doy
  write(unit=fsite, fmt='(i1.1)') site
  
  name = trim(ndir)//'/pbmrsm'//trim(fdoy)//'_site'//trim(fsite)//'.1gd4r'

end subroutine wgPBMRsm_filename



