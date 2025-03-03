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
! !ROUTINE: read_MaconLSObsdata
! \label{read_MaconLSObsdata}
!
! !REVISION HISTORY:
!  09 Jul 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_MaconLSObsdata(Obj_Space)
! !USES: 
  use ESMF
  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod,     only : LIS_logunit, LIS_verify, &
       LIS_getNextUnitNumber, LIS_releaseUnitNumber
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN
  use LIS_fileIOMod,      only : LIS_readData
  use MaconLSDataMod, only : maconlsobs_struc

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)    :: Obj_Space
!
! !DESCRIPTION:
!  
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[Obj\_State] observations state
!  \end{description}
!
!EOP
  real,    allocatable    :: lsobs(:,:)
  real,    pointer    :: obsl(:)
  type(ESMF_Field)    :: lsField
  character(len=LIS_CONST_PATH_LEN) :: lsobsdir, name
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
       lsobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(Obj_Space,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

!  call landslideobs_filename(name, lsobsdir, maconlsobs_struc(n)%site,&
!       LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn, LIS_rc%ss)
  name = trim(lsobsdir)
  
  inquire(file=name,exist=file_exists)
  
!assume that the obs is to be read at 10/01/2004, at 0z. 
  if(LIS_rc%yr.eq.2004.and.LIS_rc%mo.eq.10.and.LIS_rc%da.eq.1.and.&
       LIS_rc%hr.eq.0.and.LIS_rc%mn.eq.0.and.LIS_rc%ss.eq.0&
       .and.file_exists) then 
     
     readflag = .true. 
  else 
     readflag = .false.
  endif
  
  if (readflag) then 
     write(LIS_logunit,*)  'Reading Landslide obs data ',trim(name)
     
     call ESMF_StateGet(Obj_Space,"Landslide data",lsField,&
          rc=status)
     call LIS_verify(status)
     
     call ESMF_FieldGet(lsField,localDE=0,farrayPtr=obsl,rc=status)
     call LIS_verify(status)
     
     ftn = LIS_getNextUnitNumber()
     open(ftn,file = trim(name), form='unformatted', status='old',&
          access='direct',recl=4, iostat = istat)
     
     gridDesc(1) = 34.855
     gridDesc(2) = -84.400
     gridDesc(3) = 35.655
     gridDesc(4) = -82.850
     gridDesc(5) = 0.01
     gridDesc(6) = 0.01
     
     allocate(lsobs(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     
     if(istat.eq.0) then 
        
        call LIS_readData(n,ftn,gridDesc,lsobs)
        
        call LIS_releaseUnitNumber(ftn)
        
     endif
     
     readflag = .false.
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              obsl(LIS_domain(n)%gindex(c,r)) = lsobs(c,r)
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
  
end subroutine read_MaconLSObsdata

  
  

  

