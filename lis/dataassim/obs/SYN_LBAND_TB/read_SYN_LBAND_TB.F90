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
! !ROUTINE: read_SYN_LBAND_TB
! \label{read_SYN_LBAND_TB}
!
! !REVISION HISTORY:
!  13 Sept 2012: Sujay Kumar; Initial specification
!
! !INTERFACE: 
subroutine read_SYN_LBAND_TB(n, OBS_State, OBS_Pert_State)
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_timeMgrMod, only : LIS_isAlarmRinging
  use LIS_logMod,     only : LIS_logunit, LIS_verify
  use LIS_logMod,     only : LIS_logunit, LIS_endrun, & 
       LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify
  use SYN_LBAND_TB_Mod, only : SYN_LBAND_TB_struc
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the synthetic L-band Tb observations and
!  into an ESMF State with certain predefined 
!  attributes
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  integer                :: status, iret, ierr
  character(len=LIS_CONST_PATH_LEN) :: lbandobsdir
  character(len=LIS_CONST_PATH_LEN) :: fname
  logical                :: alarmCheck, file_exists
  integer                :: t,c,r,i,j,p,c1,r1,nc,grid_id
  integer                :: ftn
  real,          pointer :: TbH(:), TbV(:)
  type(ESMF_Field)       :: TbHfield,TbVfield

  integer             :: gid(LIS_rc%ngrid(n))
  integer             :: assimflag(LIS_rc%ngrid(n))

  logical                :: data_update
  logical                :: data_upd_flag(LIS_npes)
  logical                :: data_upd
  integer                :: ios
  real                   :: dt
  

  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       lbandobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  data_upd = .false. 
!-------------------------------------------------------------------------
!   Read the data at 0z
!-------------------------------------------------------------------------
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "SYN LBAND TB read alarm")

  if(alarmCheck) then 
     SYN_LBAND_TB_struc(n)%TbH = LIS_rc%udef
     SYN_LBAND_TB_struc(n)%TbV = LIS_rc%udef
     
     call create_SYN_LBAND_TB_filename(lbandobsdir,&
          LIS_rc%yr, LIS_rc%mo, &
          LIS_rc%da, fname)

     inquire(file=fname,exist=file_exists)

     if(file_exists) then 
        
        write(LIS_logunit,*) 'Reading ',trim(fname)

        ftn = LIS_getNextUnitNumber()
        open(ftn,file=fname,form='unformatted')
        read(ftn) SYN_LBAND_TB_struc(n)%TbH
        read(ftn) SYN_LBAND_TB_struc(n)%TbV
        call LIS_releaseUnitNumber(ftn)
        
        call ESMF_StateGet(OBS_State,"Observation01",TbHfield,&
             rc=status)
        call LIS_verify(status, 'Error: StateGet Observation01')
        
        call ESMF_StateGet(OBS_State,"Observation02",TbVfield,&
             rc=status)
        call LIS_verify(status, 'Error: StateGet Observation02')
        
        call ESMF_FieldGet(TbHfield,localDE=0,farrayPtr=TbH,rc=status)
        call LIS_verify(status, 'Error: FieldGet for TbH')
        
        call ESMF_FieldGet(TbVfield,localDE=0,farrayPtr=TbV,rc=status)
        call LIS_verify(status, 'Error: FieldGet for TbV')
        
        TbH = LIS_rc%udef
        TbV = LIS_rc%udef

        nc = (LIS_ewe_halo_ind(n,LIS_localPet+1)-LIS_ews_halo_ind(n,LIS_localPet+1))+1
        
        do r=LIS_nss_halo_ind(n,LIS_localPet+1),LIS_nse_halo_ind(n,LIS_localPet+1)
           do c=LIS_ews_halo_ind(n,LIS_localPet+1),LIS_ewe_halo_ind(n,LIS_localPet+1)
              c1 = c-LIS_ews_halo_ind(n,LIS_localPet+1)+1
              r1 = r-LIS_nss_halo_ind(n,LIS_localPet+1)+1
              grid_id = LIS_domain(n)%gindex(c1,r1)
              if(grid_id.ne.-1) then
                 TbH(grid_id) = SYN_LBAND_TB_struc(n)%TbH(c,r)
                 TbV(grid_id) = SYN_LBAND_TB_struc(n)%TbV(c,r)
              endif
           enddo
        enddo

        do t=1,LIS_rc%ngrid(n)
           gid(t) = t
           if(TbH(t).ne.-9999.0) then 
              assimflag(t) = 1
           else
              assimflag(t) = 0
           endif
        enddo
        
        call ESMF_AttributeSet(OBS_State,"Data Update Status",&
             .true. , rc=status)
        call LIS_verify(status)
        
        if(LIS_rc%ngrid(n).gt.0) then 
           call ESMF_AttributeSet(TbHfield,"Grid Number",&
                gid,itemCount=LIS_rc%ngrid(n),rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeSet(TbHfield,"Assimilation Flag",&
                assimflag,itemCount=LIS_rc%ngrid(n),rc=status)
           call LIS_verify(status)

           call ESMF_AttributeSet(TbVfield,"Grid Number",&
                gid,itemCount=LIS_rc%ngrid(n),rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeSet(TbVfield,"Assimilation Flag",&
                assimflag,itemCount=LIS_rc%ngrid(n),rc=status)
           call LIS_verify(status)
        endif

        call ESMF_AttributeSet(OBS_State,"Data Update Status",&
             .true., rc=status)
        call LIS_verify(status)     

     else
        call ESMF_AttributeSet(OBS_State,"Data Update Status",&
             .false., rc=status)
        call LIS_verify(status)     
     endif
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)          
  endif
end subroutine read_SYN_LBAND_TB

!BOP
! !ROUTINE: create_SYN_LBAND_TB_filename
! \label{create_SYN_LBAND_TB_filename}
! 
! !INTERFACE: 
subroutine create_SYN_LBAND_TB_filename(ndir,yr, mo,da, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the synthetic L-band Tb filename based on 
!  the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[name] name of the L-band Tb file
!  \item[ndir] name of the L-band Tb directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
 
  filename = trim(ndir)//'/LBAND_TB_' &
       //trim(fyr)//trim(fmo)//trim(fda)//'0000.bin'

end subroutine create_SYN_LBAND_TB_filename




