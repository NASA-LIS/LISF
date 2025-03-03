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
! !ROUTINE: read_MODISsportLAI
! \label{read_MODISsportLAI}
!
! !REVISION HISTORY:
!  21 Dec 2017    Sujay Kumar; initial specification
!
! !INTERFACE: 
subroutine read_MODISsportLAI(n, k, OBS_State, OBS_Pert_State)
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_dataAssimMod
  use LIS_DAobservationsMod
  use map_utils
  use LIS_pluginIndices
  use LIS_constantsMod, only  : LIS_CONST_PATH_LEN
  use MODISsportLAI_Mod, only : MODISsportLAI_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the AMSRE soil moisture observations 
!  from NETCDF files and applies the spatial masking for dense
!  vegetation, rain and RFI. The data is then rescaled
!  to the land surface model's climatology using rescaling
!  algorithms. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  integer                :: status
  integer                :: grid_index
  character(len=LIS_CONST_PATH_LEN) :: laiobsdir
  character(len=LIS_CONST_PATH_LEN) :: fname
  logical                :: alarmCheck, file_exists
  integer                :: t,c,r,i,j,p,jj
  real,          pointer :: obsl(:)
  type(ESMF_Field)       :: laifield, pertField
  integer                :: gid(LIS_rc%obs_ngrid(k))
  integer                :: assimflag(LIS_rc%obs_ngrid(k))
  logical                :: data_update
  logical                :: data_upd_flag(LIS_npes)
  logical                :: data_upd_flag_local
  logical                :: data_upd
  real                   :: laiobs(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  integer                :: fnd
  
  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       laiobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  data_upd = .false. 

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "MODIS SPoRT LAI read alarm")

  if(alarmCheck.or.MODISsportLAI_struc(n)%startMode) then 
     MODISsportLAI_struc(n)%startMode = .false.

     call create_MODISsportLAI_filename(laiobsdir, &
          LIS_rc%yr, LIS_rc%doy, fname)
     
     inquire(file=fname,exist=file_exists)          
     if(file_exists) then 
        write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
        call read_MODISsport_LAI_data(n,k, fname,laiobs)
        fnd = 1
     else
        fnd = 0 
        write(LIS_logunit,*) '[WARN] Missing LAI file: ',trim(fname)
     endif
     
  else
     fnd = 0 
     laiobs = LIS_rc%udef
  endif
  
  if(fnd.ne.0) then 
     call ESMF_StateGet(OBS_State,"Observation01",laifield,&
          rc=status)
     call LIS_verify(status, 'Error: StateGet Observation01')
     
     call ESMF_FieldGet(laifield,localDE=0,farrayPtr=obsl,rc=status)
     call LIS_verify(status, 'Error: FieldGet')
     
     obsl = LIS_rc%udef 
     do r=1, LIS_rc%obs_lnr(k)
        do c=1, LIS_rc%obs_lnc(k)
           if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
              obsl(LIS_obs_domain(n,k)%gindex(c,r))=&
                   laiobs(c+(r-1)*LIS_rc%obs_lnc(k))
           endif
        enddo
     enddo
     
     if(fnd.eq.0) then 
        data_upd_flag_local = .false. 
     else
        data_upd_flag_local = .true. 
     endif
     
#if (defined SPMD)
     call MPI_ALLGATHER(data_upd_flag_local,1, &
          MPI_LOGICAL, data_upd_flag(:),&
          1, MPI_LOGICAL, LIS_mpi_comm, status)
#endif
     data_upd = .false.
     do p=1,LIS_npes
        data_upd = data_upd.or.data_upd_flag(p)
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
             .true. , rc=status)
        call LIS_verify(status)
        
        if(LIS_rc%obs_ngrid(k).gt.0) then 
           call ESMF_AttributeSet(laifield,"Grid Number",&
                gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeSet(laifield,"Assimilation Flag",&
                assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
           call LIS_verify(status)
           
        endif
        
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
end subroutine read_MODISsportLAI

!BOP
! 
! !ROUTINE: read_MODISsport_LAI_data
! \label{read_MODISsport_LAI_data}
!
! !INTERFACE:
subroutine read_MODISsport_LAI_data(n, k, fname, laiobs_ip)
! 
! !USES:   

  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod
  use LIS_timeMgrMod
  use MODISsportLAI_Mod, only : MODISsportLAI_struc

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  integer                       :: k
  character (len=*)             :: fname
  real                          :: laiobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real*8                        :: cornerlat(2), cornerlon(2)


! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: 
!  This subroutine reads the GLASS LAI file and applies the data
!  quality flags to filter the data. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the RTGLASS AMSR-E file
!  \item[laiobs\_ip]    soil moisture data processed to the LIS domain
! \end{description}
!
!
!EOP
  integer                 :: ftn
  integer                 :: c,r,iret
  real                    :: lai_in(MODISsportLAI_struc(n)%nc*MODISsportLAI_struc(n)%nr)
  logical*1               :: lai_data_b(MODISsportLAI_struc(n)%nc*MODISsportLAI_struc(n)%nr)
  logical*1               :: laiobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))

  ftn = LIS_getNextUnitNumber()
  open(ftn,file=trim(fname),form='unformatted',access='direct',&
       convert="little_endian",&
       recl=MODISsportLAI_struc(n)%nc*MODISsportLAI_struc(n)%nc*4)
  read(ftn,rec=1) lai_in
  call LIS_releaseUnitNumber(ftn)
  
  do r=1,MODISsportLAI_struc(n)%nr
     do c=1,MODISsportLAI_struc(n)%nc
        if(lai_in(c+(r-1)*MODISsportLAI_struc(n)%nc).gt.0) then 
           lai_data_b(c+(r-1)*MODISsportLAI_struc(n)%nc) =  .true. 
        else
           lai_data_b(c+(r-1)*MODISsportLAI_struc(n)%nc) = .false. 
        endif
     enddo
  enddo

  if(LIS_rc%obs_gridDesc(k,10).le.0.04) then 
!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!-------------------------------------------------------------------------- 
     call bilinear_interp(LIS_rc%obs_gridDesc(k,:),&
          lai_data_b, lai_in, laiobs_b_ip, laiobs_ip, &
          MODISsportLAI_struc(n)%nc*MODISsportLAI_struc(n)%nr, &
          LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
          MODISsportLAI_struc(n)%rlat, MODISsportLAI_struc(n)%rlon,&
          MODISsportLAI_struc(n)%w11,MODISsportLAI_struc(n)%w12,&
          MODISsportLAI_struc(n)%w21,MODISsportLAI_struc(n)%w22,&
          MODISsportLAI_struc(n)%n11,MODISsportLAI_struc(n)%n12,&
          MODISsportLAI_struc(n)%n21,MODISsportLAI_struc(n)%n22,LIS_rc%udef,iret)
  else
     call upscaleByAveraging(MODISsportLAI_struc(n)%nc*MODISsportLAI_struc(n)%nr,&
          LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
          LIS_rc%udef, MODISsportLAI_struc(n)%n11,&
          lai_data_b,lai_in, laiobs_b_ip, laiobs_ip)
  endif
  
!  open(100,file='test_out.bin',form='unformatted')
!  write(100) laiobs_ip
!  close(100)
!  stop

end subroutine read_MODISsport_LAI_data




!BOP
! !ROUTINE: create_MODISsportLAI_filename
! \label{create_MODISsportLAI_filename}
! 
! !INTERFACE: 
subroutine create_MODISsportLAI_filename(ndir,  yr, doy, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr, doy
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the GLASS LAI filename
!  based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the GLASS LAI data directory
!  \item[yr]  current year
!  \item[mo]  current doy
!  \item[filename] Generated GLASS LAI filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=3) :: fdoy
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fdoy, fmt='(i3.3)') doy

  filename = trim(ndir)//'/'//trim(fyr)//'/MLAI_FLT_'//&
       trim(fyr)//trim(fdoy)//'.dat'

end subroutine create_MODISsportLAI_filename





