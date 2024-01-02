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
!
! !ROUTINE: readcrd_gefs
!  \label{readcrd_gefs}
!
! !REVISION HISTORY:
! 7 Mar 2013: Sujay Kumar, initial specification
! 1 Jul 2019: K. Arsenault, expand support for GEFS forecasts
!
! !INTERFACE:    
subroutine readcrd_gefs()
! !USES:
  use ESMF
  use LIS_coreMod,     only : LIS_rc, LIS_config
  use gefs_forcingMod, only : gefs_struc
  use LIS_logMod

! !DESCRIPTION:
!
!  This routine reads the options specific to GEFS forecast forcing 
!   from the LIS configuration file. 
!  
!EOP

  implicit none
  integer :: n,rc

  write(LIS_logunit,*) ' --- '
  write(unit=LIS_logunit,fmt=*)'[INFO] Using GEFS forecast forcing'

  call ESMF_ConfigFindLabel(LIS_config,"GEFS forecast directory:",rc=rc)
  call LIS_verify(rc, 'GEFS forecast directory: not defined ')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gefs_struc(n)%gefs_dir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GEFS forecast type:",rc=rc)
  if(rc.ne.0) then 
     write(LIS_logunit,*) '[ERR] GEFS forecast type: not defined '
     write(LIS_logunit,*) '[ERR] supported options are ..'
     write(LIS_logunit,*) "[ERR] 'Operational', 'Reforecast2'"
     call LIS_endrun()
  endif

  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gefs_struc(n)%gefs_fcsttype,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GEFS forecast run mode:",rc=rc)
  call LIS_verify(rc, 'GEFS forecast run mode: not defined ')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gefs_struc(n)%gefs_runmode,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GEFS forecast grid projection:",rc=rc)
  call LIS_verify(rc, 'GEFS forecast grid projection: not defined ')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gefs_struc(n)%gefs_proj,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GEFS forecast grid resolution:",rc=rc)
  !call LIS_verify(rc, 'GEFS forecast grid resolution: not defined ')
  if(rc.ne.0) then
     write(LIS_logunit,*) '[WARN] GEFS forecast grid resolution: not defined '
     write(LIS_logunit,*) '[WARN] GEFS forecast grid resolution set to 0.25 as default. '
     do n=1,LIS_rc%nnest
        gefs_struc(n)%gefs_res = 0.25
     enddo
  else
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config,gefs_struc(n)%gefs_res,rc=rc)
        
        if( (gefs_struc(n)%gefs_fcsttype.eq.'Reforecast2' .and. (gefs_struc(n)%gefs_res.ne.0.25 .and. gefs_struc(n)%gefs_res.ne.1.0)) &
            .or. (gefs_struc(n)%gefs_fcsttype.eq.'Operational' .and. (gefs_struc(n)%gefs_res.ne.0.25 .and. gefs_struc(n)%gefs_res.ne.0.50))) then

           write(LIS_logunit,*) '[ERR] GEFS forecast grid resolution error:'
           write(LIS_logunit,*) '[ERR] The GEFS reader currently only supports the following:'
           write(LIS_logunit,*) '[ERR] -- Reforecast2 mode with 1.0 degree or 0.25 degree data.'
           write(LIS_logunit,*) '[ERR] -- Operational mode with 0.5 degree or 0.25 degree data.'
           write(LIS_logunit,*) '[ERR] Please change GEFS mode, data source, and/or data resolution.'
           write(LIS_logunit,*) '[ERR] Stopping...'
           call LIS_endrun()
        endif
     enddo
  endif

  call ESMF_ConfigFindLabel(LIS_config,"GEFS forecast number of ensemble members:",rc=rc)
  call LIS_verify(rc, 'GEFS forecast number of ensemble members: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gefs_struc(n)%max_ens_members,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GEFS pressure level field:",rc=rc)
  call LIS_verify(rc, 'GEFS pressure level field: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gefs_struc(n)%gefs_preslevel,rc=rc)
  enddo

  do n=1,LIS_rc%nnest

     write(LIS_logunit,*) '[INFO] GEFS forecast directory:  ',&
          trim(gefs_struc(n)%gefs_dir)
     write(LIS_logunit,*) '[INFO] GEFS forecast type:  ',&
          gefs_struc(n)%gefs_fcsttype
     write(LIS_logunit,*) '[INFO] GEFS forecast run mode:  ',&
          gefs_struc(n)%gefs_runmode
     write(LIS_logunit,*) '[INFO] GEFS forecast forecast projection:  ',&
          gefs_struc(n)%gefs_proj
     write(LIS_logunit,*) '[INFO] GEFS forecast number of ensemble members:',&
          gefs_struc(n)%max_ens_members
     write(LIS_logunit,*) '[INFO] GEFS pressure level field:  ',&
          gefs_struc(n)%gefs_preslevel
     write(LIS_logunit,*) '[INFO] GEFS forecast grid resolution:  ',&
          gefs_struc(n)%gefs_res

     gefs_struc(n)%fcsttime1 = 3000.0
     gefs_struc(n)%fcsttime2 = 0.0
  enddo

end subroutine readcrd_gefs


