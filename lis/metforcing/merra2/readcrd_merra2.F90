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
! !ROUTINE: readcrd_merra2
! \label{readcrd_merra2}
!
! !REVISION HISTORY:
! 18 Mar 2015: James Geiger, initial code (based on merra-land)
!
! !INTERFACE:    
subroutine readcrd_merra2()
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_logMod
  use merra2_forcingMod, only : merra2_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to MERRA2 forcing
!  from the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: n,t,rc

  call ESMF_ConfigFindLabel(LIS_config,"MERRA2 forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,merra2_struc(n)%merra2dir,&
          rc=rc)
     call LIS_verify(rc,&
          'MERRA2 forcing directory: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"MERRA2 use lowest model level forcing:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,merra2_struc(n)%uselml,rc=rc)
     call LIS_verify(rc,&
          'MERRA2 use lowest model level forcing: not defined')
     if(merra2_struc(n)%uselml.eq.0) then
        call ESMF_ConfigFindLabel(LIS_config,"MERRA2 use 2m wind fields:",rc=rc) !Added for reference ET calculations
        call ESMF_ConfigGetAttribute(LIS_config,merra2_struc(n)%use2mwind,rc=rc)
     endif
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"MERRA2 use corrected total precipitation:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,merra2_struc(n)%usecorr,rc=rc)
     call LIS_verify(rc,&
          'MERRA2 use corrected total precipitation: not defined')
  enddo

  do n=1,LIS_rc%nnest
     merra2_struc(n)%usescalef = 0 
     merra2_struc(n)%nbins = 50 !currently hardcoded
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"MERRA2 apply precip scaling factors:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,merra2_struc(n)%usescalef,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"MERRA2 apply precip sampling from input distribution:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,merra2_struc(n)%usepcpsampling,rc=rc)
  enddo

  do n=1,LIS_rc%nnest
     if(merra2_struc(n)%usescalef.eq.1) then 
        allocate(merra2_struc(n)%refxrange(&
             merra2_struc(n)%ncold,merra2_struc(n)%nrold,1,&
             merra2_struc(n)%nbins))
        allocate(merra2_struc(n)%refcdf(&
             merra2_struc(n)%ncold,merra2_struc(n)%nrold,1,&
             merra2_struc(n)%nbins))
        if(merra2_struc(n)%usepcpsampling.eq.1) then 
           allocate(merra2_struc(n)%refmean(&
             merra2_struc(n)%ncold,merra2_struc(n)%nrold,1))
           allocate(merra2_struc(n)%refstdev(&
             merra2_struc(n)%ncold,merra2_struc(n)%nrold,1))
           allocate(merra2_struc(n)%refmean_ip(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
           allocate(merra2_struc(n)%refstdev_ip(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
           allocate(merra2_struc(n)%rseed(35,LIS_rc%ntiles(n)))

           do t=1,LIS_rc%ntiles(n)
              merra2_struc(n)%rseed(:,t) = -10000-t
           enddo

        endif

        allocate(merra2_struc(n)%merraxrange(&
             merra2_struc(n)%ncold,merra2_struc(n)%nrold,1,&
             merra2_struc(n)%nbins))
        allocate(merra2_struc(n)%merracdf(&
             merra2_struc(n)%ncold,merra2_struc(n)%nrold,1,&
             merra2_struc(n)%nbins))
     endif
     merra2_struc(n)%pcpscal_cmo = -1
  enddo
  
  call ESMF_ConfigFindLabel(LIS_Config, &
       label='MERRA2 precip scaling factor input file:',rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_Config, &
          merra2_struc(n)%scaleffile, rc=rc)
     if(rc.ne.0.and.merra2_struc(n)%usescalef.eq.1) then 
        call LIS_verify(rc,"MERRA2 precip scaling factor input file: not defined")
     endif
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*) '[INFO] Using MERRA2 forcing'
     write(LIS_logunit,*) '[INFO] MERRA2 forcing directory: ',&
          trim(merra2_struc(n)%merra2DIR)

     merra2_struc(n)%merra2time1 = 3000.0
     merra2_struc(n)%merra2time2 = 0.0

  enddo
end subroutine readcrd_merra2
