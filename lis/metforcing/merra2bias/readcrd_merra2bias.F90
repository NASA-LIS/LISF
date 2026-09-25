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
! !ROUTINE: readcrd_merra2bias
! \label{readcrd_merra2bias}
!
! !REVISION HISTORY:
! 29 Jun 2026: Kristen Whitney, initial code (based on geos-itbias)
!
! !INTERFACE:
subroutine readcrd_merra2bias()
! !USES:
  use ESMF
  use LIS_logMod
  use LIS_coreMod, only       : LIS_rc,LIS_config
  use merra2bias_forcingMod, only : merra2bias_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to MERRA2bias forcing
!  from the LIS configuration file.
!
!EOP
  implicit none

  integer :: n,rc
  logical :: usedynlapserate

  call ESMF_ConfigFindLabel(LIS_config,                            &
       "MERRA2bias forcing directory:",rc=rc)
  do n = 1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,                      &
          merra2bias_struc(n)%merra2biasdir,rc=rc)
     call LIS_verify(rc,'MERRA2bias forcing directory: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config, &
       "MERRA2bias apply dynamic lapse rates:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, &
          merra2bias_struc(n)%usedynlapserate,&
          default=0, rc=rc)
     call LIS_verify(rc,&
          'MERRA2bias apply dynamic lapse rates: not defined')
  enddo

  usedynlapserate = .true.

  do n=1,LIS_rc%nnest
     if(merra2bias_struc(n)%usedynlapserate.eq.0) then
        usedynlapserate = .false.
     endif
  enddo

  if(usedynlapserate) then
     call ESMF_ConfigFindLabel(LIS_config, &
          "MERRA2bias dynamic lapse rate data directory:",rc=rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, &
             merra2bias_struc(n)%dynlapseratedir,&
             rc=rc)
        call LIS_verify(rc,&
             'MERRA2bias dynamic lapse rate data directory: not defined')
     enddo

     call ESMF_ConfigFindLabel(LIS_config, &
          "MERRA2bias dynamic lapse rate filename prefix:",rc=rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, &
             merra2bias_struc(n)%dynlapseratepfx,&
             default='/MERRA2bias.lapse_rate.hourly.', rc=rc)
        write(LIS_logunit,*) &
             '[INFO] MERRA2bias dynamic lapse rate filename prefix:', &
             trim(merra2bias_struc(n)%dynlapseratepfx)
     enddo

     call ESMF_ConfigFindLabel(LIS_config, &
          "MERRA2bias dynamic lapse rate filename suffix: ",rc=rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, &
             merra2bias_struc(n)%dynlapseratesfx,&
             default='.global.nc', rc=rc)
        write(LIS_logunit,*) &
             '[INFO] MERRA2bias dynamic lapse rate filename suffix: ', &
             trim(merra2bias_struc(n)%dynlapseratesfx)
     enddo

     call ESMF_ConfigFindLabel(LIS_config, &
          "MERRA2bias apply double-sided dynamic lapse rate cutoff:", &
          rc=rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, &
             merra2bias_struc(n)%applydynlapseratecutoff,&
             default=0, rc=rc)
        call LIS_verify(rc,&
             'MERRA2bias apply double-sided dynamic lapse rate cutoff: not defined')
     enddo

     do n=1,LIS_rc%nnest
        if(merra2bias_struc(n)%applydynlapseratecutoff.eq.1) then
           call ESMF_ConfigFindLabel(LIS_config, &
                "MERRA2bias minimum lapse rate cutoff (K/m):",rc=rc)
           call ESMF_ConfigGetAttribute(LIS_config, &
                merra2bias_struc(n)%dynlapseratemincutoff,&
                default=-0.01, rc=rc)
           call LIS_verify(rc,&
                'MERRA2bias minimum lapse rate cutoff (K/m): not defined')
           call ESMF_ConfigFindLabel(LIS_config, &
                "MERRA2bias maximum lapse rate cutoff (K/m):",rc=rc)
           call ESMF_ConfigGetAttribute(LIS_config, &
                merra2bias_struc(n)%dynlapseratemaxcutoff,&
                   default=0.01, rc=rc)
           call LIS_verify(rc,&
                'MERRA2bias maximum lapse rate cutoff (K/m): not defined')

           ! Sanity check
           if(merra2bias_struc(n)%dynlapseratemincutoff.gt. &
                merra2bias_struc(n)%dynlapseratemaxcutoff) then
              write(LIS_logunit,*) &
                   '[ERR] MERRA2bias minimum lapse rate cutoff value should be'
              write(LIS_logunit,*) &
                   '[ERR] less than the MERRA2bias maximum lapse rate cutoff value.'
              write(LIS_logunit,*) &
                   '[ERR] Note the default value is -0.01 K/m for the minimum cutoff,'
              write(LIS_logunit,*) &
                   '[ERR] and 0.01 K/m for the maximum cutoff.'
              write(LIS_logunit,*) &
                   '[ERR] Please ensure if specifying just the minimum (maximum)'
              write(LIS_logunit,*) &
                   '[ERR] cutoff value, that it is less (greater) than the'
              write(LIS_logunit,*) &
                   '[ERR] maximum (minimum) default value.'
              write(LIS_logunit,*) '[ERR] STOPPING ....'
              call LIS_endrun()
           endif
        endif
     enddo
  endif

  do n = 1,LIS_rc%nnest
     write(LIS_logunit,*) '[INFO] Using MERRA2bias forcing'
     write(LIS_logunit,*) '[INFO] MERRA2bias forcing directory: ',    &
          trim(merra2bias_struc(n)%merra2biasdir)

     merra2bias_struc(n)%merra2biastime1 = 3000.0
     merra2bias_struc(n)%merra2biastime2 = 0.0
  enddo

end subroutine readcrd_merra2bias

