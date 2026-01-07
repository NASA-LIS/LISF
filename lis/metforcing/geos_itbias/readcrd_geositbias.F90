!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_geositbias
! \label{readcrd_geositbias}
!
! !REVISION HISTORY:
! 02 Oct 2025: Fadji Maina, initial code (based on geos-it)
! 07 Jan 2026: Kristen Whitney, initial code for using dynamic lapse rate
!
! !INTERFACE:
  subroutine readcrd_geositbias()
! !USES:
  use ESMF
  use LIS_logMod
  use LIS_coreMod, only       : LIS_rc,LIS_config
  use geositbias_forcingMod, only : geositbias_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to GEOS-ITbias forcing
!  from the LIS configuration file.
!
!EOP
  implicit none

  integer :: n,t,rc,m
  logical :: usedynlapserate
  
  call ESMF_ConfigFindLabel(LIS_config,                            &
                                 "GEOS-ITbias forcing directory:",rc=rc)
  do n = 1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,                      &
                                  geositbias_struc(n)%geositbiasdir,rc=rc)
     call LIS_verify(rc,'GEOS-ITbias forcing directory: not defined')
  enddo
 
  call ESMF_ConfigFindLabel(LIS_config,"GEOS-ITbias apply dynamic lapse rates:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,geositbias_struc(n)%usedynlapserate,&
          default=0, rc=rc)
     call LIS_verify(rc,&
          'GEOS-ITbias apply dynamic lapse rates: not defined')
  enddo

  usedynlapserate = .true.

  do n=1,LIS_rc%nnest
     if(geositbias_struc(n)%usedynlapserate.eq.0) then
        usedynlapserate = .false.
     endif
  enddo
  
  if(usedynlapserate) then
     call ESMF_ConfigFindLabel(LIS_config,"GEOS-ITbias dynamic lapse rate data directory:",rc=rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config,geositbias_struc(n)%dynlapseratedir,&
             rc=rc)
        call LIS_verify(rc,&
             'GEOS-ITbias dynamic lapse rate data directory: not defined')
     enddo

     call ESMF_ConfigFindLabel(LIS_config,"GEOS-ITbias dynamic lapse rate filename prefix:",rc=rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config,geositbias_struc(n)%dynlapseratepfx,&
                      default='/GEOS-ITbias.lapse_rate.hourly.', rc=rc)
        write(LIS_logunit,*) '[INFO] GEOS-ITbias dynamic lapse rate filename prefix:', trim(geositbias_struc(n)%dynlapseratepfx)
     enddo

     call ESMF_ConfigFindLabel(LIS_config,"GEOS-ITbias dynamic lapse rate filename suffix: ",rc=rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config,geositbias_struc(n)%dynlapseratesfx,&
                      default='.global.nc', rc=rc)
        write(LIS_logunit,*) '[INFO] GEOS-ITbias dynamic lapse rate filename suffix: ', trim(geositbias_struc(n)%dynlapseratesfx)
     enddo

     call ESMF_ConfigFindLabel(LIS_config,"GEOS-ITbias apply double-sided dynamic lapse rate cutoff:",rc=rc)
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config,geositbias_struc(n)%applydynlapseratecutoff,&
                default=0, rc=rc)
        call LIS_verify(rc,&
                'GEOS-ITbias apply double-sided dynamic lapse rate cutoff: not defined')
     enddo

     do n=1,LIS_rc%nnest
        if(geositbias_struc(n)%applydynlapseratecutoff.eq.1) then
           call ESMF_ConfigFindLabel(LIS_config,"GEOS-ITbias minimum lapse rate cutoff (K/m):",rc=rc)
           call ESMF_ConfigGetAttribute(LIS_config,geositbias_struc(n)%dynlapseratemincutoff,&
                   default=-0.01, rc=rc)
           call LIS_verify(rc,&
                   'GEOS-ITbias minimum lapse rate cutoff (K/m): not defined')
           call ESMF_ConfigFindLabel(LIS_config,"GEOS-ITbias maximum lapse rate cutoff (K/m):",rc=rc)
           call ESMF_ConfigGetAttribute(LIS_config,geositbias_struc(n)%dynlapseratemaxcutoff,&
                   default=0.01, rc=rc)
           call LIS_verify(rc,&
                   'GEOS-ITbias maximum lapse rate cutoff (K/m): not defined')

           ! Sanity check
           if(geositbias_struc(n)%dynlapseratemincutoff.gt.geositbias_struc(n)%dynlapseratemaxcutoff) then
              write(LIS_logunit,*) '[ERR] GEOS-ITbias minimum lapse rate cutoff value should be'
              write(LIS_logunit,*) '[ERR] less than the GEOS-ITbias maximum lapse rate cutoff value.'
              write(LIS_logunit,*) '[ERR] Note the default value is -0.01 K/m for the minimum cutoff,'
              write(LIS_logunit,*) '[ERR] and 0.01 K/m for the maximum cutoff.'
              write(LIS_logunit,*) '[ERR] Please ensure if specifying just the minimum (maximum)'
              write(LIS_logunit,*) '[ERR] cutoff value, that it is less (greater) than the'
              write(LIS_logunit,*) '[ERR] maximum (minimum) default value.'
              write(LIS_logunit,*) '[ERR] STOPPING ....'
              call LIS_endrun()
           endif
        endif
     enddo
  endif

  do n = 1,LIS_rc%nnest
     write(LIS_logunit,*) '[INFO] Using GEOS-ITbias forcing'
     write(LIS_logunit,*) '[INFO] GEOS-ITbias forcing directory: ',    &
                           trim(geositbias_struc(n)%geositbiasdir)

     geositbias_struc(n)%geositbiastime1 = 3000.0
     geositbias_struc(n)%geositbiastime2 = 0.0
  enddo

 end subroutine readcrd_geositbias

