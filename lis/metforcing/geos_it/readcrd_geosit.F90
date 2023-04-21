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
! !ROUTINE: readcrd_geosit
! \label{readcrd_geosit}
!
! !REVISION HISTORY:
! 18 Mar 2015: James Geiger, initial code (based on merra-land)
! 20 Apr 2023: David Mocko,  initial code (based on merra2)
!
! !INTERFACE:
      subroutine readcrd_geosit()
! !USES:
      use ESMF
      use LIS_logMod
      use LIS_coreMod, only       : LIS_rc,LIS_config
      use geosit_forcingMod, only : geosit_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to GEOS-IT forcing
!  from the LIS configuration file.
!
!EOP
      implicit none

      integer :: n,t,rc,m

      call ESMF_ConfigFindLabel(LIS_config,                            &
                                     "GEOS-IT forcing directory:",rc=rc)
      do n = 1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,                      &
                                      geosit_struc(n)%geositdir,rc=rc)
         call LIS_verify(rc,'GEOS-IT forcing directory: not defined')
      enddo

      call ESMF_ConfigFindLabel(LIS_config,                            &
                        "GEOS-IT use lowest model level forcing:",rc=rc)
      do n = 1,LIS_rc%nnest
         call ESMF_ConfigGetAttribute(LIS_config,                      &
                                      geosit_struc(n)%uselml,rc=rc)
         call LIS_verify(rc,                                           &
               'GEOS-IT use lowest model level forcing: not defined')
      enddo

      do m = 1, LIS_rc%nnest
        if (geosit_struc(m)%uselml.eq.0) then
           call ESMF_ConfigFindLabel(LIS_config,                         &
                                    "GEOS-IT use 2m wind fields:",rc=rc)
           do n = 1,LIS_rc%nnest
              call ESMF_ConfigGetAttribute(LIS_config,                   &
                                        geosit_struc(n)%use2mwind,rc=rc)
              call LIS_verify(rc,                                        &
                              'GEOS-IT use 2m wind fields: not defined')
           enddo
        endif
      enddo

      do n = 1,LIS_rc%nnest
         write(LIS_logunit,*) '[INFO] Using GEOS-IT forcing'
         write(LIS_logunit,*) '[INFO] GEOS-IT forcing directory: ',    &
                               trim(geosit_struc(n)%geositdir)

         geosit_struc(n)%geosittime1 = 3000.0
         geosit_struc(n)%geosittime2 = 0.0
      enddo

      end subroutine readcrd_geosit

