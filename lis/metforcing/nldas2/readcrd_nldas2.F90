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
! !ROUTINE: readcrd_nldas2
!  \label{readcrd_nldas2}
!
! !REVISION HISTORY:
! 02Feb2004; Sujay Kumar, Initial Code
! 24 Aug 2007: Chuck Alonge; Modified for use with NLDAS2 data
! 22 Jan 2012: K. Arsenault; Accommodate GES DISC, NCEP filename conventions
! 14 Mar 2014: David Mocko: Removed elevation file line in config file,
!                           as this data is now in the LDT input file.
!
! !INTERFACE:    
subroutine readcrd_nldas2()
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_logMod,  only : LIS_logunit, LIS_verify, LIS_endrun
  use nldas2_forcingMod, only : nldas2_struc

! !DESCRIPTION:
!
!  This routine reads the options specific to NLDAS-2 forcing from 
!  the LIS configuration file. 
!  
!EOP

  implicit none
  integer :: n,rc

  
  do n=1,LIS_rc%nnest
     nldas2_struc(n)%model_level_data  = 0
     nldas2_struc(n)%model_level_press = 0
     nldas2_struc(n)%model_pcp_data    = 0
     nldas2_struc(n)%model_dswrf_data  = 0
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"NLDAS2 forcing directory:",rc=rc)
  call LIS_verify(rc, 'NLDAS2 forcing directory: not defined')
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,nldas2_struc(n)%nldas2dir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"NLDAS2 data center source:",rc=rc)
  call LIS_verify(rc, 'NLDAS2 data center source: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,nldas2_struc(n)%nldas2_filesrc,&
          default="NCEP",rc=rc)
     if ((trim(nldas2_struc(n)%nldas2_filesrc).ne."GES-DISC").and.   &
         (trim(nldas2_struc(n)%nldas2_filesrc).ne."NCEP")) then
        write(LIS_logunit,*) '[ERR] NLDAS2 data center source:'
        write(LIS_logunit,*) '[ERR] must be equal to either "GES-DISC"'
        write(LIS_logunit,*) '[ERR] or to "NCEP".  Please check your'
        write(LIS_logunit,*) '[ERR] config file.  Stopping....'
        call LIS_endrun()
     endif
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"NLDAS2 use model level data:",rc=rc)
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,nldas2_struc(n)%model_level_data,rc=rc)
     call LIS_verify(rc,"NLDAS2 use model level data: is not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"NLDAS2 use model based swdown:",rc=rc)
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,nldas2_struc(n)%model_dswrf_data,rc=rc)
     call LIS_verify(rc,"NLDAS2 use model based swdown: not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"NLDAS2 use model based precip:",rc=rc)
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,nldas2_struc(n)%model_pcp_data,rc=rc)
     call LIS_verify(rc,"NLDAS2 use model based precip: not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"NLDAS2 use model based pressure:",rc=rc)
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,nldas2_struc(n)%model_level_press,rc=rc)
     call LIS_verify(rc,"NLDAS2 use model based pressure: not defined")
  enddo

  write(unit=LIS_logunit,fmt=*)'[INFO] Using NLDAS-2 forcing'

  do n=1,LIS_rc%nnest
     write(unit=LIS_logunit,fmt=*) '[INFO] NLDAS-2 forcing directory : ',trim(nldas2_struc(n)%nldas2dir)

     nldas2_struc(n)%ncold = 464
     nldas2_struc(n)%nrold = 224
     nldas2_struc(n)%nldas2time1 = 3000.0
     nldas2_struc(n)%nldas2time2 = 0.0
  enddo

end subroutine readcrd_nldas2
