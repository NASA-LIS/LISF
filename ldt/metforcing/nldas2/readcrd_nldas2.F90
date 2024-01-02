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
! 02 Feb 2004; Sujay Kumar, Initial Code
! 24 Aug 2007: Chuck Alonge; Modified for use with NLDAS2 data
! 22 Jan 2012: K. Arsenault; Accommodate GES DISC, NCEP filename conventions
! 14 Mar 2014: David Mocko: Removed elevation file line in config file,
!                           as this data is now in the LDT input file.
!
! !INTERFACE:    
subroutine readcrd_nldas2(findex)

! !USES:
  use ESMF
  use LDT_coreMod, only : LDT_rc, LDT_config
  use LDT_logMod,  only : LDT_logunit, LDT_verify, LDT_endrun
  use nldas2_forcingMod, only : nldas2_struc

! !DESCRIPTION:
!
!  This routine reads the options specific to NLDAS-2 forcing from 
!  the LDT configuration file. 
!  
!EOP

  implicit none
  integer, intent(in) :: findex
  integer :: n,rc

  do n=1,LDT_rc%nnest
     nldas2_struc(n)%model_level_data  = 0
     nldas2_struc(n)%model_level_press = 0
     nldas2_struc(n)%model_pcp_data    = 0
     nldas2_struc(n)%model_dswrf_data  = 0
  enddo

  if( LDT_rc%met_ecor_parms(findex) .ne."none" .and. &
      LDT_rc%runmode == "LSM parameter processing" ) then

    call ESMF_ConfigFindLabel(LDT_config,"NLDAS2 elevation difference map:",rc=rc)
    call LDT_verify(rc,"NLDAS2 elevation difference map: not defined")
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,nldas2_struc(n)%file_elevdiff,rc=rc)
    enddo

    call ESMF_ConfigFindLabel(LDT_config,"NARR terrain height map:",rc=rc)
    call LDT_verify(rc,"NARR terrain height map: not defined")
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,nldas2_struc(n)%file_narrelev,rc=rc)
    enddo
    return
  endif

  call ESMF_ConfigFindLabel(LDT_config,"NLDAS2 forcing directory:",rc=rc)
  call LDT_verify(rc, 'NLDAS2 forcing directory: not defined')
  do n=1,LDT_rc%nnest    
     call ESMF_ConfigGetAttribute(LDT_config,nldas2_struc(n)%nldas2dir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"NLDAS2 data center source:",rc=rc)
  call LDT_verify(rc, 'NLDAS2 data center source: not defined')
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,nldas2_struc(n)%nldas2_filesrc,&
          default="NCEP",rc=rc)
     if ((trim(nldas2_struc(n)%nldas2_filesrc).ne."GES-DISC").and.   &
         (trim(nldas2_struc(n)%nldas2_filesrc).ne."NCEP")) then
        write(LDT_logunit,*) 'NLDAS2 data center source:'
        write(LDT_logunit,*) ' must be equal to either "GES-DISC"'
        write(LDT_logunit,*) ' or to "NCEP".  Please check your'
        write(LDT_logunit,*) ' config file.  Stopping....'
        call LDT_endrun()
     endif
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"NLDAS2 use model level data:",rc=rc)
  do n=1,LDT_rc%nnest    
     call ESMF_ConfigGetAttribute(LDT_config,nldas2_struc(n)%model_level_data,rc=rc)
     call LDT_verify(rc,"NLDAS2 use model level data: is not defined")
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"NLDAS2 use model based swdown:",rc=rc)
  do n=1,LDT_rc%nnest    
     call ESMF_ConfigGetAttribute(LDT_config,nldas2_struc(n)%model_dswrf_data,rc=rc)
     call LDT_verify(rc,"NLDAS2 use model based swdown: not defined")
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"NLDAS2 use model based precip:",rc=rc)
  do n=1,LDT_rc%nnest    
     call ESMF_ConfigGetAttribute(LDT_config,nldas2_struc(n)%model_pcp_data,rc=rc)
     call LDT_verify(rc,"NLDAS2 use model based precip: not defined")
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"NLDAS2 use model based pressure:",rc=rc)
  do n=1,LDT_rc%nnest    
     call ESMF_ConfigGetAttribute(LDT_config,nldas2_struc(n)%model_level_press,rc=rc)
     call LDT_verify(rc,"NLDAS2 use model based pressure: not defined")
  enddo

  do n=1,LDT_rc%nnest
     write(unit=LDT_logunit,fmt=*) 'NLDAS-2 forcing directory :',&
           trim(nldas2_struc(n)%nldas2dir)

     nldas2_struc(n)%nldas2time1 = 3000.0
     nldas2_struc(n)%nldas2time2 = 0.0
  enddo

end subroutine readcrd_nldas2
