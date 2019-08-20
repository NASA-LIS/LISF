!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_nldas1
!  \label{readcrd_nldas1}
!
! !REVISION HISTORY:
!  02 Feb 2004; Sujay Kumar, Initial Code
!  04 Mar 2007; Kristi Arsenault, Implemented EDAS Elevation Correction
!  15 Feb 2012: K. Arsenault; Accommodate GES DISC, NCEP filename conventions
!  14 Mar 2014: David Mocko: Added NLDAS-1 precipitation and radiation
!                               field choices to the ldt.config file.
!                            Changed flag from "1"/"2" to "NCEP"/"GES-DISC".
!
! !INTERFACE:    
subroutine readcrd_nldas1(findex)
! !USES:
  use ESMF
  use LDT_coreMod, only : LDT_rc, LDT_config
  use LDT_logMod,  only : LDT_logunit, LDT_verify, LDT_endrun
  use nldas1_forcingMod, only : nldas1_struc

! !DESCRIPTION:
!
!  This routine reads the options specific to NLDAS-1 forcing from 
!  the LDT configuration file. 
!  
!EOP

  implicit none
  integer, intent(in) :: findex
  integer :: n,rc
  
!  if( LDT_rc%met_ecor(findex) .ne. "none" .and. &
  if( LDT_rc%met_ecor_parms(findex) .ne. "none" .and. &
      LDT_rc%runmode == "LSM parameter processing" ) then

    call ESMF_ConfigFindLabel(LDT_config,"NLDAS1 elevation difference map:",rc=rc)
    call LDT_verify(rc,"NLDAS1 elevation difference map: not defined")
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,nldas1_struc(n)%file_elevdiff,rc=rc)
    enddo

    call ESMF_ConfigFindLabel(LDT_config,"EDAS terrain height map:",rc=rc)
    call LDT_verify(rc,"EDAS terrain height map: not defined")
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,nldas1_struc(n)%file_edaselev,rc=rc)
    enddo

    return
  endif

  call ESMF_ConfigFindLabel(LDT_config,"NLDAS1 forcing directory:",rc=rc)
  call LDT_verify(rc, 'NLDAS1 forcing directory: not defined')
  do n=1,LDT_rc%nnest    
     call ESMF_ConfigGetAttribute(LDT_config,nldas1_struc(n)%nldas1dir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"NLDAS1 data center source:",rc=rc)
  call LDT_verify(rc, 'NLDAS1 data center source: not defined')
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,nldas1_struc(n)%nldas1_filesrc,rc=rc)
     if ((trim(nldas1_struc(n)%nldas1_filesrc).ne."GES-DISC").and.   &
         (trim(nldas1_struc(n)%nldas1_filesrc).ne."NCEP")) then
        write(LDT_logunit,*) 'NLDAS1 data center source:'
        write(LDT_logunit,*) ' must be equal to either "GES-DISC"'
        write(LDT_logunit,*) ' or to "NCEP".  Please check your'
        write(LDT_logunit,*) ' config file.  Stopping....'
        call LDT_endrun()
     endif
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"NLDAS1 precipitation field:",rc=rc)
  call LDT_verify(rc, 'NLDAS1 precipitation field: not defined')
  do n=1,LDT_rc%nnest    
     call ESMF_ConfigGetAttribute(LDT_config,nldas1_struc(n)%prec_field,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"NLDAS1 shortwave radiation field:",rc=rc)
  call LDT_verify(rc, 'NLDAS1 shortwave radiation field: not defined')
  do n=1,LDT_rc%nnest    
     call ESMF_ConfigGetAttribute(LDT_config,nldas1_struc(n)%swdn_field,rc=rc)
  enddo

#if 0
  call ESMF_ConfigFindLabel(LDT_config,"NLDAS1 apply CONUS mask:",rc=rc)
  call LDT_verify(rc, 'NLDAS1 apply CONUS mask: not defined')
  do n=1,LDT_rc%nnest    
     call ESMF_ConfigGetAttribute(LDT_config,nldas1_struc(n)%applymask,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LDT_config,"NLDAS1 CONUS mask file:",rc=rc)
  call LDT_verify(rc, 'NLDAS1 CONUS mask file: not defined')
  do n=1,LDT_rc%nnest    
     call ESMF_ConfigGetAttribute(LDT_config,nldas1_struc(n)%conusmaskfile,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LDT_config,"NLDAS1 mask lower left lat:",rc=rc)
  do n=1,LDT_rc%nnest    
     call ESMF_ConfigGetAttribute(LDT_config,nldas1_struc(n)%mask_gd(1),rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LDT_config,"NLDAS1 mask lower left lon:",rc=rc)
  do n=1,LDT_rc%nnest    
     call ESMF_ConfigGetAttribute(LDT_config,nldas1_struc(n)%mask_gd(2),rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LDT_config,"NLDAS1 mask upper right lat:",rc=rc)
  do n=1,LDT_rc%nnest    
     call ESMF_ConfigGetAttribute(LDT_config,nldas1_struc(n)%mask_gd(3),rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LDT_config,"NLDAS1 mask upper right lon:",rc=rc)
  do n=1,LDT_rc%nnest    
     call ESMF_ConfigGetAttribute(LDT_config,nldas1_struc(n)%mask_gd(4),rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LDT_config,"NLDAS1 mask resolution (dx):",rc=rc)
  do n=1,LDT_rc%nnest    
     call ESMF_ConfigGetAttribute(LDT_config,nldas1_struc(n)%mask_gd(5),rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LDT_config,"NLDAS1 mask resolution (dy):",rc=rc)
  do n=1,LDT_rc%nnest    
     call ESMF_ConfigGetAttribute(LDT_config,nldas1_struc(n)%mask_gd(6),rc=rc)
  enddo
#endif

  do n=1,LDT_rc%nnest
     write(unit=LDT_logunit,fmt=*) 'NLDAS-1 forcing directory :', &
           trim(nldas1_struc(n)%nldas1dir)
     nldas1_struc(n)%nldas1time1 = 3000.0
     nldas1_struc(n)%nldas1time2 = 0.0
  enddo

end subroutine readcrd_nldas1
