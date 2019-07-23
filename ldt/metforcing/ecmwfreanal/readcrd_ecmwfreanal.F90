!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_ecmwfreanal
!  \label{readcrd_ecmwfreanal}
!
! !REVISION HISTORY:
! 26Jan2004; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readcrd_ecmwfreanal(findex)
! !USES:
  use ESMF
  use LDT_logMod,  only : LDT_logunit
  use LDT_coreMod, only : LDT_rc, LDT_config
  use ecmwfreanal_forcingMod, only : ecmwfreanal_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to ECMWF reanalysis forcing from 
!  the LDT configuration file. 
!  
!EOP
  implicit none
  integer, intent(in) :: findex
  integer :: n, rc

  if( LDT_rc%met_ecor_parms(findex) .ne. "none" .and. &
      LDT_rc%runmode == "LSM parameter processing" ) then

    call ESMF_ConfigFindLabel(LDT_config,"ECMWF Reanalysis elevation map:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,ecmwfreanal_struc(n)%elevfile,rc=rc)
    enddo
    call ESMF_ConfigFindLabel(LDT_config,"ECMWF Reanalysis elevation spatial transform:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,ecmwfreanal_struc(n)%elevtransform,rc=rc)
    enddo

    return
  endif

  call ESMF_ConfigFindLabel(LDT_config,"ECMWF Reanalysis forcing directory:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,ecmwfreanal_struc(n)%ecmwfreanaldir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"ECMWF Reanalysis maskfile:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,ecmwfreanal_struc(n)%emaskfile,rc=rc)
  enddo

  do n=1,LDT_rc%nnest
     write(LDT_logunit,*) 'ECMWF Reanalysis forcing directory: ', &
           trim(ecmwfreanal_struc(n)%ecmwfreanaldir)
     ecmwfreanal_struc(n)%fmodeltime1 = 3000.0
     ecmwfreanal_struc(n)%fmodeltime2 = 0.0
  enddo

end subroutine readcrd_ecmwfreanal
