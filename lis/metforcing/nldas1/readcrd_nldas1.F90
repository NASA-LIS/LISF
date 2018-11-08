!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
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
!                               field choices to the lis.config file.
!                            Changed flag from "1"/"2" to "NCEP"/"GES-DISC".
!
! !INTERFACE:    
subroutine readcrd_nldas1()
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_logMod,  only : LIS_logunit, LIS_verify, LIS_endrun
  use nldas1_forcingMod, only : nldas1_struc

! !DESCRIPTION:
!
!  This routine reads the options specific to NLDAS-1 forcing from 
!  the LIS configuration file. 
!  
!EOP

  implicit none
  integer :: n,rc
  
  call ESMF_ConfigFindLabel(LIS_config,"NLDAS1 forcing directory:",rc=rc)
  call LIS_verify(rc, 'NLDAS1 forcing directory: not defined')
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,nldas1_struc(n)%nldas1dir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"NLDAS1 data center source:",rc=rc)
  call LIS_verify(rc, 'NLDAS1 data center source: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,nldas1_struc(n)%nldas1_filesrc,rc=rc)
     if ((trim(nldas1_struc(n)%nldas1_filesrc).ne."GES-DISC").and.   &
         (trim(nldas1_struc(n)%nldas1_filesrc).ne."NCEP")) then
        write(LIS_logunit,*) 'NLDAS1 data center source:'
        write(LIS_logunit,*) ' must be equal to either "GES-DISC"'
        write(LIS_logunit,*) ' or to "NCEP".  Please check your'
        write(LIS_logunit,*) ' config file.  Stopping....'
        call LIS_endrun()
     endif
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"NLDAS1 precipitation field:",rc=rc)
  call LIS_verify(rc, 'NLDAS1 precipitation field: not defined')
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,nldas1_struc(n)%prec_field,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"NLDAS1 shortwave radiation field:",rc=rc)
  call LIS_verify(rc, 'NLDAS1 shortwave radiation field: not defined')
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,nldas1_struc(n)%swdn_field,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"NLDAS1 apply CONUS mask:",rc=rc)
  call LIS_verify(rc, 'NLDAS1 apply CONUS mask: not defined')
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,nldas1_struc(n)%applymask,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"NLDAS1 CONUS mask file:",rc=rc)
  call LIS_verify(rc, 'NLDAS1 CONUS mask file: not defined')
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,nldas1_struc(n)%conusmaskfile,rc=rc)
  enddo
  
  call ESMF_ConfigFindLabel(LIS_config,"NLDAS1 mask lower left lat:",rc=rc)
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,nldas1_struc(n)%mask_gd(1),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"NLDAS1 mask lower left lon:",rc=rc)
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,nldas1_struc(n)%mask_gd(2),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"NLDAS1 mask upper right lat:",rc=rc)
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,nldas1_struc(n)%mask_gd(3),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"NLDAS1 mask upper right lon:",rc=rc)
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,nldas1_struc(n)%mask_gd(4),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"NLDAS1 mask resolution (dx):",rc=rc)
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,nldas1_struc(n)%mask_gd(5),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"NLDAS1 mask resolution (dy):",rc=rc)
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,nldas1_struc(n)%mask_gd(6),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"NLDAS1 elevation difference map:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,nldas1_struc(n)%ediff_file,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"EDAS height map:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,nldas1_struc(n)%elevfile,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"EDAS height lower left lat:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,nldas1_struc(n)%edas_gridDesc(1),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"EDAS height lower left lon:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,nldas1_struc(n)%edas_gridDesc(2),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"EDAS height upper right lat:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,nldas1_struc(n)%edas_gridDesc(3),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"EDAS height upper right lon:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,nldas1_struc(n)%edas_gridDesc(4),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"EDAS height resolution (dx):",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,nldas1_struc(n)%edas_gridDesc(5),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"EDAS height resolution (dy):",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,nldas1_struc(n)%edas_gridDesc(6),rc=rc)
  enddo

  write(unit=LIS_logunit,fmt=*)'Using NLDAS-1 forcing'

  do n=1,LIS_rc%nnest
     write(unit=LIS_logunit,fmt=*) 'NLDAS-1 forcing directory :',nldas1_struc(n)%nldas1dir

     nldas1_struc(n)%nldas1time1 = 3000.0
     nldas1_struc(n)%nldas1time2 = 0.0
     
  enddo

end subroutine readcrd_nldas1
