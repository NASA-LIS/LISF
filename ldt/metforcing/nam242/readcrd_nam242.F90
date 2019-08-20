!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT), v1.0
!
! See RELEASE_NOTES.txt for more information.
!
! The LDT source code and documentation are not in the public domain
! and may not be freely distributed.  Only qualified entities may receive 
! the source code and documentation. 
!
! Qualified entities must be covered by a Software Usage Agreement. 
! The Software Usage Agreement contains all the terms and conditions
! regarding the release of the LDT software.
!
! NASA GSFC MAKES NO REPRESENTATIONS ABOUT THE SUITABILITY OF THE
! SOFTWARE FOR ANY PURPOSE.  IT IS PROVIDED AS IS WITHOUT EXPRESS OR
! IMPLIED WARRANTY.  NEITHER NASA GSFC NOR THE US GOVERNMENT SHALL BE
! LIABLE FOR ANY DAMAGES SUFFERED BY THE USER OF THIS SOFTWARE.
!
! See the Software Usage Agreement for the full disclaimer of warranty.
!
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_nam242
! \label{readcrd_nam242}
!
!
! !REVISION HISTORY:
!     Sep 2012: NOHRSC/NOAA: Initial specification
!
! !INTERFACE:    
subroutine readcrd_nam242(findex)

! !USES:
  use ESMF
  use LDT_coreMod,       only : LDT_rc, LDT_config
  use LDT_logMod,        only : LDT_logunit
  use nam242_forcingMod, only : nam242_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to NAM242 forcing from 
!  the LDT configuration file. 
!  
!EOP
  implicit none
  integer, intent(in) :: findex

  integer :: n, rc

  if( LDT_rc%met_ecor_parms(findex) .ne. "none" .and. &
      LDT_rc%runmode == "LSM parameter processing" ) then

    call ESMF_ConfigFindLabel(LDT_config,"NAM242 elevation map:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,nam242_struc(n)%elevfile,rc=rc)
    enddo

    return
  endif

  call ESMF_ConfigFindLabel(LDT_config,"NAM242 forcing directory:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,nam242_struc(n)%namdir,rc=rc)
  enddo

  do n=1,LDT_rc%nnest
     write(LDT_logunit,*) 'NAM242 forcing directory : ',trim(nam242_struc(n)%namdir)
     nam242_struc(n)%namtime1  = 3000.0
     nam242_struc(n)%namtime2  = 0.0
  enddo

end subroutine readcrd_nam242
