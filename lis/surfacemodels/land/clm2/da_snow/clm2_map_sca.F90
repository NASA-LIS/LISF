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
! !ROUTINE: clm2_map_sca
! \label{clm2_map_sca}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar; Updated for the ESMF design
! 03May2007: K. Arseanult; Added MODIS SCA (MOD10A1) Obs Option
! 16Feb2009: K. Arsenault; Updated for MODIS SCF Obs 
!
! !INTERFACE:
subroutine clm2_map_sca( n, OBS_State, LSM_Incr_State )
! !USES:
  use ESMF
  use LIS_coreMod,  only : LIS_rc
  use LIS_logMod,   only : LIS_logunit, LIS_verify
  use clm2_varcon,   only  : bdsno, tfrz
  use clm2_lsmMod,   only : clm2_struc

  implicit none

! !ARGUMENTS: 
  integer, intent(in)      :: n
  type(ESMF_State)         :: OBS_State
  type(ESMF_State)         :: LSM_Incr_State

! !DESCRIPTION:
!
!  This subroutine directly maps the observation state to the corresponding 
!  variables in the LSM state for SCF data assimilation.
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF State for observations \newline
!  \item[LSM\_State] ESMF State for LSM state variables \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)         :: obs_sca_field
  type(ESMF_Field)         :: sweincrField
  type(ESMF_Field)         :: snodincrField

  real, pointer            :: scaobs(:)
  real, pointer            :: sweincr(:)
  real, pointer            :: snodincr(:)

  real, allocatable        :: lsm_swe(:)
  real, allocatable        :: lsm_snod(:)

  integer                  :: t, index
  integer                  :: status
  integer                  :: obs_state_count
  character*100,allocatable    :: obs_state_objs(:)

! -------------------------------------------------------------------------------
  allocate(lsm_swe(LIS_rc%ntiles(n)))
  allocate(lsm_snod(LIS_rc%ntiles(n)))

  !- SWE State Initiated::
  call ESMF_StateGet( LSM_Incr_State, "Total Snow Water", sweincrField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet( sweincrField, localDE=0,farrayPtr=sweincr, rc=status)
  call LIS_verify(status)

  !- Snow Depth State Initiated::
  call ESMF_StateGet( LSM_Incr_State,"Total Snow Depth", snodincrField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet( snodincrField, localDE=0,farrayPtr=snodincr, rc=status)
  call LIS_verify(status)

  !- Number of Obs Count::
  call ESMF_StateGet( OBS_State, itemCount=obs_state_count,rc=status)
  call LIS_verify(status)
  allocate(obs_state_objs(obs_state_count))
  call ESMF_StateGet( OBS_State, itemNameList=obs_state_objs,rc=status)
  call LIS_verify(status)

  !- Obtain the Obs State Field for the MODIS SCF Data::
  call ESMF_StateGet( OBS_State, obs_state_objs(1), obs_sca_field,&
       rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet( obs_sca_field,localDE=0,farrayPtr= scaobs,rc=status)
  call LIS_verify(status)

  ! ------------------------------------------------------------------------------

  do t=1,LIS_rc%ntiles(n)
     lsm_swe(t)  = clm2_struc(n)%clm(t)%h2osno
     lsm_snod(t) = clm2_struc(n)%clm(t)%snowdp
  end do

  !- Algorithm Option Based on Rodell and Houser (2004) ::
  write(LIS_logunit,*) " .... Performing Rodell-Houser (2004) DI .... "
  call di_method_rh04 ( LIS_rc%ntiles(n), scaobs, lsm_swe, lsm_snod )

  !else
  !      write(*,*) " ### NO Observation Operator Allowed At This Time With RH04 DI Option ... " 
  !      write(*,*) " ### Please change 'Obs Operator Option' in lis.config file to:  0  "
  !      write(*,*) " (stopping) ...";  stop
  !   end if

  !- Algorithm Option Based on K. Arsenault et al. (2007; 2009)::
  !not supported currently
#if 0
elseif ( disca_option == 1 ) then

  if ( LIS_rc%H_scf == 0 ) then
     write(LIS_logunit,*) " .... Performing DI with CLM2 Snow-Temp. Adjustment .... "
     !         call di_method_snwtmp_adj ( LIS_rc%ntiles(n), scaobs, lsm_swe, lsm_snod,   &
     !                                     lsm_snotmp1, lsm_snotmp2, lsm_snotmp3,   &
     !                                     lsm_snotmp4, lsm_snotmp5, clm2_struc(n)%clm%forc_t )

  else
     write(*,*) " ### NO Observation Operator Allowed At This Time With SnowTemp_Adj DI Option ... "
     write(*,*) " ### Please change 'Obs Operator Option' in lis.config file to:  0  "
     write(*,*) " (stopping) ...";  stop
  end if

end if
#endif

do t=1,LIS_rc%ntiles(n)
  sweincr(t)  = lsm_swe(t)  - clm2_struc(n)%clm(t)%h2osno
  snodincr(t) = lsm_snod(t) - clm2_struc(n)%clm(t)%snowdp
end do
deallocate(lsm_swe)
deallocate(lsm_snod)
deallocate(obs_state_objs)

end subroutine clm2_map_sca

!BOP
! !ROUTINE: di_method_rh04
! \label{di_method_rh04}
!
! !REVISION HISTORY:
! 16Feb2009: K. Arsenault -- Updated for MODIS SCF Obs
!
! !INTERFACE:
subroutine di_method_rh04 ( num_tiles, obs_scf, lsm_swe, lsm_snod )

  ! !USES:
  use LIS_coreMod,  only : LIS_rc
  use LIS_logMod,   only : LIS_logunit
  use clm2_varcon, only    : bdsno, tfrz

  implicit none

! !DESCRIPTION: 
!
!  This routine is based on the simple direction insertion method 
!   of Rodell and Houser (2004). 
!  
!  The arguments are:
!  \begin{description}
!  \item[num\_tiles] ?
!  \item[obs\_scf] ?
!  \item[lsm\_swe] ?
!  \item[lsm\_snod] ?
!  \end{description}
! 
!EOP

  integer, intent(in) :: num_tiles
  real,    intent(in) :: obs_scf(num_tiles)
  real,    intent(out):: lsm_swe(num_tiles)
  real,    intent(out):: lsm_snod(num_tiles)
  
  integer :: t

! -----------------------------------------------------------------------

!- Based on SCF, we update the SWE based on the rule-based approach::
  do t = 1, num_tiles 

     != Only include SCF % Values::
     if ( obs_scf(t) >= 0.    .and. &   
          obs_scf(t) <= 100.  .and. &
          obs_scf(t) .ne. LIS_rc%udef ) then

        != 1) MODIS SCF > 40 %; LSM == NO SNOW:: ADD NOMINAL 5 mm SWE LAYER --
        if ( ( obs_scf(t) > 40. ) .and. &         ! MODIS SCF Present
             ( lsm_swe(t) == 0. ) ) then          ! NO LSM SWE

           lsm_swe(t) =  5.0         ! In mm
           lsm_snod(t)= lsm_swe(t) / bdsno

           !              write(LIS_logunit,*) 'updated snow(clm2_map_sca) ... ',t,&
           !                    obs_scf(t), lsm_swe(t)


           != 2) MODIS SCF < 10%;  LSM SWE (mm) > 0. ::
        elseif( ( obs_scf(t) < 10.) .and. &       ! NO MODIS SCF
             ( lsm_swe(t) > 0. ) ) then        ! LSM SWE Present

           lsm_swe(t) = 0.0 
           lsm_snod(t)= 0.0

           !               write(LIS_logunit,*) 'scf obs < 10%; model swe > 0., ',t, &
           !                     obs_scf(t), lsm_swe(t)

        end if
     end if
  end do

end subroutine di_method_rh04

