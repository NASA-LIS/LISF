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
! !ROUTINE: get_metForcGenerated
!  \label{get_metForcGenerated}
!
! !REVISION HISTORY:
! 10Oct2014: KR Arsenault: Initial specification
!
! !INTERFACE:
subroutine get_metForcGenerated(n, findex)

! !USES:
  use LIS_coreMod,      only : LIS_rc, LIS_domain
  use LIS_logMod,       only : LIS_logunit, LIS_endrun
  use LIS_timeMgrMod,   only : LIS_tick, LIS_get_nstep
  use metForcGenerated_forcingMod, only : metForcGen_struc
  use metForcGen_VariablesMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens and reads in needed met forcing file, which
!   is from LDT-forcing processing routines.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    index of the met forcing dataset
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    Determines the LDT-generated data times
!  \item[get\_metForcGen\_filename](\ref{get_metForcGen_filename}) \newline
!    Puts together appropriate file name for data intervals
!  \end{description}
!
!EOP
  integer        :: i, gindex
  integer        :: c, r
  integer        :: kk      ! Forecast index
  integer        :: metforc_hrts
  integer        :: metforc_mnts
  character(len=LIS_CONST_PATH_LEN) :: fullfilename
  logical        :: file_exists

! Date/time parameters for file get/read:
! - LDT-generated forcing 1:
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1
! - LDT-generated forcing 2:
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2
  real    :: gmt1, gmt2, ts1, ts2

! Logical flag for determining if need to retrieve next file
  logical :: retrieve_file2

! _________________________________________________________

  retrieve_file2 = .false.
  fullfilename   = ""

! Determine metforcing timestep equivalent in hours/minutes:
  metforc_hrts = nint(metForcGen_struc%ts)/3600.
  if( nint(metForcGen_struc%ts) < 3600 ) then
     metforc_mnts = nint(metForcGen_struc%ts)
     write(*,*) "[WARN] LDT-Generated Forcing Reader HAS NOT"
     write(*,*) "[WARN]  been tested for forcing data timesteps < 3600s ..."
     write(*,*) "[WARN] Please proceed with caution and check files being opened and when."
  else
     metforc_mnts = 0
  endif

!----------------------------------------------------------
! Determine LDT-generated Forcing 1 Time 
!----------------------------------------------------------
   yr1 = LIS_rc%yr
   mo1 = LIS_rc%mo
   da1 = LIS_rc%da
   hr1 = metforc_hrts*(LIS_rc%hr/metforc_hrts)
   mn1 = metforc_mnts 
   ss1 = 0  
   ts1 = 0
   call LIS_tick( metForcGen_struc%metforc_time1, doy1, gmt1, &
        yr1, mo1, da1, hr1, mn1, ss1, ts1 )
!----------------------------------------------------------
! Determine LDT-generated Forcing 2 Time 
!----------------------------------------------------------
   yr2 = LIS_rc%yr
   mo2 = LIS_rc%mo
   da2 = LIS_rc%da
   hr2 = metforc_hrts*(LIS_rc%hr/metforc_hrts)
   mn2 = metforc_mnts  
   ss2 = 0             
   ts2 = metForcGen_struc%ts
   call LIS_tick( metForcGen_struc%metforc_time2, doy2, gmt2, &
        yr2, mo2, da2, hr2, mn2, ss2, ts2 )
!----------------------------------------------------------

#if 0
   print *, " -----------------"
   print *, " - Current LIS timestep: ",LIS_get_nstep(LIS_rc,n)
   write(*,999) "LIS   time : ", LIS_rc%time, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn
   write(*,999) "Force time1: ", metForcGen_struc%metforc_time1, mo1, da1, hr1, mn1
   write(*,999) "Force time2: ", metForcGen_struc%metforc_time2, mo2, da2, hr2, mn2
 999 format(a14,1x,f12.7,1x,4(i2,2x))
#endif

!- Determine if Forcing File 2 Needs to be Retrieved:
 ! LIS_time > Force_time1 (e.g., LIS_ts < Force_ts ):
   if( LIS_rc%time > metForcGen_struc%metforc_time1 ) then
      retrieve_file2 = .true.
    ! Assign Metforcing Data 2 to Metforcing Data 1:
      metForcGen_struc%metdata1(:,:,:)=metForcGen_struc%metdata2(:,:,:)

 ! LIS_time == Force_time1 (e.g., LIS_ts == Force_ts ):
   elseif( LIS_rc%time == metForcGen_struc%metforc_time1 ) then
     if( LIS_rc%ts == metForcGen_struc%ts ) then
       retrieve_file2 = .true.
     else
       retrieve_file2 = .false.
     endif
   endif

!- Initial LIS run timestep:
   if( LIS_get_nstep(LIS_rc,n) == 1 .or. &
        LIS_rc%rstflag(n) == 1 ) then
      LIS_rc%rstflag(n) = 0
      retrieve_file2 = .true.
      write(LIS_logunit,*)" [ERR] Currently Forcing 1 File is NOT opened in first timestep ..."
   endif

!- Retrieve LDT-generated Forcing File 2:
   if( retrieve_file2 ) then

     do kk=metForcGen_struc%st_iterId, metForcGen_struc%en_iterId

       ! Obtain LDT-generated Forcing Full Path/Filename:
       call get_metForcGen_filename( n, kk, findex, &
                yr2, mo2, da2, hr2, mn2, &
                metForcGen_struc%directory, fullfilename )

       inquire( file=trim(fullfilename), exist=file_exists)
       if( file_exists ) then
          write(LIS_logunit,*) "[INFO] Reading in forcing file: "
          write(LIS_logunit,*) trim(fullfilename)
       else
          write(LIS_logunit,*) "[ERR] -- LDT-generated forcing file: "
          write(LIS_logunit,*) "[ERR]",trim(fullfilename)//", does not exist -- "
          write(LIS_logunit,*) "[ERR] Calling End run "
          call LIS_endrun
       endif

       ! Read in LDT-generated Forcing File 2:
       !  Also, spatially reproject/reinterpolate LDT-Forcing file.

       call metForcGen_Variables_read( kk, findex, fullfilename,&
               metForcGen_struc%nc, metForcGen_struc%nr )

       do r = 1,LIS_rc%lnr(n)
          do c = 1,LIS_rc%lnc(n)
            if( LIS_domain(n)%gindex(c,r).ne.-1 ) then
               gindex = LIS_domain(n)%gindex(c,r)
 
               if( forcopts%read_airtmp ) then
                 metForcGen_struc%metdata2(kk,forcopts%index_airtmp,gindex) = &
                     forcvars_struc(n)%airtmp(c,r)
               endif

               if( forcopts%read_spechum ) then
                 metForcGen_struc%metdata2(kk,forcopts%index_spechum,gindex) = &
                     forcvars_struc(n)%spechum(c,r)
               endif

               if( forcopts%read_psurf ) then
                 metForcGen_struc%metdata2(kk,forcopts%index_psurf,gindex) = &
                     forcvars_struc(n)%psurf(c,r)
               endif

               if( forcopts%read_swdown ) then
                 metForcGen_struc%metdata2(kk,forcopts%index_swdown,gindex) = &
                     forcvars_struc(n)%swdown(c,r)
               endif

               if( forcopts%read_lwdown ) then
                 metForcGen_struc%metdata2(kk,forcopts%index_lwdown,gindex) = &
                     forcvars_struc(n)%lwdown(c,r)
               endif
 
               if( forcopts%read_uwind ) then
                 metForcGen_struc%metdata2(kk,forcopts%index_uwind,gindex) = &
                     forcvars_struc(n)%uwind(c,r)
               endif
 
               if( forcopts%read_vwind ) then
                 metForcGen_struc%metdata2(kk,forcopts%index_vwind,gindex) = &
                     forcvars_struc(n)%vwind(c,r)
               endif

               if( forcopts%read_rainf ) then
                 metForcGen_struc%metdata2(kk,forcopts%index_rainf,gindex) = &
                     forcvars_struc(n)%rainf(c,r)
!                 if(r==100 .and. c==87) print *, "metdata2 rainf: ",forcvars_struc(n)%rainf(c,r)
               endif

               if( forcopts%read_cpcp ) then
                 metForcGen_struc%metdata2(kk,forcopts%index_cpcp,gindex) = &
                     forcvars_struc(n)%cpcp(c,r)
               endif

            endif
         enddo  ! c
       enddo    ! r

       ! IF first timestep:
       if( LIS_get_nstep(LIS_rc,n) == 1 ) then
          ! Assign Metforcing Data 2 to Metforcing Data 1:
          metForcGen_struc%metdata1(:,:,:) = metForcGen_struc%metdata2(:,:,:)
       endif

     end do  ! End forecast index loop
   endif     ! Check if need to retrieve forcing file 2

end subroutine get_metForcGenerated

