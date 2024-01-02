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
! !ROUTINE: get_climatology
!  \label{get_climatology}
!
! !REVISION HISTORY:
! 27Sep2016: KR Arsenault: Initial specification
!
! !INTERFACE:
subroutine get_climatology(n, findex)

! !USES:
  use LIS_coreMod,      only : LIS_rc, LIS_domain
  use LIS_logMod,       only : LIS_logunit, LIS_endrun
  use LIS_timeMgrMod,   only : LIS_tick, LIS_get_nstep
  use climatology_forcingMod,  only : clim_struc
  use climatology_VariablesMod
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
!    Determines the LDT-climatology data times
!  \item[get\_climatology\_filename](\ref{get_climatology_filename}) \newline
!    Puts together appropriate file name for data intervals
!  \end{description}
!
!EOP
  integer        :: i, gindex
  integer        :: c, r
  integer        :: metforc_hrts
  integer        :: metforc_mnts
  character(len=LIS_CONST_PATH_LEN) :: fullfilename
  logical        :: file_exists

! Date/time parameters for file get/read:
! - LDT-climatology forcing 1:
  integer  :: doy1, yr1, mo1, da1, hr1, mn1, ss1
! - LDT-climatology forcing 2:
!  integer  :: doy2, yr2, mo2, da2, hr2, mn2, ss2
  real     :: gmt1, gmt2, ts1, ts2
  integer  :: tdel, tindex, doy
  real*8   :: LIS_time

! Logical flag for determining if need to retrieve next file
  logical  :: retrieve_file
  integer  :: hr_int1, hr_int2
! _________________________________________________________

  retrieve_file = .false.
  clim_struc%findtime1 = 0
  clim_struc%findtime2 = 0

!----------------------------------------------------------
! Determine LDT-climatology Forcing File 1 Date/Time
!----------------------------------------------------------
   ! Identify when time is exactly at 00:00Z ==
   if( LIS_rc%hr == 0 .and. &
       LIS_rc%mn == 0 .and. &
       LIS_rc%ss == 0 ) then
!     print *, " -----------------"
!     print *, 'Time @ 00:00Z = ',LIS_rc%hr, LIS_rc%mn, LIS_rc%ss
!     print *, ' Date @ 00Z   = ',LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%doy 
     clim_struc%findtime1 = 1
   endif
   ! Determine time in year.seconds at these 00Z points:
   yr1 = LIS_rc%yr
   mo1 = LIS_rc%mo
   da1 = LIS_rc%da
   hr1 = LIS_rc%hr
   mn1 = 0
   ss1 = 0
   ts1 = 0
!   ts1 = clim_struc%ts
   call LIS_tick( clim_struc%metforc_time1, doy1, gmt1, &
        yr1, mo1, da1, hr1, mn1, ss1, ts1 )
!__________________________________________________________

   ! Determine time in year.seconds at these 00Z points:
   yr1 = LIS_rc%yr
   mo1 = LIS_rc%mo
   da1 = LIS_rc%da
   hr1 = LIS_rc%hr
   mn1 = 0
   ss1 = 0
   ts1 = clim_struc%ts
   call LIS_tick( clim_struc%metforc_time2, doy1, gmt1, &
        yr1, mo1, da1, hr1, mn1, ss1, ts1 )
!__________________________________________________________

!#if 0
!   print *, " "
!   print *, " -- Current LIS timestep: ",LIS_get_nstep(LIS_rc,n)
!   write(*,999) "LIS   time : ", LIS_rc%time, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn
!   write(*,999) "Force time1: ", clim_struc%metforc_time1, mo1, da1, hr1, doy1
!   write(*,999) "Force time2: ", clim_struc%metforc_time2, mo2, da2, hr2, doy2
!   print *, " --"
! 999 format(a14,1x,f12.7,1x,4(i2,2x))
!#endif

   ! Flag date/time for when to open and read next climatology file:
   if( LIS_get_nstep(LIS_rc,n) == 1 &
      .or. (clim_struc%reset_flag .eqv. .true.) &
      .or. (clim_struc%findtime1 == 1 ) ) then
      retrieve_file = .true.
      clim_struc%findtime1 = 0
   endif

   ! Retrieve and check if next climatology file is present: 
   if( retrieve_file ) then

    ! Assign Metforcing Data 2 to Metforcing Data 1:
      clim_struc%metdata1(:,:)=clim_struc%metdata2(:,:)

    ! Obtain LDT-climatology Forcing Full Path/Filename:
      fullfilename = "none"

      doy = LIS_rc%doy
      
      if((mod(LIS_rc%yr,4) .eq. 0 .and. mod(LIS_rc%yr, 100).ne.0) &!leap year
           .or.(mod(LIS_rc%yr,400) .eq.0)) then 
         ! decrease the doy by a day after Feb 28       
         if(LIS_rc%doy.gt.59) then 
            doy = doy - 1
         endif
      endif

      call get_climatology_filename( doy, &
               clim_struc%directory, fullfilename  )

!      print *, " RETRIEVING NEXT FILE ... 00:00Z !!" 
!      print *, "Current LIS DOY :: ", LIS_rc%doy
!      print *, " Getting Climatology-file: ",&
!             trim(fullfilename)

      ! Check if climatology file is present
      inquire( file=trim(fullfilename), exist=file_exists)
      if( file_exists ) then
         write(LIS_logunit,*) "[INFO] Reading in Climatology Forcing File: ",&
               trim(fullfilename)
      else  ! File missing
          write(LIS_logunit,*) "[ERR] LDT-climatology forcing file, "
          write(LIS_logunit,*) "[ERR]",trim(fullfilename)//", does not exist -- "
          write(LIS_logunit,*) "[ERR] Calling End run "
          call LIS_endrun
      endif
      clim_struc%findtime1 = 0
      retrieve_file = .false.

      ! Read in LDT-climatology Forcing File:
      !  Also, spatially reproject/reinterpolate LDT-Forcing file.
      call climatology_variables_read( findex, fullfilename, &
                       clim_struc%nc, clim_struc%nr, clim_struc%ntimes )

   endif  ! End file/data retrieval

   ! Assign forcing fields to final metdata1 and 2 forcing arrays:
   if( LIS_get_nstep(LIS_rc,n) == 1  .or. &
       LIS_rc%mn == 0 .and. LIS_rc%ss == 0 ) then
     do r = 1,LIS_rc%lnr(n)
        do c = 1,LIS_rc%lnc(n)
          if( LIS_domain(n)%gindex(c,r).ne.-1 ) then
             gindex = LIS_domain(n)%gindex(c,r)

             tdel = (24/clim_struc%ntimes) 
             tindex = LIS_rc%hr/tdel+1

             if( forcopts%read_airtmp ) then
               clim_struc%metdata2(forcopts%index_airtmp,gindex) = &
                   climvars_struc(n)%airtmp(c,r,tindex)
             endif

             if( forcopts%read_spechum ) then
               clim_struc%metdata2(forcopts%index_spechum,gindex) = &
                   climvars_struc(n)%spechum(c,r,tindex)
             endif

             if( forcopts%read_psurf ) then
               clim_struc%metdata2(forcopts%index_psurf,gindex) = &
                   climvars_struc(n)%psurf(c,r,tindex)
             endif

             if( forcopts%read_swdown ) then
               clim_struc%metdata2(forcopts%index_swdown,gindex) = &
                   climvars_struc(n)%swdown(c,r,tindex)
             endif

             if( forcopts%read_lwdown ) then
               clim_struc%metdata2(forcopts%index_lwdown,gindex) = &
                   climvars_struc(n)%lwdown(c,r,tindex)
             endif

             if( forcopts%read_uwind ) then
               clim_struc%metdata2(forcopts%index_uwind,gindex) = &
                   climvars_struc(n)%uwind(c,r,tindex)
             endif

             if( forcopts%read_vwind ) then
               clim_struc%metdata2(forcopts%index_vwind,gindex) = &
                   climvars_struc(n)%vwind(c,r,tindex)
             endif

             if( forcopts%read_rainf ) then
               clim_struc%metdata2(forcopts%index_rainf,gindex) = &
                   climvars_struc(n)%rainf(c,r,tindex)
             endif

             if( forcopts%read_cpcp ) then
               clim_struc%metdata2(forcopts%index_cpcp,gindex) = &
                   climvars_struc(n)%cpcp(c,r,tindex)
             endif

          endif ! End tile mask check
       enddo    ! end col loop
     enddo      ! end row loop
   endif

 ! IF first timestep:
   if( LIS_get_nstep(LIS_rc,n) == 1 ) then
    ! Assign Metforcing Data 2 to Metforcing Data 1:
      clim_struc%metdata1(:,:) = clim_struc%metdata2(:,:)
   endif

end subroutine get_climatology

