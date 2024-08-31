!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: get_pptEnsFcst
!  \label{get_pptEnsFcst}
!
! !REVISION HISTORY:
! 17Dec2016: KR Arsenault: Initial specification
!
! !INTERFACE:
subroutine get_pptEnsFcst(n, findex)

! !USES:
  use LIS_coreMod,      only : LIS_rc, LIS_domain
  use LIS_logMod,       only : LIS_logunit, LIS_verify, LIS_endrun
  use LIS_timeMgrMod,   only : LIS_tick, LIS_get_nstep
  use LIS_metforcingMod,only : LIS_forc
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use pptEnsFcst_forcingMod,  only : pptensfcst_struc
  use pptEnsFcst_VariablesMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens and reads in needed ensemble forecast met forcing 
!   files, which are derived from outside processing 
!   routines/scripts.
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
!    Determines the precipitation ensemble forecast data times
!  \item[get\_pptEnsFcst\_filename](\ref{get_pptEnsFcst_filename}) \newline
!    Puts together appropriate file name for data intervals
!  \end{description}
!
!EOP
  integer        :: i, gindex, ios
  integer        :: c, r, f, m
  integer        :: metforc_hrts
  integer        :: metforc_mnts
  character(len=LIS_CONST_PATH_LEN) :: fullfilename
  logical        :: file_exists

! Date/time parameters for file get/read:
! - PPTEnsFcst forcing 1:
  integer  :: doy1, yr1, mo1, da1, hr1, mn1, ss1
! - PPTEnsFcst forcing 2:
!  integer  :: doy2, yr2, mo2, da2, hr2, mn2, ss2
  real     :: gmt1, gmt2, ts1, ts2

! Logical flag for determining if need to retrieve next file
  logical  :: retrieve_file
  integer  :: hr_int1, hr_int2

  integer  :: ensnum

  integer  :: nid, ncId, nrId, ntimesId
  integer  :: tdel, hindex, tindex

! _________________________________________________________

   retrieve_file = .false.
   pptensfcst_struc%findtime1 = 0
   pptensfcst_struc%findtime2 = 0

   ! Determine timestep in hours:
   tdel = (24/pptensfcst_struc%time_intv)

   ! Determine final time index (tindex) of the point in the file to be read:
   hindex = (LIS_rc%hr/tdel)+1                                   ! Hour index
   tindex = pptensfcst_struc%time_intv * (LIS_rc%da-1) + hindex  ! File time index

   ! Identify when time is exactly at valid time step intervals ==
   ! To open files ... for all ensemble members:
   if( (mod(LIS_rc%hr,tdel)==0) .and. &
       (LIS_rc%mn==0 .and. LIS_rc%ss==0) ) then
!     print *, " -----------------"
!     print *, ' Date, time == ',LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, tindex
     write(LIS_logunit,*)" [INFO] pptEns - yr,mo,da,hr,tindex:",&
           LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, tindex
     pptensfcst_struc%findtime1 = 1
   endif

   ! Flag date/time for when to open and read next pptEnsFcst file:
   if( LIS_get_nstep(LIS_rc,n) == 1 &
      .or. (pptensfcst_struc%reset_flag .eqv. .true.) &
      .or. (pptensfcst_struc%findtime1 == 1 ) ) then
      retrieve_file = .true.
      pptensfcst_struc%findtime1 = 0
   endif

   ! Retrieve and check if next pptEnsFcst file is present: 
   if( retrieve_file ) then

    ! Assign Metforcing Data 2 to Metforcing Data 1:
      do f=1,LIS_rc%met_nf(findex)
         pptensfcst_struc%metdata1(f,:,:)=pptensfcst_struc%metdata2(f,:,:)
      enddo

    ! Obtain PPTEnsFcst Forcing Full Path/Filename:
      fullfilename = "none"

    ! Assemble the forecast filename:

      do m = 1, pptensfcst_struc%max_ens_members
         ensnum = m
         call get_pptEnsFcst_filename( pptensfcst_struc%fcst_type,&
!               LIS_rc%syr, LIS_rc%smo, &   ! Original code (prior to 09-02-2022)
               pptensfcst_struc%fcst_inityr, pptensfcst_struc%fcst_initmo, &  ! New config file entry
               ensnum, LIS_rc%yr, LIS_rc%mo, &
               pptensfcst_struc%directory, fullfilename  )

        ! Check if pptEnsFcst file is present
        inquire( file=trim(fullfilename), exist=file_exists)
        if( file_exists ) then
           write(LIS_logunit,*) "[INFO] Reading in PPTEnsFcst Forcing File: ",&
                 trim(fullfilename)
        else  ! File missing
           write(LIS_logunit,*) "[ERR] PPTEnsFcst forcing file,"
           write(LIS_logunit,*) "[ERR] ",trim(fullfilename)," does not exist --"
           write(LIS_logunit,*) "[ERR] Calling LIS_endrun"
           call LIS_endrun
        endif
        pptensfcst_struc%findtime1 = 0
        retrieve_file = .false.

        ! Read in PPTEnsFcst Forcing File:
        !  Also, spatially reproject/reinterpolate pptEnsFcst file.
        call pptEnsFcst_Variables_read( findex, fullfilename, &
                       pptensfcst_struc%nc, pptensfcst_struc%nr, tindex )


        ! Assign final forcing fields to metdata2:
        do r = 1,LIS_rc%lnr(n)
           do c = 1,LIS_rc%lnc(n)
             if( LIS_domain(n)%gindex(c,r).ne.-1 ) then
                gindex = LIS_domain(n)%gindex(c,r)
 
                if( forcopts%read_airtmp ) then
                  pptensfcst_struc%metdata2(forcopts%index_airtmp,m,gindex) = &
                     ensfcstppt_struc(n)%airtmp(c,r)
                endif

                if( forcopts%read_spechum ) then
                  pptensfcst_struc%metdata2(forcopts%index_spechum,m,gindex) = &
                      ensfcstppt_struc(n)%spechum(c,r)
                endif

                if( forcopts%read_psurf ) then
                  pptensfcst_struc%metdata2(forcopts%index_psurf,m,gindex) = &
                      ensfcstppt_struc(n)%psurf(c,r)
                endif

                if( forcopts%read_swdown ) then
                  pptensfcst_struc%metdata2(forcopts%index_swdown,m,gindex) = &
                      ensfcstppt_struc(n)%swdown(c,r)
                endif

                if( forcopts%read_lwdown ) then
                  pptensfcst_struc%metdata2(forcopts%index_lwdown,m,gindex) = &
                      ensfcstppt_struc(n)%lwdown(c,r)
                endif

                if( forcopts%read_uwind ) then
                  pptensfcst_struc%metdata2(forcopts%index_uwind,m,gindex) = &
                      ensfcstppt_struc(n)%uwind(c,r)
                endif

                if( forcopts%read_vwind ) then
                  pptensfcst_struc%metdata2(forcopts%index_vwind,m,gindex) = &
                      ensfcstppt_struc(n)%vwind(c,r)
                endif

                if( forcopts%read_rainf ) then
                  pptensfcst_struc%metdata2(forcopts%index_rainf,m,gindex) = &
                      ensfcstppt_struc(n)%rainf(c,r)
                endif

                if( forcopts%read_cpcp ) then
                  pptensfcst_struc%metdata2(forcopts%index_cpcp,m,gindex) = &
                      ensfcstppt_struc(n)%cpcp(c,r)
                endif

            endif ! End tile mask check
         enddo    ! end col loop
       enddo      ! end row loop

     enddo  ! End ensemble member loop
   endif    ! End file/data retrieval

 ! IF first timestep:
   if( LIS_get_nstep(LIS_rc,n) == 1 ) then
    ! Assign Metforcing Data 2 to Metforcing Data 1:
     do f=1,LIS_rc%met_nf(findex)
        pptensfcst_struc%metdata1(f,:,:) = pptensfcst_struc%metdata2(f,:,:)
     end do
   endif

!----------------------------------------------------------
! Determine PPTEnsFcst Forcing Files 1 and 2 Date/Time
!----------------------------------------------------------
   ! Determine time in year.seconds at these 00Z points:
   yr1 = LIS_rc%yr
   mo1 = LIS_rc%mo
   da1 = LIS_rc%da
   hr1 = LIS_rc%hr
   mn1 = 0
   ss1 = 0
   ts1 = 0
   call LIS_tick( pptensfcst_struc%metforc_time1, doy1, gmt1, &
        yr1, mo1, da1, hr1, mn1, ss1, ts1 )

   ! Determine time in year.seconds at these 00Z points:
   yr1 = LIS_rc%yr
   mo1 = LIS_rc%mo
   da1 = LIS_rc%da
   hr1 = LIS_rc%hr
   mn1 = 0
   ss1 = 0
   ts1 = pptensfcst_struc%ts
   call LIS_tick( pptensfcst_struc%metforc_time2, doy1, gmt1, &
        yr1, mo1, da1, hr1, mn1, ss1, ts1 )
!__________________________________________________________

#if 0
   print *, " -- Current LIS timestep: ",LIS_get_nstep(LIS_rc,n)
   write(*,999) "LIS   time : ", LIS_rc%time, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn
   write(*,999) "Force time1: ", pptensfcst_struc%metforc_time1, mo1, da1, hr1, doy1
   write(*,999) "Force time2: ", pptensfcst_struc%metforc_time2, mo2, da2, hr2, doy2
   print *, " --"
 999 format(a14,1x,f12.7,1x,4(i2,2x))
#endif

end subroutine get_pptEnsFcst

