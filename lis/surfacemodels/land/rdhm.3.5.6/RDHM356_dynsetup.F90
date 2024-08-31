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
! !ROUTINE: RDHM356_dynsetup
! \label{RDHM356_dynsetup}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   11/5/13: Shugong Wang; initial implementation for LIS 7 and RDHM356
!
! !INTERFACE:
subroutine RDHM356_dynsetup(n)
! !USES:
    use LIS_logMod, only     : LIS_logunit
    use LIS_coreMod, only    : LIS_rc, LIS_domain, LIS_surface
    use LIS_timeMgrMod, only : LIS_date2time, LIS_tick
    use LIS_fileIOMod, only  : LIS_read_param
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN
    use RDHM356_lsmMod, only : RDHM356_struc
   !use any other modules 
!
! !DESCRIPTION:
!  This routine sets up the time-dependent variables in RDHM356
!
!  The routines invoked are:
!  \begin{description}
!  \item[LIS\_read\_param](\ref{LIS_read_param}) \newline
!    retrieves LIS parameter data from NetCDF file
!  \item[RDHM356\_gen\_date\_str](\ref{RDHM356_gen_date_str}) \newline
!    generate date-time string for reading dynamic spatial parameters \newline
!    format: YYYY, YYYY\_MM, YYYY\_MM\_DD, YYYY\_MM\_DD\_HH
!  \end{description}
!
!EOP
    implicit none
    integer, intent(in) :: n
    
    integer   :: tid
    integer   :: t, gid, change, local_hour
    integer   :: locdoy, locyr, locmo, locda, lochr, locmn, locss
    real*8    :: loctime
    real      :: interp_fraction
    real      :: locgmt
    integer   :: col, row, ncount(LIS_rc%npatch(n, LIS_rc%lsm_index))
 
    integer           :: mtype 
    character(len=128):: ncvar_name  
    character(len=64) :: dt_str
    character(len=LIS_CONST_PATH_LEN):: tmp_paramfile
    real, allocatable :: placeholder(:,:)
    
    ! if tmxmn_file is not specified, do not read daily maximun and mimimun temperature 
    if(trim(RDHM356_struc(n)%tmxmn_dir) .eq. "none")  return 

    ! only read the daily data at the first time step of the first time step of a day
    if(LIS_rc%tscount(n) .eq. 1) then
      RDHM356_struc(n)%day = LIS_rc%da
    endif
    ! OHD time has one hour difference than LIS time  
    if((LIS_rc%tscount(n) .eq. 1) .or. &
       ((LIS_rc%da .ne. RDHM356_struc(n)%day) .and. (LIS_rc%hr .ge. 1))) then
        ! by default, LIS_read_param reads LIS_rc%paramfile(n). 
        ! the following line temporarily assign tmxmn_file to LIS_rc%paramfile(n)
        ! and restore its value in the end of dynsetup 
        call RDHM356_gen_date_str(dt_str, 'day')
        tmp_paramfile = LIS_rc%paramfile(n)
        LIS_rc%paramfile(n) = trim(RDHM356_struc(n)%tmxmn_dir)//"/TAIR_MAX_MIN_"//trim(dt_str)//".nc"
        mtype = LIS_rc%lsm_index

        ! allocate memory for place holder
        allocate(placeholder(LIS_rc%lnc(n), LIS_rc%lnr(n)))

        !------------------------------------!
        ! reading dynamic spatial spatial parameters !
        !------------------------------------!
        ! read dynamic spatial parameter: Tair_min
        !call RDHM356_gen_date_str(dt_str, 'day')
        !ncvar_name = trim(RDHM356_struc(n)%LDT_ncvar_Tair_min)//'_'//trim(dt_str)
        ncvar_name = "TAIR_MIN"
        write(LIS_logunit,*) 'RDHM356: reading parameter TAIR_MIN from ', trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(ncvar_name), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%tair_min = placeholder(col, row)
        enddo 

        ! read dynamic spatial parameter: Tair_max
        !call RDHM356_gen_date_str(dt_str, 'day')
        !ncvar_name = trim(RDHM356_struc(n)%LDT_ncvar_Tair_max)//'_'//trim(dt_str)
        ncvar_name = "TAIR_MAX" 
        write(LIS_logunit,*) 'RDHM356: reading parameter TAIR_MAX from ', trim(LIS_rc%paramfile(n))
        call LIS_read_param(n, trim(ncvar_name), placeholder)
        do t = 1, LIS_rc%npatch(n, mtype)
            col = LIS_surface(n, mtype)%tile(t)%col
            row = LIS_surface(n, mtype)%tile(t)%row
            RDHM356_struc(n)%rdhm356(t)%tair_max = placeholder(col, row)
        enddo 

        ! free memory for place holder
        deallocate(placeholder)
        !TODO: add code here if needed.
        LIS_rc%paramfile(n) = tmp_paramfile

        ! set the RDHM day to the new day 
        RDHM356_struc(n)%day = LIS_rc%da
    endif

end subroutine RDHM356_dynsetup
 
!BOP
!
! !ROUTINE: RDHM356_gen_date_str
! \label{RDHM356_gen_date_str}
! !INTERFACE:
subroutine RDHM356_gen_date_str(dt_str, time_option)
! !USES:
    use ESMF 
    use LIS_coreMod, only : LIS_rc
! !DESCRIPTION:
! generate date/time string for reading time-dependent variables 
!EOP
    implicit none
    character(len=*), intent(out) :: dt_str
    character(len=*), intent(in)  :: time_option
        
    integer :: year, month, day, hour, minute

    ! ESMF variables are used to adjust LIS time into OHD time 
    type(ESMF_Time) :: lis_time, ohd_time
    type(ESMF_TimeInterval) :: timeinterval1h
    integer :: rc 

    ! YEAR, MONTH, DAY and HOUR are used in solar radiation approximation if OHD forcing 
    ! data are used. OHD RDHM time is one hour behind LIS time. The following code does
    ! conversion from LIS time into RDHM time.  04/07/2014 Shugong Wang 
    
    ! initialize time 1 
    call ESMF_TimeSet(lis_time, yy = LIS_rc%yr, &
                                mm = LIS_rc%mo, &
                                dd = LIS_rc%da, &
                                h  = LIS_rc%hr, &
                                m  = LIS_rc%mn, s=0)
    
    ! set time interval of 1 hour
    call ESMF_TimeIntervalSet(timeinterval1h, h=1, rc=rc)
    ohd_time = lis_time - timeinterval1h
    call ESMF_TimeGet(ohd_time, yy=year, mm=month, dd=day, h=hour, m=minute, rc=rc);
    
!    year  = LIS_rc%yr
!    month = LIS_rc%mo
!    day   = LIS_rc%da
!    hour  = LIS_rc%hr

    if (trim(time_option) .eq. 'year') then
        write(dt_str, '(I4)') year
    elseif (trim(time_option) .eq. 'month') then
        write(dt_str, '(I4,I2.2)') year, month
    elseif (trim(time_option) .eq. 'day') then
        write(dt_str, '(I4,I2.2,I2.2)') year, month, day
    elseif (trim(time_option) .eq. 'hour') then
        write(dt_str, '(I4,I2.2,I2.2,I2.2)') year, month, day, hour
    endif
end subroutine RDHM356_gen_date_str
