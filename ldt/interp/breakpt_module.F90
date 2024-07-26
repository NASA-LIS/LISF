!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! Program Breakpt.F90
!
! Original Copyright October 2004
!
! Current Copyright May 2005
!
! by Matthew Garcia
!    Research Associate, GEST Center 
!    NASA Hydrological Sciences Branch
!    GSFC, Code 614.3
!    Greenbelt, Maryland 20771
!    mgarcia@hsb.gsfc.nasa.gov
!
!----------------------------------------------------------------------------------------
!
! This source code and documentation are in the public domain,
! available without fee for educational, research, non-commercial and
! commercial purposes.  Users may distribute the binary or source
! code to third parties provided this statement appears on all copies and
! that no charge is made for such copies.
!
! NASA/GSFC AND THE AUTHOR(S) MAKE NO REPRESENTATIONS ABOUT THE SUITABILITY 
! OF THIS SOFTWARE FOR ANY PURPOSE.  IT IS PROVIDED AS IS WITHOUT EXPRESS 
! OR IMPLIED WARRANTY.  NEITHER NASA/GSFC, THE U.S. GOVERNMENT, NOR THE 
! INDIVIDUAL AUTHOR(S) SHALL BE LIABLE FOR ANY DAMAGES SUFFERED BY THE USER 
! OF THIS SOFTWARE.
!
!----------------------------------------------------------------------------------------
!
! Version History
!
!  -Date-   -Programmer-    -Comments-
!  20Sep04  Matthew Garcia  Initial specification completed
!  13Oct04  Matthew Garcia  Revisions to reduce memory usage
!  13Oct04  Matthew Garcia  Additional documentation for release
!  14Oct04  Matthew Garcia  Final pre-release specification
!  14Oct04  Matthew Garcia  Delivered to USDA ARS (SWRC, Tucson AZ)
!  14Oct04  Matthew Garcia  Delivered to USACE TEC for ARMS preprocessor use  
!  01Mar05  Matthew Garcia  Revised for integration into NASA/GSFC LISv4.1
!  10May05  Matthew Garcia  Bug fixes for long inter-break periods
!  11May05  Matthew Garcia  Extracted from LISv4.1 for re-release to USDA ARS
!  19May05  Matthew Garcia  One (1) bug fix for intervals over month-end in leapyear
!
!----------------------------------------------------------------------------------------
!
! Processing of Precipitation Breakpoint Data Files
!
! - module breakpoint       -- driver, calls functions/subroutines
!  - function checkintfile  -- find preprocessed data, if it exists
!   - function filedate     -- extract file timestamps <-- CONTAINS SYSTEM DEPENDENCY
!    - module dfport        -- Compaq Fortran library dependency <-- SYSTEM MODULE
!  - subroutine readcard    -- cardfile information (ARS file, dates, interval, output)   
!   - function days         -- calculate number of days between two dates
!    - function leapyr      -- logical check:  is this a leap year?
!  - subroutine readheader  -- header information (# lines, including gauge numbers)             
!  - subroutine procstndata -- process observations (precipitation breakpoint) data       
!   - function leapyr       -- logical check:  is this a leap year?
!  - subroutine breakagg    -- aggregation/disaggregation of precipitation breakpoints    
!   - subroutine intindex   -- interval index corresponding to a particular date
!    - function days        -- calculate number of days between two dates
!     - function leapyr     -- logical check:  is this a leap year?
!  - subroutine metric      -- checking calculated event total precip against records  
!  - subroutine writestnrec -- writing of station records to file(s)                      
!   - function datestr      -- format of date/time as string for output
!    - function leapyr      -- logical check:  is this a leap year?
!
! *****USER NOTE:  FUNCTION "FILEDATE" CONTAINS SYSTEM DEPENDENCY*****
!
! This dependency has been written specifically for the corresponding library routine 
!  "dfport" that is included with Compaq Visual Fortram v6.6a for MS Windows.  For 
!  porting to another compiler/system, you will need to edit the dependency so that 
!  file timestamps are obtained according to the target system's procedures/syntax. 
!
! ***USER NOTE:  THIS PROGRAM WAS WRITTEN SPECIFICALLY FOR USDA/ARS/SWRC DATAFILES***
!
! That is, for precipitation breakpoint datafiles from the U.S. Department of 
!  Agriculture (USDA), Agricultural Research Service (ARS), Southwest Watershed 
!  Research Center (SWRC).  You may need to modify some aspects of this program 
!  related to dataline format and the number/contents of datafile header lines for 
!  use with datafiles from other data sources, including other USDA/ARS watershed 
!  research centers.
!
! ***USER NOTE:  THIS PROGRAM REQUIRES A "CARDFILE" FOR EXECUTION***
!
! You should have received a sample cardfile with this program and/or executable.
!  If not, an example of the cardfile contents (with explanations) is given here:
!
!  WG_20020101-20040916.dat  <-- name of raw breakpoint data file in current directory
!  01/01/2002                <-- starting date of processing
!  09/16/2004                <-- ending date of processing
!  85                        <-- number of gauges (confirmed with data header information)
!  20                        <-- desired regular interval, in minutes
!  1                         <-- desired output format 
!                                1 = single file including all gauges
!                                2 = individual file for each gauge
!
! The name of the required cardfile is currently hardcoded in subroutine "readcard" as 
!  "Breakpt.crd"
!
! ***USER NOTE:  THIS PROGRAM CREATES NEW FILES***
! 
! This program creates output files according to user specification (see cardfile 
!  description above):
!   <datafile>.out -- the results of breakpoint disaggregation to regular intervals
! *NOTE* Depending on the user specification of output format in the card file, for 
!  N gauges there will be either one output file for all gauges or N output files, 
!  one for each gauge.
!
! This program creates one other file of interest to the user named:
!   <datafile>.log -- execution diagnostic messages and results, including QC results
! This file may be *large*, depending on the size of the raw breakpoint data file and 
!  the specified interval length.
!
! This program also creates two new intermediate processing files named:
!   <datafile>.ppd -- preprocessed breakpoints, time-adjusted for raw data "offset"
!   <datafile>.idx -- gauge-specific total numbers of breakpoints and events
! These files are currently *not* removed at the end of program execution.  Both of 
!  these files aid in subsequent program executions on the same data file by providing 
!  data that needs to be calculated anyway, so we just read that data from the inter-
!  mediate files instead.  The idea here is that reprocessing of the same input data-
!  file for a new interval of interest (e.g. 15 minutes instead of 60 minutes) will thus 
!  require changing only the card file--full recalculation of preprocessing fields (the 
!  contents of these intermediate files) is not necessary because you would get the same 
!  intermediate file contents as for the previous execution.
!
!----------------------------------------------------------------------------------------
!
!program Breakpt
module breakpoint_module
!
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  implicit none
!
!  logical :: lerr
!  character(80) :: brkfile
!  character(10) :: datebeg,dateend
!  integer :: ngauges,intmins,outfmt
!
!  print *,'Program Breakpt.F90 -- starting precipitation breakpoint disaggregation'
!  lerr = .FALSE.
!  call readcard(brkfile,datebeg,dateend,ngauges,intmins,outfmt,lerr)
!  if (.not.lerr) then
!    call brkprecip(brkfile,datebeg,dateend,ngauges,intmins,outfmt)
!  end if
!
!  print *,'Program Breakpt.F90 -- done'
!
!
  contains
!
!
! begin internal functions 
!
!
#if ( defined PC_NT )
integer function filedate(log1,filename)
!
! dependencies
  use dfport
!
! argument variables
  integer, intent(INOUT) :: log1
  character(len=LDT_CONST_PATH_LEN) :: filename
!
! local variables
  integer, dimension(12) :: fileinfo
  character*80 :: wrtline
  character*3 :: errnum
  integer :: istat
!
  filedate = 0
  istat = stat(filename,fileinfo)
  if (istat.eq.0) then ! successful inquiry of file
    filedate = fileinfo(10)
  else
    write(errnum,'(I3)') istat
    wrtline = 'Error '//trim(adjustl(errnum))//' accessing file '&
                  //trim(adjustl(filename))
        print *,wrtline
    write(log1,'(A80)') wrtline
    stop
  end if
!
  return
end function filedate
#endif
!
!
logical function checkintfile(log1,dataf1,dataf2)
!
! argument variables
  integer, intent(INOUT) :: log1
  character*80, intent(INOUT) :: dataf1,dataf2
!
! local variables
  logical :: fileexists = .FALSE.
  integer :: datef1,datef2,dateprod
!
  checkintfile = .FALSE.
  datef1 = 0
  datef2 = 0
  inquire(file=trim(adjustl(dataf2)), exist=fileexists)
  if (fileexists) then
#if ( defined PC_NT )
    datef1 = filedate(log1,dataf1)
    datef2 = filedate(log1,dataf2)
#endif
    dateprod = datef1 * datef2
    if (dateprod.ne.0.and.datef2.gt.datef1) checkintfile = .TRUE.
  end if
!
  return
end function checkintfile
!
!
logical function leapyr(yy,mm)
!
! argument variables
  integer, intent(IN) :: yy,mm
!
  leapyr = .FALSE.
  if (mm.eq.2.and.mod(yy,4).eq.0) then
        if (mod(yy,100).ne.0.or.mod(yy,1000).eq.0) then 
          leapyr = .TRUE.
        end if
  end if
!
  return
end function leapyr
!
!
integer function days(dts)
!
! argument variables
  integer, intent(IN) :: dts(:,:)
!
! local variables
  integer :: months
  integer :: curryy,currmm
  integer, dimension(12) :: mmdays = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
  integer :: i
!
  days = 0
  months = 12 * (dts(2,1) - dts(1,1)) + dts(2,2) - dts(1,2) + 1
  curryy = dts(1,1)
  currmm = dts(1,2)
  do i = 1,months
    if (currmm.gt.12) then
          curryy = curryy + 1
      currmm = 1
        end if
        if (i.eq.1.and.months.eq.1) then ! dts(1,:) and dts(2,:) occur in the same month
      days = dts(2,3) - dts(1,3) + 1
        else if (i.eq.1.and.months.gt.1) then ! dts(1,:) and dts(2,:) occur in different months
          days = mmdays(currmm) - dts(1,3) + 1
        else if (i.eq.months.and.months.gt.1) then
          days = days + dts(2,3)
        else
      days = days + mmdays(currmm)
        end if
        if (leapyr(curryy,currmm) .and. currmm.gt.2) then
          days = days + 1
        end if
    currmm = currmm + 1
  end do
!
  return
end function days
!
!
character(16) function datestr(k,intm,dts)
!
! argument variables
  integer, intent(IN) :: k,intm
  integer, intent(IN) :: dts(:,:)
!
! local variables
  integer, dimension(12) :: mmdays = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
  integer :: intdd,inthh
  integer :: currmn,currhh,currdd,currmm,curryy
  integer :: dys
  integer :: j
  character(2) :: mn,hh,dd,mm
  character(4) :: yyyy
!
  intdd = 1440 / intm
  inthh = int(intdd / 24)
  currmn = mod(k,inthh) * intm
  currhh = int(mod(k,intdd) / inthh)
  currdd = dts(1,3)
  currmm = dts(1,2)
  curryy = dts(1,1)
  dys = int(k / intdd)
  do j = 1,dys
    currdd = currdd + 1
        if (currdd.gt.mmdays(currmm)) then
          if (leapyr(curryy,currmm)) then
            currmm = currmm + 0
          else
            currmm = currmm + 1
                currdd = 1
      end if
        end if
    if (currmm.gt.12) then
          curryy = curryy + 1
          currmm = 1
    end if
  end do
  select case (currmn)
    case (10:)
      write(mn,'(I2)') currmn
        case (:9)
      write(mn,'(I1)') currmn
          write(mn,'(A2)') '0'//trim(adjustl(mn))
  end select
  select case (currhh)
    case (10:)
      write(hh,'(I2)') currhh
        case (:9)
      write(hh,'(I1)') currhh
          write(hh,'(A2)') '0'//trim(adjustl(hh))
  end select
  select case (currdd)
    case (10:)
      write(dd,'(I2)') currdd
        case (:9)
      write(dd,'(I1)') currdd
          write(dd,'(A2)') '0'//trim(adjustl(dd))
  end select
  select case (currmm)
    case (10:)
      write(mm,'(I2)') currmm
        case (:9)
      write(mm,'(I1)') currmm
          write(mm,'(A2)') '0'//trim(adjustl(mm))
  end select
  write(yyyy,'(I4)') curryy
  write(datestr,'(A16)') mm//'/'//dd//'/'//yyyy//' '//hh//':'//mn
!
  return
end function datestr
!
!
! end internal functions 
!
! begin internal subroutines
!
!
! Duplicate declaration... check if this is needed
#if 0 
subroutine readcard(brkfile,datebeg,dateend,ng,intmins,outfmt,error)
!
! argument variables
  character(80), intent(OUT) :: brkfile
  character(10), intent(OUT) :: datebeg,dateend
  integer, intent(OUT) :: ng,intmins,outfmt
  logical, intent(OUT) :: error
!
! local variables
  logical :: fileexists
  character(11) :: cardfile
  integer :: funit = 10
!
  cardfile = 'Breakpt.crd'
  inquire(file=trim(cardfile), exist=fileexists)
  error = .not.fileexists
  if (error) then
    print *,'Card file '//trim(cardfile)//' does not exist'
    return
  else
    print *,'Opening descriptor file:  '//cardfile
    open(funit,file=cardfile,status='old')
    read(funit,'(A80)') brkfile
    read(funit,'(A10)') datebeg
    read(funit,'(A10)') dateend
    read(funit,*) ng
    read(funit,*) intmins
    read(funit,*) outfmt
  end if
!
  return
end subroutine readcard
#endif
!
!
subroutine brkprecip(brkfile,datebeg,dateend,ngauges,intmins,outfmt)
!
! argument variables
  character(len=*), intent(INOUT) :: brkfile
  character*10, intent(INOUT) :: datebeg,dateend
  integer, intent(INOUT) :: ngauges,intmins,outfmt
!
! local variables
  logical :: brkexists,ppdexists,idxexists,intfiles
  character*80 :: wrtline
  character(len=LDT_CONST_PATH_LEN) :: logfile,ppdfile,idxfile,outfile
  integer :: brkfunit = 50
  integer :: log1 = 60
  integer :: ppdfunit,idxfunit
  integer :: nintervals
  integer :: totdays
  integer :: hlines
  integer :: i,j
  integer, allocatable :: dates(:,:)
  integer, allocatable :: gaugeindxs(:,:)
  integer, allocatable :: accuracy(:,:)
  real, allocatable :: precipitation(:,:)
!
  allocate(dates(2,3))
  i = index(brkfile,'.',.TRUE.)
  logfile = trim(adjustl(brkfile(1:i-1)))//'.log'
  ppdfile = trim(adjustl(brkfile(1:i-1)))//'.ppd'
  ppdfunit = brkfunit + 1
  idxfile = trim(adjustl(brkfile(1:i-1)))//'.idx'
  idxfunit = brkfunit + 2
  outfile = trim(adjustl(brkfile(1:i-1)))//'.out'
  open(log1,file=logfile,status='replace')
  inquire(file=trim(adjustl(brkfile)), exist=brkexists)
  if (.not.brkexists) then
    wrtline = 'ERR: breakpoint -- data file '//trim(adjustl(brkfile))//' does not exist'
    print *, wrtline
    write(log1,'(A80)') wrtline
        stop
  end if
!
  if (mod(1440,intmins).ne.0) then
    wrtline = 'ERR: breakpoint -- output interval must evenly divide one day (1440 minutes)'
    print *, wrtline
    write(log1,'(A80)') wrtline
    stop
  end if
!
  if (outfmt.gt.2) then
    wrtline = 'ERR: breakpoint -- choose a supported output format'
    print *, wrtline
    write(log1,'(A80)') wrtline
    wrtline = '                   1 = single file including all stations'
    print *, wrtline
    write(log1,'(A80)') wrtline
    wrtline = '                   2 = individual file for each station'
    print *, wrtline
    write(log1,'(A80)') wrtline
    stop
  end if
!
  i = index(datebeg,'/',.FALSE.)
  j = index(datebeg,'/',.TRUE.)
  read(datebeg(1:i-1),'(I2)') dates(1,2) !recbegmm
  read(datebeg(i+1:j-1),'(I2)') dates(1,3) !recbegdd
  read(datebeg(j+1:len_trim(adjustl(datebeg))),'(I4)') dates(1,1) !recbegyy
  wrtline = 'Begin date: '//trim(adjustl(datebeg))
!  print *, wrtline
  write(log1,'(A80)') wrtline
  i = index(dateend,'/',.FALSE.)
  j = index(dateend,'/',.TRUE.)
  read(dateend(1:i-1),'(I2)') dates(2,2) !recendmm
  read(dateend(i+1:j-1),'(I2)') dates(2,3) !recenddd
  read(dateend(j+1:len_trim(adjustl(dateend))),'(I4)') dates(2,1) !recendyy
  wrtline = 'End date: '//trim(adjustl(dateend))
!  print *, wrtline
  write(log1,'(A80)') wrtline
!
  totdays = days(dates)
  nintervals = totdays * 1440 / intmins
  write(wrtline,'(I6)') nintervals
  wrtline = trim(adjustl(wrtline))//' intervals requested at'
!  print *, wrtline
  write(log1,'(A80)') wrtline
  write(wrtline,'(I6)') intmins
  wrtline = trim(adjustl(wrtline))//' minutes each over'
!  print *, wrtline
  write(log1,'(A80)') wrtline
  write(wrtline,'(I6)') totdays
  wrtline = trim(adjustl(wrtline))//' total days for'
!  print *, wrtline
  write(log1,'(A80)') wrtline
  write(wrtline,'(I6)') ngauges
  wrtline = trim(adjustl(wrtline))//' precipitation gauges'
!  print *, wrtline
  write(log1,'(A80)') wrtline
!
  allocate(gaugeindxs(ngauges,3))
  gaugeindxs = 0
!
  ppdexists = checkintfile(log1,brkfile,ppdfile)
  idxexists = checkintfile(log1,brkfile,idxfile)
  if (ppdexists .and. idxexists) then
    intfiles = .TRUE.
    wrtline = 'MSG: breakpoint -- using pre-processed data file: '&
              //trim(adjustl(ppdfile))
    print *, wrtline
    write(log1,'(A80)') wrtline
    wrtline = 'MSG: breakpoint -- using pre-processed index file: '&
              //trim(adjustl(idxfile))
    print *, wrtline
    write(log1,'(A80)') wrtline
  else
    intfiles = .FALSE.
    wrtline = 'MSG: breakpoint -- processing data file: '//trim(adjustl(brkfile))
    print *, wrtline
    write(log1,'(A80)') wrtline
    call readheader(log1,brkfunit,brkfile,intfiles,gaugeindxs,hlines)
    call procstndata(log1,brkfile,hlines,ngauges,gaugeindxs)
    intfiles = .TRUE.
  end if
!
  call readheader(log1,ppdfunit,ppdfile,intfiles,gaugeindxs,hlines)
  allocate(precipitation(ngauges,nintervals))
  precipitation = 0.0
  allocate(accuracy(ngauges,2))
  accuracy = 0
  call breakagg(log1,ppdfunit,ppdfile,hlines,ngauges,nintervals,intmins,dates,&
                gaugeindxs,precipitation,accuracy)
  call metric(log1,ngauges,gaugeindxs,accuracy)
  call writestnrec(log1,ngauges,nintervals,intmins,outfmt,dates,&
                   gaugeindxs,precipitation,outfile)
!
  wrtline = 'MSG: breakpoint -- precipitation disaggregation completed'
  print *, wrtline
  write(log1,'(A80)') wrtline
  close(log1)
!
  return
end subroutine brkprecip
!
!
subroutine intindex(dates,tcurr,intm,idx,dm)
!
! argument variables
  integer, intent(INOUT) :: dates(:,:)
  integer, intent(IN) :: tcurr(:)
  integer, intent(IN) :: intm
  integer, intent(OUT) :: idx,dm
!
! local variables
  integer :: elapsdays,elapsints,dints
!
  dates(2,1) = tcurr(1)
  dates(2,2) = tcurr(2)
  dates(2,3) = tcurr(3)
  elapsdays = days(dates) - 1
  elapsints = elapsdays * 1440 / intm
  dm = tcurr(4) * 60 + tcurr(5) ! time of day in minutes
  dints = int(dm/intm) + 1 ! intervals as a function of time of day
  idx = elapsints + dints
  dm = mod(dm,intm)
!
  return
end subroutine intindex
!
!
subroutine readheader(log1,funit,datafile,intf,gidxs,hlines)
!
! argument variables
  integer, intent(IN) :: log1,funit
  character*80, intent(IN) :: datafile
  logical, intent(IN) :: intf
  integer, intent (INOUT) :: gidxs(:,:)
  integer, intent(OUT) :: hlines
!
! local variables
  character*80 :: wrtline
  character*120 :: hline
  character*80 :: dataf2,dataf3
  integer :: ng
  integer :: ppdfunit,idxfunit
  integer :: i,j,k,m,n
!
  wrtline = ' '
  write(log1,'(A80)') wrtline
  wrtline = '<<<<< begin data file header information >>>>>'
  write(log1,'(A80)') wrtline
  open(funit,file=datafile,status='old')
  if (.not.intf) then
    ppdfunit = funit + 1
    idxfunit = funit + 2
        i = index(datafile,'.',.TRUE.)
    dataf2 =  trim(adjustl(datafile(1:i-1)))//'.ppd'
    dataf3 =  trim(adjustl(datafile(1:i-1)))//'.idx'
    open(ppdfunit,file=dataf2,status='replace')
    open(idxfunit,file=dataf3,status='replace')
  else
    idxfunit = funit + 1
        i = index(datafile,'.',.TRUE.)
    dataf3 =  trim(adjustl(datafile(1:i-1)))//'.idx'
    open(idxfunit,file=dataf3,status='old')    
  end if
  hlines = 0
  ng = 0
  do
    read(funit,'(A120)') hline
    if (hline(1:1).ne.'#') then
          exit ! header lines are done
    else
      hlines = hlines + 1
	  write(log1,'(A120)') hline
	  if (.not.intf) then
	    write(ppdfunit,'(A120)') hline
	    write(idxfunit,'(A120)') hline
	  end if
	end if
    if (hline(1:6).eq.'#Gages' .or. hline(1:6).eq.'#gages') then
      j = index(hline,':',.FALSE.)
      read(hline(j+2:len_trim(adjustl(hline))),'(A120)') hline
      do
        j = index(hline,'-',.FALSE.)
        k = index(hline,',',.FALSE.)
        if (k.eq.0) then    ! last gauge or series of gauges in line
                  if (j.eq.0) then
                    exit            ! end of line
                  else
                    k = j + 4       
                  end if
                end if
        if (j.lt.k) then    ! line ends with series of gauges, e.g. 'Xx-Yy'
                  if (j.eq.0) then
            read(hline(j+1:k-1),'(I3)') m
                  else
            read(hline(1:j-1),'(I3)') m
                  end if
              read(hline(j+1:k-1),'(I3)') n
              do i = m,n
                ng = ng + 1
            gidxs(ng,1) = i
              end do
          read(hline(k+1:len_trim(adjustl(hline))),'(A120)') hline
        else                ! line ends with single gauge, e.g. '...,Xx'
              read(hline(1:k-1),'(I3)') n
              ng = ng + 1
          gidxs(ng,1) = n
          read(hline(k+1:len_trim(adjustl(hline))),'(A120)') hline
        end if
            j = 0
            k = 0
      end do
      if (hline.ne.' ') then
        read(hline(1:len_trim(adjustl(hline))),'(I3)') n
            ng = ng + 1
        gidxs(ng,1) = n
      end if
        end if
  end do
  wrtline = '<<<<<< end data file header information >>>>>>'
  write(log1,'(A80)') wrtline
  wrtline = ' '
  write(log1,'(A80)') wrtline
!
  close(funit)
  if (intf) close(ppdfunit)
  if (intf) close(idxfunit)
!
  return
end subroutine readheader
!
!
subroutine procstndata(log1,datafile,hlines,ng,gidxs)
!
! argument variables
  integer, intent(IN) :: log1,hlines,ng
  character*80, intent(IN) :: datafile
  integer, intent(INOUT) :: gidxs(:,:)
!
! local variables
  character*80 :: dataf2,dataf3
  character*80 :: wrtline
  character*80 :: hline
  character*40 :: tprec,next_tprec
  character*7 :: nrecs
  integer :: fileend = 0
  integer :: gnum,month,day,year,hour,minute,offset,eoeflag
  integer :: next_offset
  integer :: funit,ppdfunit,idxfunit
  integer :: i,j,k,m
  integer, dimension(12) :: mmdays = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
  real :: strmp,prate
!
  funit = 50
  open(funit,file=datafile,status='old')
  ppdfunit = funit + 1
  i = index(datafile,'.',.TRUE.)
  dataf2 =  trim(adjustl(datafile(1:i-1)))//'.ppd'
  open(ppdfunit,file=dataf2,status='replace')
  idxfunit = funit + 2
  dataf3 =  trim(adjustl(datafile(1:i-1)))//'.idx'
  open(idxfunit,file=dataf3,status='replace')    
  wrtline = 'MSG: breakpoint -- reading station, time and breakpoint data records'
!  print *, wrtline
  write(log1,'(A80)') wrtline
  do i = 1,hlines
    read(funit,'(A80)') hline
  end do
  m = 0
  read(funit,'(A40)',iostat=fileend) tprec
  do 
    read(funit,'(A40)',iostat=fileend) next_tprec 
        if (fileend.ge.0 .and. trim(adjustl(next_tprec)).eq.' ') cycle ! blank line encountered
        m = m + 1
! extract gauge number
    i = index(tprec,',',.FALSE.)
    read(tprec(1:i-1),'(I3)') gnum
        tprec = tprec(i+1:len_trim(adjustl(tprec)))
! extract date (mm/dd/yyyy)
    i = index(tprec,',',.FALSE.)
    j = index(tprec,'/',.FALSE.)
    k = index(tprec,'/',.TRUE.)
    read(tprec(1:j-1),'(I2)') month
        read(tprec(j+1:k-1),'(I2)') day
        read(tprec(k+1:i-1),'(I4)') year
        tprec = tprec(i+1:len_trim(adjustl(tprec)))
! extract time (hh:mm)
    i = index(tprec,',',.FALSE.)
    j = index(tprec,':',.FALSE.)
    read(tprec(1:j-1),'(I2)') hour
        read(tprec(j+1:i-1),'(I2)') minute
        tprec = tprec(i+1:len_trim(adjustl(tprec)))
! extract record offset (minutes)
    i = index(tprec,',',.FALSE.)
        read(tprec(1:i-1),'(I3)') offset
        tprec = tprec(i+1:len_trim(adjustl(tprec)))
! extract storm cumulative precipitation (units [L])
    i = index(tprec,',',.FALSE.)
        read(tprec(1:i-1),'(F6.2)') strmp
        tprec = tprec(i+1:len_trim(adjustl(tprec)))
! extract precipitation rate (units [L/T])
        read(tprec(1:len_trim(adjustl(tprec))),'(F7.3)') prate
! adjust record time for offset
    minute = minute + offset
        do 
          if (minute.ge.60) then
            minute = minute - 60
        hour = hour + 1
          else
            exit
          end if
          if (hour.ge.24) then
                hour = hour - 24
        day = day + 1
          end if
          if (leapyr(year,month)) then
        if (day.gt.mmdays(month)+1) then
          day = day - mmdays(month) - 1
                  month = month + 1
                end if  
        else if (day.gt.mmdays(month)) then
                day = day - mmdays(month)
        month = month + 1
      end if
          if (month.gt.12) then
        month = month - 12
        year = year + 1
          end if
    end do
        eoeflag = 0
! evaluate next breakpoint record to identify end of event
    if (prate.eq.0) then
          tprec = next_tprec
! remove gauge number
      i = index(tprec,',',.FALSE.)
          tprec = tprec(i+1:len_trim(adjustl(tprec)))
! remove date (mm/dd/yyyy)
      i = index(tprec,',',.FALSE.)
          tprec = tprec(i+1:len_trim(adjustl(tprec)))
! remove time (hh:mm)
      i = index(tprec,',',.FALSE.)
          tprec = tprec(i+1:len_trim(adjustl(tprec)))
! extract record offset (minutes)
      i = index(tprec,',',.FALSE.)
          read(tprec(1:i-1),'(I3)') next_offset
      if (next_offset.eq.0) eoeflag = 1
        end if
! write line data to intermediate processing file
    do i = 1,ng
      if (gnum.eq.gidxs(i,1)) then
                gidxs(i,2) = gidxs(i,2) + 1 ! increment gauge-record counter
        if (eoeflag.eq.1) then
                  gidxs(i,3) = gidxs(i,3) + 1 ! increment gauge-event counter
                end if
        write(ppdfunit,'(I3,1X,I4,1X,4(I2,1X),F6.2,1X,F7.3,1X,I1)') &
                  gnum,year,month,day,hour,minute,strmp,prate,eoeflag
                exit
          end if
        end do
        if (fileend.lt.0) exit ! end of file
        tprec = next_tprec
  end do
! write gauge summary values to intermediate processing file
  do i = 1,ng
    write(idxfunit,'(I3,1X,I5,1X,I5)') gidxs(i,1),gidxs(i,2),gidxs(i,3)
  end do
!
  write(nrecs,'(I7)') m
  wrtline = 'MSG: breakpoint -- done reading '//trim(adjustl(nrecs))//' breakpoint records'
!  print *, wrtline
  write(log1,'(A80)') wrtline
  wrtline = 'MSG: breakpoint -- closing data and intermediate processing files'
!  print *, wrtline
  write(log1,'(A80)') wrtline
  close(funit)
  close(ppdfunit)
  close(idxfunit)
!
  return
end subroutine procstndata
!
!
subroutine breakagg(log1,funit,datafile,hlines,ng,nint,intmins,dates,gidxs,precip,accy)
!
! argument variables
  integer, intent(IN) :: log1,funit
  character*80, intent(IN) :: datafile
  integer, intent(IN) :: hlines,ng,nint
  integer, intent(INOUT) :: intmins
  integer, intent(INOUT) :: dates(:,:)
  integer, intent(INOUT) :: gidxs(:,:)
  real, intent(OUT) :: precip(:,:)
  integer, intent(OUT) :: accy(:,:)
!
! local variables
  character*80 :: wrtline
  character*80 :: hline
  character*80 :: idxdataf
  character*3 :: gage
  character*6 :: brecs,events
  integer, dimension(2,5) :: tdata = 0 
  real, dimension(2,2) :: pdata = 0.0
  real :: errtol = 0.3 ! according to David Goodrich, ARS-Tucson
  real :: ppart,lwrmarg,uprmarg,calcptot
  integer :: gnum,dmins1,dmins2,mins
  integer :: eoeflag
  integer :: idxfunit
  integer :: fileend = 0
  integer :: i,j,k,l,m,n,mm
!
! Array tdata content description
!   tdata(:,1) = year
!   tdata(:,2) = month
!   tdata(:,3) = day
!   tdata(:,4) = hour
!   tdata(:,5) = minute
!
! Array pdata content description
!   pdata(:,1) = strmp
!   pdata(:,2) = prate
!   
  wrtline = 'MSG: breakpoint -- reading intermediate index file'
!  print *, wrtline
  write(log1,'(A80)') wrtline
  idxfunit = funit + 1
  i = index(datafile,'.',.TRUE.)
  idxdataf =  trim(adjustl(datafile(1:i-1)))//'.idx'
  open(idxfunit,file=idxdataf,status='old')
  do i = 1,hlines
    read(idxfunit,'(A80)') hline
  end do
  do i = 1,ng
    read(idxfunit,'(I3,1X,I5,1X,I5)') gidxs(i,1),gidxs(i,2),gidxs(i,3)
  end do
  close(idxfunit)
!
  wrtline = 'MSG: breakpoint -- beginning disaggregation to regular intervals'
!  print *, wrtline
  write(log1,'(A80)') wrtline
  do i = 1,ng
    wrtline = ' '
    write(log1,'(A80)') wrtline
    write(gage,'(I3)') gidxs(i,1)
    write(brecs,'(I6)') gidxs(i,2)
    write(events,'(I6)') gidxs(i,3)
    wrtline = 'MSG: breakpoint -- processing gauge no. '//trim(adjustl(gage))//&
              ': '//trim(adjustl(events))//' events, '//&
              trim(adjustl(brecs))//' breakpoint records'
!    print *, wrtline
    write(log1,'(A80)') wrtline
    open(funit,file=datafile,status='old')
    do j = 1,hlines
      read(funit,'(A80)') hline
    end do
!
    calcptot = 0.0
    j = 1 ! breakpoint record counter
       l = 1 ! gauge-event counter
    do
       if (j.gt.gidxs(i,2)) exit ! no need to go further in file
          if (l.gt.gidxs(i,3)) exit ! no need to go further in file
      read(funit,'(I3,1X,I4,1X,4(I2,1X),F6.2,1X,F7.3,1X,I1)',iostat=fileend) &
           gnum,tdata(2,1),tdata(2,2),tdata(2,3),tdata(2,4),tdata(2,5),&
           pdata(2,1),pdata(2,2),eoeflag
      if (fileend.lt.0) exit ! end of file
! pick out only breakpoint records corresponding to the current gauge
      if (gnum.eq.gidxs(i,1)) then
! don't need to do any processing on the first breakpoint in a new event
        if (pdata(2,1).gt.0) then 
! find interval index m at previous breakpoint of precipitation event
          if (j.eq.1) then ! the first breakpoint in the gauge record
             dmins1 = 0
          else
             call intindex(dates,tdata(1,:),intmins,m,dmins1)
          end if
! find interval index n at current breakpoint of precipitation event
          call intindex(dates,tdata(2,:),intmins,n,dmins2)
          if (j.eq.1) then ! the first breakpoint in the gauge record
             m = n - 1
          end if
! interval precipitation calculations
          if (n.gt.m) then ! split the breakpoint over two or more intervals
             mins = intmins - dmins1
             ppart = mins * pdata(1,2) / 60
             precip(i,m) = precip(i,m) + ppart
             calcptot = calcptot + ppart
             write(wrtline,'(A3,6I6,4F8.3)') ' 1:',l,j,m,dmins1,intmins,mins,pdata(1,2),ppart,precip(i,m),calcptot
            write(log1,'(A80)') wrtline
            dmins1 = 0
            mm = m + 1
            do
               if (mm.eq.n) then
                  mins = dmins2 - dmins1
                  ppart = mins * pdata(1,2) / 60
                  precip(i,mm) = precip(i,mm) + ppart
                  calcptot = calcptot + ppart
                write(wrtline,'(A3,6I6,4F8.3)') ' 3:',l,j,mm,dmins1,dmins2,mins,pdata(1,2),ppart,precip(i,mm),calcptot
               else
                  mins = intmins - dmins1
                  ppart = mins * pdata(1,2) / 60
                  precip(i,mm) = precip(i,mm) + ppart
                  calcptot = calcptot + ppart
                  write(wrtline,'(A3,6I6,4F8.3)') ' 2:',l,j,mm,dmins1,intmins,mins,pdata(1,2),ppart,precip(i,mm),calcptot
               end if
              write(log1,'(A80)') wrtline
                    mm = mm + 1
              if (mm.gt.n) exit
              end do
              else if (n.eq.m) then ! current breakpoint occurs in same interval as previous
                 mins = dmins2 - dmins1
                 ppart = mins * pdata(1,2) / 60
                 precip(i,n) = precip(i,n) + ppart
                 calcptot = calcptot + ppart
                 write(wrtline,'(A3,6I6,4F8.3)') ' 4:',l,j,n,dmins1,dmins2,mins,pdata(1,2),ppart,precip(i,n),calcptot
            write(log1,'(A80)') wrtline
          else ! (n.lt.m) ERROR in breakpoint sorting or calculation of indices
		    wrtline = 'ERROR in calculation of interval indices -- stopping'
			print *,wrtline
			write(log1,'(A80)') wrtline
			write(wrtline,'(A5,I3,A7,I5)') 'gauge',gidxs(i,1),'  event',l
			print *,wrtline
			write(log1,'(A80)') wrtline
			write(wrtline,'(A5,I4,4(1X,I2),A4,I6)') 'date ',&
			      tdata(1,1),tdata(1,2),tdata(1,3),tdata(1,4),tdata(1,5),', m:',m
			print *,wrtline
            write(log1,'(A80)') wrtline
			write(wrtline,'(A5,I4,4(1X,I2),A4,I6)') 'date ',&
			      tdata(2,1),tdata(2,2),tdata(2,3),tdata(2,4),tdata(2,5),', n:',n
			print *,wrtline
            write(log1,'(A80)') wrtline
			stop
	      end if
		end if
        j = j + 1 ! increment breakpoint record counter
        if (eoeflag.eq.1) then ! breakpoint record at end of event
           lwrmarg = pdata(2,1) - errtol ! error tolerance, lower bound
           uprmarg = pdata(2,1) + errtol ! error tolerance, upper bound
           if (lwrmarg.le.calcptot.and.calcptot.le.uprmarg) then ! within error tolerances
            accy(i,1) = accy(i,1) + 1 ! increment possible score
            accy(i,2) = accy(i,2) + 1 ! increment actual score
            write(wrtline,'(A10,I3,A11,F6.2,A13,F6.2,A5)') ' -- Event ',l,&
                          '  recorded:',pdata(2,1),'  calculated:',calcptot,&
                          '  HIT'
           else
            accy(i,1) = accy(i,1) + 1 ! increment possible score
            write(wrtline,'(A10,I3,A11,F6.2,A13,F6.2,A6)') ' -- Event ',l,&
                          '  recorded:',pdata(2,1),'  calculated:',calcptot,&
                          '  MISS'
           end if
          write(log1,'(A80)') wrtline
          l = l + 1 ! increment gauge-event counter
          calcptot = 0.0
        end if
! move tdata/pdata values in row 2 (current breakpoint) to row 1 (previous breakpoint)
        tdata(1,:) = tdata(2,:)
        pdata(1,:) = pdata(2,:)
      end if
      end do ! j = 1,gidxs(i,2) or l = 1,gidxs(i,3) by internal counters
    write(wrtline,'(A17,I3,A7,I3,A7)') ' Gauge summary:  ',accy(i,2),&
                  ' hits  ',accy(i,1),' events'
    write(log1,'(A80)') wrtline
    close(funit)
  end do ! i = 1,ng by loop-end counter
!
  return
end subroutine breakagg
!
!
subroutine metric(log1,ng,gidxs,accy)
!
! argument variables
  integer, intent(IN) :: log1,ng
  integer, intent(IN) :: gidxs(:,:)
  integer, intent(IN) :: accy(:,:)
!
! local variables
  character*3 :: gage,hits
  character*6 :: evnts,acc
  character*80 :: wrtline
  real :: accr
  integer :: totgevents = 0
  integer :: totcorr = 0
  integer :: i
!
  wrtline = ' '
  write(log1,'(A80)') wrtline
  wrtline = 'MSG: breakpoint -- quality control (QC) processing on gauge-event basis'
!  print *, wrtline
  write(log1,'(A80)') wrtline
!
  do i = 1,ng
	write(gage,'(I3)') gidxs(i,1)
	write(evnts,'(I6)') gidxs(i,3)
	write(hits,'(I3)') accy(i,2)
	if (accy(i,1).gt.0) then
	  accr = (real(accy(i,2)) / real(accy(i,1))) * 100.0
	  totgevents = totgevents + accy(i,1)
	  totcorr = totcorr + accy(i,2)
	else
	  accr = 0.0
	end if
	write(acc,'(F6.2)') accr
	wrtline = 'QC -- gauge '//trim(adjustl(gage))//': '//&
		            trim(adjustl(hits))//' hits, '//&
		            trim(adjustl(evnts))//' events, '//&
					trim(adjustl(acc))//'% accuracy'
	write(log1,'(A80)') wrtline

  end do
  wrtline = ' '
  write(log1,'(A80)') wrtline
  write(evnts,'(I6)') totgevents
  accr = (real(totcorr) / real(totgevents)) * 100.0
  write(acc,'(F6.2)') accr
  wrtline = 'MSG: breakpoint -- QC overall results: '//trim(adjustl(evnts))//&
            ' gauge-events, '//trim(adjustl(acc))//'% accuracy'
  print *,wrtline
  write(log1,'(A80)') wrtline
  wrtline = ' '
  write(log1,'(A80)') wrtline
!
  return
end subroutine metric
!
!
subroutine writestnrec(log1,ng,nint,intmins,outfmt,dates,gidxs,precip,outfile)
!
! argument variables
  integer, intent(IN) :: log1,ng,nint,intmins,outfmt
  character*80, intent(INOUT) :: outfile
  integer, intent(INOUT) :: dates(:,:)
  integer, intent(IN) :: gidxs(:,:)
  real, intent(IN) :: precip(:,:)
!
! local variables
  character*80 :: wrtline,datafile
  character*16 :: datetime ! MM/DD/YYYY hh:mm 
  character*3 :: gnum
  integer :: funit = 10
  integer :: i,j
!
  select case (outfmt)
    case (1)
       wrtline = 'MSG: breakpoint -- writing regular interval records to single output file'
!      print *, wrtline
      write(log1,'(A80)') wrtline
      open(funit,file=outfile,status='replace')
           datetime = 'MM/DD/YYYY hh:mm'
      write(funit,'(A16)', advance='NO') datetime
      do j = 1,ng-1
        write(funit,'(I7)', advance='NO') gidxs(j,1)
      end do
      write(funit,'(I7)', advance='YES') gidxs(ng,1)
      do i = 1,nint
         datetime = datestr(i,intmins,dates)
        write(funit,'(A16)', advance='NO') datetime
        do j = 1,ng-1
          write(funit,'(F7.2)', advance='NO') precip(j,i)
        end do
        write(funit,'(F7.2)', advance='YES') precip(ng,i)
        end do
    case (2)
       wrtline = 'MSG:  breakpoint -- writing regular interval records to individual output files'
!      print *, wrtline
      write(log1,'(A80)') wrtline
      do i = 1,ng
         select case (gidxs(i,1))
          case (100:)
             write(gnum,'(I3)') gidxs(i,1)
          case (10:99)
             write(gnum,'(I2)') gidxs(i,1)
             write(gnum,'(A3)') '0'//trim(adjustl(gnum))
          case (:9)
             write(gnum,'(I1)') gidxs(i,1)
             write(gnum,'(A3)') '00'//trim(adjustl(gnum))
          end select
          datafile = outfile
        j = index(datafile,'.',.TRUE.)
        outfile = trim(adjustl(datafile(1:j-1)))//&
                       'WG'//trim(adjustl(gnum))//'.out'
        open(funit,file=outfile,status='replace')
        write(funit,'(A12)') outfile
        do j = 1,nint
           datetime = datestr(j,intmins,dates)
           write(funit,'(A16,F7.2)') datetime,precip(i,j)
        end do
        close(funit)
        write(wrtline,'(I3)') gidxs(i,1)
        wrtline = 'MSG: breakpoint -- wrote record file for gauge '//trim(adjustl(wrtline))
!	    print *, wrtline
		write(log1,'(A80)') wrtline
      end do
  end select
!
  return
end subroutine writestnrec
!
!
! end internal subroutines
!
!
end module breakpoint_module
!end program Breakpt
