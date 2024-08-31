!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
subroutine readLISlsmSMANNdata(n, iomode, p_s, p_e)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
#if (defined USE_GRIBAPI)
  use grib_api
#endif
  use LDT_coreMod,      only : LDT_rc
  use LDT_ANNMod
  use LDT_historyMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_timeMgrMod
  use LISlsmSM_ANNdataMod

  implicit none

  integer,   intent(in) :: n
  integer,   intent(in) :: iomode
  integer,   intent(in) :: p_s
  integer,   intent(in) :: p_e

  character(len=LDT_CONST_PATH_LEN) :: fname 
  logical          :: file_exists
  real             :: sm_data(LDT_rc%lnc(n), LDT_rc%lnr(n))
  real             :: gvf_data(LDT_rc%lnc(n), LDT_rc%lnr(n))
  real             :: swe_data(LDT_rc%lnc(n), LDT_rc%lnr(n))
  real             :: st_data(LDT_rc%lnc(n), LDT_rc%lnr(n))
  real             :: snowf_data(LDT_rc%lnc(n), LDT_rc%lnr(n))
  real             :: pcp_data(LDT_rc%lnc(n), LDT_rc%lnr(n))

  character*50     :: sm_units, pcp_units, snowf_units, st_units
  character*50     :: gvf_units, swe_units

  integer          :: t, index
  integer          :: ftn
  integer          :: iret
  real             :: topLev,botLev,param_num,trange
  integer          :: igrib,nvars
  real             :: smvalue(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  integer          :: c,r
  character*3      :: fnest 
  logical          :: alarmCheck

#if (defined USE_GRIBAPI)
  write(fnest,'(i3.3)') n
  alarmCheck = LDT_isAlarmRinging(LDT_rc,&
       "LIS LSM SM alarm "//trim(fnest))
  
  if(alarmCheck) then

     sm_data = LDT_rc%udef
     gvf_data = LDT_rc%udef
     swe_data = LDT_rc%udef
     st_data = LDT_rc%udef
     pcp_data = LDT_rc%udef

     call create_ANNlsm_output_filename(lsmsmANNdata%nest, lsmsmANNdata%format,&
          fname, lsmsmANNdata%odir, lsmsmANNdata%wstyle, &
          lsmsmANNdata%wopt)
     
     inquire(file=trim(fname),exist=file_exists)
     
     if(file_exists) then 
        write(LDT_logunit,*) 'reading LSM output ',trim(fname)
        if(lsmsmANNdata%format.eq."binary") then 
           write(LDT_logunit,*) 'DA preprocessing on the binary format is not '
           write(LDT_logunit,*) 'currently supported. Program stopping....'
           call LDT_endrun()
           
        elseif(lsmsmANNdata%format.eq."grib1") then 
           if(lsmsmANNdata%wstyle.ne."WMO convention") then 
              write(LDT_logunit,*) 'LDT currently does not support this style of grib output'
              call LDT_endrun
           endif
           
           call grib_open_file(ftn,trim(fname),'r',iret)
           if(iret.ne.0) then 
              write(LDT_logunit,*) 'Could not open file: ',trim(fname)
              call LDT_endrun()
           endif
           call grib_multi_support_on
           
           do 
              call grib_new_from_file(ftn,igrib,iret)
              if(iret==GRIB_END_OF_FILE) then 
                 exit
              endif
              
              call grib_get(igrib, 'indicatorOfParameter',param_num, iret)
              call LDT_verify(iret, &
                   'grib_get: indicatorOfParameter failed in readLISlsmsmANNdata')
              
              call grib_get(igrib, 'topLevel',topLev, iret)
              call LDT_verify(iret, &
                   'grib_get: topLevel failed in readLISlsmsmANNdata')
              
              call grib_get(igrib, 'bottomLevel',botLev, iret)
              call LDT_verify(iret, &
                   'grib_get: bottomLevel failed in readLISlsmsmANNdata')
              
              call grib_get(igrib, 'timeRangeIndicator',trange, iret)
              call LDT_verify(iret, &
                   'grib_get: timeRangeIndicator failed in readLISlsmsmANNdata')
              
              !right now specifically geared for AFWA outputs.            
              if(param_num.eq.201.and.topLev.eq.0.and.botlev.eq.10.and.&
                   trange.eq.1) then 
                 
                 call grib_get(igrib,'values',smvalue,iret)
                 call LDT_verify(iret,&
                      'grib_get: values failed in readLISlsmSMANNdata')
                 
                 do r=1,LDT_rc%lnr(n)
                    do c=1,LDT_rc%lnc(n)
                       if(smvalue(c+(r-1)*LDT_rc%lnc(n)).gt.0.and.&
                            smvalue(c+(r-1)*LDT_rc%lnc(n)).le.0.5) then 
                          sm_data(c,r) = smvalue(c+(r-1)*LDT_rc%lnc(n))
                       else
                          sm_data(c,r) = LDT_rc%udef
                       endif
                    enddo
                 enddo
              endif
              
              call grib_release(igrib,iret)
              call LDT_verify(iret, 'error in grib_release in readLISlsmsmANNdata')
           enddo
           
           call grib_close_file(ftn)
           
        elseif(lsmsmANNdata%format.eq."netcdf") then 
#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
           
           iret = nf90_open(path=trim(fname),mode=nf90_nowrite, ncid=ftn)
           call LDT_verify(iret, 'Error opening file '//trim(fname))
           
           call LDT_readLISSingleNetcdfVar(n,ftn, "SoilMoist_tavg",&
                1,sm_data, units=sm_units)

           call LDT_readLISSingleNetcdfVar(n,ftn, "Rainf_tavg",&
                1,pcp_data, units=pcp_units)

           call LDT_readLISSingleNetcdfVar(n,ftn, "Snowf_tavg",&
                1,snowf_data, units=snowf_units)

           call LDT_readLISSingleNetcdfVar(n,ftn, "Greenness_inst",&
                1,gvf_data, units=gvf_units)

           call LDT_readLISSingleNetcdfVar(n,ftn, "SoilTemp_tavg",&
                1,st_data, units=st_units)

           call LDT_readLISSingleNetcdfVar(n,ftn, "SWE_tavg",&
                1,swe_data, units=swe_units)

           iret = nf90_close(ftn)
           call LDT_verify(iret,'Error in nf90_close')
#endif
        endif
     else
        write(LDT_logunit,*) 'LIS file '//trim(fname)
        write(LDT_logunit,*) 'not found. Program stopping....'
        call LDT_endrun()
     endif
     
     call LDT_logSingleANNdata(n,&
          sm_data,  &
          pindex=p_s, &
          iomode = iomode, &
          name = "SoilMoist",&
          units = sm_units)

     call LDT_logSingleANNdata(n,&
          pcp_data,  &
          pindex=p_s+1, &
          iomode = iomode, &
          name = "Rainf",&
          units=pcp_units)

     call LDT_logSingleANNdata(n,&
          snowf_data,  &
          pindex=p_s+2, &
          iomode = iomode, &
          name = "Snowf",&
          units=snowf_units)

     call LDT_logSingleANNdata(n,&
          gvf_data,  &
          pindex=p_s+3, &
          iomode = iomode,&
          name = "Greenness",&
          units=gvf_units)

     call LDT_logSingleANNdata(n,&
          st_data,  &
          pindex=p_s+4, &
          iomode = iomode,&
          name = "SoilTemp",&
          units=st_units)

     call LDT_logSingleANNdata(n,&
          swe_data,  &
          pindex=p_s+5, &
          iomode = iomode,&
          name = "SWE",&
          units=swe_units)
  endif
#endif
end subroutine readLISlsmSMANNdata
   

!BOP
!
! !ROUTINE: create_ANNlsm_output_filename
! \label{create_ANNlsm_output_filename}
!
! !INTERFACE:
subroutine create_ANNlsm_output_filename(n,form,fname, odir,&
     wstyle, wopt)
! !USES:
   use LDT_coreMod,  only : LDT_rc
   use LDT_logMod
   use LDT_constantsMod, only : LDT_CONST_PATH_LEN

   implicit none 
! !ARGUMENTS:
   integer,   intent(IN)   :: n 
   character(len=*)        :: fname
   character(len=*)        :: form
   character(len=*)        :: odir
   character(len=*)        :: wstyle
   character(len=*)        :: wopt
! 
! !DESCRIPTION:  
!  Create the file name for the output data files. It creates both the GSWP
!  style of output filenames and the standard LDT style. The convention used
!  in LDT creates a filename in a hierarchical style (output directory, 
!  model name, date, file extention)
!
!  2 level hierarchy
!   <output directory>/<model name>/LIS_HIST_<yyyymmddhhmnss>.<extension>
!  3 level hierarchy
!   <output directory>/<model name>/<yyyymm>/LIS_HIST_<yyyymmddhhmnss>.<extension>
!  4 level hierarchy
!   <output directory>/<model name>/<yyyy>/<yyyymm>/LIS_HIST_<yyyymmddhhmnss>.<extension>
!  WMO convention
!   <output directory>/<AFWA Weather product style>
!   A filename in the convention of weather products (such as):  
!   PS.AFWA\_SC.U\_DI.C\_DC.ANLYS\_GP.LDT\_GR.C0P25DEG\_AR.GLOBAL\_PA.03-HR-SUM\_DD.YYYYMMDD\_DT.HH00\_DF.GR1
!   where
!    PS = Product source
!    SC = security classification
!    DI = distribution classification
!    DC = data category 
!    GP = generating process
!    GR = grid 
!    AR = area of data
!    PA = parameter 
!    DD = date
!    DT = data time
!    DF = data format
! 
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest
!   \item [fname]
!     the created file name. 
!   \item [model\_name]
!    string describing the name of the model 
!   \item [writeint]
!    output writing interval  of the model
!   \item [style]
!    style option as described above
!  \end{description}
!EOP
   character(len=8)        :: date
   character(len=10)       :: time
   character(len=5)        :: zone
   integer, dimension(8)   :: values
   character(len=20)       :: mname 
   character(len=10)       :: cdate
   character(len=14)       :: cdate1
   character(len=2)        :: fint
   character(len=10)       :: fres
   character(len=10)       :: fres2
   character(len=10)       :: fres3
   character*1             :: fres1(10)
   character(len=1)        :: fproj
   integer                 :: curr_mo = 0
   character(len=LDT_CONST_PATH_LEN)       :: dname
   character(len=LDT_CONST_PATH_LEN), save :: out_fname
   integer                  :: i, c

   mname = 'SURFACEMODEL'
   if(wstyle.eq."4 level hierarchy") then 
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LDT_rc%yr, LDT_rc%mo, &
           LDT_rc%da, LDT_rc%hr,LDT_rc%mn
      
      dname = trim(odir)//'/'
      dname = trim(dname)//trim(mname)//'/'
      
      write(unit=cdate, fmt='(i4.4)') LDT_rc%yr
      dname = trim(dname)//trim(cdate)//'/'
      
      write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') LDT_rc%yr, LDT_rc%mo, LDT_rc%da
      dname = trim(dname)//trim(cdate)
      
      out_fname = trim(dname)//'/LIS_HIST_'//trim(cdate1)
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      select case ( form )
      case ( "binary" )
         if(wopt.eq."1d tilespace") then 
            out_fname = trim(out_fname)//'.ts4r'
         elseif(wopt.eq."2d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         elseif(wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         endif
      case ("grib1")
         out_fname = trim(out_fname)//'.grb'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.gr2'
      case default
         call ldt_log_msg('ERR: create_output_filename -- '// &
              'Unrecognized output format')
         call LDT_endrun 
      endselect
   elseif(wstyle.eq."3 level hierarchy") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LDT_rc%yr, LDT_rc%mo, &
           LDT_rc%da, LDT_rc%hr,LDT_rc%mn
      
      dname = trim(odir)//'/'
      dname = trim(dname)//trim(mname)//'/'
      
      write(unit=cdate, fmt='(i4.4, i2.2)') LDT_rc%yr, LDT_rc%mo
      dname = trim(dname)//trim(cdate)//'/'

      out_fname = trim(dname)//'LIS_HIST_'//trim(cdate1)
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      select case ( form )
      case ("binary")
         if(wopt.eq."1d tilespace") then 
            out_fname = trim(out_fname)//'.ts4r'
         elseif(wopt.eq."2d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         elseif(wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         endif
      case ("grib1")
         out_fname = trim(out_fname)//'.grb'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.gr2'
      case default
         call ldt_log_msg('ERR: create_output_filename -- '// &
              'Unrecognized form value')
         call LDT_endrun 
      endselect
   elseif(wstyle.eq."2 level hierarchy") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LDT_rc%yr, LDT_rc%mo, &
           LDT_rc%da, LDT_rc%hr,LDT_rc%mn
      
      dname = trim(odir)//'/'
      dname = trim(dname)//trim(mname)//'/'
      
      out_fname = trim(dname)//'LIS_HIST_'//trim(cdate1)
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      select case ( form )
      case ("binary")
         if(wopt.eq."1d tilespace") then 
            out_fname = trim(out_fname)//'.ts4r'
         elseif(wopt.eq."2d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         elseif(wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         endif
      case ("grib1")
         out_fname = trim(out_fname)//'.grb'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.gr2'
      case default
         call ldt_log_msg('ERR: create_output_filename -- '// &
              'Unrecognized form value')
         call LDT_endrun 
      endselect
   elseif(wstyle.eq."WMO convention") then 

      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2)') &
           LDT_rc%yr, LDT_rc%mo, LDT_rc%da
      
      write(unit=cdate, fmt='(i2.2, i2.2)') LDT_rc%hr, LDT_rc%mn
      
      out_fname = trim(odir)//'/'//trim(cdate1)//&
           '/PS.AFWA_SC.U_DI.C_DC.ANLYS_GP.LIS_GR.C0P25DEG_AR.GLOBAL_PA.03-HR-SUM_DD.'//&
           trim(cdate1)//'_DT.'//trim(cdate)//'_DF.GR1'
   endif
   fname = out_fname
 end subroutine create_ANNlsm_output_filename


