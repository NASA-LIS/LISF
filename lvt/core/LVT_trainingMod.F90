!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
#include "LVT_NetCDF_inc.h"
!BOP
! 
! !MODULE: LVT_trainingMod
! \label(LVT_trainingMod)
!
! !INTERFACE:
module LVT_trainingMod
! 
! !USES:   
  use LVT_trainingAlg_pluginMod
  use LVT_histDataMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   The code in this file provides interfaces to manage the operation
!   of different training algorithms for benchmarking. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  14 Jul 2015    Sujay Kumar  Initial Specification
! 
!EOP
!BOP
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_trainingInit          !initialize specified domains
  public :: LVT_runTraining
  public :: LVT_writeTrainingOutput
!EOP
  public :: LVT_trainingctl

  type, public ::  trainingctl
     character*10           :: algname
  end type trainingctl

  type(trainingctl)                   :: LVT_trainingctl

contains
!BOP
! 
! !ROUTINE: LVT_trainingInit
! \label{LVT_trainingInit}
!
! !INTERFACE: 
  subroutine LVT_trainingInit()
! 
! !USES:   
    use ESMF
    use LVT_coreMod
    use LVT_logMod
    use map_utils

    implicit none
!
! !DESCRIPTION: 
!  
! 
!EOP

    integer :: status

!read decision space file, number of generations, 
    call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%trainingAlg,&
         label="Training algorithm for benchmarking:", rc=status)
    call LVT_verify(status, 'Training algorithm for benchmarking: not defined')

    call LVT_trainingAlg_plugin  

    call traininginit(trim(LVT_rc%trainingAlg)//char(0))

  end subroutine LVT_trainingInit

!BOP
! 
! !ROUTINE: LVT_runTraining
!  \label{LVT_runTraining}
!
! !INTERFACE: 
  subroutine LVT_runTraining(pass)
! 
! !USES:   
    use LVT_coreMod
    use LVT_logMod

    implicit none    
!
! !INPUT PARAMETERS: 

    integer, intent(in)    :: pass
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine reads the training output for the specified iteration number. 
!  For the optimization algorithms (GA), the output file is read, and for the 
!  uncertainty estimation algorithms (MCMC, DEMC), the restart file is read. 
!
!  The arguments are: 
!  \begin{description}
!   \item[iterNo] generation number 
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 

!EOP

    call trainingrun(trim(LVT_rc%trainingAlg)//char(0), pass)
  end subroutine LVT_runTraining


!BOP
! 
! !ROUTINE: LVT_writeTrainingOutput
!  \label{LVT_writeTrainingOutput}
!
! !INTERFACE: 
  subroutine LVT_writeTrainingOutput(pass)
! 
! !USES:   
    use LVT_coreMod
    use LVT_logMod

    implicit none    
!
! !INPUT PARAMETERS: 

    integer, intent(in)    :: pass
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine reads the training output for the specified iteration number. 
!  For the optimization algorithms (GA), the output file is read, and for the 
!  uncertainty estimation algorithms (MCMC, DEMC), the restart file is read. 
!
!  The arguments are: 
!  \begin{description}
!   \item[iterNo] generation number 
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 

!EOP
    character*100           :: filename
    character(len=14)       :: cdate
    integer                 :: ftn
    integer                 :: varId
    integer                 :: shuffle, deflate, deflate_level
    character(len=8)        :: date
    character(len=10)       :: time
    character(len=5)        :: zone
    integer, dimension(8)   :: values
    integer                 :: dimID(3), dimID_t(2)
    integer                 :: tdimID
    character*8             :: xtime_begin_date
    character*6             :: xtime_begin_time
    character*50            :: xtime_units
    character*50            :: xtime_timeInc
    integer                 :: c,r,t,gid,iret,col,row
    real                    :: gtmp(LVT_rc%gnc,LVT_rc%gnr)

    type(LVT_metadataEntry), pointer :: ds2

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    if(pass.eq.2) then 
       if(LVT_rc%computeFlag) then 
          
          call LVT_getDataStream2Ptr(ds2)
          
          shuffle = NETCDF_shuffle
          deflate = NETCDF_deflate
          deflate_level =NETCDF_deflate_level
!---------------------------------------------------------------------------
! write outputs at every time averaging interval, only after the 
! training is completed. 
!---------------------------------------------------------------------------

          call system("mkdir -p "//trim(LVT_rc%statsodir)//"/TRAINING")
          write(unit=cdate,fmt='(i4.4,i2.2,i2.2,i2.2,i2.2,i2.2)') &
               LVT_rc%dyr(2), LVT_rc%dmo(2), LVT_rc%dda(2), LVT_rc%dhr(2), &
               LVT_rc%dmn(2), LVT_rc%dss(2)
          filename = trim(LVT_rc%statsodir)//"/TRAINING/"//&
               "LVT_HIST_OUT_"//trim(cdate)//".nc"
          

#if (defined USE_NETCDF4)
          iret = nf90_create(path=filename,cmode=nf90_hdf5,&
               ncid = ftn)
          call LVT_verify(iret,'creating netcdf file failed in LVT_trainingMod')
#endif          
#if (defined USE_NETCDF3)
          iret = nf90_create(path=filename,cmode=nf90_clobber,&
               ncid = ftn)
          call LVT_verify(iret,'creating netcdf file failed in LVT_trainingMod')
#endif      
          
!write lat/lon fields first

          if(allocated(LVT_histData%xlat%value)) then 
             deallocate(LVT_histData%xlat%value)
          endif
          if(allocated(LVT_histData%xlat%unittypes)) then 
             deallocate(LVT_histData%xlat%unittypes)
          endif
          if(allocated(LVT_histData%xlon%value)) then 
             deallocate(LVT_histData%xlon%value)
          endif
          if(allocated(LVT_histData%xlon%unittypes)) then 
             deallocate(LVT_histData%xlon%unittypes)
          endif
          call date_and_time(date,time,zone,values)
          LVT_histData%xlat%short_name = "latitude"
          LVT_histData%xlat%long_name = "latitude"
          LVT_histData%xlat%standard_name = "latitude"
          LVT_histData%xlat%units = "degree_north"
          LVT_histData%xlat%nunits = 1
          LVT_histData%xlat%format = 'F'
          LVT_histData%xlat%vlevels = 1
          LVT_histData%xlat%timeAvgOpt = 0 
          LVT_histData%xlat%selectNlevs = 1
          if(LVT_rc%lvt_wopt.eq."1d tilespace") then 
             allocate(LVT_histData%xlat%value(LVT_LIS_rc(1)%ntiles,&
                  1,&
                  LVT_histData%xlat%vlevels))
          else
             allocate(LVT_histData%xlat%value(LVT_rc%ngrid,&
                  1,&
                  LVT_histData%xlat%vlevels))
          endif
          allocate(LVT_histData%xlat%unittypes(1))
          LVT_histData%xlat%unittypes(1) = "degree_north"
          LVT_histData%xlon%short_name = "longitude"
          LVT_histData%xlon%long_name = "longitude"
          LVT_histData%xlon%standard_name = "longitude"
          LVT_histData%xlon%units = "degree_east"
          LVT_histData%xlon%nunits = 1
          LVT_histData%xlon%format = 'F'
          LVT_histData%xlon%vlevels = 1
          LVT_histData%xlon%timeAvgOpt = 0 
          LVT_histData%xlon%selectNlevs = 1
          if(LVT_rc%lvt_wopt.eq."1d tilespace") then 
             allocate(LVT_histData%xlon%value(LVT_LIS_rc(1)%ntiles,&
                  1,&
                  LVT_histData%xlat%vlevels))
          else
             allocate(LVT_histData%xlon%value(LVT_rc%ngrid,&
                  1,&
                  LVT_histData%xlon%vlevels))
          endif
          allocate(LVT_histData%xlon%unittypes(1))
          LVT_histData%xlon%unittypes(1) = "degree_north"
          
          if(LVT_rc%lvt_wopt.eq."1d tilespace") then 
             call LVT_verify(nf90_def_dim(ftn,'ntiles',LVT_LIS_rc(1)%ntiles,&
                  dimID(1)))
          elseif(LVT_rc%lvt_wopt.eq."2d gridspace") then 
             call LVT_verify(nf90_def_dim(ftn,'east_west',LVT_rc%gnc,dimID(1)))
             call LVT_verify(nf90_def_dim(ftn,'north_south',LVT_rc%gnr,dimID(2)))
          endif
!LVT output is writing output for a single time
          call LVT_verify(nf90_def_dim(ftn,'time',1,tdimID))
          call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"missing_value",&
               LVT_rc%udef))
!write lat/lon information
!write lat/lon information
          if(trim(LVT_rc%lvt_wopt).eq."1d tilespace") then
             call LVT_verify(nf90_def_var(ftn,&
                  trim(LVT_histData%xlat%short_name),&
                  nf90_float,&
                  dimids = dimID(1:1), varID=LVT_histData%xlat%varid_def))
#if(defined USE_NETCDF4) 
             call LVT_verify(nf90_def_var_deflate(ftn,&
                  LVT_histData%xlat%varid_def,&
                  shuffle,deflate,deflate_level))
#endif
             call LVT_verify(nf90_def_var(ftn,&
                  trim(LVT_histData%xlon%short_name),&
                  nf90_float,&
                  dimids = dimID(1:1), varID=LVT_histData%xlon%varid_def))
#if(defined USE_NETCDF4) 
             call LVT_verify(nf90_def_var_deflate(ftn,&
                  LVT_histData%xlon%varid_def,&
                  shuffle,deflate,deflate_level))
#endif             

          elseif(trim(LVT_rc%lvt_wopt).eq."2d gridspace") then 
             call LVT_verify(nf90_def_var(ftn,&
                  trim(LVT_histData%xlat%short_name),&
                  nf90_float,&
                  dimids = dimID(1:2), varID=LVT_histData%xlat%varid_def))
#if(defined USE_NETCDF4) 
             call LVT_verify(nf90_def_var_deflate(ftn,&
                  LVT_histData%xlat%varid_def,&
                  shuffle,deflate,deflate_level))
#endif
             call LVT_verify(nf90_def_var(ftn,&
                  trim(LVT_histData%xlon%short_name),&
                  nf90_float,&
                  dimids = dimID(1:2), varID=LVT_histData%xlon%varid_def))
#if(defined USE_NETCDF4) 
             call LVT_verify(nf90_def_var_deflate(ftn,&
                  LVT_histData%xlon%varid_def,&
                  shuffle,deflate,deflate_level))
#endif
          endif

          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlat%varid_def,&
               "units",trim(LVT_histData%xlat%units)))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlat%varid_def,&
               "standard_name",trim(LVT_histData%xlat%standard_name)))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlat%varid_def,&
               "long_name",trim(LVT_histData%xlat%long_name)))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlat%varid_def,&
               "scale_factor",1.0))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlat%varid_def,&
               "add_offset",0.0))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlat%varid_def,&
               "missing_value",LVT_rc%udef))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlat%varid_def,&
               "_FillValue",LVT_rc%udef))

          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlon%varid_def,&
               "units",trim(LVT_histData%xlon%units)))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlon%varid_def,&
               "standard_name",trim(LVT_histData%xlon%standard_name)))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlon%varid_def,&
               "long_name",trim(LVT_histData%xlon%long_name)))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlon%varid_def,&
               "scale_factor",1.0))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlon%varid_def,&
               "add_offset",0.0))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlon%varid_def,&
               "missing_value",LVT_rc%udef))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlon%varid_def,&
               "_FillValue",LVT_rc%udef))


!define time field
          call LVT_verify(nf90_def_var(ftn,'time',&
               nf90_float,dimids = tdimID,varID=LVT_histData%xtimeID))
          write(xtime_units,200) LVT_rc%yr, LVT_rc%mo, LVT_rc%da, &
               LVT_rc%hr, LVT_rc%mn, LVT_rc%ss
200       format ('minutes since ',I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':', &
               I2.2,':',I2.2)
          write(xtime_begin_date, fmt='(I4.4,I2.2,I2.2)') &
               LVT_rc%yr, LVT_rc%mo, LVT_rc%da
          write(xtime_begin_time, fmt='(I2.2,I2.2,I2.2)') &
               LVT_rc%hr, LVT_rc%mn, LVT_rc%ss
          write(xtime_timeInc, fmt='(I20)') &
               LVT_rc%ts
          
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xtimeID,&
               "units",trim(xtime_units)))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xtimeID,&
               "long_name","time"))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xtimeID,&
               "time_increment",trim(adjustl(xtime_timeInc))))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xtimeID,&
               "begin_date",xtime_begin_date))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xtimeID,&
               "begin_time",xtime_begin_time))

          call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"title", &
               "LVT land surface analysis output"))
          call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"institution", &
               trim(LVT_rc%institution)))
          call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"history", &
               "created on date: "//date(1:4)//"-"//date(5:6)//"-"//&
               date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10)))
          call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"references", &
               "Kumar_etal_GMD_2012"))
          call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"comment", &
               "website: http://lis.gsfc.nasa.gov/"))

          !grid information
          if(trim(LVT_rc%domain).eq."latlon") then !latlon
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
                  "EQUIDISTANT CYLINDRICAL"))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
                  LVT_rc%gridDesc(4)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LVT_rc%gridDesc(5)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
                  LVT_rc%gridDesc(9)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
                  LVT_rc%gridDesc(10)))       
          elseif(trim(LVT_rc%domain).eq."mercator") then 
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
                  "MERCATOR"))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
                  LVT_rc%gridDesc(4)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
                  LVT_rc%gridDesc(5)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
                  LVT_rc%gridDesc(10)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
                  LVT_rc%gridDesc(11)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
                  LVT_rc%gridDesc(8)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
                  LVT_rc%gridDesc(9)))
          elseif(trim(LVT_rc%domain).eq."lambert") then !lambert conformal
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
                  "LAMBERT CONFORMAL"))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
                  LVT_rc%gridDesc(4)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
                  LVT_rc%gridDesc(5)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
                  LVT_rc%gridDesc(10)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT2", &
                  LVT_rc%gridDesc(7)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
                  LVT_rc%gridDesc(11)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
                  LVT_rc%gridDesc(8)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
                  LVT_rc%gridDesc(9)))
             
          elseif(trim(LVT_rc%domain).eq."polar") then ! polar stereographic
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
                  "POLAR STEREOGRAPHIC"))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
                  LVT_rc%gridDesc(4)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
                  LVT_rc%gridDesc(5)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
                  LVT_rc%gridDesc(10)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"ORIENT", &
                  LVT_rc%gridDesc(7)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
                  LVT_rc%gridDesc(11)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
                  LVT_rc%gridDesc(8)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
                  LVT_rc%gridDesc(9)))
          endif

          if(LVT_rc%lvt_wopt.eq."1d tilespace") then
             if(ds2%vlevels.gt.1) then 
                call LVT_verify(nf90_def_dim(ftn,trim(ds2%short_name)//'_profiles',&
                     ds2%vlevels, dimID(3)))
             endif
             if(ds2%vlevels.eq.1) then 
                call LVT_verify(nf90_def_var(ftn,trim(ds2%short_name),nf90_float,&
                     dimids = dimID(1:1), varID=varId))
#if(defined USE_NETCDF4) 
                call LVT_verify(nf90_def_var_deflate(ftn,varID,&
                     shuffle,deflate,deflate_level))
#endif
             else
                dimID_t(1) = dimID(1)
                dimID_t(2) = dimID(3)
                call LVT_verify(nf90_def_var(ftn,trim(ds2%short_name),nf90_float,&
                     dimids = dimID_t(1:2), varID=varId),&
                     'nf90_def_var failed in LVT_historyMod')
#if(defined USE_NETCDF4) 
                call LVT_verify(nf90_def_var_deflate(ftn,varID,&
                     shuffle,deflate,deflate_level),&
                     'nf90_def_var_deflate failed in LVT_historyMod')
#endif
             endif
             
          elseif(LVT_rc%lvt_wopt.eq."2d gridspace") then 
             if(ds2%vlevels.gt.1) then 
                call LVT_verify(nf90_def_dim(ftn,trim(ds2%short_name)//'_profiles',&
                     ds2%vlevels, dimID(3)))
             endif
             if(ds2%vlevels.eq.1) then 
                call LVT_verify(nf90_def_var(ftn,trim(ds2%short_name),nf90_float,&
                     dimids = dimID(1:2), varID=varId))
#if(defined USE_NETCDF4) 
                call LVT_verify(nf90_def_var_deflate(ftn,varID,&
                  shuffle,deflate,deflate_level))
#endif
             else
                call LVT_verify(nf90_def_var(ftn,trim(ds2%short_name),nf90_float,&
                     dimids = dimID, varID=varId),&
                     'nf90_def_var failed in LVT_historyMod')
#if(defined USE_NETCDF4) 
                call LVT_verify(nf90_def_var_deflate(ftn,varID,&
                     shuffle,deflate,deflate_level),&
                     'nf90_def_var_deflate failed in LVT_historyMod')
#endif
             endif
          endif
          
          call LVT_verify(nf90_put_att(ftn,varID,&
               "units",trim(ds2%units)))
          call LVT_verify(nf90_put_att(ftn,varID,&
               "standard_name",trim(ds2%standard_name)))
          call LVT_verify(nf90_put_att(ftn,varID,&
               "long_name",trim(ds2%long_name)))
          call LVT_verify(nf90_put_att(ftn,varID,&
               "scale_factor",1.0))
          call LVT_verify(nf90_put_att(ftn,varID,&
               "add_offset",0.0))
          call LVT_verify(nf90_put_att(ftn,varID,&
               "missing_value",LVT_rc%udef))
          call LVT_verify(nf90_put_att(ftn,varID,&
               "_FillValue",LVT_rc%udef))
          
          call LVT_verify(nf90_enddef(ftn))
          
          call LVT_verify(nf90_put_var(ftn,LVT_histData%xtimeID,0.0))
          
          if(LVT_rc%lvt_wopt.eq."1d tilespace") then 
             do t=1,LVT_LIS_rc(1)%ntiles
                col = LVT_LIS_domain(1)%tile(t)%col
                row = LVT_LIS_domain(1)%tile(t)%row
                gid = LVT_domain%gindex(col,row)
                LVT_histData%xlat%value(t,1,1) = LVT_domain%grid(gid)%lat
                LVT_histData%xlon%value(t,1,1) = LVT_domain%grid(gid)%lon
             enddo
          else
             ! set lat/lons 
             do t=1,LVT_rc%ngrid
                LVT_histData%xlat%value(t,1,1) = LVT_domain%grid(t)%lat
                LVT_histData%xlon%value(t,1,1) = LVT_domain%grid(t)%lon
             enddo
          endif
          
          gtmp = LVT_rc%udef
          do r=1,LVT_rc%gnr
             do c=1,LVT_rc%gnc
                gid = LVT_domain%gindex(c,r)
                if(gid.ne.-1) then 
                   gtmp(c,r) = LVT_histData%xlat%value(gid,1,1)
                endif
             enddo
          enddo
          
          if(LVT_rc%lvt_wopt.eq."1d tilespace") then 
             iret = nf90_put_var(ftn,LVT_histData%xlat%varid_def, &
                  LVT_histData%xlat%value(:,1,1), (/1/),&
                  (/LVT_LIS_rc(1)%ntiles/))
          else
             iret = nf90_put_var(ftn,LVT_histData%xlat%varid_def, &
                  gtmp, (/1,1/),&
                  (/LVT_rc%gnc,LVT_rc%gnr/))
          endif
          
          gtmp = LVT_rc%udef
          do r=1,LVT_rc%gnr
             do c=1,LVT_rc%gnc
                gid = LVT_domain%gindex(c,r)
                if(gid.ne.-1) then 
                   gtmp(c,r) = LVT_histData%xlon%value(gid,1,1)
                endif
             enddo
          enddo
          
          if(LVT_rc%lvt_wopt.eq."1d tilespace") then 
             iret = nf90_put_var(ftn,LVT_histData%xlon%varid_def, &
                  LVT_histData%xlon%value(:,1,1), (/1/),&
                  (/LVT_LIS_rc(1)%ntiles/))
          else
             iret = nf90_put_var(ftn,LVT_histData%xlon%varid_def, &
                  gtmp, (/1,1/),&
                  (/LVT_rc%gnc,LVT_rc%gnr/))
          endif
          
          
          !finally write data: 
          if(trim(LVT_rc%lvt_wopt).eq."2d gridspace") then 
             gtmp = LVT_rc%udef
             do r=1,LVT_rc%gnr
                do c=1,LVT_rc%gnc
                   gid = LVT_domain%gindex(c,r)
                   if(gid.ne.-1) then 
                      gtmp(c,r) = ds2%value(gid,1,1)
                   endif
                enddo
             enddo
             
             iret = nf90_put_var(ftn,varid, gtmp, (/1,1/),&
                  (/LVT_rc%gnc,LVT_rc%gnr/))
          elseif(LVT_rc%lvt_wopt.eq."1d tilespace") then 
             iret = nf90_put_var(ftn,varid,ds2%value(:,1,1), (/1/),&
                  (/LVT_LIS_rc(1)%ntiles/))
          endif
          
          iret = nf90_close(ftn)
       endif
    endif
#endif
    
  end subroutine LVT_writeTrainingOutput
  
end module LVT_trainingMod

