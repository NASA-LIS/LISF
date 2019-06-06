!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
#include "LDT_NetCDF_inc.h"
module LDT_rstProcMod
!BOP
!
! !MODULE: LDT_rstProcMod
! 
! !DESCRIPTION: 
!   The code in this file provides interfaces to manage 
!   the processing of climatological restart files for land 
!   surface models and routing models. 
!
! !REVISION HISTORY: 
!  26 Jan 2016    Sujay Kumar  Initial Specification
!  18 Oct 2017    Hiroko Beaudoing   Added binary format and VIC4.1.2.
! 
  use ESMF
  use LDT_coreMod
  use LDT_logMod
  use LDT_fileIOMod
  use LDT_timeMgrMod
  use VIC_parmsMod
  use LDT_vic412rstMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LDT_rstProcInit
  public :: LDT_diagnoseRstData
  
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: LDT_rst_struc

  type, public :: rstoutdec
     integer                   :: nDims
     integer     , allocatable :: dims(:)
     real        , allocatable :: outVar(:,:,:)
     integer     , allocatable :: noutVar(:,:,:)     
  end type rstoutdec

  type, public ::  ldtrstdec
     integer,     allocatable     :: ftn(:)
     integer                      :: nVars
     logical                      :: startFlag
     character*100                :: oDir
     real                         :: outInterval
     character*50                 :: wstyle
     character*50                 :: wformat
     real                         :: modelTS

     character*50                 :: outMode
     character*50                 :: outIntervalType
     integer                      :: nfiles
     character*100, allocatable   :: outfname(:)
     type(rstoutdec), allocatable :: rstout(:)

  end type ldtrstdec
  
  type (ldtrstdec)  :: LDT_rst_struc

  contains

!BOP
! !ROUTINE: LDT_rstProcInit
! label{LDT_rstProcInit}
! 
! !INTERFACE: 
    subroutine LDT_rstProcInit
! 
! !DESCRIPTION: 
! 
!  This routine performs initialization in the restart processing mode by 
!  reading the relevant config file entries. 
!  
!  Example entries in the config file are as follows: 
!
! \begin{verabtim}
!   LIS restart source: "LSM"  # "LSM" "Routing"
!   Input restart file directory:        ../../OL_NLDAS/OUTPUT
!   Input restart file output interval:  "1mo"
!   Input restart file naming style:      "3 level hierarchy"
!   Input restart model timestep used:    "15mn"
!
!   Output restart file generation mode:   "climatological average"
!   Output restart file averaging interval type:  "monthly" #"3 monthly", "daily"
!  \end{verbatim}
!
!EOP  
  
      integer                   :: n,i 
      integer                   :: status
      character*20              :: stime
      character*3               :: month_name(12)
      character*100             :: model_name      

      n = 1
      
      month_name = (/"JAN","FEB","MAR","APR","MAY","JUN",&
           "JUL","AUG","SEP","OCT","NOV","DEC"/)

      call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%rstsource,&
           label="LIS restart source:",&
           rc=status)
      call LDT_verify(status,'LIS restart source: not defined')

      call ESMF_ConfigGetAttribute(LDT_config,LDT_rst_struc%odir,&
           label="Input restart file directory:",&
           rc=status)
      call LDT_verify(status,'Input restart file directory: not defined')

      call ESMF_ConfigGetAttribute(LDT_config,stime, &
           label="Input restart file output interval:",&
           rc=status)
      call LDT_verify(status,'Input restart file output interval: not defined')
      call LDT_parseTimeString(stime, LDT_rst_struc%outInterval)

      call ESMF_ConfigGetAttribute(LDT_config,LDT_rst_struc%wstyle,&
           label="Input restart file naming style:",&
           rc=status)
      call LDT_verify(status,'Input restart file naming style: not defined')
      
      call ESMF_ConfigGetAttribute(LDT_config,stime, &
           label="Input restart model timestep used:",&
           rc=status)
      call LDT_verify(status,'Input restart model timestep used: not defined')
      call LDT_parseTimeString(stime, LDT_rst_struc%modelTs)

!hkb add input file format option. if not specified in config, assume netcdf
      call ESMF_ConfigGetAttribute(LDT_config,LDT_rst_struc%wformat, &
           label="Input restart file format:",&
           default="netcdf",rc=status)

      LDT_rc%pass = 1

      call LDT_registerAlarm("LIS restart alarm",&
           LDT_rst_struc%modelTs, &
           LDT_rst_struc%outInterval)

      call LDT_update_timestep(LDT_rc, n, LDT_rst_struc%modelTs)


      call ESMF_ConfigGetAttribute(LDT_config,LDT_rst_struc%outMode,&
           label="Output restart file generation mode:",&
           rc=status)
      call LDT_verify(status,'Output restart file generation mode: not defined')

      call ESMF_ConfigGetAttribute(LDT_config,LDT_rst_struc%outIntervalType,&
           label="Output restart file averaging interval type:",&
           rc=status)
      call LDT_verify(status,'Output restart file averaging interval type: not defined')

      if(LDT_rc%rstsource.eq."LSM") then 
         if(LDT_rc%lsm.eq."Noah.3.2") then 
            model_name = "NOAH32"
         elseif(LDT_rc%lsm.eq."Noah.3.3") then 
            model_name = "NOAH33"
         elseif(LDT_rc%lsm.eq."Noah.3.6") then 
            model_name = "NOAH36"
         elseif(LDT_rc%lsm.eq."Noah.3.9") then 
            model_name = "NOAH39"
         elseif(LDT_rc%lsm.eq."Noah.2.7.1") then 
            model_name = "NOAH271"
         elseif(LDT_rc%lsm.eq."Noah-MP.3.6") then 
            model_name = "NOAHMP36"
         elseif(LDT_rc%lsm.eq."Noah-MP.4.0.1") then 
            model_name = "NOAHMP401"
         elseif(LDT_rc%lsm.eq."CLSMF2.5") then 
            model_name = "CLSMF25"
         elseif(LDT_rc%lsm.eq."RUC.3.7") then 
            model_name = "RUC37"
         elseif(LDT_rc%lsm.eq."VIC.4.1.1") then 
            model_name = "VIC411"
         elseif(LDT_rc%lsm.eq."VIC.4.1.2") then 
            model_name = "VIC412"
         else
            write(LDT_logunit,*) "[INFO] Climatological Restart File Generation - LSMs supported: "
            write(LDT_logunit,*) "  -- CLSMF2.5, Noah.3.2, Noah.3.3, Noah.3.6, Noah.3.9, "
            write(LDT_logunit,*) "  -- Noah-MP.3.6, Noah-MP.4.0.1, "
            write(LDT_logunit,*) "     Noah.2.7.1, RUC.3.7, VIC.4.1.1, VIC.4.1.2 "
            write(LDT_logunit,*) "[ERR] No other LSMs supported at this time ... stopping."
            call LDT_endrun() 
         endif

         if(LDT_rst_struc%outMode.eq."climatological average") then 
            if(LDT_rst_struc%outIntervalType.eq."monthly") then 
               
               LDT_rst_struc%nfiles = 12 
               
               allocate(LDT_rst_struc%outfname(LDT_rst_struc%nfiles))
               allocate(LDT_rst_struc%ftn(LDT_rst_struc%nfiles))
               
               call system('mkdir -p '//(LDT_rc%odir))
               do i=1,12
!hkb add netcdf and binary format selections below.
                 if(LDT_rst_struc%wformat.eq."netcdf") then
                  LDT_rst_struc%outfname(i) = &
                       trim(LDT_rc%odir)//"/LDT_CLIM_RST_"//trim(month_name(i))//&
                       "_"//trim(model_name)//".nc"
                  
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
               
#if (defined USE_NETCDF4)
                  call LDT_verify(nf90_create(path=LDT_rst_struc%outfname(i),&
                       cmode=nf90_netcdf4,&
                       ncid = LDT_rst_struc%ftn(i)),&
                       'creating netcdf file failed in LDT_rstProcMod')
#endif
#if (defined USE_NETCDF3)
                  call LDT_verify(nf90_create(path=LDT_rst_struc%outfname(i),&
                       cmode=nf90_clobber,&
                       ncid = LDT_rst_struc%ftn(i)),&
                       'creating netcdf file failed in LDT_rstProcMod')
#endif
#endif
                 elseif (LDT_rst_struc%wformat.eq."binary") then
                  LDT_rst_struc%outfname(i) = &
                       trim(LDT_rc%odir)//"/LDT_CLIM_RST_"//trim(month_name(i))//&
                       "_"//trim(model_name)//".bin"
                  LDT_rst_struc%ftn(i) = LDT_getNextUnitNumber()
                  open(LDT_rst_struc%ftn(i),&
                       file=trim(LDT_rst_struc%outfname(i)), &
                       form='unformatted')
                 endif  !wformat
               enddo
            else
               write(LDT_logunit,*) '[ERR] Averaging interval type ',&
                    trim(LDT_rst_struc%outIntervalType)
               write(LDT_logunit,*) '[ERR] is not currently supported...'
               call LDT_endrun()
               
            endif
         else
            write(LDT_logunit,*) '[ERR] Restart processing mode ',&
                 trim(LDT_rst_struc%outMode)
            write(LDT_logunit,*) '[ERR] is not currently supported...'
            call LDT_endrun()
         endif
      elseif(LDT_rc%rstsource.eq."Routing") then 
         model_name = LDT_rc%routingmodel

         if(LDT_rst_struc%outMode.eq."climatological average") then 
            if(LDT_rst_struc%outIntervalType.eq."monthly") then 
               
               LDT_rst_struc%nfiles = 12 
               
               allocate(LDT_rst_struc%outfname(LDT_rst_struc%nfiles))
               allocate(LDT_rst_struc%ftn(LDT_rst_struc%nfiles))
               
               call system('mkdir -p '//(LDT_rc%odir))
               do i=1,12
                  LDT_rst_struc%outfname(i) = &
                       trim(LDT_rc%odir)//"/LDT_CLIM_RST_"//trim(month_name(i))//&
                       "_"//trim(model_name)//".bin"
                  
                   LDT_rst_struc%ftn(i) = LDT_getNextUnitNumber()
                   open(LDT_rst_struc%ftn(i),&
                        file=trim(LDT_rst_struc%outfname(i)), &
                        form='unformatted')
                enddo
            else
               write(LDT_logunit,*) '[ERR] Averaging interval type ',&
                    trim(LDT_rst_struc%outIntervalType)
               write(LDT_logunit,*) '[ERR] is not currently supported...'
               call LDT_endrun()
               
            endif
         else
            write(LDT_logunit,*) '[ERR] Restart processing mode ',&
                 trim(LDT_rst_struc%outMode)
            write(LDT_logunit,*) '[ERR] is not currently supported...'
            call LDT_endrun()
         endif

      endif

      LDT_rst_struc%startFlag = .true. 

    end subroutine LDT_rstProcInit

!BOP
! 
! !ROUTINE: LDT_diagnoseRstData
!  \label{LDT_diagnoseRstData}
! 
! !ROUTINE: 
    subroutine LDT_diagnoseRstData(n)
!
! !DESCRIPTION: 
! 
!   This routine reads each LIS restart file, and then stores the variables 
!   for the processing of climatological restart files. 
! 
!   The arguments are:
!   \begin{description}
!   \item [n]
!     index of the nest
!   \end{description}
!EOP
      integer,      intent(in)  :: n 

      integer                   :: k
      integer                   :: c,r
      integer                   :: i
      integer                   :: kk
      integer                   :: t
      integer                   :: tindex
      integer                   :: ftn
      character*100             :: dir_name
      character*100             :: model_name
      character*100             :: fname
      integer                   :: nDims
      integer                   :: nVars
      integer                   :: nGlobalAtts
      integer                   :: unlimdimID
      integer,      allocatable :: dimID(:)
      integer,      allocatable :: dims(:)
      integer,      allocatable :: n_dimids(:)
      character*50, allocatable :: dimName(:)
      character*50              :: varName
      character*1               :: fd
      integer                   :: tdimId
      character*33              :: units
      character*50              :: tIncr
      character*8               :: beg_date
      character*6               :: beg_time
      integer                   :: nvardims
      character*100             :: standard_name
      real                      :: scale_factor
      real                      :: offset
      real                      :: vmin
      real                      :: vmax
      real                      :: gmt
      real*8                    :: time
      real,        allocatable  :: var(:,:)
      integer                   :: yr, mo, da, hr, mn, ss,doy
      logical                   :: alarmCheck
      logical                   :: file_exists
      real,        allocatable  :: tmptilen(:)
      integer                   :: nc, nr, npatch, state_chunk_size

      alarmCheck = LDT_isAlarmRinging(LDT_rc, &
           "LIS restart alarm")
      
      if(alarmCheck) then 
         !read each file, save to variables. 
         
         if(LDT_rc%rstsource.eq."LSM") then 
           dir_name = 'SURFACEMODEL'
           if(LDT_rc%lsm.eq."Noah.3.2") then
              model_name = "NOAH32"
           elseif(LDT_rc%lsm.eq."Noah.3.3") then
              model_name = "NOAH33"
           elseif(LDT_rc%lsm.eq."Noah.3.6") then
              model_name = "NOAH36"
           elseif(LDT_rc%lsm.eq."Noah.3.9") then
              model_name = "NOAH39"
           elseif(LDT_rc%lsm.eq."Noah.2.7.1") then
              model_name = "NOAH271"
           elseif(LDT_rc%lsm.eq."Noah-MP.3.6") then
              model_name = "NOAHMP36"
           elseif(LDT_rc%lsm.eq."Noah-MP.4.0.1") then
              model_name = "NOAHMP401"
           elseif(LDT_rc%lsm.eq."CLSMF2.5") then 
              model_name = "CLSMF25"
           elseif(LDT_rc%lsm.eq."RUC.3.7") then    
              model_name = "RUC37"
           elseif(LDT_rc%lsm.eq."VIC.4.1.1") then
              model_name = "VIC411"
           elseif(LDT_rc%lsm.eq."VIC.4.1.2") then
              model_name = "VIC412"
           else
              write(LDT_logunit,*) "[INFO] Climatological Restart File Generation - LSMs supported: "
              write(LDT_logunit,*) "  -- CLSMF2.5, Noah.3.2, Noah.3.3, Noah.3.6, Noah.3.9, "
              write(LDT_logunit,*) "  -- Noah-MP.3.6, Noah-MP.4.0.1, "
              write(LDT_logunit,*) "     Noah.2.7.1, RUC.3.7, VIC.4.1.1, VIC.4.1.2 "
              write(LDT_logunit,*) "[ERR] No other LSMs supported at this time ... stopping."
              call LDT_endrun() 
            endif
         elseif(LDT_rc%rstsource.eq."Routing") then 
            dir_name = 'ROUTING'
            model_name = LDT_rc%routingmodel
         endif

         yr = LDT_rc%yr
         mo = LDT_rc%mo
         da = LDT_rc%da
         hr = LDT_rc%hr
         mn = LDT_rc%mn
         ss = LDT_rc%ss

         if(LDT_rst_struc%outMode.eq."climatological average") then 
            if(LDT_rst_struc%outIntervalType.eq."monthly") then 
               tindex = LDT_rc%mo
            endif
         endif
         
         if(LDT_rc%rstsource.eq."LSM") then 
!hkb add netcdf and binary format selections below.
          if(LDT_rst_struc%wformat.eq."netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
           call LDT_create_restart_filename(n,&
              LDT_rst_struc%odir,   &
              dir_name,        &   
              model_name,      &
              yr,       &
              mo,       &
              da,       &
              hr,       &
              mn,       &
              ss,       &
              LDT_rst_struc%wstyle, &
              "netcdf",&
              fname)

           inquire(file=trim(fname),exist=file_exists) 

           if(.not.file_exists) then 
              write(LDT_logunit,*) '[ERR] restart file ',trim(fname)
              write(LDT_logunit,*) '[ERR] does not exist'
              call LDT_endrun()
           endif

           write(LDT_logunit,*) '[INFO] Reading restart file ',trim(fname)

            call LDT_verify(nf90_open(path=fname,&
                 mode=nf90_NOWRITE,ncid=ftn),&
                 'failed to open '//trim(fname))
            
            call LDT_verify(nf90_inquire(ftn,nDims,nVars,nGlobalAtts,unlimdimId),&
                 'nf90_inquire failed in LDT_rstProcMod')
            
            LDT_rst_struc%nVars = nVars

            allocate(dimID(nDims))
            allocate(dims(nDims))
            allocate(dimName(nDims))
            
            call LDT_verify(nf90_inq_dimId(ftn,"ntiles",dimId(1)),&
                 'nf90_inq_dimId 1 failed in LDT_rstProcMod')
            call LDT_verify(nf90_inquire_dimension(ftn,dimId(1),len=dims(1)),&
                 'nf90_inquire_dimension failed in LDT_rstProcMod')

            if ( nDims .gt. 2 ) then
             do k=2,nDims - 1
                write(unit=fd,fmt='(I1)') k-1
                dimName(k) = 'dim'//trim(fd)
                call LDT_verify(nf90_inq_dimId(ftn,dimName(k),dimId(k)),&
                     'nf90_inq_dimId 2 failed in LDT_rstProcMod')
                call LDT_verify(nf90_inquire_dimension(ftn,dimID(k),len=dims(k)),&
                     'nf90_inquire failed in LDT_rstProcMod')
             enddo

             call LDT_verify(nf90_inq_dimId(ftn,"time",tdimId),&
                  'nf90_inq_dimId 3 failed for LDT_rstProcMod')
            else  ! nDims = 2 
             do k=2,nDims 
                write(unit=fd,fmt='(I1)') k-1
                dimName(k) = 'dim'//trim(fd)
                call LDT_verify(nf90_inq_dimId(ftn,dimName(k),dimId(k)),&
                     'nf90_inq_dimId 2 failed in LDT_rstProcMod')
                call LDT_verify(nf90_inquire_dimension(ftn,dimID(k),len=dims(k)),&
                     'nf90_inquire failed in LDT_rstProcMod')
             enddo
            endif

            if(LDT_rst_struc%startFlag) then 
               
              if ( nVars .gt. 1 ) then
               allocate(LDT_rst_struc%rstOut(nVars-1))
              else  ! nVars = 1
               allocate(LDT_rst_struc%rstOut(nVars))
              endif
               
               do i=1,LDT_rst_struc%nfiles
                  call writeglobalheader(LDT_rst_struc%ftn(i),model_name, &
                       nDims, dims, dimID)
                  call LDT_verify(nf90_enddef(LDT_rst_struc%ftn(i)))
               enddo

            endif
            
            do k=1,nVars
               call LDT_verify(nf90_inquire_variable(ftn,k,varName,nDims=nvardims),&
                 'nf90_inquire_variable failed in LDT_rstProcMod')
               allocate(n_dimids(nvardims))
               call LDT_verify(nf90_inquire_variable(ftn,k,varName,dimids=n_dimids),&
                    'nf90_inquire_variable failed in LDT_rstProcMod')
               if(varname.eq."time") then 
                  call LDT_verify(nf90_get_att(ftn,k,"units",units),&
                       'nf90_get_att failed for time in LDT_rstProcMod')
                  call LDT_verify(nf90_get_att(ftn,k,"time_increment",tincr),&
                       'nf90_get_att failed for time_increment in LDT_rstProcMod')
                  call LDT_verify(nf90_get_att(ftn,k,"begin_date",beg_date),&
                       'nf90_get_att failed for begin_date in LDT_rstProcMod')
                  call LDT_verify(nf90_get_att(ftn,k,"begin_time",beg_time),&
                       'nf90_get_att failed for begin_time in LDT_rstProcMod')
               else
                  call LDT_verify(nf90_get_att(ftn,k,"units",units),&
                       'nf90_get_att failed in LDT_rstProcMod')
                  call LDT_verify(nf90_get_att(ftn,k,"standard_name",standard_name),&
                       'nf90_get_att failed in LDT_rstProcMod')
                  call LDT_verify(nf90_get_att(ftn,k,"scale_factor",scale_factor),&
                       'nf90_get_att failed in LDT_rstProcMod')
                  call LDT_verify(nf90_get_att(ftn,k,"add_offset",offset),&
                       'nf90_get_att failed in LDT_rstProcMod')
                  call LDT_verify(nf90_get_att(ftn,k,"vmin",vmin),&
                       'nf90_get_att failed in LDT_rstProcMod')
                  call LDT_verify(nf90_get_att(ftn,k,"vmax",vmax),&
                       'nf90_get_att failed in LDT_rstProcMod')

                  if(LDT_rst_struc%startFlag) then 
                     if(nvardims.gt.1) then 
                        LDT_rst_struc%rstOut(k)%nDims = 2
                        allocate(LDT_rst_struc%rstOut(k)%dims(LDT_rst_struc%rstOut(k)%nDims))
                        LDT_rst_struc%rstOut(k)%dims(1) = dims(1)
                        LDT_rst_struc%rstOut(k)%dims(2) = dims(n_dimids(2))

                        allocate(LDT_rst_struc%rstOut(k)%outVar(&
                             dims(1),dims(n_dimids(2)),LDT_rst_struc%nfiles))
                        allocate(LDT_rst_struc%rstOut(k)%noutVar(&
                             dims(1),dims(n_dimids(2)),LDT_rst_struc%nfiles))
                        LDT_rst_struc%rstOut(k)%outVar = 0.0
                        LDT_rst_struc%rstOut(k)%noutVar = 0.0
!hkb allocate vic arrays
                        if(LDT_rc%lsm.eq."VIC.4.1.2") then
                           call LDT_vic412rstInit(n,LDT_rst_struc%nfiles, &
                                                  dims(1),dims(n_dimids(2)))
                           allocate(VIC_struc(n)%state_chunk(dims(1),dims(n_dimids(2))))
                        endif
                     else
                        LDT_rst_struc%rstOut(k)%nDims = 1
                        allocate(LDT_rst_struc%rstOut(k)%dims(LDT_rst_struc%rstOut(k)%nDims))
                        LDT_rst_struc%rstOut(k)%dims(1) = dims(1)

                        allocate(LDT_rst_struc%rstOut(k)%outVar(dims(1),1,LDT_rst_struc%nfiles))
                        allocate(LDT_rst_struc%rstOut(k)%noutVar(dims(1),1,LDT_rst_struc%nfiles))
                        
                        LDT_rst_struc%rstOut(k)%outVar = 0.0
                        LDT_rst_struc%rstOut(k)%noutVar = 0.0
                     endif

                     do i=1,LDT_rst_struc%nfiles
                        call writeheader_restart(LDT_rst_struc%ftn(i),&
                             nvardims,            &
                             n_dimIds,            &
                             k,                   &
                             standard_name,       &
                             units,               &
                             scale_factor,        &
                             offset,              &
                             vmin,                &
                             vmax)
                     enddo
                  endif
                  
                  if(nvardims.gt.1) then 

                     allocate(var(dims(1),dims(n_dimids(2))))
                     call LDT_verify(nf90_get_var(ftn,k,var),&
                          'nf90_get_var failed in LDT_rstProcMod')

!hkb vic needs unpacking and accumulate here !!!
                    if(LDT_rc%lsm.eq."VIC.4.1.2") then
                     call LDT_vic412rstDiagnose(n,tindex, &
                          dims(1),dims(n_dimids(2)),var)
                     do t=1,dims(1)
                        LDT_rst_struc%rstOut(k)%noutVar(t,1,tindex) = & 
                             LDT_rst_struc%rstOut(k)%noutVar(t,1,tindex) + 1
                     enddo
                    else  ! all other LSMs
                     do kk=1,dims(n_dimids(2))
                        do t=1,dims(1)
                           LDT_rst_struc%rstOut(k)%outVar(t,kk,tindex) = & 
                                LDT_rst_struc%rstOut(k)%outVar(t,kk,tindex) + & 
                                var(t,kk) 
                           LDT_rst_struc%rstOut(k)%noutVar(t,kk,tindex) = & 
                                LDT_rst_struc%rstOut(k)%noutVar(t,kk,tindex) + 1
                        enddo
                     enddo
                    endif   ! VIC

                  else

                     allocate(var(dims(1),1))
                     call LDT_verify(nf90_get_var(ftn,k,var),&
                          'nf90_get_var failed in LDT_rstProcMod')

                     do t=1,dims(1)
                        LDT_rst_struc%rstOut(k)%outVar(t,1,tindex) = & 
                             LDT_rst_struc%rstOut(k)%outVar(t,1,tindex) + & 
                             var(t,1) 
                        LDT_rst_struc%rstOut(k)%noutVar(t,1,tindex) = & 
                             LDT_rst_struc%rstOut(k)%noutVar(t,1,tindex) + 1
                     enddo
                  endif

                  deallocate(var)
               endif
               
               deallocate( n_dimids )
            enddo

            LDT_rst_struc%startFlag = .false. 
            
            deallocate(dimID)
            deallocate(dims)
            deallocate(dimName)
#endif
          elseif(LDT_rst_struc%wformat.eq."binary") then

           call LDT_create_restart_filename(n,&
              LDT_rst_struc%odir,   &
              dir_name,        &   
              model_name,      &
              yr,       &
              mo,       &
              da,       &
              hr,       &
              mn,       &
              ss,       &
              LDT_rst_struc%wstyle, &
              "binary",&
              fname)

           inquire(file=trim(fname),exist=file_exists) 

           if(.not.file_exists) then 
              write(LDT_logunit,*) '[ERR] restart file ',trim(fname)
              write(LDT_logunit,*) '[ERR] does not exist'
              call LDT_endrun()
           endif

           write(LDT_logunit,*) '[INFO] Reading restart file ',trim(fname)
              ftn = LDT_getNextUnitNumber()
              open(ftn,file=trim(fname),form='unformatted')
              read(ftn) nc,nr,npatch  !time, veg class, no. tiles

              if (nc.ne.LDT_rc%gnc(n) .or. nr.ne.LDT_rc%gnr(n)) then
                 write(LDT_logunit,*) '[ERR] grid space mismatch - VIC 4.1.2 halted'
                 call LDT_endrun()
              endif

           if(LDT_rc%lsm.eq."VIC.4.1.2") then
              call count_model_state_412(n,npatch,state_chunk_size)
              LDT_rst_struc%nVars = 1

             if(LDT_rst_struc%startFlag) then
               allocate(LDT_rst_struc%rstOut(LDT_rst_struc%nVars))

               do k=1,LDT_rst_struc%nVars
                  LDT_rst_struc%rstOut(k)%nDims = 2
                  allocate(LDT_rst_struc%rstOut(k)%dims(LDT_rst_struc%rstOut(k)%nDims))
                  
                  LDT_rst_struc%rstOut(k)%dims(1) = npatch
                  LDT_rst_struc%rstOut(k)%dims(2) = state_chunk_size
                  
                  allocate(LDT_rst_struc%rstOut(k)%outVar(&
                       LDT_rst_struc%rstOut(k)%dims(1),&
                       LDT_rst_struc%rstOut(k)%dims(2),&
                       LDT_rst_struc%nfiles))
                  
                  allocate(LDT_rst_struc%rstOut(k)%noutVar(&
                       LDT_rst_struc%rstOut(k)%dims(1),&
                       LDT_rst_struc%rstOut(k)%dims(2),&
                       LDT_rst_struc%nfiles))
                  
                  LDT_rst_struc%rstOut(k)%outVar = 0.0
                  LDT_rst_struc%rstOut(k)%noutVar = 0.0
               enddo

               call LDT_vic412rstInit(n,LDT_rst_struc%nfiles, &
                                                  npatch,state_chunk_size)
               allocate(VIC_struc(n)%state_chunk(npatch,state_chunk_size))

               LDT_rst_struc%startFlag = .false. 
             endif

             allocate(var(npatch,state_chunk_size))
             allocate(tmptilen(npatch))

             do k=1,LDT_rst_struc%nVars
                do kk=1,state_chunk_size
                  read(ftn) tmptilen
                  do t=1,npatch
                     var(t,kk) = tmptilen(t)
                  enddo
                end do
                call LDT_vic412rstDiagnose(n,tindex,npatch,state_chunk_size,var)
                do t=1,npatch
                   LDT_rst_struc%rstOut(k)%noutVar(t,1,tindex) = & 
                          LDT_rst_struc%rstOut(k)%noutVar(t,1,tindex) + 1
                enddo
             enddo

             call LDT_releaseUnitNumber(ftn)

             deallocate(var)
             deallocate(tmptilen)

           else ! all other models
              write(LDT_logunit,*) '[ERR] binary format restart file is not'
              write(LDT_logunit,*) '[ERR] supported for ',trim(LDT_rc%lsm)
              write(LDT_logunit,*) '[ERR] ...stopping.'
              call LDT_endrun()
           endif  ! VIC

          endif ! binary or netcdf format

         elseif(LDT_rc%rstsource.eq."Routing") then 
            model_name = trim(model_name)//"_router"
            call LDT_create_restart_filename(n,&
                 LDT_rst_struc%odir,   &
                 dir_name,        &   
                 model_name,      &
                 yr,       &
                 mo,       &
                 da,       &
                 hr,       &
                 mn,       &
                 ss,       &
                 LDT_rst_struc%wstyle, &
                 "binary", & 
                 fname)
            
            inquire(file=trim(fname),exist=file_exists) 
            
            if(.not.file_exists) then 
               write(LDT_logunit,*) '[ERR] restart file ',trim(fname)
               write(LDT_logunit,*) '[ERR] does not exist'
               call LDT_endrun()
            endif

            if(LDT_rst_struc%startFlag) then 
               LDT_rst_struc%nVars = 4    ! HYMAP has 4 variables in the restart file
               allocate(LDT_rst_struc%rstOut(LDT_rst_struc%nVars))

               do k=1,LDT_rst_struc%nVars
                  LDT_rst_struc%rstOut(k)%nDims = 2
                  allocate(LDT_rst_struc%rstOut(k)%dims(LDT_rst_struc%rstOut(k)%nDims))
                  
                  LDT_rst_struc%rstOut(k)%dims(1) = LDT_rc%gnc(n)
                  LDT_rst_struc%rstOut(k)%dims(2) = LDT_rc%gnr(n)
                  
                  allocate(LDT_rst_struc%rstOut(k)%outVar(&
                       LDT_rst_struc%rstOut(k)%dims(1),&
                       LDT_rst_struc%rstOut(k)%dims(2),&
                       LDT_rst_struc%nfiles))
                  
                  allocate(LDT_rst_struc%rstOut(k)%noutVar(&
                       LDT_rst_struc%rstOut(k)%dims(1),&
                       LDT_rst_struc%rstOut(k)%dims(2),&
                       LDT_rst_struc%nfiles))
                  
                  LDT_rst_struc%rstOut(k)%outVar = 0.0
                  LDT_rst_struc%rstOut(k)%noutVar = 0.0
               enddo

               LDT_rst_struc%startFlag = .false. 
            endif

            write(LDT_logunit,*) '[INFO] Reading restart file ',trim(fname)

            allocate(var(LDT_rc%gnc(n), LDT_rc%gnr(n)))
            
            ftn = LDT_getNextUnitNumber()
            
            open(ftn,file=trim(fname), form='unformatted')

            do k=1,LDT_rst_struc%nVars
               
               read(ftn) var
               do r=1,LDT_rc%gnr(n)
                  do c=1,LDT_rc%gnc(n)
                     if(var(c,r).ne.LDT_rc%udef) then 
                        LDT_rst_struc%rstOut(k)%outVar(c,r,tindex) = & 
                             LDT_rst_struc%rstOut(k)%outVar(c,r,tindex) + & 
                             var(c,r) 
                        
                        LDT_rst_struc%rstOut(k)%noutVar(c,r,tindex) = & 
                             LDT_rst_struc%rstOut(k)%noutVar(c,r,tindex) + 1
                     endif
                  enddo
               enddo
            enddo
            call LDT_releaseUnitNumber(ftn)

            deallocate(var)
         endif
      endif

      if(LDT_rc%endtime.eq.1) then 
         call writeRstData(n,LDT_rst_struc%nVars,LDT_rst_struc%wformat)
!hkb deallocate vic arrays
         if((LDT_rc%lsm.eq."VIC.4.1.2").and.(LDT_rc%rstsource.eq."LSM")) then
            call LDT_vic412rstFinalize(LDT_rst_struc%nfiles, &
                                       LDT_rst_struc%rstOut(1)%dims(1))

            deallocate(VIC_struc(n)%state_chunk)
         endif  ! VIC
      end if
      

    end subroutine LDT_diagnoseRstData

!BOP
! !ROUTINE: writeRstData
! \label{writeRstData}
!
! !INTERFACE: 
    subroutine writeRstData(n,nVars,wformat)
! !ARGUMENTS:       
      integer           :: n      ! nest index
      integer           :: nVars
      character*50      :: wformat
!
! !DESCRIPTION:
!  This routine writes the time averaged restart data into the
!  climatological restart files. The climatology restart files
!  for LSMs and routing models are in NetCDF and binary files, 
!  respectively. 
!
!   The arguments are:
!   \begin{description}
!   \item [nVars]
!     number of variables in the climatological restart file
!   \item [wformat]
!     format of the climatological restart file
!   \end{description}
! 
!EOP
      integer           :: k
      integer           :: kk
      integer           :: l
      integer           :: t
      integer           :: dims(2)
      integer           :: nCount
      real, allocatable :: tmptilen(:)

      if(LDT_rc%rstsource.eq."LSM") then 
        if ( nVars .gt. 1 ) then
         nCount = nVars-1
        else   ! nVars = 1
         nCount = nVars
        endif
      elseif(LDT_rc%rstsource.eq."Routing") then 
         nCount = nVars
      endif

      do k=1,nCount
         if(LDT_rst_struc%rstOut(k)%nDims.gt.1) then 
            dims(1) = LDT_rst_struc%rstOut(k)%dims(1)
            dims(2) = LDT_rst_struc%rstOut(k)%dims(2)
         else
            dims(1) = LDT_rst_struc%rstOut(k)%dims(1)
            dims(2) = 1
         endif
!hkb vic needs averaging & packing here !!!
        if((LDT_rc%lsm.eq."VIC.4.1.2").and.(LDT_rc%rstsource.eq."LSM")) then
         do l=1,LDT_rst_struc%nfiles
            call LDT_vic412rstAvePack(n,l,dims(1),dims(2),LDT_rst_struc%rstOut(k)%noutVar(:,1,l))
            do t=1,dims(1)
               do kk=1,dims(2)
                  LDT_rst_struc%rstOut(k)%outVar(t,kk,l) = VIC_struc(n)%state_chunk(t,kk)
               enddo
            enddo
         enddo
        else  ! all other LSMs
         do l=1,LDT_rst_struc%nfiles
            do t=1,dims(1)
               do kk=1,dims(2)
                  LDT_rst_struc%rstOut(k)%outVar(t,kk,l) = & 
                       LDT_rst_struc%rstOut(k)%outVar(t,kk,l) /&
                       LDT_rst_struc%rstOut(k)%noutVar(t,kk,l) 
                  
               enddo
            enddo
         enddo
        endif  ! VIC
      enddo

      if(LDT_rc%rstsource.eq."LSM") then 
!hkb add netcdf and binary format selections below.
        if(LDT_rst_struc%wformat.eq."netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)       
      
         do l=1,LDT_rst_struc%nfiles
!hkb            do k=1, LDT_rst_struc%nVars-1
            do k=1, nCount
               if(LDT_rst_struc%rstOut(k)%nDims.gt.1) then 
                  call LDT_verify(nf90_put_var(LDT_rst_struc%ftn(l),&
                       k,&
                       LDT_rst_struc%rstOut(k)%outVar(:,:,l),&
                       (/1,1/),&
                       (/LDT_rst_struc%rstOut(k)%dims(1),LDT_rst_struc%rstOut(k)%dims(2)/)),&
                       'nf90_put_var failed in LDT_rstProcMod ')
               else
                  call LDT_verify(nf90_put_var(LDT_rst_struc%ftn(l),&
                       k,&
                       LDT_rst_struc%rstOut(k)%outVar(:,:,l),&
                       (/1,1/),&
                       (/LDT_rst_struc%rstOut(k)%dims(1),1/)),&
                       'nf90_put_var failed in LDT_rstProcMod ')
               endif
            enddo
            call LDT_verify(nf90_close(LDT_rst_struc%ftn(l)))
         enddo
#endif
        elseif(LDT_rst_struc%wformat.eq."binary") then
         allocate(tmptilen(dims(1)))

         do l=1,LDT_rst_struc%nfiles
            write(LDT_rst_struc%ftn(l)) LDT_rc%gnc(n),LDT_rc%gnr(n),LDT_rst_struc%rstOut(1)%dims(1)
            do k=1, nCount
               do kk=1,dims(2)
                  tmptilen = 0
                  do t=1,dims(1)
                     tmptilen(t) = LDT_rst_struc%rstOut(k)%outVar(t,kk,l)
                  enddo
                  write(LDT_rst_struc%ftn(l)) tmptilen
               enddo
            enddo
            call LDT_releaseUnitNumber(LDT_rst_struc%ftn(l))
         enddo
         deallocate(tmptilen)
        endif

      elseif(LDT_rc%rstsource.eq."Routing") then 

         do l=1,LDT_rst_struc%nfiles
            do k=1, LDT_rst_struc%nVars
               write(LDT_rst_struc%ftn(l))  LDT_rst_struc%rstOut(k)%outVar(:,:,l)
            enddo
         
            call LDT_releaseUnitNumber(LDT_rst_struc%ftn(l))
         enddo
                     
      endif
      
    end subroutine writeRstData

!BOP
! 
! \label{writeglboalheader}
! 
! !INTERFACE: 
    subroutine writeglobalheader(ftn,model_name, &
         nDims, dims, dimID)
! 
! !DESCRIPTION: 
!  This routine writes the global NetCDF header information
!  in the climatological restart files. 
!
!   The arguments are:
!   \begin{description}
!   \item [ftn]
!     unit number of the file being written
!   \item [model\_name]
!     name of the model 
!   \item [nDims]
!     number of dimensions in the NetCDF file
!   \item [dims]
!     values of the dimension in the NetCDF file
!   \item [dimID]
!     NetCDF dimension ID values 
!   \end{description}
!EOP
      integer               :: n 
      integer               :: ftn
      character(len=*)      :: model_name
      integer               :: nDims
      integer               :: dims(nDims)
      integer               :: dimID(nDims)
      integer               :: k 
      character*50          :: dimName
      character*1           :: fd
      character(len=8)      :: date
      character(len=10)     :: time
      character(len=5)      :: zone
      integer, dimension(8) :: values

      n = 1
#if (defined USE_NETCDF3 || defined USE_NETCDF4)       
      call date_and_time(date,time,zone,values)
      call LDT_verify(nf90_def_dim(ftn,'ntiles',dims(1),&
           dimID(1)),&
           'nf90_def_dim failed for ntiles in LDT_rstProcMod')

      if ( nDims .gt. 2 ) then
       do k=2,nDims-1 
         write(unit=fd,fmt='(I1)') k-1
         dimName = 'dim'//trim(fd)
         call LDT_verify(nf90_def_dim(ftn,dimName,dims(k),&
              dimID(k)),&
              'nf90_def_dim failed for ntiles in LDT_rstProcMod')
       enddo
       call LDT_verify(nf90_def_dim(ftn,'time',1,&
            dimID(nDims)),&
            'nf90_def_dim failed for ntiles in LDT_rstProcMod')
      else   ! nDims = 2
       do k=2,nDims
         write(unit=fd,fmt='(I1)') k-1
         dimName = 'dim'//trim(fd)
         call LDT_verify(nf90_def_dim(ftn,dimName,dims(k),&
              dimID(k)),&
              'nf90_def_dim failed for ntiles in LDT_rstProcMod')
       enddo
      endif

      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"missing_value", &
           LDT_rc%udef),'nf90_put_att failed for missing_value')
      
      
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"title", &
           "LIS land surface model restart"),&
           'nf90_put_att failed for title')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"institution", &
           trim(LDT_rc%institution)),&
           'nf90_put_att failed for institution')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"source",&
           trim(model_name)),&
           'nf90_put_att failed for source')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"history", &
           "created on date: "//date(1:4)//"-"//date(5:6)//"-"//&
           date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10)),&
           'nf90_put_att failed for history')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"references", &
           "Kumar_etal_EMS_2006, Peters-Lidard_etal_ISSE_2007"),&
           'nf90_put_att failed for references')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"conventions", &
           "CF-1.6"),'nf90_put_att failed for conventions') !CF version 1.6
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"comment", &
           "website: http://lis.gsfc.nasa.gov/"),&
           'nf90_put_att failed for comment')
      
!grid information
      if(LDT_rc%lis_map_proj.eq."latlon") then !latlon
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
              "EQUIDISTANT CYLINDRICAL"),&
              'nf90_put_att failed for MAP_PROJECTION')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
              "SOUTH_WEST_CORNER_LAT", &
              LDT_rc%gridDesc(n,4)),&
              'nf90_put_att failed for SOUTH_WEST_CORNER_LAT')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
              "SOUTH_WEST_CORNER_LON", &
              LDT_rc%gridDesc(n,5)),&
              'nf90_put_att failed for SOUTH_WEST_CORNER_LON')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
              LDT_rc%gridDesc(n,9)),&
              'nf90_put_att failed for DX')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
              LDT_rc%gridDesc(n,10)),&
                  'nf90_put_att failed for DY')

      elseif(LDT_rc%lis_map_proj.eq."mercator") then 
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
              "MERCATOR"),&
              'nf90_put_att failed for MAP_PROJECTION')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
              "SOUTH_WEST_CORNER_LAT", &
              LDT_rc%gridDesc(n,4)),&
              'nf90_put_att failed for SOUTH_WEST_CORNER_LAT')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
              "SOUTH_WEST_CORNER_LON", &
              LDT_rc%gridDesc(n,5)),&
              'nf90_put_att failed for SOUTH_WEST_CORNER_LON') 
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
              LDT_rc%gridDesc(n,10)),&
              'nf90_put_att failed for TRUELAT1')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
              LDT_rc%gridDesc(n,11)),&
              'nf90_put_att failed for STANDARD_LON')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
              LDT_rc%gridDesc(n,8)),&
              'nf90_put_att failed for DX')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
              LDT_rc%gridDesc(n,9)),&
              'nf90_put_att failed for DY')

      elseif(LDT_rc%lis_map_proj.eq."lambert") then !lambert conformal
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
              "LAMBERT CONFORMAL"),&
              'nf90_put_att failed for MAP_PROJECTION')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
              "SOUTH_WEST_CORNER_LAT", &
              LDT_rc%gridDesc(n,4)),&
              'nf90_put_att failed for SOUTH_WEST_CORNER_LAT')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
              "SOUTH_WEST_CORNER_LON", &
              LDT_rc%gridDesc(n,5)),&
              'nf90_put_att failed for SOUTH_WEST_CORNER_LON')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
              LDT_rc%gridDesc(n,10)),&
              'nf90_put_att failed for TRUELAT1')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT2", &
              LDT_rc%gridDesc(n,7)),&
              'nf90_put_att failed for TRUELAT2')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
              LDT_rc%gridDesc(n,11)),&
              'nf90_put_att failed for STANDARD_LON')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
              LDT_rc%gridDesc(n,8)),&
              'nf90_put_att failed for DX')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
              LDT_rc%gridDesc(n,9)),&
              'nf90_put_att failed for DY')
         
      elseif(LDT_rc%lis_map_proj.eq."polar") then ! polar stereographic
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
              "POLAR STEREOGRAPHIC"),&
              'nf90_put_att failed for MAP_PROJECTION')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
              "SOUTH_WEST_CORNER_LAT", &
              LDT_rc%gridDesc(n,4)),&
              'nf90_put_att failed for SOUTH_WEST_CORNER_LAT')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
              "SOUTH_WEST_CORNER_LON", &
              LDT_rc%gridDesc(n,5)),&
              'nf90_put_att failed for SOUTH_WEST_CORNER_LON')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
              LDT_rc%gridDesc(n,10)),&
              'nf90_put_att failed for TRUELAT1')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"ORIENT", &
              LDT_rc%gridDesc(n,7)),&
              'nf90_put_att failed for ORIENT')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
              LDT_rc%gridDesc(n,11)),&
              'nf90_put_att failed for STANDARD_LON')                  
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
              LDT_rc%gridDesc(n,8)),&
              'nf90_put_att failed for DX')                  
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
              LDT_rc%gridDesc(n,9)),&
              'nf90_put_att failed for DY')                  
      endif
      
#endif
    end subroutine writeglobalheader

!BOP
! 
! !ROUTINE: writeheader_restart
! \label{writeheader_restart}
!
! !INTERFACE: 
    subroutine writeheader_restart(ftn,nvardims,&
         dimID,varId,&
         standardName,units,&
         scale_factor, offset,&
         vmin,vmax)
!
! !DESCRIPTION: 
!
!  This routine writes the NetCDF header information
!  in the climatological restart files.    
!
!   \begin{description}
!   \item [ftn]
!     unit number of the file being written
!   \item [nvardims]
!     number of dimensions of the variable
!   \item [dimID]
!     NetCDF dimension ID values 
!   \item [varID]
!     variable Id to be encoded in the NetCDF file
!   \item [standardName]
!     standard name of the variable
!   \item [units]
!     units of the variable
!   \item [scale\_factor]
!     scale factor of the variable
!   \item [offset]
!     offset value of the variable
!   \item [vmin]
!    minimum value of the variable
!   \item [vmax]
!    maximum value of the variable
!   \end{description}
!EOP

      integer            :: ftn
      integer            :: nvardims
      integer            :: dimID(nvardims)
      integer            :: varID
      character(len=*)   :: standardName
      character(len=*)   :: units
      real               :: scale_factor
      real               :: offset
      real               :: vmin
      real               :: vmax

      integer :: shuffle, deflate, deflate_level
      integer :: dimID_t(2)

      shuffle = NETCDF_shuffle
      deflate = NETCDF_deflate
      deflate_level =NETCDF_deflate_level

      dimID_t(1) = dimID(1)
      
      if(nvardims.gt.1) then 
         dimID_t(2) = dimID(2)

         call LDT_verify(nf90_def_var(ftn,trim(standardName),&
              nf90_float, dimids = dimID_t(1:2), varID=varID),&
              'nf90_def_var(2d) failed in LDT_ensRstMod')
         
#if(defined USE_NETCDF4)
         call LDT_verify(nf90_def_var_deflate(ftn,&
               varID,shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate(2d) failed in LDT_ensRstMod')
#endif
      else
         call LDT_verify(nf90_def_var(ftn,trim(standardName),&
              nf90_float,dimids = dimID_t(1:1), varID=varID),&
              'nf90_def_var(1d) failed in LDT_ensRstMod')
#if(defined USE_NETCDF4)                
         call LDT_verify(nf90_def_var_deflate(ftn,&
              varID, shuffle, deflate, deflate_level),&
              'nf90_def_var_deflate(1d) failed in LDT_ensRstMod')
#endif
      endif
      
      call LDT_verify(nf90_put_att(ftn,varID,&
           "units",trim(units)),&
           'nf90_put_att failed for units')
      call LDT_verify(nf90_put_att(ftn,varID,&
           "standard_name",trim(standardName)))
      call LDT_verify(nf90_put_att(ftn,varID,&
           "long_name",trim(standardName)),&
           'nf90_put_att failed for long_name')
      call LDT_verify(nf90_put_att(ftn,varID,&
           "scale_factor",scale_factor),&
           'nf90_put_att failed for scale_factor')
      call LDT_verify(nf90_put_att(ftn,varID,&
           "add_offset",offset),&
           'nf90_put_att failed for add_offset')
      call LDT_verify(nf90_put_att(ftn,varID,&
           "missing_value",LDT_rc%udef),&
           'nf90_put_att failed for missing_value')
      call LDT_verify(nf90_put_att(ftn,varID,&
           "_FillValue",LDT_rc%udef),&
           'nf90_put_att failed for _FillValue')
      call LDT_verify(nf90_put_att(ftn,varID,&
           "vmin",vmin),&
           'nf90_put_att failed for vmin')
      call LDT_verify(nf90_put_att(ftn,varID,&
           "vmax",vmax),&
           'nf90_put_att failed for vmax')
    end subroutine writeheader_restart

  end module LDT_rstProcMod
