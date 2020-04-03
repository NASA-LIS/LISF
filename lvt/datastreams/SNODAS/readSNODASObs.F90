!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !ROUTINE: readSNODASObs
! \label{readSNODASObs}
!
! !INTERFACE:
subroutine readSNODASObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_logMod
  use LVT_timeMgrMod
  use SNODAS_obsMod

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for SNODAS SWE and snowdepth data.
! LVT expects the data to be provided in a timestamped manner. The raw
! SNODAS binary files are read, and the data is interpolated to the
! LIS model output. The data is interpolated using the bilinear
! averaging technique.
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  6 May 2010: Sujay Kumar, Initial Specification
!  1 Aug 2019: David Mocko, Added/fixed log messages
! 
!EOP
  integer                :: source
  integer                :: yr, mo, da, hr
  integer                :: i,j,t,c,r
  integer                :: stn_col, stn_row
  real                   :: col,row
  character*100          :: snodas_swefilename
  character*100          :: snodas_snwdfilename
  logical                :: file_exists
  logical                :: readflag
  integer                :: ftn, ios
  integer                :: status
  type(ESMF_Time)        :: snodastime, snodastime1
  integer                :: stnindex,tind
  real                   :: offset
  integer                :: iret
  logical*1              :: lb(snodasobs(source)%nc*snodasobs(source)%nr)
  logical*1              :: lo(LVT_rc%lnc,LVT_rc%lnr)
  integer*2              :: var(snodasobs(source)%nc, snodasobs(source)%nr)
  real                   :: var1d(snodasobs(source)%nc*snodasobs(source)%nr)
  real                   :: var_ip(LVT_rc%lnc,LVT_rc%lnr)
  
  var_ip  = LVT_rc%udef

  call ESMF_TimeSet(snodastime1,yy=LVT_rc%dyr(source), &
       mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), h=LVT_rc%dhr(source), &
       m=LVT_rc%dmn(source), s=LVT_rc%dss(source), calendar=LVT_calendar, &
       rc=status)
  call LVT_verify(status, 'snodastime1 set failed')

  offset = (snodastime1-snodasobs(source)%starttime)/snodasobs(source)%timestep

  if ((nint(offset)-offset).eq.0) then !only when LIS time matches the observation
  
     call create_SNODAS_swefilename(snodasobs(source)%odir, &
          LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source), &
          snodas_swefilename)
     
     inquire(file=trim(snodas_swefilename),exist=file_exists)

     if (file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading SNODAS SWE data ', &
                             trim(snodas_swefilename)
        
        ftn=LVT_getNextUnitNumber()
        open(ftn,file=trim(snodas_swefilename),form='unformatted', &
             access='direct',convert='big_endian',                 &
             recl=snodasobs(source)%nc*snodasobs(source)%nr*2)
        
        read(ftn,rec=1) var
        lb = .false. 
        do r=1,snodasobs(source)%nr
           do c=1,snodasobs(source)%nc
              var1d(c+(r-1)*snodasobs(source)%nc) = var(c,(snodasobs(source)%nr-r+1))
              if(var(c,(snodasobs(source)%nr-r+1)).ge.0) then 
                 var1d(c+(r-1)*snodasobs(source)%nc) = var1d(c+(r-1)*snodasobs(source)%nc)/1000.0
                 lb(c+(r-1)*snodasobs(source)%nc) = .true. 
              endif
           enddo
        enddo
        
        call bilinear_interp(LVT_rc%gridDesc,lb,var1d,lo,var_ip,&
             snodasobs(source)%nc*snodasobs(source)%nr,LVT_rc%lnc*LVT_rc%lnr,&
             snodasobs(source)%rlat,snodasobs(source)%rlon,&
             snodasobs(source)%w11, snodasobs(source)%w12, &
             snodasobs(source)%w21, snodasobs(source)%w22, &
             snodasobs(source)%n11, snodasobs(source)%n12, &
             snodasobs(source)%n21, snodasobs(source)%n22, LVT_rc%udef,iret)
        readflag = .false. 
        call LVT_releaseUnitNumber(ftn)

        write(LVT_logunit,*) '[INFO] Successfully processed SWE from ', &
                             trim(snodas_swefilename)
     else
        write(LVT_logunit,*) '[WARN] Could not find SNODAS SWE file ', &
                             trim(snodas_swefilename)
     endif
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_SWE,source,var_ip,vlevel=1,units="m")

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(var_ip(c,r).ne.-9999.0) then 
           var_ip(c,r) = var_ip(c,r)*1000.0
        endif
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_SWE,source,var_ip,vlevel=1,units="kg/m2")

!process snow depth data
  var_ip  = LVT_rc%udef

  call ESMF_TimeSet(snodastime1,yy=LVT_rc%dyr(source), &
       mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), h=LVT_rc%dhr(source), &
       m=LVT_rc%dmn(source), s=LVT_rc%dss(source), calendar=LVT_calendar, &
       rc=status)
  call LVT_verify(status, 'snodastime1 set failed')

  offset = (snodastime1-snodasobs(source)%starttime)/snodasobs(source)%timestep

  if ((nint(offset)-offset).eq.0) then !only when LIS time matches the observation
  
     call create_SNODAS_snwdfilename(snodasobs(source)%odir, &
          LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source), &
          snodas_snwdfilename)
     
     inquire(file=trim(snodas_snwdfilename),exist=file_exists)

     if (file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading SNODAS SNWD data ', &
                             trim(snodas_snwdfilename)
        
        ftn=LVT_getNextUnitNumber()
        open(ftn,file=trim(snodas_snwdfilename),form='unformatted', &
             access='direct',convert='big_endian',                  &
             recl=snodasobs(source)%nc*snodasobs(source)%nr*2)
        
        read(ftn,rec=1) var
        lb = .false. 
        do r=1,snodasobs(source)%nr
           do c=1,snodasobs(source)%nc
              var1d(c+(r-1)*snodasobs(source)%nc) = var(c,(snodasobs(source)%nr-r+1))
              if(var(c,(snodasobs(source)%nr-r+1)).ge.0) then 
                 var1d(c+(r-1)*snodasobs(source)%nc) = var1d(c+(r-1)*snodasobs(source)%nc)/1000.0
                 lb(c+(r-1)*snodasobs(source)%nc) = .true. 
              endif
           enddo
        enddo
        
        call bilinear_interp(LVT_rc%gridDesc,lb,var1d,lo,var_ip,&
             snodasobs(source)%nc*snodasobs(source)%nr,LVT_rc%lnc*LVT_rc%lnr,&
             snodasobs(source)%rlat,snodasobs(source)%rlon,&
             snodasobs(source)%w11, snodasobs(source)%w12, &
             snodasobs(source)%w21, snodasobs(source)%w22, &
             snodasobs(source)%n11, snodasobs(source)%n12, &
             snodasobs(source)%n21, snodasobs(source)%n22, LVT_rc%udef,iret)
        readflag = .false. 
        call LVT_releaseUnitNumber(ftn)

        write(LVT_logunit,*) '[INFO] Successfully processed SNWD from ', &
                             trim(snodas_snwdfilename)
     else
        write(LVT_logunit,*) '[WARN] Could not find SNODAS SNWD file ', &
                             trim(snodas_snwdfilename)
     endif
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_SNOWDEPTH,source,var_ip,vlevel=1,units="m")

end subroutine readSNODASObs

!BOP
! 
! !ROUTINE: create_SNODAS_swefilename
! \label{create_SNODAS_swefilename}
!
! !INTERFACE:
subroutine create_SNODAS_swefilename(odir, yr,mo,da, snodasname)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
  character(len=*), intent(in)  :: odir
  integer,          intent(in)  :: yr
  integer,          intent(in)  :: mo
  integer,          intent(in)  :: da
! 
! !OUTPUT PARAMETERS:
  character(len=*), intent(out) :: snodasname
!
! !DESCRIPTION: 
! 
! This routine creates a timestamped SNODAS filename.
! The filename convention is as follows: 
!   us : region (united states)
!   ssm : simple snow model
!   v1  : operational snow model output
!   1034 : Snow water equivalent
!   ts__ : integral through all layers of the snow pack
!   T0001 : A 1-hour snapshot
! 
!  The arguments are: 
!  \begin{description}
!   \item[odir] SNODAS base directory
!   \item[yr]   year of data 
!   \item[mo]   month of data
!   \item[da]   day of data
!   \item[snodasname]  Name of the SNODAS file  
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

  character*4             :: fyr
  character*2             :: fmo
  character*2             :: fda
  
  write(fyr, '(i4.4)' ) yr
  write(fmo, '(i2.2)' ) mo
  write(fda, '(i2.2)' ) da

  snodasname = trim(odir)//'/'//trim(fyr)//trim(fmo)//trim(fda)// &
       '/us_ssmv11034tS__T0001TTNATS' &
       //trim(fyr)//trim(fmo)//trim(fda)//'05HP001.dat'
  
end subroutine create_SNODAS_swefilename

!BOP
! 
! !ROUTINE: create_SNODAS_snwdfilename
! \label{create_SNODAS_snwdfilename}
!
! !INTERFACE:
subroutine create_SNODAS_snwdfilename(odir, yr,mo,da, snodasname)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
  character(len=*), intent(in)  :: odir
  integer,          intent(in)  :: yr
  integer,          intent(in)  :: mo
  integer,          intent(in)  :: da
! 
! !OUTPUT PARAMETERS:
  character(len=*), intent(out) :: snodasname
!
! !DESCRIPTION: 
! This routine creates a timestamped SNODAS filename.
! The filename convention is as follows: 
!   us : region (united states)
!   ssm : simple snow model
!   v1  : operational snow model output
!   1034 : Snow water equivalent
!   ts__ : integral through all layers of the snow pack
!   T0001 : A 1-hour snapshot
! 
!  The arguments are: 
!  \begin{description}
!   \item[odir] SNODAS base directory
!   \item[yr]   year of data 
!   \item[mo]   month of data
!   \item[da]   day of data
!   \item[snodasname]  Name of the SNODAS file  
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

  character*4             :: fyr
  character*2             :: fmo
  character*2             :: fda
  
  write(fyr, '(i4.4)' ) yr
  write(fmo, '(i2.2)' ) mo
  write(fda, '(i2.2)' ) da

  snodasname = trim(odir)//'/'//trim(fyr)//trim(fmo)//trim(fda)// &
       '/us_ssmv11036tS__T0001TTNATS' &
       //trim(fyr)//trim(fmo)//trim(fda)//'05HP001.dat'
  
end subroutine create_SNODAS_snwdfilename
