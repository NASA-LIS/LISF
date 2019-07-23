!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !ROUTINE: readISMNObs
! \label{readISMNObs}
!
! !INTERFACE: 
subroutine readISMNObs(source)
! 
! !USES: 
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_timeMgrMod
  use LVT_logMod
  use ISMN_obsMod
  use map_utils

  implicit none
!
! !INPUT PARAMETERS: 
  integer,   intent(in) :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  16 May 2011: Sujay Kumar, Initial Specification
! 
!EOP
  character*400           :: cline,cline2
  integer                 :: iloc
  integer                 :: num_files
  integer                 :: status
  integer                 :: k,m,v,i,t,st,et,c,r
  logical                 :: file_exists
  integer                 :: ftn
  integer                 :: ios
  integer                 :: yr,mo,da,hr,mn,ss
  character*500           :: name1, name2
  real                    :: lat, lon, elev
  real                    :: depthfrom, depthto
  real                    :: sm_value
  real                    :: col,row
  integer                 :: stn_col,stn_row
  integer                 :: tind
  real,       allocatable     :: depth(:)
  real,       allocatable     :: sf_wt(:),rz_wt(:)
  real                    :: sfsm, rzsm
  character*1             :: sm_flag
  logical                 :: dataCheck
  type(ESMF_Time)         :: obstime, obstime1, obstime2
  real*8                  :: lis_prevtime
  integer                 :: doy
  real                    :: gmt
  real                    :: sfsm_grid(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: rzsm_grid(LVT_rc%lnc,LVT_rc%lnr)
  integer                 :: fnd


  if((LVT_rc%dyr(source).ne.ISMNobs(source)%yr).or.&
       LVT_rc%resetFlag(source)) then 
     LVT_rc%resetFlag(source) = .false. 
     deallocate(ISMNobs(source)%stn)
     
     ISMNobs(source)%yr = LVT_rc%dyr(source)
     
     call getNumberOfISMNfiles(ISMNobs(source)%odir, &
          LVT_rc%dyr(source), num_files)
     
     call generateISMNstationInfo(source, &
          ISMNobs(source)%odir, LVT_rc%dyr(source), num_files)
     
     call ESMF_TimeSet(ISMNobs(source)%startTime,&
          yy=LVT_rc%dyr(source),&
          mm=LVT_rc%dmo(source),&
          dd=1,&
          h=0,&
          m=0,&
          calendar=LVT_calendar,&
          rc=status)
     call LVT_verify(status,'error in setting ISMN start time')
     
     do k = 1, ISMNobs(source)%n_stns
        
        allocate(depth(ISMNobs(source)%stn(k)%vlevels))
        allocate(sf_wt(ISMNobs(source)%stn(k)%vlevels))
        allocate(rz_wt(ISMNobs(source)%stn(k)%vlevels))

        depth = LVT_rc%udef
        
        do v=1, ISMNobs(source)%stn(k)%vlevels
           
           inquire(file=trim(ISMNobs(source)%stn(k)%fname(v)),&
                exist=file_exists)

           if(file_exists) then 
              write(LVT_logunit,*) '[INFO] Reading ',&
                   trim(ISMNobs(source)%stn(k)%fname(v))
              ftn = LVT_getNextUnitNumber()
              open(ftn,file=trim(ISMNobs(source)%stn(k)%fname(v)),&
                   form='formatted')
              ios = 0 

              do while(ios.eq.0) 
                 read(ftn,'(a)',iostat=ios) cline
                 if(ios.ne.0.or.len(trim(adjustl(cline))).eq.0) exit

                 iloc = index(cline," ")
                 read(cline(1:iloc-1),&
                      '(I4.4,1X,I2.2,1X,I2.2)') yr,mo,da
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 read(cline(1:iloc-1),&
                      '(I2.2,1X,I2.2)') hr,mn
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 read(cline(1:iloc-1),*) lat
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 read(cline(1:iloc-1),*) lon
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 read(cline(1:iloc-1),*) elev
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 read(cline(1:iloc-1),*) depthfrom
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 read(cline(1:iloc-1),*) depthto
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 read(cline(1:iloc-1),*) sm_value
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 read(cline(1:iloc-1),*) sm_flag

                 call ESMF_TimeSet(obstime, yy=yr,mm=mo,dd=da,h=hr,m=mn,&
                      calendar=LVT_calendar,rc=status)
                 call LVT_verify(status, 'ESMF_TimeSet in readISMNObs')
                 
                 
                 tind = nint((obstime-ISMNobs(source)%startTime)/&
                      ISMNobs(source)%timestep)+1
                 ISMNobs(source)%stn(k)%lat = lat
                 ISMNobs(source)%stn(k)%lon = lon

                 if(tind.gt.0.and.tind.le.8784) then 
                    
                    !if we choose only good ('G') data, then we have no data. So including the
                    !unchecked ones as well. 
                    
                    if(sm_flag.eq.'G'.or.sm_flag.eq.'U') then 

                       ISMNobs(source)%stn(k)%sm(tind,v) = sm_value                 
                       depth(v) = (depthfrom+depthto)/2.0
!                       print*, yr,mo,da,hr,sm_value
                    endif
                 endif
              enddo
              call LVT_releaseUnitNumber(ftn)
           endif
        enddo

        call compute_vinterp_weights(&
             ISMNobs(source)%stn(k)%vlevels,LVT_rc%lis_sf_d,&
             LVT_rc%lis_rz_d,&
             depth,sf_wt,rz_wt)
        
        sfsm = 0 
        rzsm = 0 

        do tind = 1, ISMNobs(source)%nts
           dataCheck = .true. 
           do v=1,ISMNobs(source)%stn(k)%vlevels
              if(ISMNobs(source)%stn(k)%sm(tind,v).eq.LVT_rc%udef) then 
                 dataCheck = .false. 
              endif
           enddo
           if(dataCheck) then 
              sfsm = 0 
              rzsm = 0         
              do v=1,ISMNobs(source)%stn(k)%vlevels
                 if(ISMNobs(source)%stn(k)%sm(tind,v).ne.LVT_rc%udef) then 
                    sfsm = sfsm + sf_wt(v)*ISMNobs(source)%stn(k)%sm(tind,v)
                    rzsm = rzsm + rz_wt(v)*ISMNobs(source)%stn(k)%sm(tind,v)
                 endif
              enddo
              if(sfsm.gt.0.001) then 
                 ISMNobs(source)%stn(k)%sfsm(tind) = sfsm
              endif
              if(rzsm.gt.0.001) then 
                 ISMNobs(source)%stn(k)%rzsm(tind) = rzsm
              endif
           endif
        enddo
        deallocate(depth)
        deallocate(sf_wt)
        deallocate(rz_wt)
     end do
  end if

  
  call ESMF_TimeSet(obstime1,yy=LVT_rc%dyr(source), &
       mm=LVT_rc%dmo(source), dd=LVT_rc%dda(source), &
       h=LVT_rc%dhr(source), m=LVT_rc%dmn(source), &
       s = LVT_rc%dss(source), calendar=LVT_calendar, rc=status)
  call LVT_verify(status, 'obstime1 set failed')
  
  t = nint((obstime1-ISMNobs(source)%starttime)/ISMNobs(source)%timestep) + 1

  fnd = 0 

  sfsm_grid = LVT_rc%udef
  rzsm_grid = LVT_rc%udef

  do i=1,ISMNobs(source)%n_stns
     call latlon_to_ij(LVT_domain%lvtproj,ISMNobs(source)%stn(i)%lat, &
          ISMNobs(source)%stn(i)%lon,col,row)
     stn_col = nint(col)
     stn_row = nint(row)
     
     if(stn_col.ge.1.and.stn_col.le.LVT_rc%lnc.and.&
          stn_row.ge.1.and.stn_row.le.LVT_rc%lnr) then 
        if(ISMNobs(source)%stn(i)%sfsm(t).ne.LVT_rc%udef) then 
           sfsm_grid(stn_col,stn_row) =&
                ISMNobs(source)%stn(i)%sfsm(t)
           fnd =1
        endif
        
        if(ISMNobs(source)%stn(i)%rzsm(t).ne.LVT_rc%udef) then 
           rzsm_grid(stn_col,stn_row) = &
                ISMNobs(source)%stn(i)%rzsm(t)
        endif
     endif
  enddo
  
  call LVT_logSingleDataStreamVar(LVT_MOC_SOILMOIST,source,sfsm_grid,&
       vlevel=1,units="m3/m3")
  call LVT_logSingleDataStreamVar(LVT_MOC_ROOTMOIST,source,rzsm_grid,&
       vlevel=1,units="m3/m3")
  
end subroutine readISMNObs

!BOP
! 
! !ROUTINE: getNumberOfISMNfiles
! \label(getNumberOfISMNfiles)
!
! !INTERFACE:
subroutine getNumberOfISMNfiles(odir, yr, num_files)
! 
! !USES:   
  use LVT_logMod

  implicit none
!
! !INPUT PARAMETERS: 
  character(len=*), intent(in) :: odir
  integer  ,        intent(in) :: yr
  integer                      :: num_files
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  character*500                :: ls_comm
  character*500                :: cmd2
  character*4                  :: fyr
  integer                      :: ftn 

  write(unit=fyr, fmt='(i4.4)') yr  
  ls_comm = 'ls '//trim(odir)//'/'//trim(fyr)//'/*sm*stm > file_list'  
  cmd2 = 'wc -w file_list > file_list.wc'

  call system(ls_comm)
  call system(cmd2)

  ftn = LVT_getNextUnitNumber()
  open(ftn,file='file_list.wc',form='formatted',action='read')
  read(ftn,*) num_files
  call LVT_releaseUnitNumber(ftn)

end subroutine getNumberOfISMNfiles


subroutine generateISMNstationInfo(source, odir, yr, num_files)

  use LVT_coreMod
  use LVT_logMod
  use ISMN_obsMod

  implicit none

  integer  ,        intent(in) :: source
  character(len=*), intent(in) :: odir
  integer  ,        intent(in) :: yr
  integer  ,        intent(in) :: num_files


  real                       :: depthfrom
  real                       :: depthto
  
  character*500, allocatable     :: filenames(:)
  character*500              :: filename2
  character*500              :: name1, name2, name3
  character*500              :: stnname
  character*500              :: checkString

  
  integer                      :: iloc
  integer                      :: n_stns
  character*500                :: ls_comm
  character*500                :: cmd2
  character*4                  :: fyr
  integer                      :: k 
  integer                      :: ftn

  write(unit=fyr, fmt='(i4.4)') yr  
  ls_comm = 'ls '//trim(odir)//'/'//trim(fyr)//'/*sm*stm > file_list'  
  cmd2 = 'wc -w file_list > file_list.wc'

  call system(ls_comm)

  allocate(filenames(num_files))

  if(num_files.gt.0) then 
     ftn = LVT_getNextUnitNumber()
     do k=1,num_files
        open(ftn,file='file_list',form='formatted',action='read')
        read(ftn,'(a)') filenames(k)
     enddo
     call LVT_releaseUnitNumber(ftn)
  endif

  n_stns = 0 
  checkString = "dummy"

  do k=1, num_files
     filename2 = filenames(k)

     iloc = index(filename2,"_")
     read(filename2(1:iloc-1),'(a)') name1
     filename2 = filename2(iloc+1:len(filename2))

     iloc = index(filename2,"_")
     read(filename2(1:iloc-1),'(a)') name2
     filename2 = filename2(iloc+1:len(filename2))

     iloc = index(filename2,"_")
     read(filename2(1:iloc-1),'(a)') name3
     filename2 = filename2(iloc+1:len(filename2))

     stnname = trim(name1)//"_"//trim(name2)//"_"//trim(name3)
     if(checkString.ne.stnname) then 
        n_stns = n_stns + 1
        checkString = stnname
     endif
  enddo

  ISMNobs(source)%n_stns = n_stns

  allocate(ISMNobs(source)%stn(ISMNobs(source)%n_stns))

  do k=1,ISMNobs(source)%n_stns

     ISMNobs(source)%stn(k)%vlevels = 0 

  enddo

  n_stns = 0 
  checkString = "dummy"
  do k=1, num_files
     filename2 = filenames(k)
     
     iloc = index(filename2,"_")
     read(filename2(1:iloc-1),'(a)') name1
     filename2 = filename2(iloc+1:len(filename2))
     
     iloc = index(filename2,"_")
     read(filename2(1:iloc-1),'(a)') name2
     filename2 = filename2(iloc+1:len(filename2))
     
     iloc = index(filename2,"_")
     read(filename2(1:iloc-1),'(a)') name3
     filename2 = filename2(iloc+1:len(filename2))
     
     stnname = trim(name1)//"_"//trim(name2)//"_"//trim(name3)
     
     if(checkString.ne.stnname) then 
        n_stns = n_stns + 1
        checkString = stnname
     endif         
   
     ISMNobs(source)%stn(n_stns)%vlevels = ISMNobs(source)%stn(n_stns)%vlevels +1

  enddo

  do k=1,ISMNobs(source)%n_stns 
     allocate(ISMNobs(source)%stn(k)%fname(ISMNobs(source)%stn(k)%vlevels))
     allocate(ISMNobs(source)%stn(k)%sm(ISMNobs(source)%nts,ISMNobs(source)%stn(k)%vlevels))
     allocate(ISMNobs(source)%stn(k)%sfsm(ISMNobs(source)%nts))
     allocate(ISMNobs(source)%stn(k)%rzsm(ISMNobs(source)%nts))
     ISMNobs(source)%stn(k)%vlevels = 0 
     ISMNobs(source)%stn(k)%sm = LVT_rc%udef
     ISMNobs(source)%stn(k)%sfsm = LVT_rc%udef
     ISMNobs(source)%stn(k)%rzsm = LVT_rc%udef
  enddo

  n_stns = 0
  checkString = "dummy"
  do k=1, num_files
     filename2 = filenames(k)
     
     iloc = index(filename2,"_")
     read(filename2(1:iloc-1),'(a)') name1
     filename2 = filename2(iloc+1:len(filename2))
     
     iloc = index(filename2,"_")
     read(filename2(1:iloc-1),'(a)') name2
     filename2 = filename2(iloc+1:len(filename2))
     
     iloc = index(filename2,"_")
     read(filename2(1:iloc-1),'(a)') name3
     filename2 = filename2(iloc+1:len(filename2))
     
     stnname = trim(name1)//"_"//trim(name2)//"_"//trim(name3)
     
     if(checkString.ne.stnname) then 
        n_stns = n_stns + 1
        checkString = stnname
     endif         

     ISMNobs(source)%stn(n_stns)%vlevels = ISMNobs(source)%stn(n_stns)%vlevels +1
     ISMNobs(source)%stn(n_stns)%fname(ISMNobs(source)%stn(n_stns)%vlevels) = filenames(k)
  enddo

end subroutine generateISMNstationInfo
   
